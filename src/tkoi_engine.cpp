// tkoi_engine.cpp
//
// Native engine for tKOI. Everything performance-critical in run_tkoi()
// lives here:
//
//   * tkoi_build_csr()    undirected graph -> symmetric normalized CSR matrix
//   * tkoi_ppr_null()     batched personalized PageRank (observed + null)
//   * tkoi_seed_reach()   number of seeds within 1 and 2 hops of every node
//   * tkoi_sample_null()  degree-matched null seed sampling (R's RNG)
//   * tkoi_total_memory() physical memory, used to size work buffers
//
// Personalized PageRank on an undirected graph solves
//
//     (I - d W) z = u,   W = A D^-1,   x = z / sum(z)
//
// which, with y = D^-1/2 z, is the symmetric positive definite system
//
//     (I - d S) y = D^-1/2 u,   S = D^-1/2 A D^-1/2.
//
// Its condition number is at most (1 + d) / (1 - d), so conjugate gradient
// converges in a few dozen sparse products. Many right-hand sides are solved
// together ("blocks"), so every pass over the 10M-edge matrix serves up to
// 16 PageRank vectors at once.
//
// Conventions match igraph::page_rank(algo = "prpack", directed = FALSE):
// undirected self-loops contribute 2 * weight to the degree, multi-edges add
// up, and probability mass that reaches a vertex without edges returns to the
// personalization vector (handled by the final normalization).
//
// Determinism: rows are split into a fixed number of chunks that does not
// depend on the thread count, and every reduction is summed chunk by chunk
// in a fixed order. Results are bit-identical for any n_cores.

#include <Rcpp.h>
#include <R_ext/Random.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <initializer_list>
#include <mutex>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace Rcpp;

namespace {

// ---------------------------------------------------------------------------
// Threading
// ---------------------------------------------------------------------------

// Run fn(chunk) for every chunk in [0, n_chunks) on up to n_threads threads.
// Chunks are handed out dynamically, so faster cores take on more work. If
// the OS refuses to create a thread, the calling thread does the remaining
// work itself instead of failing.
void parallel_chunks(int n_chunks, int n_threads, const std::function<void(int)>& fn) {
  if (n_threads <= 1 || n_chunks <= 1) {
    for (int c = 0; c < n_chunks; ++c) fn(c);
    return;
  }
  const int nt = std::min(n_threads, n_chunks);
  std::atomic<int> next(0);
  std::exception_ptr error = nullptr;
  std::mutex error_mutex;
  auto worker = [&]() {
    try {
      for (int c = next++; c < n_chunks; c = next++) fn(c);
    } catch (...) {
      std::lock_guard<std::mutex> lock(error_mutex);
      if (!error) error = std::current_exception();
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nt - 1);
  for (int t = 0; t < nt - 1; ++t) {
    try {
      pool.emplace_back(worker);
    } catch (...) {
      break;  // could not spawn: the calling thread picks up the slack
    }
  }
  worker();
  for (auto& th : pool) th.join();
  if (error) std::rethrow_exception(error);
}

// Split rows [0, n) into at most max_chunks contiguous chunks that carry
// roughly equal numbers of stored entries (plus rows, for sparse regions).
std::vector<int> make_chunks(const int* row_ptr, int n, int max_chunks) {
  std::vector<int> bounds;
  bounds.push_back(0);
  if (n <= 0) {
    bounds.push_back(0);
    return bounds;
  }
  const int n_chunks = std::max(1, std::min(max_chunks, n));
  const double total = static_cast<double>(row_ptr[n]) + n;
  for (int ch = 1; ch < n_chunks; ++ch) {
    const double target = total * ch / n_chunks;
    int lo = bounds.back();
    int hi = n;
    while (lo < hi) {
      const int mid = lo + (hi - lo) / 2;
      if (static_cast<double>(row_ptr[mid]) + mid < target) lo = mid + 1;
      else hi = mid;
    }
    if (lo > bounds.back() && lo < n) bounds.push_back(lo);
  }
  bounds.push_back(n);
  return bounds;
}

const int kChunks = 1024;

// ---------------------------------------------------------------------------
// Block conjugate gradient
// ---------------------------------------------------------------------------

struct Csr {
  const int* row_ptr;
  const int* col;
  const double* val;
  const double* scale;
  int n;
};

struct Seed {
  int node;
  double weight;
};

// Q = P - d * S * P for rows [i0, i1); also returns the chunk's P.Q partials.
template <int BB>
inline void spmm_chunk(const Csr& g, const double* P, double* Q, double d,
                       int i0, int i1, double* pq_out) {
  double pq[BB];
  for (int c = 0; c < BB; ++c) pq[c] = 0.0;
  for (int i = i0; i < i1; ++i) {
    double acc[BB];
    for (int c = 0; c < BB; ++c) acc[c] = 0.0;
    for (int k = g.row_ptr[i]; k < g.row_ptr[i + 1]; ++k) {
      const double v = g.val[k];
      const double* pj = P + static_cast<std::size_t>(g.col[k]) * BB;
      for (int c = 0; c < BB; ++c) acc[c] += v * pj[c];
    }
    const double* pi = P + static_cast<std::size_t>(i) * BB;
    double* qi = Q + static_cast<std::size_t>(i) * BB;
    for (int c = 0; c < BB; ++c) {
      qi[c] = pi[c] - d * acc[c];
      pq[c] += pi[c] * qi[c];
    }
  }
  for (int c = 0; c < BB; ++c) pq_out[c] = pq[c];
}

struct BlockResult {
  std::vector<int> iterations;
  std::vector<double> rel_residual;
};

// Solve up to BB right-hand sides at once. Column c of the block uses
// seeds[c]; unused padding columns stay zero. On return, X holds the
// normalized PageRank vectors in row-major (n x BB) layout.
//
// A column stops when both relative residuals are at most `tol`:
//   * the 2-norm residual of the symmetric system, which CG minimizes, and
//   * the 1-norm residual of the original PageRank system, ||D^1/2 r||_1 /
//     ||u||_1, which bounds the 1-norm error of the normalized PageRank
//     vector by about 2 * tol however widely degrees and weights vary.
template <int BB>
BlockResult solve_block(const Csr& g, const std::vector<std::vector<Seed>>& seeds,
                        double d, double tol, int max_iter, int n_threads,
                        const std::vector<int>& chunks, std::vector<double>& X) {
  const int n = g.n;
  const std::size_t len = static_cast<std::size_t>(n) * BB;
  const int n_chunks = static_cast<int>(chunks.size()) - 1;
  const int n_cols = static_cast<int>(seeds.size());

  X.assign(len, 0.0);
  std::vector<double> R(len, 0.0), P(len, 0.0), Q(len, 0.0);
  std::vector<double> part(static_cast<std::size_t>(n_chunks) * 2 * BB, 0.0);
  std::vector<double> rr(BB, 0.0), l1(BB, 0.0), b_norm(BB, 0.0), u_norm(BB, 0.0);
  std::vector<double> alpha(BB, 0.0), beta(BB, 0.0);
  std::vector<char> active(BB, 0);
  std::vector<int> column_iterations(BB, 0);

  // Right-hand side b = D^-1/2 u. Seeds are applied by assignment, so a
  // vertex listed twice keeps its last weight (same as prob[idx] = w in R).
  for (int c = 0; c < n_cols; ++c) {
    for (const Seed& s : seeds[c]) {
      R[static_cast<std::size_t>(s.node) * BB + c] = g.scale[s.node] * s.weight;
    }
  }
  P = R;

  // Column sums of two per-row quantities, reduced chunk by chunk in a fixed
  // order so the result does not depend on the thread count.
  auto column_sums = [&](const std::function<void(int, double*, double*)>& partial,
                         double* out1, double* out2) {
    parallel_chunks(n_chunks, n_threads, [&](int ch) {
      double* slot = &part[static_cast<std::size_t>(ch) * 2 * BB];
      for (int c = 0; c < 2 * BB; ++c) slot[c] = 0.0;
      partial(ch, slot, slot + BB);
    });
    for (int c = 0; c < BB; ++c) {
      double s1 = 0.0, s2 = 0.0;
      for (int ch = 0; ch < n_chunks; ++ch) {
        s1 += part[static_cast<std::size_t>(ch) * 2 * BB + c];
        s2 += part[static_cast<std::size_t>(ch) * 2 * BB + BB + c];
      }
      out1[c] = s1;
      if (out2) out2[c] = s2;
    }
  };

  // ||r||_2^2 and ||D^1/2 r||_1 of the current residual.
  auto residual_norms = [&](int ch, double* sq, double* abs1) {
    for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
      const double* ri = &R[static_cast<std::size_t>(i) * BB];
      const double inv_scale = 1.0 / g.scale[i];
      for (int c = 0; c < BB; ++c) {
        sq[c] += ri[c] * ri[c];
        abs1[c] += std::fabs(ri[c]) * inv_scale;
      }
    }
  };

  column_sums(residual_norms, rr.data(), u_norm.data());
  for (int c = 0; c < n_cols; ++c) {
    if (!(rr[c] > 0.0)) stop("A PageRank vector has no seed with positive weight.");
  }
  auto converged = [&](int c) {
    return std::sqrt(rr[c]) <= tol * b_norm[c] && l1[c] <= tol * u_norm[c];
  };
  for (int c = 0; c < BB; ++c) {
    b_norm[c] = std::sqrt(rr[c]);
    l1[c] = u_norm[c];
    active[c] = (c < n_cols && rr[c] > 0.0) ? 1 : 0;
  }

  int it = 0;
  for (; it < max_iter; ++it) {
    bool any_active = false;
    for (int c = 0; c < BB; ++c) any_active = any_active || active[c];
    if (!any_active) break;

    // q = M p and p.q
    column_sums([&](int ch, double* pq, double*) {
      spmm_chunk<BB>(g, P.data(), Q.data(), d, chunks[ch], chunks[ch + 1], pq);
    }, alpha.data(), nullptr);
    for (int c = 0; c < BB; ++c) {
      const double pq = alpha[c];
      alpha[c] = (active[c] && pq > 0.0) ? rr[c] / pq : 0.0;
    }

    // x += alpha p, r -= alpha q, then the residual norms.
    std::vector<double> rr_new(BB, 0.0);
    column_sums([&](int ch, double* sq, double* abs1) {
      for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
        const std::size_t o = static_cast<std::size_t>(i) * BB;
        const double inv_scale = 1.0 / g.scale[i];
        for (int c = 0; c < BB; ++c) {
          X[o + c] += alpha[c] * P[o + c];
          R[o + c] -= alpha[c] * Q[o + c];
          sq[c] += R[o + c] * R[o + c];
          abs1[c] += std::fabs(R[o + c]) * inv_scale;
        }
      }
    }, rr_new.data(), l1.data());

    for (int c = 0; c < BB; ++c) {
      if (!active[c]) {
        beta[c] = 0.0;
        continue;
      }
      beta[c] = rr_new[c] / rr[c];
      rr[c] = rr_new[c];
      column_iterations[c] = it + 1;
      if (converged(c)) active[c] = 0;
    }

    // p = r + beta p
    parallel_chunks(n_chunks, n_threads, [&](int ch) {
      for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
        const std::size_t o = static_cast<std::size_t>(i) * BB;
        for (int c = 0; c < BB; ++c) P[o + c] = R[o + c] + beta[c] * P[o + c];
      }
    });

    Rcpp::checkUserInterrupt();
  }

  BlockResult res;
  res.iterations.assign(column_iterations.begin(), column_iterations.begin() + n_cols);
  res.rel_residual.assign(n_cols, 0.0);
  for (int c = 0; c < n_cols; ++c) {
    const double r2 = b_norm[c] > 0.0 ? std::sqrt(rr[c]) / b_norm[c] : 0.0;
    const double r1 = u_norm[c] > 0.0 ? l1[c] / u_norm[c] : 0.0;
    res.rel_residual[c] = std::max(r1, r2);
  }

  // Back-transform z = D^1/2 y, then normalize every column to sum to one.
  std::vector<double> totals(BB, 0.0);
  column_sums([&](int ch, double* acc, double*) {
    for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
      const std::size_t o = static_cast<std::size_t>(i) * BB;
      for (int c = 0; c < BB; ++c) {
        X[o + c] /= g.scale[i];
        acc[c] += X[o + c];
      }
    }
  }, totals.data(), nullptr);
  parallel_chunks(n_chunks, n_threads, [&](int ch) {
    for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
      const std::size_t o = static_cast<std::size_t>(i) * BB;
      for (int c = 0; c < BB; ++c) {
        if (totals[c] != 0.0) X[o + c] /= totals[c];
      }
    }
  });
  return res;
}

}  // namespace

// ---------------------------------------------------------------------------
// Exported functions
// ---------------------------------------------------------------------------

// [[Rcpp::export(.tkoi_build_csr)]]
List tkoi_build_csr(IntegerVector from, IntegerVector to, NumericVector weight, int n) {
  const R_xlen_t m = from.size();
  if (to.size() != m) stop("'from' and 'to' must have the same length.");
  if (n < 0) stop("'n' must be non-negative.");
  const bool weighted = weight.size() > 0;
  if (weighted && weight.size() != m) stop("'weight' must be empty or match the number of edges.");

  std::vector<double> degree(n, 0.0);
  std::vector<std::int64_t> count(n, 0);
  for (R_xlen_t e = 0; e < m; ++e) {
    if (from[e] == NA_INTEGER || to[e] == NA_INTEGER) stop("Edge endpoint out of range.");
    const int a = from[e] - 1;
    const int b = to[e] - 1;
    if (a < 0 || a >= n || b < 0 || b >= n) stop("Edge endpoint out of range.");
    const double w = weighted ? weight[e] : 1.0;
    if (!std::isfinite(w) || w < 0.0) stop("Edge weights must be finite and non-negative.");
    if (a == b) {
      degree[a] += 2.0 * w;
      count[a] += 1;
    } else {
      degree[a] += w;
      degree[b] += w;
      count[a] += 1;
      count[b] += 1;
    }
  }

  std::int64_t nnz = 0;
  for (int i = 0; i < n; ++i) nnz += count[i];
  if (nnz > 2147483647LL) stop("The network is too large for the native engine (more than 2^31 matrix entries).");

  IntegerVector row_ptr(n + 1);
  row_ptr[0] = 0;
  for (int i = 0; i < n; ++i) row_ptr[i + 1] = row_ptr[i] + static_cast<int>(count[i]);
  IntegerVector col(static_cast<R_xlen_t>(nnz));
  NumericVector val(static_cast<R_xlen_t>(nnz));
  NumericVector scale(n);
  for (int i = 0; i < n; ++i) scale[i] = degree[i] > 0.0 ? 1.0 / std::sqrt(degree[i]) : 1.0;

  std::vector<int> pos(row_ptr.begin(), row_ptr.end() - 1);
  for (R_xlen_t e = 0; e < m; ++e) {
    const int a = from[e] - 1;
    const int b = to[e] - 1;
    const double w = weighted ? weight[e] : 1.0;
    if (a == b) {
      const int p = pos[a]++;
      col[p] = a;
      val[p] = 2.0 * w * scale[a] * scale[a];
    } else {
      const double v = w * scale[a] * scale[b];
      int p = pos[a]++;
      col[p] = b;
      val[p] = v;
      p = pos[b]++;
      col[p] = a;
      val[p] = v;
    }
  }
  return List::create(_["row_ptr"] = row_ptr, _["col"] = col, _["val"] = val,
                      _["scale"] = scale, _["n"] = n);
}

// Personalized PageRank for the observed seeds (column 0) and for every
// null seed set (columns 1..n_perm). Null set p uses the vertices in column
// p of `null_node` (1-based) with the weights in `null_weight`. Returns the
// observed vector, optionally every null vector, and the null mean and
// sample standard deviation per vertex.
//
// [[Rcpp::export(.tkoi_ppr_null)]]
List tkoi_ppr_null(List csr, IntegerVector observed_node, NumericVector observed_weight,
                   IntegerMatrix null_node, NumericVector null_weight, double damping,
                   double tol, int max_iter, int n_threads, int max_block, bool keep_null,
                   Nullable<Function> progress = R_NilValue) {
  IntegerVector row_ptr = csr["row_ptr"];
  IntegerVector col = csr["col"];
  NumericVector val = csr["val"];
  NumericVector scale = csr["scale"];
  const int n = scale.size();
  const int n_perm = null_node.ncol();
  const int n_null_seeds = null_node.nrow();
  const int K = 1 + n_perm;

  if (!(damping > 0.0 && damping < 1.0)) stop("'damping' must be in (0, 1).");
  if (!(tol > 0.0)) stop("'tol' must be positive.");
  if (max_iter < 1) stop("'max_iter' must be a positive integer.");
  if (observed_node.size() != observed_weight.size()) stop("Observed seed vectors differ in length.");
  if (null_weight.size() != n_null_seeds) stop("'null_weight' must have one weight per row of 'null_node'.");

  auto check_nodes = [&](const int* nodes, R_xlen_t count) {
    for (R_xlen_t j = 0; j < count; ++j) {
      if (nodes[j] == NA_INTEGER || nodes[j] < 1 || nodes[j] > n) stop("Seed vertex out of range.");
    }
  };
  auto check_weights = [&](const NumericVector& weights) {
    double total = 0.0;
    for (R_xlen_t j = 0; j < weights.size(); ++j) {
      if (!std::isfinite(weights[j]) || weights[j] < 0.0) stop("Seed weights must be finite and non-negative.");
      total += weights[j];
    }
    if (!(total > 0.0)) stop("Seed weights must have a positive total.");
  };
  check_nodes(observed_node.begin(), observed_node.size());
  check_weights(observed_weight);
  if (n_perm > 0) {
    check_nodes(null_node.begin(), static_cast<R_xlen_t>(n_null_seeds) * n_perm);
    check_weights(null_weight);
  }

  Csr g{row_ptr.begin(), col.begin(), val.begin(), scale.begin(), n};
  const std::vector<int> chunks = make_chunks(g.row_ptr, n, kChunks);
  const int n_chunks = static_cast<int>(chunks.size()) - 1;

  auto seeds_of = [&](int k) {
    std::vector<Seed> out;
    if (k == 0) {
      out.reserve(observed_node.size());
      for (R_xlen_t j = 0; j < observed_node.size(); ++j) out.push_back({observed_node[j] - 1, observed_weight[j]});
    } else {
      out.reserve(n_null_seeds);
      for (int j = 0; j < n_null_seeds; ++j) out.push_back({null_node(j, k - 1) - 1, null_weight[j]});
    }
    return out;
  };

  NumericVector observed(n);
  List null_columns(keep_null ? n_perm : 0);
  NumericVector null_mean(n, 0.0), null_m2(n, 0.0);
  IntegerVector iterations(K);
  NumericVector rel_residual(K);

  std::vector<double> X;
  int done = 0;
  while (done < K) {
    // Full blocks use the widest allowed width. The last block is padded up
    // to the next supported width: a padded column costs far less than an
    // extra pass over the matrix with a narrower block.
    const int remaining = K - done;
    int width = 1;
    for (int w : {16, 8, 4, 2, 1}) {
      if (w <= max_block) {
        width = w;
        break;
      }
    }
    if (remaining < width) {
      int padded = 1;
      while (padded < remaining) padded *= 2;
      width = std::min(width, padded);
    }
    const int used = std::min(width, remaining);
    std::vector<std::vector<Seed>> block_seeds;
    block_seeds.reserve(used);
    for (int c = 0; c < used; ++c) block_seeds.push_back(seeds_of(done + c));

    BlockResult res;
    switch (width) {
      case 16: res = solve_block<16>(g, block_seeds, damping, tol, max_iter, n_threads, chunks, X); break;
      case 8:  res = solve_block<8>(g, block_seeds, damping, tol, max_iter, n_threads, chunks, X); break;
      case 4:  res = solve_block<4>(g, block_seeds, damping, tol, max_iter, n_threads, chunks, X); break;
      case 2:  res = solve_block<2>(g, block_seeds, damping, tol, max_iter, n_threads, chunks, X); break;
      default: res = solve_block<1>(g, block_seeds, damping, tol, max_iter, n_threads, chunks, X); break;
    }
    const int stride = width;

    // Copy results out; allocate R vectors on the main thread only.
    std::vector<double*> dest(used, nullptr);
    for (int c = 0; c < used; ++c) {
      const int k = done + c;
      iterations[k] = res.iterations[c];
      rel_residual[k] = res.rel_residual[c];
      if (k == 0) {
        dest[c] = observed.begin();
      } else if (keep_null) {
        NumericVector column(n);
        null_columns[k - 1] = column;
        dest[c] = column.begin();
      }
    }

    // Null statistics: Welford updates one null column at a time, in column
    // order, so the result does not depend on how columns fall into blocks.
    const int first_null = (done == 0) ? 1 : 0;
    const int nulls_before = (done == 0) ? 0 : done - 1;
    parallel_chunks(n_chunks, n_threads, [&](int ch) {
      for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
        const double* xi = &X[static_cast<std::size_t>(i) * stride];
        for (int c = 0; c < used; ++c) {
          if (dest[c]) dest[c][i] = xi[c];
        }
        double mean = null_mean[i];
        double m2 = null_m2[i];
        for (int c = first_null; c < used; ++c) {
          const double count = nulls_before + (c - first_null) + 1.0;
          const double delta = xi[c] - mean;
          mean += delta / count;
          m2 += delta * (xi[c] - mean);
        }
        null_mean[i] = mean;
        null_m2[i] = m2;
      }
    });
    done += used;

    if (progress.isNotNull()) {
      Function f(progress.get());
      f(done, K);
    }
  }

  NumericVector null_sd(n);
  for (int i = 0; i < n; ++i) {
    null_sd[i] = n_perm > 1 ? std::sqrt(null_m2[i] / (n_perm - 1.0)) : NA_REAL;
    if (n_perm == 0) null_mean[i] = NA_REAL;
  }

  return List::create(_["observed"] = observed, _["null_columns"] = null_columns,
                      _["null_mean"] = null_mean, _["null_sd"] = null_sd,
                      _["iterations"] = iterations, _["rel_residual"] = rel_residual);
}

// For every vertex, count how many seeds lie within 1 hop (closed
// neighborhood) and within 2 hops. Equivalent to tabulating
// igraph::ego(order = 1 / 2, mode = "all") over the seeds: a seed listed
// twice is counted twice. Seeds are processed 64 at a time as bitsets.
//
// [[Rcpp::export(.tkoi_seed_reach)]]
List tkoi_seed_reach(IntegerVector row_ptr, IntegerVector col, IntegerVector seeds, int n_threads) {
  const int n = row_ptr.size() - 1;
  const int n_seeds = seeds.size();
  for (int j = 0; j < n_seeds; ++j) {
    if (seeds[j] == NA_INTEGER || seeds[j] < 1 || seeds[j] > n) stop("Seed vertex out of range.");
  }
  const int* rp = row_ptr.begin();
  const int* ci = col.begin();
  const std::vector<int> chunks = make_chunks(rp, n, kChunks);
  const int n_chunks = static_cast<int>(chunks.size()) - 1;

  const int kWords = 8;  // 512 seeds per pass over the matrix
  std::vector<int> direct(n, 0), indirect(n, 0);
  std::vector<std::uint64_t> b1, b2;
  for (int start = 0; start < n_seeds; start += 64 * kWords) {
    const int batch = std::min(64 * kWords, n_seeds - start);
    const int words = (batch + 63) / 64;
    b1.assign(static_cast<std::size_t>(n) * words, 0);
    b2.assign(static_cast<std::size_t>(n) * words, 0);
    for (int j = 0; j < batch; ++j) {
      const int s = seeds[start + j] - 1;
      const int w = j / 64;
      const std::uint64_t bit = std::uint64_t(1) << (j % 64);
      b1[static_cast<std::size_t>(s) * words + w] |= bit;
      for (int k = rp[s]; k < rp[s + 1]; ++k) b1[static_cast<std::size_t>(ci[k]) * words + w] |= bit;
    }
    parallel_chunks(n_chunks, n_threads, [&](int ch) {
      for (int i = chunks[ch]; i < chunks[ch + 1]; ++i) {
        const std::size_t oi = static_cast<std::size_t>(i) * words;
        for (int w = 0; w < words; ++w) b2[oi + w] = b1[oi + w];
        for (int k = rp[i]; k < rp[i + 1]; ++k) {
          const std::size_t oj = static_cast<std::size_t>(ci[k]) * words;
          for (int w = 0; w < words; ++w) b2[oi + w] |= b1[oj + w];
        }
        int c1 = 0, c2 = 0;
        for (int w = 0; w < words; ++w) {
          std::uint64_t x = b1[oi + w], y = b2[oi + w];
          while (x) { x &= x - 1; ++c1; }
          while (y) { y &= y - 1; ++c2; }
        }
        direct[i] += c1;
        indirect[i] += c2;
      }
    });
    Rcpp::checkUserInterrupt();
  }
  return List::create(_["direct"] = wrap(direct), _["indirect"] = wrap(indirect));
}

// Degree-matched null sampling. For permutation p and seed gene i (in
// order), draw uniformly from the gene's candidate pool minus the genes
// already drawn in this permutation; if nothing is left, draw from the full
// pool. Uses R's RNG exactly like
//
//   pool = setdiff(candidates[[i]], drawn); drawn[i] = sample(pool, 1)
//
// so results follow set.seed() and match the legacy implementation.
// Candidates are 1-based indices into a universe of size n_universe and
// must be unique within each pool. Returns an n_genes x n_perm matrix.
//
// [[Rcpp::export(.tkoi_sample_null)]]
IntegerMatrix tkoi_sample_null(IntegerVector cand_ptr, IntegerVector cand_idx,
                               int n_universe, int n_perm) {
  const int n_genes = cand_ptr.size() - 1;
  if (n_genes < 0 || n_perm < 0 || n_universe < 0) stop("Invalid dimensions.");
  if (cand_ptr[0] != 0 || cand_ptr[n_genes] != cand_idx.size()) stop("Malformed 'cand_ptr'.");
  for (int i = 0; i < n_genes; ++i) {
    if (cand_ptr[i + 1] <= cand_ptr[i]) stop("Every candidate pool must be non-empty.");
  }
  for (R_xlen_t j = 0; j < cand_idx.size(); ++j) {
    if (cand_idx[j] == NA_INTEGER || cand_idx[j] < 1 || cand_idx[j] > n_universe) {
      stop("Candidate index out of range.");
    }
  }
  IntegerMatrix out(n_genes, n_perm);
  std::vector<char> used(n_universe + 1, 0);
  std::vector<int> drawn(n_genes);

  GetRNGstate();
  for (int p = 0; p < n_perm; ++p) {
    for (int i = 0; i < n_genes; ++i) {
      const int a = cand_ptr[i];
      const int b = cand_ptr[i + 1];
      int available = 0;
      for (int j = a; j < b; ++j) available += !used[cand_idx[j]];
      int pick;
      if (available == 0) {
        const int k = static_cast<int>(R_unif_index(static_cast<double>(b - a)));
        pick = cand_idx[a + k];
      } else {
        int k = static_cast<int>(R_unif_index(static_cast<double>(available)));
        pick = -1;
        for (int j = a; j < b; ++j) {
          if (used[cand_idx[j]]) continue;
          if (k == 0) {
            pick = cand_idx[j];
            break;
          }
          --k;
        }
      }
      used[pick] = 1;
      drawn[i] = pick;
      out(i, p) = pick;
    }
    for (int i = 0; i < n_genes; ++i) used[drawn[i]] = 0;
    if ((p & 15) == 15) {
      PutRNGstate();
      Rcpp::checkUserInterrupt();
      GetRNGstate();
    }
  }
  PutRNGstate();
  return out;
}

// Total physical memory in bytes (NA if unknown).
//
// [[Rcpp::export(.tkoi_total_memory)]]
double tkoi_total_memory() {
#if defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status)) return static_cast<double>(status.ullTotalPhys);
  return NA_REAL;
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGESIZE);
  if (pages > 0 && page_size > 0) return static_cast<double>(pages) * static_cast<double>(page_size);
  return NA_REAL;
#else
  return NA_REAL;
#endif
}

// Number of hardware threads reported by the C++ runtime (0 if unknown).
//
// [[Rcpp::export(.tkoi_hardware_threads)]]
int tkoi_hardware_threads() {
  return static_cast<int>(std::thread::hardware_concurrency());
}
