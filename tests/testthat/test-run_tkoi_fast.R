###########################
# Author: Wanjun Gu
# Email: wanjun.gu@ucsf.edu
# Date: 2026-04-20
###########################

library(tkoi)

expression_data = data.table::fread(testthat::test_path("ra.csv"))

# Run both with identical seed and permutation count so outputs are comparable
n_perm = 10

set.seed(42)
t_orig = system.time({
  result_orig = run_tkoi(
    expression_data = expression_data,
    subnetwork = tkoi::tkoi_net,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.25,
    indirect_link_threshold = 3,
    topology_similarity = 0.9,
    n_permutation = n_perm,
    damping_factor = 0.85,
    maximum_iteration = 500
  )
})

set.seed(42)
t_fast = system.time({
  result_fast = run_tkoi_fast(
    expression_data = expression_data,
    subnetwork = tkoi::tkoi_net,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.25,
    indirect_link_threshold = 3,
    topology_similarity = 0.9,
    n_permutation = n_perm,
    damping_factor = 0.85,
    maximum_iteration = 500,
    n_cores = 1
  )
})

# ── Timing comparison ─────────────────────────────────────────────────────────
test_that("run_tkoi_fast is at least as fast as run_tkoi", {
  message(sprintf("\nrun_tkoi       total: %.1fs", t_orig["elapsed"]))
  message(sprintf("run_tkoi_fast  total: %.1fs", t_fast["elapsed"]))
  message(sprintf("Speedup:              %.2fx", t_orig["elapsed"] / t_fast["elapsed"]))

  step_t = attr(result_fast, "step_timings")
  message("\nrun_tkoi_fast per-step breakdown:")
  for(step in names(step_t)){
    message(sprintf("  %-35s %6.2fs  (%4.1f%%)",
                    step,
                    step_t[[step]],
                    100 * step_t[[step]] / step_t$total))
  }

  expect_true(is.list(step_t))
  expect_true(all(c("gene_filtering", "metapath_filters", "observed_pagerank",
                    "candidate_pool_precompute", "permutation",
                    "statistics", "annotation", "total") %in% names(step_t)))
})

# ── Output structure ──────────────────────────────────────────────────────────
test_that("run_tkoi_fast returns a tKOIList with identical structure", {
  expect_s4_class(result_fast, "tKOIList")
  expect_equal(names(result_fast@network_summary_statistics),
               names(result_orig@network_summary_statistics))
  expect_equal(nrow(result_fast@pagerank_data), nrow(result_orig@pagerank_data))
  expect_equal(ncol(result_fast@pagerank_data), ncol(result_orig@pagerank_data))
})

# ── Numerical correctness (identical output with same seed, n_cores=1) ────────
test_that("observed PageRank scores are identical to run_tkoi", {
  expect_equal(result_fast@pagerank_data$node_id,
               result_orig@pagerank_data$node_id)
  expect_equal(result_fast@pagerank_data$pagerank,
               result_orig@pagerank_data$pagerank)
})

test_that("permutation columns are identical to run_tkoi", {
  perm_cols = grep("^perm\\.", names(result_orig@pagerank_data), value = TRUE)
  for(col in perm_cols){
    expect_equal(result_fast@pagerank_data[[col]],
                 result_orig@pagerank_data[[col]],
                 label = col)
  }
})

test_that("beta and p_value are identical across all node types", {
  for(ntype in names(result_orig@network_summary_statistics)){
    orig_df = result_orig@network_summary_statistics[[ntype]]
    fast_df = result_fast@network_summary_statistics[[ntype]]
    if(nrow(orig_df) == 0) next

    # match rows by node_id since sort order could differ
    orig_sorted = orig_df[order(orig_df$node_id), ]
    fast_sorted = fast_df[order(fast_df$node_id), ]

    expect_equal(orig_sorted$beta, fast_sorted$beta,
                 label = paste0(ntype, "$beta"))
    expect_equal(orig_sorted$p_value, fast_sorted$p_value,
                 label = paste0(ntype, "$p_value"))
    expect_equal(orig_sorted$fdr, fast_sorted$fdr,
                 label = paste0(ntype, "$fdr"))
  }
})

test_that("stored parameters are correct in run_tkoi_fast result", {
  expect_equal(result_fast@n_permutation, n_perm)
  expect_equal(result_fast@pvalue_threshold, 0.05)
  expect_equal(result_fast@logfc_threshold, 0.25)
  expect_equal(result_fast@damping_factor, 0.85)
})
