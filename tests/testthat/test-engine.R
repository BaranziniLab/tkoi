# Unit tests of the native engine (src/tkoi_engine.cpp) and the internal
# helpers in R/engine.R. Everything runs on tiny igraph graphs or on the toy
# network, so the whole file takes a few seconds.

# Helpers ----------------------------------------------------------------------

# Normalized CSR matrix of an igraph object, weighted by its `weight` edge
# attribute when it has one. Directed graphs are treated as undirected.
graph_csr = function(graph, weight = igraph::E(graph)$weight) {
  edges = igraph::as_edgelist(graph, names = FALSE)
  .tkoi_build_csr(
    as.integer(edges[, 1]),
    as.integer(edges[, 2]),
    if (is.null(weight)) numeric(0) else as.numeric(weight),
    igraph::vcount(graph)
  )
}

# Dense version of a CSR matrix (duplicate entries add up).
csr_to_dense = function(csr) {
  dense = matrix(0, csr$n, csr$n)
  rows = rep(seq_len(csr$n), diff(csr$row_ptr))
  for (k in seq_along(csr$col)) {
    dense[rows[k], csr$col[k] + 1] = dense[rows[k], csr$col[k] + 1] + csr$val[k]
  }
  dense
}

# Dense adjacency matrix with undirected self-loops counted twice.
dense_adjacency = function(from, to, weight, n) {
  adjacency = matrix(0, n, n)
  for (e in seq_along(from)) {
    i = from[e]
    j = to[e]
    if (i == j) {
      adjacency[i, i] = adjacency[i, i] + 2 * weight[e]
    } else {
      adjacency[i, j] = adjacency[i, j] + weight[e]
      adjacency[j, i] = adjacency[j, i] + weight[e]
    }
  }
  adjacency
}

# Reference D^-1/2 A D^-1/2 built with plain R.
dense_normalized = function(from, to, weight, n) {
  adjacency = dense_adjacency(from, to, weight, n)
  degree = rowSums(adjacency)
  scale = ifelse(degree > 0, 1 / sqrt(degree), 1)
  list(matrix = adjacency * outer(scale, scale), scale = scale)
}

# Reference personalized PageRank from a dense solve of (I - d A D^-1) z = u,
# x = z / sum(z). A vertex without edges keeps its restart mass, as in prpack.
# The system has a 1-norm condition number of at most (1 + d) / (1 - d)
# whatever the edge weights, so the solve is accurate to a few ulps.
dense_ppr = function(graph, nodes, weights, damping = 0.85) {
  edges = igraph::as_edgelist(graph, names = FALSE)
  weight = igraph::E(graph)$weight
  if (is.null(weight)) {
    weight = rep(1, nrow(edges))
  }
  n = igraph::vcount(graph)
  adjacency = dense_adjacency(edges[, 1], edges[, 2], weight, n)
  degree = colSums(adjacency)
  transition = sweep(adjacency, 2, ifelse(degree > 0, degree, 1), "/")
  personalized = numeric(n)
  personalized[nodes] = weights
  z = solve(diag(n) - damping * transition, personalized)
  z / sum(z)
}

# Seed sets in the form .tkoi_ppr_null() takes: the observed seeds, and null
# seed sets that share one weight vector. `null_node` is a list with one
# vector of vertices per null set (or a matrix with one column per set).
seed_sets = function(observed_node, observed_weight, null_node = list(), null_weight = numeric(0)) {
  if (is.list(null_node)) {
    stopifnot(all(lengths(null_node) == length(null_weight)))
    null_node = matrix(as.integer(unlist(null_node)), nrow = length(null_weight), ncol = length(null_node))
  }
  storage.mode(null_node) = "integer"
  list(
    observed_node = as.integer(observed_node),
    observed_weight = as.numeric(observed_weight),
    null_node = null_node,
    null_weight = as.numeric(null_weight)
  )
}

# Every seed set as vertices and weights: the observed set, then the null sets.
all_sets = function(seeds) {
  columns = lapply(seq_len(ncol(seeds$null_node)), function(p) seeds$null_node[, p])
  list(
    nodes = c(list(seeds$observed_node), columns),
    weights = c(list(seeds$observed_weight), rep(list(seeds$null_weight), length(columns)))
  )
}

run_ppr = function(
  csr,
  seeds,
  damping = 0.85,
  tol = 1e-14,
  max_iter = 1000L,
  n_threads = 1L,
  max_block = 16L,
  keep_null = TRUE,
  progress = NULL
) {
  .tkoi_ppr_null(
    csr,
    seeds$observed_node,
    seeds$observed_weight,
    seeds$null_node,
    seeds$null_weight,
    damping,
    tol,
    as.integer(max_iter),
    as.integer(n_threads),
    as.integer(max_block),
    keep_null,
    progress
  )
}

# igraph reference: one personalized PageRank vector (column) per seed set,
# observed set first.
igraph_ppr = function(graph, seeds, damping = 0.85) {
  sets = all_sets(seeds)
  columns = vapply(seq_along(sets$nodes), function(k) {
    personalized = numeric(igraph::vcount(graph))
    personalized[sets$nodes[[k]]] = sets$weights[[k]]
    igraph::page_rank(
      graph,
      algo = "prpack",
      personalized = personalized,
      damping = damping,
      directed = FALSE
    )$vector
  }, numeric(igraph::vcount(graph)))
  unname(columns)
}

null_matrix = function(result) {
  do.call(cbind, result$null_columns)
}

# Checks the engine against igraph for the seed sets on `graph`.
expect_ppr_matches_igraph = function(graph, seeds, damping = 0.85, csr = graph_csr(graph), tol = 1e-14) {
  result = run_ppr(csr, seeds, damping = damping, tol = tol)
  reference = igraph_ppr(graph, seeds, damping = damping)
  expect_equal(result$observed, reference[, 1], tolerance = 1e-10)
  expect_lt(max(abs(result$observed - reference[, 1])), 1e-10)
  if (ncol(reference) > 1) {
    expect_equal(null_matrix(result), reference[, -1, drop = FALSE], tolerance = 1e-10)
    expect_lt(max(abs(null_matrix(result) - reference[, -1])), 1e-10)
  }
  expect_true(all(result$rel_residual <= tol))
  invisible(result)
}

# A small undirected graph with self-loops, multi-edges, and isolated vertices.
messy_graph = function(n = 40, m = 90, seed = 1) {
  withr::with_seed(seed, {
    graph = igraph::sample_gnm(n, m)
    graph = igraph::add_edges(graph, c(1, 1, 2, 2, 5, 5, 3, 4, 3, 4, 3, 4, 6, 7, 6, 7))
    graph = igraph::add_vertices(graph, 3)
    igraph::E(graph)$weight = round(stats::runif(igraph::ecount(graph), 0.1, 3), 3)
  })
  graph
}

# Number of seeds within `order` hops of every vertex, from igraph::ego().
ego_counts = function(graph, seeds, order) {
  reached = unlist(lapply(igraph::ego(graph, order = order, nodes = seeds, mode = "all"), as.integer))
  tabulate(reached, nbins = igraph::vcount(graph))
}

# The released null sampler: for each seed gene in turn, draw one gene from
# its pool minus the genes already drawn in this permutation, falling back
# to the full pool once that is empty. Pools are character IDs, as in the
# released code, so sample() never sees a length-one numeric vector.
legacy_sample_null = function(pools, n_perm) {
  pools = lapply(pools, function(pool) paste0("u", pool))
  n_genes = length(pools)
  draws = matrix(0L, n_genes, n_perm)
  for (p in seq_len(n_perm)) {
    substitute_genes = character(n_genes)
    for (i in seq_len(n_genes)) {
      pool = setdiff(pools[[i]], substitute_genes)
      if (length(pool) == 0L) pool = pools[[i]]
      substitute_genes[i] = sample(pool, 1)
    }
    draws[, p] = as.integer(substring(substitute_genes, 2))
  }
  draws
}

unpack_pools = function(pools) {
  sizes = diff(pools$ptr)
  unname(split(pools$idx, factor(rep(seq_along(sizes), sizes), levels = seq_along(sizes))))
}

brute_force_pools = function(seed_degree, universe_degree, topology_similarity) {
  lapply(seed_degree, function(degree) {
    which(
      universe_degree >= topology_similarity * degree &
        universe_degree <= (2 - topology_similarity) * degree
    )
  })
}

detected_cores = function() {
  cores = suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(cores) || cores < 1) {
    cores = max(1L, .tkoi_hardware_threads())
  }
  as.integer(cores)
}

# Serves fake first lines of system files (named by path) to the engine
# helpers; every other file reads as missing.
with_system_files = function(files, code) {
  testthat::with_mocked_bindings(
    code,
    .tkoi_read_first_line = function(path) if (path %in% names(files)) files[[path]] else NULL,
    .package = "tkoi"
  )
}

# Removes the Slurm, cgroup, and R CMD check core limits for the calling
# test, so the expected core count is the hardware count on any machine.
local_unlimited_cores = function(env = parent.frame()) {
  withr::local_envvar(c(SLURM_CPUS_PER_TASK = NA, "_R_CHECK_LIMIT_CORES_" = NA), .local_envir = env)
  testthat::local_mocked_bindings(.tkoi_read_first_line = function(path) NULL, .package = "tkoi", .env = env)
}

local_memory_limit = function(bytes, env = parent.frame()) {
  testthat::local_mocked_bindings(.tkoi_memory_limit = function() bytes, .package = "tkoi", .env = env)
}

# Counts the matrix builds made by .tkoi_prepare_network() in the calling test.
local_build_counter = function(env = parent.frame()) {
  counter = new.env(parent = emptyenv())
  counter$calls = 0L
  build = .tkoi_build_csr
  testthat::local_mocked_bindings(
    .tkoi_build_csr = function(...) {
      counter$calls = counter$calls + 1L
      build(...)
    },
    .package = "tkoi",
    .env = env
  )
  counter
}

small_network = function() {
  graph = igraph::make_graph(c(1, 2, 2, 3, 3, 3, 3, 4, 1, 4, 1, 2), n = 5, directed = FALSE)
  igraph::V(graph)$name = paste0("node", 1:5)
  igraph::V(graph)$labels = c("['Gene']", "['Gene']", "['Disease']", "['Pathway']", "['Gene']")
  igraph::V(graph)$identifier = paste0("ID:", 1:5)
  igraph::V(graph)$degree = c(10, 20, 30, 40, 50)
  graph
}

small_network_csr = function(weight = numeric(0)) {
  edges = igraph::as_edgelist(small_network(), names = FALSE)
  .tkoi_build_csr(edges[, 1], edges[, 2], weight, 5L)
}

# .tkoi_build_csr() ------------------------------------------------------------

test_that(".tkoi_build_csr() builds the normalized matrix of a hand-worked graph", {
  # 1-2 (w = 1), 2-3 (w = 2), loop 3-3 (w = 0.5), vertex 4 isolated.
  # Degrees: 1, 3, 2 + 2 * 0.5 = 3, 0.
  csr = .tkoi_build_csr(c(1L, 2L, 3L), c(2L, 3L, 3L), c(1, 2, 0.5), 4L)

  expect_named(csr, c("row_ptr", "col", "val", "scale", "n"))
  expect_type(csr$row_ptr, "integer")
  expect_type(csr$col, "integer")
  expect_type(csr$val, "double")
  expect_type(csr$scale, "double")
  expect_identical(csr$n, 4L)
  expect_identical(csr$row_ptr, c(0L, 1L, 3L, 5L, 5L))
  expect_length(csr$col, 5)
  expect_length(csr$val, 5)
  expect_equal(csr$scale, c(1, 1 / sqrt(3), 1 / sqrt(3), 1))

  expected = matrix(0, 4, 4)
  expected[1, 2] = expected[2, 1] = 1 / sqrt(3)
  expected[2, 3] = expected[3, 2] = 2 / 3
  expected[3, 3] = 2 * 0.5 / 3
  expect_equal(csr_to_dense(csr), expected)

  # Every stored column index is 0-based and in range.
  expect_true(all(csr$col >= 0 & csr$col < 4))
})

test_that(".tkoi_build_csr() counts an undirected self-loop as 2 * weight", {
  csr = .tkoi_build_csr(1L, 1L, 3, 1L)
  expect_identical(csr$row_ptr, c(0L, 1L))
  expect_identical(csr$col, 0L)
  expect_equal(csr$scale, 1 / sqrt(6))
  expect_equal(csr$val, 1)

  # A loop next to an ordinary edge: degree(1) = 2 * 2 + 1 = 5.
  csr = .tkoi_build_csr(c(1L, 1L), c(1L, 2L), c(2, 1), 2L)
  expect_equal(csr$scale, c(1 / sqrt(5), 1))
  expect_equal(csr_to_dense(csr), matrix(c(4 / 5, 1 / sqrt(5), 1 / sqrt(5), 0), 2, 2))
})

test_that(".tkoi_build_csr() matches a dense reference on a random graph", {
  graph = messy_graph()
  edges = igraph::as_edgelist(graph, names = FALSE)
  weight = igraph::E(graph)$weight
  csr = graph_csr(graph)
  reference = dense_normalized(edges[, 1], edges[, 2], weight, igraph::vcount(graph))

  dense = csr_to_dense(csr)
  expect_equal(dense, reference$matrix, tolerance = 1e-14)
  expect_true(isSymmetric(dense))
  expect_equal(csr$scale, reference$scale, tolerance = 1e-14)

  # One stored entry per edge end (a loop is stored once).
  loops = edges[, 1] == edges[, 2]
  entries = tabulate(c(edges[, 1], edges[!loops, 2]), nbins = igraph::vcount(graph))
  expect_identical(diff(csr$row_ptr), entries)
  expect_identical(csr$row_ptr[length(csr$row_ptr)], as.integer(2 * sum(!loops) + sum(loops)))

  # Isolated vertices have no entries and a scale of 1.
  isolated = which(igraph::degree(graph) == 0)
  expect_gt(length(isolated), 0)
  expect_true(all(diff(csr$row_ptr)[isolated] == 0))
  expect_true(all(csr$scale[isolated] == 1))
})

test_that(".tkoi_build_csr() treats an empty weight vector as unit weights", {
  graph = messy_graph()
  edges = igraph::as_edgelist(graph, names = FALSE)
  unweighted = .tkoi_build_csr(edges[, 1], edges[, 2], numeric(0), igraph::vcount(graph))
  unit = .tkoi_build_csr(edges[, 1], edges[, 2], rep(1, nrow(edges)), igraph::vcount(graph))
  expect_identical(unweighted, unit)

  weighted = graph_csr(graph)
  expect_identical(weighted$row_ptr, unweighted$row_ptr)
  expect_identical(weighted$col, unweighted$col)
  expect_false(isTRUE(all.equal(weighted$val, unweighted$val)))

  # Scaling every weight by a constant leaves the normalized matrix unchanged.
  scaled = graph_csr(graph, weight = 7 * igraph::E(graph)$weight)
  expect_equal(scaled$val, weighted$val, tolerance = 1e-14)

  # Zero weights are allowed; a vertex whose edges all weigh zero behaves
  # like an isolated vertex.
  zero = .tkoi_build_csr(c(1L, 2L), c(2L, 3L), c(0, 1), 3L)
  expect_equal(zero$scale, c(1, 1, 1))
  expect_equal(csr_to_dense(zero), matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3))
})

test_that(".tkoi_build_csr() handles graphs without edges or vertices", {
  empty = .tkoi_build_csr(integer(0), integer(0), numeric(0), 3L)
  expect_identical(empty$row_ptr, c(0L, 0L, 0L, 0L))
  expect_identical(empty$col, integer(0))
  expect_identical(empty$val, numeric(0))
  expect_identical(empty$scale, c(1, 1, 1))

  none = .tkoi_build_csr(integer(0), integer(0), numeric(0), 0L)
  expect_identical(none$row_ptr, 0L)
  expect_identical(none$n, 0L)
})

test_that(".tkoi_build_csr() validates its inputs", {
  expect_error(.tkoi_build_csr(1:2, 1L, numeric(0), 3L), "same length")
  expect_error(.tkoi_build_csr(1L, 2L, numeric(0), -1L), "non-negative")
  expect_error(.tkoi_build_csr(1L, 2L, numeric(0), NA_integer_), "non-negative")
  expect_error(.tkoi_build_csr(1L, 2L, c(1, 2), 3L), "match the number of edges")
  expect_error(.tkoi_build_csr(1L, 4L, numeric(0), 3L), "out of range")
  expect_error(.tkoi_build_csr(4L, 1L, numeric(0), 3L), "out of range")
  expect_error(.tkoi_build_csr(0L, 2L, numeric(0), 3L), "out of range")
  expect_error(.tkoi_build_csr(2L, 0L, numeric(0), 3L), "out of range")
  expect_error(.tkoi_build_csr(1L, 1L, numeric(0), 0L), "out of range")
  expect_error(.tkoi_build_csr(1L, 2L, -1, 3L), "finite and non-negative")
  expect_error(.tkoi_build_csr(1L, 2L, Inf, 3L), "finite and non-negative")
  expect_error(.tkoi_build_csr(1L, 2L, NaN, 3L), "finite and non-negative")
  expect_error(.tkoi_build_csr(1L, 2L, NA_real_, 3L), "finite and non-negative")
  expect_error(.tkoi_build_csr(c(1L, 2L), c(2L, 3L), c(1, -1e-300), 3L), "finite and non-negative")
})

test_that(".tkoi_build_csr() rejects missing edge endpoints", {
  message = "Edge endpoint out of range."
  expect_error(.tkoi_build_csr(NA_integer_, 2L, numeric(0), 3L), message, fixed = TRUE)
  expect_error(.tkoi_build_csr(1L, NA_integer_, numeric(0), 3L), message, fixed = TRUE)
  expect_error(.tkoi_build_csr(c(1L, NA), c(2L, 3L), c(1, 1), 3L), message, fixed = TRUE)
  expect_error(.tkoi_build_csr(c(1L, 2L), c(2L, NA), numeric(0), 3L), message, fixed = TRUE)
  expect_error(.tkoi_build_csr(NA_integer_, NA_integer_, numeric(0), 3L), message, fixed = TRUE)
})

# .tkoi_ppr_null(): agreement with igraph -------------------------------------

test_that(".tkoi_ppr_null() matches igraph on a weighted graph with loops, multi-edges, and isolated vertices", {
  graph = messy_graph()
  n = igraph::vcount(graph)
  isolated = which(igraph::degree(graph) == 0)
  expect_length(isolated, 3)

  seeds = seed_sets(
    observed_node = c(1, 5, isolated[1]),
    observed_weight = withr::with_seed(4, stats::runif(3, 0.1, 1)),
    null_node = list(
      c(2, 3, 10),
      rep(isolated[2], 3),
      withr::with_seed(2, sample(n, 3)),
      withr::with_seed(3, sample(n, 3))
    ),
    null_weight = c(0.3, 0.6, 0.9)
  )
  result = expect_ppr_matches_igraph(graph, seeds)

  # Unseeded isolated vertices get nothing; a seeded one keeps its restart mass.
  expect_identical(result$observed[isolated[2:3]], c(0, 0))
  expect_gt(result$observed[isolated[1]], 0)
  expect_equal(result$null_columns[[2]], replace(numeric(n), isolated[2], 1))

  columns = cbind(result$observed, null_matrix(result))
  expect_equal(unname(colSums(columns)), rep(1, 5), tolerance = 1e-12)
  expect_true(all(columns >= 0))
})

test_that(".tkoi_ppr_null() matches igraph on disconnected components", {
  graph = igraph::disjoint_union(igraph::make_ring(12), igraph::make_star(9, mode = "undirected"))
  graph = igraph::add_edges(graph, c(13, 13))
  second = 13:21

  # The observed set and the null sets may have different sizes.
  seeds = seed_sets(c(1, 4, 7), c(1, 3, 2), list(c(14, 20), c(2, 13), c(13, 13)), c(0.2, 0.8))
  result = expect_ppr_matches_igraph(graph, seeds)

  # Components without seeds are exactly zero.
  expect_true(all(result$observed[second] == 0))
  expect_true(all(result$null_columns[[1]][1:12] == 0))
  expect_true(all(result$null_columns[[3]][1:12] == 0))
})

test_that(".tkoi_ppr_null() treats a directed igraph as undirected", {
  graph = withr::with_seed(11, igraph::sample_gnm(35, 110, directed = TRUE, loops = TRUE))
  graph = igraph::add_edges(graph, c(1, 2, 2, 1, 3, 3, 4, 5, 4, 5))
  igraph::E(graph)$weight = withr::with_seed(12, stats::runif(igraph::ecount(graph), 0.5, 2))
  expect_true(igraph::is_directed(graph))

  seeds = withr::with_seed(13, seed_sets(
    sample(35, 5),
    stats::runif(5),
    replicate(3, sample(35, 6), simplify = FALSE),
    stats::runif(6)
  ))
  expect_ppr_matches_igraph(graph, seeds)

  # The same graph without weights.
  unweighted = igraph::delete_edge_attr(graph, "weight")
  expect_ppr_matches_igraph(unweighted, seeds)
})

test_that(".tkoi_ppr_null() matches igraph for other damping factors", {
  graph = messy_graph(seed = 21)
  seeds = withr::with_seed(22, seed_sets(
    sample(40, 4),
    1:4,
    replicate(2, sample(40, 3), simplify = FALSE),
    1:3
  ))
  for (damping in c(0.3, 0.5, 0.95)) {
    expect_ppr_matches_igraph(graph, seeds, damping = damping)
  }
})

test_that(".tkoi_ppr_null() matches igraph on many small graph shapes", {
  graphs = list(
    star = igraph::make_star(15, mode = "undirected"),
    complete = igraph::make_full_graph(8),
    path = igraph::make_ring(20, circular = FALSE),
    ring = igraph::make_ring(17),
    lattice = igraph::make_lattice(c(5, 5)),
    tree = igraph::make_tree(31, 2, mode = "undirected"),
    bipartite = igraph::make_full_bipartite_graph(4, 6),
    sparse = withr::with_seed(23, igraph::sample_gnm(30, 25)),
    dense = withr::with_seed(24, igraph::sample_gnp(25, 0.6))
  )
  for (name in names(graphs)) {
    graph = graphs[[name]]
    n = igraph::vcount(graph)
    seeds = withr::with_seed(25, seed_sets(
      sample(n, 2),
      stats::runif(2),
      replicate(3, sample(n, 2), simplify = FALSE),
      stats::runif(2)
    ))
    result = run_ppr(graph_csr(graph), seeds)
    reference = igraph_ppr(graph, seeds)
    expect_lt(max(abs(cbind(result$observed, null_matrix(result)) - reference)), 1e-10, label = name)
    expect_true(all(result$rel_residual <= 1e-14), label = name)
  }
})

test_that(".tkoi_ppr_null() matches igraph on the toy network", {
  graph = toy_network()
  expect_true(any(igraph::which_loop(graph)))
  n = igraph::vcount(graph)
  loop_vertex = igraph::ends(graph, which(igraph::which_loop(graph))[1], names = FALSE)[1]
  seeds = withr::with_seed(31, seed_sets(
    c(loop_vertex, sample(n, 5)),
    stats::runif(6),
    list(sample(n, 6), sample(n, 6)),
    stats::runif(6)
  ))
  expect_ppr_matches_igraph(graph, seeds, csr = .tkoi_prepare_network(graph)$csr)
  .tkoi_clear_cache()
})

test_that(".tkoi_ppr_null() handles a graph without edges", {
  csr = .tkoi_build_csr(integer(0), integer(0), numeric(0), 4L)
  result = run_ppr(csr, seed_sets(c(1, 3), c(1, 3), list(2, 4), 5))
  expect_equal(result$observed, c(0.25, 0, 0.75, 0))
  expect_equal(result$null_columns, list(c(0, 1, 0, 0), c(0, 0, 0, 1)))
  expect_equal(result$null_mean, c(0, 0.5, 0, 0.5))
  expect_equal(result$null_sd, c(0, sqrt(0.5), 0, sqrt(0.5)))
})

# .tkoi_ppr_null(): accuracy ---------------------------------------------------

test_that(".tkoi_ppr_null() is accurate when edge weights span twelve orders of magnitude", {
  graph = withr::with_seed(91, igraph::sample_gnm(80, 240))
  # A path through every vertex keeps the graph connected; two self-loops.
  graph = igraph::add_edges(graph, c(rbind(1:79, 2:80), 5, 5, 9, 9))
  # Log-uniform weights between 1e-6 and 1e6.
  igraph::E(graph)$weight = withr::with_seed(92, 10^stats::runif(igraph::ecount(graph), -6, 6))
  expect_gt(max(igraph::E(graph)$weight) / min(igraph::E(graph)$weight), 1e10)

  seeds = seed_sets(
    c(1, 40, 77),
    c(1, 1e-3, 5),
    withr::with_seed(93, replicate(5, sample(80, 4), simplify = FALSE)),
    10^c(-3, 0, 1, 3)
  )
  result = run_ppr(graph_csr(graph), seeds, tol = 1e-14)
  expect_true(all(result$rel_residual <= 1e-14))

  sets = all_sets(seeds)
  columns = c(list(result$observed), result$null_columns)
  for (k in seq_along(columns)) {
    reference = dense_ppr(graph, sets$nodes[[k]], sets$weights[[k]])
    expect_lt(sum(abs(columns[[k]] - reference)), 1e-10)
  }
})

test_that(".tkoi_ppr_null() leaves vertices beyond its Krylov horizon at exactly zero", {
  # After k conjugate gradient iterations the solution lies in the Krylov
  # space of the seeds, so vertices farther than k hops from every seed are
  # exactly zero. On a long path the true PageRank there is positive but far
  # below the tolerance, so the 1-norm error stays tiny; elementwise equality
  # with prpack's tiny values is therefore not expected.
  path = igraph::make_ring(300, circular = FALSE)
  seeds = seed_sets(1, 1, list(150), 1)
  result = run_ppr(graph_csr(path), seeds)
  iterations = result$iterations[1]
  expect_true(all(result$rel_residual <= 1e-14))
  expect_gt(iterations, 10)
  expect_lt(iterations, 100)

  from_end = 0:299
  from_middle = abs(1:300 - 150)
  expect_true(all(result$observed[from_end > iterations] == 0))
  expect_true(all(result$null_columns[[1]][from_middle > iterations] == 0))
  expect_true(all(result$observed[from_end < iterations] > 0))
  expect_true(all(result$null_columns[[1]][from_middle < iterations] > 0))

  # The exact PageRank is positive everywhere, but what the engine sets to
  # zero sums to less than the tolerance.
  reference = dense_ppr(path, 1, 1)
  expect_true(all(reference > 0))
  expect_lt(sum(reference[from_end > iterations]), 1e-14)
  expect_lt(sum(abs(result$observed - reference)), 1e-12)
})

# .tkoi_ppr_null(): null statistics and options ---------------------------------

test_that(".tkoi_ppr_null() null mean and sd equal rowMeans() and sd() of the null columns", {
  graph = messy_graph(seed = 5)
  n = igraph::vcount(graph)
  seeds = seed_sets(
    1:3,
    c(1, 1, 1),
    withr::with_seed(6, replicate(8, sample(n, 5), simplify = FALSE)),
    c(0.4, 0.3, 0.1, 0.1, 0.1)
  )
  result = run_ppr(graph_csr(graph), seeds, max_block = 4L)

  expect_named(result, c("observed", "null_columns", "null_mean", "null_sd", "iterations", "rel_residual"))
  expect_length(result$null_columns, 8)
  null = null_matrix(result)
  expect_equal(result$null_mean, rowMeans(null), tolerance = 1e-13)
  expect_equal(result$null_sd, apply(null, 1, stats::sd), tolerance = 1e-12)
  expect_type(result$iterations, "integer")
  expect_length(result$iterations, 9)
  expect_length(result$rel_residual, 9)
})

test_that(".tkoi_ppr_null() with keep_null = FALSE returns the same statistics without the columns", {
  graph = messy_graph(seed = 7)
  seeds = seed_sets(1:4, 4:1, withr::with_seed(8, replicate(5, sample(40, 4), simplify = FALSE)), 1:4)
  csr = graph_csr(graph)

  kept = run_ppr(csr, seeds, keep_null = TRUE, max_block = 2L)
  dropped = run_ppr(csr, seeds, keep_null = FALSE, max_block = 2L)
  expect_identical(dropped$null_columns, list())
  expect_length(kept$null_columns, 5)
  expect_identical(dropped$observed, kept$observed)
  expect_identical(dropped$null_mean, kept$null_mean)
  expect_identical(dropped$null_sd, kept$null_sd)
  expect_identical(dropped$rel_residual, kept$rel_residual)
  expect_identical(dropped$iterations, kept$iterations)
})

test_that(".tkoi_ppr_null() runs the observed vector alone when there are no null sets", {
  graph = messy_graph(seed = 13)
  n = igraph::vcount(graph)
  csr = graph_csr(graph)
  seeds = seed_sets(c(2, 8), c(1, 3))
  expect_identical(dim(seeds$null_node), c(0L, 0L))

  result = run_ppr(csr, seeds)
  expect_named(result, c("observed", "null_columns", "null_mean", "null_sd", "iterations", "rel_residual"))
  expect_identical(result$null_columns, list())
  expect_identical(result$null_mean, rep(NA_real_, n))
  expect_identical(result$null_sd, rep(NA_real_, n))
  expect_length(result$iterations, 1)
  expect_length(result$rel_residual, 1)
  expect_lte(result$rel_residual, 1e-14)
  expect_equal(result$observed, igraph_ppr(graph, seeds)[, 1], tolerance = 1e-10)

  # A null matrix with rows but no columns works the same, with or without
  # keeping null columns.
  empty_columns = seed_sets(c(2, 8), c(1, 3), matrix(integer(0), 3, 0), c(1, 1, 1))
  expect_identical(run_ppr(csr, empty_columns), result)
  expect_identical(run_ppr(csr, empty_columns, keep_null = FALSE), result)

  # The observed vector does not depend on the null sets solved next to it.
  with_nulls = run_ppr(csr, seed_sets(c(2, 8), c(1, 3), list(c(1, 2), c(3, 4)), c(1, 1)))
  expect_identical(with_nulls$observed, result$observed)
  expect_identical(with_nulls$rel_residual[1], result$rel_residual)
})

test_that(".tkoi_ppr_null() handles a single null vector", {
  csr = graph_csr(igraph::make_ring(10))
  # With one null vector, the mean is that vector and the sd is undefined.
  one_null = run_ppr(csr, seed_sets(1, 1, list(2), 1))
  expect_length(one_null$null_columns, 1)
  expect_identical(one_null$null_mean, one_null$null_columns[[1]])
  expect_true(all(is.na(one_null$null_sd)))
  expect_length(one_null$iterations, 2)
})

test_that(".tkoi_ppr_null() applies duplicate seed vertices by assignment (last weight wins)", {
  graph = messy_graph(seed = 9)
  csr = graph_csr(graph)
  duplicated_seeds = run_ppr(
    csr,
    seed_sets(c(3, 7, 3), c(0.2, 0.5, 0.9), list(c(4, 9, 9), c(9, 4, 4)), c(1, 5, 2))
  )
  last_weight = run_ppr(csr, seed_sets(c(7, 3), c(0.5, 0.9), list(c(4, 9), c(9, 4)), c(1, 2)))
  expect_identical(duplicated_seeds$observed, last_weight$observed)
  expect_identical(duplicated_seeds$null_columns, last_weight$null_columns)
  expect_identical(duplicated_seeds$null_mean, last_weight$null_mean)
  expect_identical(duplicated_seeds$null_sd, last_weight$null_sd)

  # Same semantics as `personalized[c(3, 7, 3)] = c(0.2, 0.5, 0.9)` in R.
  reference = igraph_ppr(graph, seed_sets(c(3, 7, 3), c(0.2, 0.5, 0.9)))
  expect_equal(duplicated_seeds$observed, reference[, 1], tolerance = 1e-10)

  # The first weight is not used.
  first_weight = run_ppr(csr, seed_sets(c(3, 7), c(0.2, 0.5)))
  expect_false(isTRUE(all.equal(duplicated_seeds$observed, first_weight$observed)))
})

test_that(".tkoi_ppr_null() reports a large residual when maximum iterations are too few", {
  graph = messy_graph(seed = 10)
  seeds = seed_sets(c(1, 2), c(1, 1), list(c(3, 4), c(5, 6)), c(1, 2))
  csr = graph_csr(graph)

  short = run_ppr(csr, seeds, tol = 1e-12, max_iter = 2L)
  expect_identical(short$iterations, rep(2L, 3))
  expect_true(all(short$rel_residual > 1e-12))

  one = run_ppr(csr, seeds, tol = 1e-12, max_iter = 1L)
  expect_identical(one$iterations, rep(1L, 3))
  expect_true(all(one$rel_residual > short$rel_residual))

  converged = run_ppr(csr, seeds, tol = 1e-12, max_iter = 500L)
  expect_true(all(converged$rel_residual <= 1e-12))
  expect_true(all(converged$iterations > 2L & converged$iterations < 500L))
})

test_that(".tkoi_ppr_null() calls progress with (done, total), ending at total", {
  csr = graph_csr(messy_graph(seed = 12))
  # One observed and six null vectors.
  seeds = seed_sets(1, 1, as.list(2:7), 1)

  record_progress = function(max_block) {
    log = new.env()
    log$calls = list()
    progress = function(done, total) {
      log$calls[[length(log$calls) + 1]] = c(done = done, total = total)
    }
    run_ppr(csr, seeds, max_block = max_block, progress = progress)
    do.call(rbind, log$calls)
  }

  for (max_block in c(1L, 2L, 4L, 16L)) {
    calls = record_progress(max_block)
    expect_true(all(calls[, "total"] == 7))
    expect_true(all(diff(calls[, "done"]) > 0))
    expect_identical(unname(calls[nrow(calls), "done"]), 7L)
  }
  expect_identical(unname(record_progress(1L)[, "done"]), 1:7)
  expect_identical(unname(record_progress(2L)[, "done"]), c(2L, 4L, 6L, 7L))
  expect_identical(nrow(record_progress(16L)), 1L)

  # No callback is fine too.
  expect_no_error(run_ppr(csr, seeds, progress = NULL))
})

# .tkoi_ppr_null(): determinism --------------------------------------------------

ppr_determinism_case = function() {
  graph = messy_graph(n = 300, m = 900, seed = 41)
  n = igraph::vcount(graph)
  # K = 1 + 10 = 11 vectors, not a multiple of any block width.
  seeds = withr::with_seed(42, seed_sets(
    sample(n, 8),
    stats::runif(8),
    replicate(10, sample(n, 8), simplify = FALSE),
    stats::runif(8)
  ))
  list(csr = graph_csr(graph), seeds = seeds)
}

test_that(".tkoi_ppr_null() is bit-identical for any number of threads", {
  case = ppr_determinism_case()
  for (max_block in c(1L, 4L, 16L)) {
    reference = run_ppr(case$csr, case$seeds, n_threads = 1L, max_block = max_block)
    expect_identical(run_ppr(case$csr, case$seeds, n_threads = 2L, max_block = max_block), reference)
  }
})

test_that(".tkoi_ppr_null() is bit-identical on four threads", {
  skip_on_cran()
  case = ppr_determinism_case()
  for (max_block in c(2L, 8L, 16L)) {
    reference = run_ppr(case$csr, case$seeds, n_threads = 1L, max_block = max_block)
    expect_identical(run_ppr(case$csr, case$seeds, n_threads = 4L, max_block = max_block), reference)
  }
})

test_that(".tkoi_ppr_null() is bit-identical for any block width", {
  # run_tkoi() picks the block width from the machine's memory, so no output,
  # including the null mean and sd, may depend on it.
  case = ppr_determinism_case()
  reference = run_ppr(case$csr, case$seeds, max_block = 16L)
  expect_length(reference$null_columns, 10)
  for (max_block in c(1L, 2L, 4L, 8L)) {
    result = run_ppr(case$csr, case$seeds, max_block = max_block)
    for (output in names(reference)) {
      expect_identical(result[[output]], reference[[output]], info = paste0("max_block = ", max_block, ": ", output))
    }
  }
})

# .tkoi_ppr_null(): validation ------------------------------------------------------

test_that(".tkoi_ppr_null() validates its inputs", {
  csr = graph_csr(igraph::make_ring(10))
  ppr = function(
    observed_node = 1L,
    observed_weight = 1,
    null_node = matrix(2L, 1, 1),
    null_weight = 1,
    damping = 0.85,
    tol = 1e-12,
    max_iter = 100L
  ) {
    .tkoi_ppr_null(
      csr, observed_node, observed_weight, null_node, null_weight,
      damping, tol, max_iter, 1L, 16L, TRUE
    )
  }
  expect_no_error(ppr())
  expect_no_error(ppr(max_iter = 1L))

  for (damping in c(0, 1, -0.5, 1.5, NaN, NA_real_)) {
    expect_error(ppr(damping = damping), "'damping' must be in (0, 1).", fixed = TRUE)
  }
  for (tol in c(0, -1e-12, NaN, NA_real_)) {
    expect_error(ppr(tol = tol), "'tol' must be positive.", fixed = TRUE)
  }
  for (max_iter in c(0L, -5L, NA_integer_)) {
    expect_error(ppr(max_iter = max_iter), "'max_iter' must be a positive integer.", fixed = TRUE)
  }

  expect_error(ppr(observed_node = 1:2), "Observed seed vectors differ in length.", fixed = TRUE)
  expect_error(ppr(observed_weight = c(1, 1)), "Observed seed vectors differ in length.", fixed = TRUE)
  weight_length = "'null_weight' must have one weight per row of 'null_node'."
  expect_error(ppr(null_weight = c(1, 1)), weight_length, fixed = TRUE)
  expect_error(ppr(null_weight = numeric(0)), weight_length, fixed = TRUE)
  expect_error(ppr(null_node = matrix(integer(0), 2, 0)), weight_length, fixed = TRUE)
  expect_error(ppr(null_node = 2L), "matrix")

  out_of_range = "Seed vertex out of range."
  for (vertex in c(0L, 11L, -3L, NA_integer_)) {
    expect_error(ppr(observed_node = vertex), out_of_range, fixed = TRUE)
    expect_error(ppr(null_node = matrix(vertex, 1, 1)), out_of_range, fixed = TRUE)
    expect_error(ppr(null_node = matrix(c(2L, 3L, vertex), 1, 3)), out_of_range, fixed = TRUE)
  }
  expect_error(
    ppr(null_node = matrix(c(1L, NA, 3L, 4L), 2, 2), null_weight = c(1, 1)),
    out_of_range,
    fixed = TRUE
  )

  not_finite = "Seed weights must be finite and non-negative."
  for (weight in c(-1, -1e-300, Inf, -Inf, NaN, NA_real_)) {
    expect_error(ppr(observed_weight = weight), not_finite, fixed = TRUE)
    expect_error(ppr(null_weight = weight), not_finite, fixed = TRUE)
  }

  zero_total = "Seed weights must have a positive total."
  expect_error(ppr(observed_node = 1:2, observed_weight = c(0, 0)), zero_total, fixed = TRUE)
  expect_error(ppr(observed_node = integer(0), observed_weight = numeric(0)), zero_total, fixed = TRUE)
  expect_error(ppr(null_weight = 0), zero_total, fixed = TRUE)
  expect_error(
    ppr(null_node = matrix(1:4, 2, 2), null_weight = c(0, 0)),
    zero_total,
    fixed = TRUE
  )
  # A zero weight next to a positive one is allowed.
  expect_no_error(ppr(observed_node = 1:2, observed_weight = c(0, 1)))
})

# .tkoi_seed_reach() ---------------------------------------------------------------

test_that(".tkoi_seed_reach() matches igraph::ego() counts", {
  graph = messy_graph(seed = 51)
  csr = graph_csr(graph)
  # Vertex 1 has a self-loop; 2 and 7 appear twice; 41 is isolated.
  seeds = c(1L, 2L, 2L, 7L, 12L, 7L, 30L, 41L)
  reach = .tkoi_seed_reach(csr$row_ptr, csr$col, seeds, 1L)

  expect_named(reach, c("direct", "indirect"))
  expect_type(reach$direct, "integer")
  expect_identical(reach$direct, ego_counts(graph, seeds, 1))
  expect_identical(reach$indirect, ego_counts(graph, seeds, 2))

  # A duplicated seed counts twice, including at itself.
  expect_identical(reach$direct[2], sum(seeds %in% igraph::ego(graph, 1, 2)[[1]]))
  expect_gte(reach$direct[7], 2L)
  # The isolated seed only reaches itself.
  expect_identical(reach$direct[41], 1L)
  expect_identical(reach$indirect[41], 1L)

  expect_identical(.tkoi_seed_reach(csr$row_ptr, csr$col, seeds, 2L), reach)
})

test_that(".tkoi_seed_reach() handles more than 512 seeds (several bitset batches)", {
  graph = withr::with_seed(52, igraph::sample_gnm(400, 700))
  graph = igraph::add_edges(graph, c(3, 3, 3, 4, 3, 4))
  csr = graph_csr(graph)
  # 1100 seeds = batches of 512, 512, and 76 (a partial last word).
  seeds = withr::with_seed(53, sample(400L, 1100, replace = TRUE))
  seeds[c(1, 600, 1100)] = 3L

  reach = .tkoi_seed_reach(csr$row_ptr, csr$col, seeds, 2L)
  expect_identical(reach$direct, ego_counts(graph, seeds, 1))
  expect_identical(reach$indirect, ego_counts(graph, seeds, 2))
  expect_identical(.tkoi_seed_reach(csr$row_ptr, csr$col, seeds, 1L), reach)

  # Exactly 512 and 513 seeds.
  for (k in c(512L, 513L)) {
    reach = .tkoi_seed_reach(csr$row_ptr, csr$col, seeds[seq_len(k)], 1L)
    expect_identical(reach$direct, ego_counts(graph, seeds[seq_len(k)], 1))
    expect_identical(reach$indirect, ego_counts(graph, seeds[seq_len(k)], 2))
  }
})

test_that(".tkoi_seed_reach() matches igraph::ego() on the toy network", {
  graph = toy_network()
  csr = .tkoi_prepare_network(graph)$csr
  .tkoi_clear_cache()
  loop_vertex = igraph::ends(graph, which(igraph::which_loop(graph))[1], names = FALSE)[1]
  genes = which(igraph::V(graph)$labels == "['Gene']")
  seeds = c(loop_vertex, withr::with_seed(54, sample(genes, 30)))
  seeds = c(seeds, seeds[2:4])

  reach = .tkoi_seed_reach(csr$row_ptr, csr$col, seeds, 2L)
  expect_identical(reach$direct, ego_counts(graph, seeds, 1))
  expect_identical(reach$indirect, ego_counts(graph, seeds, 2))
})

test_that(".tkoi_seed_reach() handles no seeds and rejects out-of-range or missing seeds", {
  csr = graph_csr(igraph::make_ring(6))
  none = .tkoi_seed_reach(csr$row_ptr, csr$col, integer(0), 1L)
  expect_identical(none, list(direct = integer(6), indirect = integer(6)))

  message = "Seed vertex out of range."
  expect_error(.tkoi_seed_reach(csr$row_ptr, csr$col, 0L, 1L), message, fixed = TRUE)
  expect_error(.tkoi_seed_reach(csr$row_ptr, csr$col, 7L, 1L), message, fixed = TRUE)
  expect_error(.tkoi_seed_reach(csr$row_ptr, csr$col, NA_integer_, 1L), message, fixed = TRUE)
  expect_error(.tkoi_seed_reach(csr$row_ptr, csr$col, c(1L, NA_integer_), 1L), message, fixed = TRUE)
  expect_error(.tkoi_seed_reach(csr$row_ptr, csr$col, c(NA_integer_, 2L, 3L), 2L), message, fixed = TRUE)
})

# .tkoi_sample_null() --------------------------------------------------------------

sample_null_case = function(topology_similarity, seed = 61, n_universe = 150, n_genes = 25) {
  withr::with_seed(seed, {
    universe_degree = sample(c(1:30, 50, 51, 52, 100), n_universe, replace = TRUE)
    seed_degree = sample(universe_degree, n_genes, replace = TRUE)
  })
  # Three genes compete for the two degree-100 candidates.
  universe_degree[universe_degree == 100] = 99
  universe_degree[c(10, 20)] = 100
  seed_degree[1:3] = 100
  pools = .tkoi_candidate_pools(seed_degree, universe_degree, topology_similarity)
  list(pools = pools, n_universe = n_universe, candidates = unpack_pools(pools))
}

test_that(".tkoi_sample_null() reproduces the released sampling loop", {
  for (topology_similarity in c(0, 0.5, 0.9, 1)) {
    case = sample_null_case(topology_similarity)
    draws = withr::with_seed(71, .tkoi_sample_null(case$pools$ptr, case$pools$idx, case$n_universe, 12L))
    legacy = withr::with_seed(71, legacy_sample_null(case$candidates, 12L))
    expect_identical(draws, legacy)
    expect_identical(dim(draws), c(25L, 12L))
  }
})

test_that(".tkoi_sample_null() falls back to the full pool once it is exhausted", {
  # Degree 5 has three candidates for five genes, degree 10 two candidates
  # for three genes, and degree 20 a single candidate for two genes.
  universe_degree = c(5, 10, 5, 20, 5, 10)
  seed_degree = c(5, 10, 5, 5, 20, 10, 5, 5, 20, 10)
  pools = .tkoi_candidate_pools(seed_degree, universe_degree, 1)
  candidates = unpack_pools(pools)
  expect_identical(candidates[[1]], c(1L, 3L, 5L))
  expect_identical(candidates[[5]], 4L)

  draws = withr::with_seed(72, .tkoi_sample_null(pools$ptr, pools$idx, 6L, 20L))
  legacy = withr::with_seed(72, legacy_sample_null(candidates, 20L))
  expect_identical(draws, legacy)

  # Every draw comes from the gene's own pool.
  for (i in seq_along(candidates)) {
    expect_true(all(draws[i, ] %in% candidates[[i]]))
  }
  # The first three degree-5 genes take all three candidates; the single
  # degree-20 candidate is drawn for both degree-20 genes.
  expect_true(all(apply(draws[c(1, 3, 4), ], 2, sort) == c(1L, 3L, 5L)))
  expect_true(all(draws[c(5, 9), ] == 4L))
})

test_that(".tkoi_sample_null() leaves the RNG stream where the released loop leaves it", {
  # 20 permutations: the engine hands the RNG state back to R every 16.
  case = sample_null_case(0.9, seed = 62)
  after_new = withr::with_seed(73, {
    draws = .tkoi_sample_null(case$pools$ptr, case$pools$idx, case$n_universe, 20L)
    list(draws = draws, next_numbers = stats::runif(5))
  })
  after_legacy = withr::with_seed(73, {
    draws = legacy_sample_null(case$candidates, 20L)
    list(draws = draws, next_numbers = stats::runif(5))
  })
  expect_identical(after_new, after_legacy)

  # No permutations: nothing is drawn and the stream does not move.
  untouched = withr::with_seed(74, {
    draws = .tkoi_sample_null(case$pools$ptr, case$pools$idx, case$n_universe, 0L)
    list(dim = dim(draws), next_numbers = stats::runif(3))
  })
  expect_identical(untouched, list(dim = c(25L, 0L), next_numbers = withr::with_seed(74, stats::runif(3))))
})

test_that(".tkoi_sample_null() is reproducible and draws distinct genes when pools allow", {
  pools = .tkoi_candidate_pools(rep(10, 12), rep(10, 60), 0.9)
  first = withr::with_seed(75, .tkoi_sample_null(pools$ptr, pools$idx, 60L, 10L))
  second = withr::with_seed(75, .tkoi_sample_null(pools$ptr, pools$idx, 60L, 10L))
  other = withr::with_seed(76, .tkoi_sample_null(pools$ptr, pools$idx, 60L, 10L))
  expect_identical(first, second)
  expect_false(identical(first, other))
  expect_type(first, "integer")
  expect_true(all(first >= 1L & first <= 60L))
  expect_true(all(apply(first, 2, function(column) !anyDuplicated(column))))
})

test_that(".tkoi_sample_null() validates its inputs", {
  expect_error(.tkoi_sample_null(c(0L, 2L, 2L), 1:2, 5L, 3L), "non-empty")
  expect_error(.tkoi_sample_null(c(0L, 2L), 1:2, 5L, -1L), "Invalid dimensions")
  expect_error(.tkoi_sample_null(c(0L, 2L), 1:2, -1L, 3L), "Invalid dimensions")
  expect_error(.tkoi_sample_null(c(0L, 2L), 1:2, 5L, NA_integer_), "Invalid dimensions")
  expect_error(.tkoi_sample_null(c(0L, 2L), 1:2, NA_integer_, 3L), "Invalid dimensions")
  expect_error(.tkoi_sample_null(integer(0), integer(0), 5L, 3L), "Invalid dimensions")

  # No genes is a valid, empty request.
  expect_identical(dim(.tkoi_sample_null(0L, integer(0), 5L, 3L)), c(0L, 3L))
})

test_that(".tkoi_sample_null() rejects candidate indices that are missing or out of range", {
  message = "Candidate index out of range."
  expect_error(.tkoi_sample_null(c(0L, 2L), c(1L, 6L), 5L, 3L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 2L), c(0L, 1L), 5L, 3L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 2L), c(-2L, 1L), 5L, 3L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 1L), NA_integer_, 5L, 3L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 3L), c(1L, NA, 2L), 5L, 3L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 1L, 3L), c(1L, 2L, NA), 5L, 3L), message, fixed = TRUE)
})

test_that(".tkoi_sample_null() rejects offsets that do not span cand_idx", {
  # Offsets must start at 0 and end at length(cand_idx), as they do for
  # .tkoi_candidate_pools(). Offsets past the end of cand_idx used to make the
  # engine read and write out of bounds; they now stop with an error.
  message = "Malformed 'cand_ptr'."
  expect_error(.tkoi_sample_null(c(2L, 3L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(1L, 3L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 2L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 4L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 1000000L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, 1L, 1000000L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(NA, 3L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(0L, NA), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(c(-1L, 3L), 1:3, 3L, 4L), message, fixed = TRUE)
  expect_error(.tkoi_sample_null(0L, 1:3, 3L, 4L), message, fixed = TRUE)

  # Offsets that start and end right but do not increase are rejected too.
  expect_error(.tkoi_sample_null(c(0L, 1000000L, 3L), 1:3, 3L, 4L), "non-empty")
  expect_error(.tkoi_sample_null(c(0L, NA, 3L), 1:3, 3L, 4L), "non-empty")
  expect_error(.tkoi_sample_null(c(0L, 2L, 1L, 3L), 1:3, 3L, 4L), "non-empty")
})

# .tkoi_candidate_pools() ----------------------------------------------------------

test_that(".tkoi_candidate_pools() matches a brute-force search", {
  # Every degree from 1 to 60 occurs, so even exact matching finds candidates.
  universe_degree = withr::with_seed(81, sample(c(1:60, sample(1:60, 240, replace = TRUE))))
  seed_degree = c(10, 10, 1, 60, 33, withr::with_seed(82, sample(1:60, 20, replace = TRUE)))
  for (topology_similarity in c(0, 0.25, 0.5, 0.9, 1)) {
    pools = .tkoi_candidate_pools(seed_degree, universe_degree, topology_similarity)
    expect_type(pools$ptr, "integer")
    expect_type(pools$idx, "integer")
    expect_identical(pools$ptr[1], 0L)
    expect_length(pools$ptr, length(seed_degree) + 1)
    expect_identical(pools$ptr[length(pools$ptr)], length(pools$idx))
    expect_identical(
      unpack_pools(pools),
      brute_force_pools(seed_degree, universe_degree, topology_similarity)
    )
  }
})

test_that(".tkoi_candidate_pools() includes the exact boundary degrees", {
  universe_degree = c(12, 8, 9, 10, 11, 8.9, 11.1, 10)
  pools = .tkoi_candidate_pools(10, universe_degree, 0.9)
  # [0.9 * 10, 1.1 * 10] = [9, 11], in universe order.
  expect_identical(pools$idx, c(3L, 4L, 5L, 8L))
  expect_identical(pools$ptr, c(0L, 4L))

  expect_identical(.tkoi_candidate_pools(10, universe_degree, 1)$idx, c(4L, 8L))
  expect_identical(.tkoi_candidate_pools(10, universe_degree, 0)$idx, 1:8)
  expect_identical(.tkoi_candidate_pools(10, universe_degree, 0.8)$idx, 1:8)
})

test_that(".tkoi_candidate_pools() stops when a seed gene has no candidates", {
  expect_error(
    .tkoi_candidate_pools(c(10, 1000), c(9, 10, 11), 0.9),
    "No degree-matched candidates"
  )
  expect_error(.tkoi_candidate_pools(5, c(9, 10, 11), 0.9), "lower `topology_similarity`", fixed = TRUE)

  empty = .tkoi_candidate_pools(numeric(0), c(9, 10, 11), 0.9)
  expect_identical(empty$ptr, 0L)
  expect_identical(empty$idx, integer(0))
})

test_that(".tkoi_candidate_pools() output is accepted by .tkoi_sample_null()", {
  pools = .tkoi_candidate_pools(c(3, 3, 7, 12), c(3, 4, 7, 7, 8, 12, 12, 3), 0.8)
  draws = withr::with_seed(83, .tkoi_sample_null(pools$ptr, pools$idx, 8L, 6L))
  expect_identical(dim(draws), c(4L, 6L))
  candidates = unpack_pools(pools)
  for (i in seq_along(candidates)) {
    expect_true(all(draws[i, ] %in% candidates[[i]]))
  }
})

# .tkoi_resolve_cores() and .tkoi_available_cores() ---------------------------------

test_that(".tkoi_resolve_cores() uses every core or the tkoi.n_cores option", {
  local_unlimited_cores()
  withr::local_options(tkoi.n_cores = NULL)
  expect_identical(.tkoi_resolve_cores(), detected_cores())
  expect_identical(.tkoi_resolve_cores(NULL), detected_cores())
  expect_identical(.tkoi_resolve_cores(1), 1L)
  expect_identical(.tkoi_resolve_cores(1L), 1L)

  withr::local_options(tkoi.n_cores = 1)
  expect_identical(.tkoi_resolve_cores(), 1L)
  expect_identical(.tkoi_resolve_cores(NULL), 1L)
  # An explicit value wins over the option.
  expect_identical(.tkoi_resolve_cores(2), min(2L, detected_cores()))

  withr::local_options(tkoi.n_cores = 1e6)
  expect_identical(.tkoi_resolve_cores(), detected_cores())
})

test_that(".tkoi_resolve_cores() caps requests at the available cores", {
  local_unlimited_cores()
  expect_identical(.tkoi_resolve_cores(detected_cores() + 5), detected_cores())
  expect_identical(.tkoi_resolve_cores(1e6), detected_cores())
  expect_identical(.tkoi_resolve_cores(detected_cores()), detected_cores())
})

test_that(".tkoi_resolve_cores() respects _R_CHECK_LIMIT_CORES_", {
  local_unlimited_cores()
  withr::local_options(tkoi.n_cores = NULL)
  for (value in c("TRUE", "true", "warn")) {
    withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = value))
    expect_identical(.tkoi_resolve_cores(), min(2L, detected_cores()))
    expect_identical(.tkoi_resolve_cores(64), min(2L, detected_cores()))
    expect_identical(.tkoi_resolve_cores(1), 1L)
  }
  for (value in c("FALSE", "false", "")) {
    withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = value))
    expect_identical(.tkoi_resolve_cores(), detected_cores())
  }
})

test_that(".tkoi_resolve_cores() rejects invalid values", {
  local_unlimited_cores()
  bad_values = list(0, -1, 1.5, NA, NA_real_, NA_integer_, NaN, "2", c(1, 2), TRUE, list(2), numeric(0))
  for (value in bad_values) {
    expect_error(.tkoi_resolve_cores(value), "single positive whole number")
  }
  withr::local_options(tkoi.n_cores = 0)
  expect_error(.tkoi_resolve_cores(), "single positive whole number")
  withr::local_options(tkoi.n_cores = "four")
  expect_error(.tkoi_resolve_cores(), "single positive whole number")
})

test_that(".tkoi_available_cores() respects a Slurm CPU allocation", {
  local_unlimited_cores()
  withr::local_options(tkoi.n_cores = NULL)
  cores = detected_cores()
  expect_identical(.tkoi_available_cores(), cores)

  withr::local_envvar(SLURM_CPUS_PER_TASK = "1")
  expect_identical(.tkoi_available_cores(), 1L)
  expect_identical(.tkoi_resolve_cores(), 1L)
  expect_identical(.tkoi_resolve_cores(8), 1L)

  for (value in c("2", "2.7", " 2 ")) {
    withr::local_envvar(SLURM_CPUS_PER_TASK = value)
    expect_identical(.tkoi_available_cores(), min(2L, cores), label = value)
  }
  withr::local_envvar(SLURM_CPUS_PER_TASK = "100000")
  expect_identical(.tkoi_available_cores(), cores)

  # Values that are not a CPU count are ignored.
  for (value in c("", "0", "0.5", "-4", "abc", "NA")) {
    withr::local_envvar(SLURM_CPUS_PER_TASK = value)
    expect_identical(.tkoi_available_cores(), cores, label = value)
  }

  # R CMD check's limit applies on top of the allocation.
  withr::local_envvar(c(SLURM_CPUS_PER_TASK = "8", "_R_CHECK_LIMIT_CORES_" = "TRUE"))
  expect_identical(.tkoi_available_cores(), min(2L, cores))
})

test_that(".tkoi_available_cores() respects a cgroup CPU quota", {
  local_unlimited_cores()
  cores = detected_cores()
  available = function(cpu_max) {
    with_system_files(list("/sys/fs/cgroup/cpu.max" = cpu_max), .tkoi_available_cores())
  }

  # "max" means no quota.
  expect_identical(available("max 100000"), cores)
  expect_identical(available("200000 100000"), min(2L, cores))
  expect_identical(available("400000 100000"), min(4L, cores))
  # Fractional quotas round down, but never below one core.
  expect_identical(available("150000 100000"), 1L)
  expect_identical(available("50000 100000"), 1L)
  expect_identical(available("1 100000"), 1L)
  expect_identical(available("100000000 100000"), cores)
  # Unreadable quotas are ignored.
  for (cpu_max in c("", "garbage", "100000", "100000 0", "100000 abc")) {
    expect_identical(available(cpu_max), cores, label = cpu_max)
  }

  # Only the cpu.max file is consulted for cores.
  other = with_system_files(list("/sys/fs/cgroup/memory.max" = "1 1"), .tkoi_available_cores())
  expect_identical(other, cores)

  # The smallest of the quota, the Slurm allocation, and the hardware wins.
  withr::local_envvar(SLURM_CPUS_PER_TASK = "3")
  expect_identical(available("400000 100000"), min(3L, cores))
  expect_identical(available("200000 100000"), min(2L, cores))
  resolved = with_system_files(list("/sys/fs/cgroup/cpu.max" = "200000 100000"), .tkoi_resolve_cores(16))
  expect_identical(resolved, min(2L, cores))
})

# .tkoi_read_first_line() and .tkoi_memory_limit() -------------------------------------

test_that(".tkoi_read_first_line() returns the trimmed first line, or NULL", {
  expect_null(.tkoi_read_first_line(file.path(tempdir(), "tkoi-no-such-file")))

  path = withr::local_tempfile()
  writeLines(c("  max 100000  ", "second line"), path)
  expect_identical(.tkoi_read_first_line(path), "max 100000")

  empty = withr::local_tempfile()
  file.create(empty)
  expect_null(.tkoi_read_first_line(empty))

  # A file without a final newline is read without a warning.
  unterminated = withr::local_tempfile()
  cat("8000000000", file = unterminated)
  expect_no_warning(expect_identical(.tkoi_read_first_line(unterminated), "8000000000"))
})

test_that(".tkoi_memory_limit() lowers physical memory to a cgroup limit", {
  local_mocked_bindings(.tkoi_total_memory = function() 16e9, .package = "tkoi")
  v2 = "/sys/fs/cgroup/memory.max"
  v1 = "/sys/fs/cgroup/memory/memory.limit_in_bytes"
  limit = function(files) with_system_files(files, .tkoi_memory_limit())

  expect_identical(limit(list()), 16e9)
  # cgroup v2: "max" means no limit.
  expect_identical(limit(stats::setNames(list("max"), v2)), 16e9)
  expect_identical(limit(stats::setNames(list("8000000000"), v2)), 8e9)
  expect_identical(limit(stats::setNames(list("32000000000"), v2)), 16e9)
  # cgroup v1 reports "no limit" as a huge number.
  expect_identical(limit(stats::setNames(list("4000000000"), v1)), 4e9)
  expect_identical(limit(stats::setNames(list("9223372036854771712"), v1)), 16e9)
  # The lower of the two limits wins.
  expect_identical(limit(stats::setNames(list("8000000000", "4000000000"), c(v2, v1))), 4e9)
  expect_identical(limit(stats::setNames(list("2000000000", "4000000000"), c(v2, v1))), 2e9)
  # Unusable values are ignored.
  for (value in c("0", "-1", "abc", "")) {
    expect_identical(limit(stats::setNames(list(value), v2)), 16e9, label = value)
  }
})

test_that(".tkoi_memory_limit() falls back to the cgroup limit when physical memory is unknown", {
  local_mocked_bindings(.tkoi_total_memory = function() NA_real_, .package = "tkoi")
  expect_identical(with_system_files(list(), .tkoi_memory_limit()), NA_real_)
  expect_identical(
    with_system_files(list("/sys/fs/cgroup/memory.max" = "8000000000"), .tkoi_memory_limit()),
    8e9
  )
  expect_identical(
    with_system_files(list("/sys/fs/cgroup/memory.max" = "max"), .tkoi_memory_limit()),
    NA_real_
  )
})

test_that(".tkoi_memory_limit() reports this machine's memory", {
  limit = .tkoi_memory_limit()
  total = .tkoi_total_memory()
  expect_type(limit, "double")
  expect_length(limit, 1)
  expect_true(is.na(limit) || limit > 0)
  if (!is.na(total)) {
    expect_lte(limit, total)
  }
})

# .tkoi_plan_memory() --------------------------------------------------------------

test_that(".tkoi_plan_memory() uses the widest block when memory is ample", {
  local_memory_limit(16e9)
  expect_identical(.tkoi_plan_memory(5260, 817900, 11, 10, TRUE), 16L)
  expect_identical(.tkoi_plan_memory(5260, 817900, 11, 10, FALSE), 16L)
  expect_identical(.tkoi_plan_memory(1, 0, 1, 2, TRUE), 16L)
  expect_identical(.tkoi_plan_memory(1, 0, 0, 0, FALSE), 16L)
})

test_that(".tkoi_plan_memory() uses the widest block when memory is unknown", {
  local_memory_limit(NA_real_)
  expect_identical(.tkoi_plan_memory(1e15, 1e15, 1e4, 1e6, TRUE), 16L)
  expect_identical(.tkoi_plan_memory(1e15, 1e15, 1e4, 1e6, FALSE), 16L)
})

test_that(".tkoi_plan_memory() fits the toy network in this machine's memory", {
  limit = .tkoi_memory_limit()
  skip_if(is.na(limit) || limit < 2e9, "less than 2 GB of memory is available")
  expect_identical(.tkoi_plan_memory(5260, 817900, 11, 10, TRUE), 16L)
})

test_that(".tkoi_plan_memory() narrows the block to fit the budget", {
  # 1,000 vertices, no edges or seeds, keep_permutations = FALSE: the plan
  # needs 56,000 + 40,000 * width bytes of a budget of 3/4 of the memory.
  plan = function(limit) {
    local_memory_limit(limit)
    .tkoi_plan_memory(1000, 0, 0, 5, FALSE)
  }
  expect_identical(plan(928000), 16L)
  expect_identical(plan(927996), 8L)
  expect_identical(plan(501336), 8L)
  expect_identical(plan(501332), 4L)
  expect_identical(plan(288000), 4L)
  expect_identical(plan(287996), 2L)
  expect_identical(plan(181336), 2L)
  expect_identical(plan(181332), 1L)
  expect_identical(plan(128000), 1L)
  expect_error(plan(127996), "smaller network")
})

test_that(".tkoi_plan_memory() counts the network, the null draws, and the kept permutations", {
  # 1,000 vertices, 5,000 matrix entries, 10 seeds, 20 permutations:
  # network 84,000 bytes, draws 2,400 bytes, outputs 168,000 bytes when every
  # permutation is kept (32,000 otherwise), 40,000 bytes per block column.
  local_memory_limit(392536)
  expect_identical(.tkoi_plan_memory(1000, 5000, 10, 20, TRUE), 1L)
  expect_identical(.tkoi_plan_memory(1000, 5000, 10, 20, FALSE), 4L)

  local_memory_limit(392532)
  expect_error(.tkoi_plan_memory(1000, 5000, 10, 20, TRUE), "keep_permutations = FALSE", fixed = TRUE)
  expect_identical(.tkoi_plan_memory(1000, 5000, 10, 20, FALSE), 4L)
})

test_that(".tkoi_plan_memory() stops with advice when the outputs cannot fit", {
  local_memory_limit(8e9)
  expect_error(
    .tkoi_plan_memory(1e8, 1e9, 1000, 100, TRUE),
    paste(
      "run_tkoi() needs about 99.2 GB but only 8.0 GB of memory is available.",
      "Set `keep_permutations = FALSE` or reduce `n_permutation`."
    ),
    fixed = TRUE
  )
  expect_error(
    .tkoi_plan_memory(1e8, 1e9, 1000, 100, FALSE),
    paste(
      "run_tkoi() needs about 21.6 GB but only 8.0 GB of memory is available.",
      "Reduce `n_permutation` or use a smaller network."
    ),
    fixed = TRUE
  )
  expect_error(.tkoi_plan_memory(1e15, 1e15, 101, 100, TRUE), "reduce `n_permutation`", fixed = TRUE)
  expect_error(.tkoi_plan_memory(1e15, 1e15, 101, 100, FALSE), "smaller network", fixed = TRUE)
})

test_that(".tkoi_plan_memory() plans within a cgroup memory limit", {
  local_mocked_bindings(.tkoi_total_memory = function() 16e9, .package = "tkoi")
  plan = function(memory_max) {
    with_system_files(
      list("/sys/fs/cgroup/memory.max" = memory_max),
      .tkoi_plan_memory(1000, 0, 0, 5, FALSE)
    )
  }
  expect_identical(plan("max"), 16L)
  expect_identical(plan("288000"), 4L)
  expect_error(plan("127996"), "smaller network")
})

test_that("the engine reports memory and hardware threads", {
  total = .tkoi_total_memory()
  expect_type(total, "double")
  expect_true(is.na(total) || total > 0)
  threads = .tkoi_hardware_threads()
  expect_type(threads, "integer")
  expect_gte(threads, 0L)
})

# .tkoi_clean_node_type() and .tkoi_vertex_types() -------------------------------

test_that(".tkoi_clean_node_type() strips Neo4j label formatting", {
  expect_identical(.tkoi_clean_node_type("['Gene']"), "Gene")
  expect_identical(
    .tkoi_clean_node_type(c("['BiologicalProcess']", "Disease", "[Pathway]", "'Compound'")),
    c("BiologicalProcess", "Disease", "Pathway", "Compound")
  )
  expect_identical(.tkoi_clean_node_type(factor(c("['Anatomy']", "['CellType']"))), c("Anatomy", "CellType"))
  expect_identical(.tkoi_clean_node_type(character(0)), character(0))
})

test_that(".tkoi_clean_node_type() names missing and empty types Unknown", {
  expect_identical(.tkoi_clean_node_type(c("['Gene']", NA)), c("Gene", "Unknown"))
  expect_identical(.tkoi_clean_node_type(c("", "[]", "['']", "''")), rep("Unknown", 4))
  expect_identical(.tkoi_clean_node_type(NA), "Unknown")
  expect_identical(.tkoi_clean_node_type(factor(c("['Anatomy']", NA))), c("Anatomy", "Unknown"))
})

test_that(".tkoi_vertex_types() reads `labels`, then `label`", {
  graph = igraph::make_ring(3)
  expect_null(.tkoi_vertex_types(graph))

  igraph::V(graph)$label = c("['Gene']", NA, "Disease")
  expect_identical(.tkoi_vertex_types(graph), c("Gene", "Unknown", "Disease"))

  igraph::V(graph)$labels = c("['Pathway']", "", "['Compound']")
  expect_identical(.tkoi_vertex_types(graph), c("Pathway", "Unknown", "Compound"))
})

# .tkoi_clean_expression() -----------------------------------------------------------

clean_expression_input = function() {
  data.frame(
    gene_name = c("ENSG1", "ENSG2", NA, "", "ENSG1", "ENSG3", "ENSG2", ""),
    logfc = c(1, -2, 3, 4, 5, 6, 7, 8),
    pvalue = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08),
    extra = letters[1:8],
    stringsAsFactors = FALSE
  )
}

clean_expression_expected = function() {
  data.frame(
    gene_name = c("ENSG1", "ENSG2", "ENSG3"),
    logfc = c(1, -2, 6),
    pvalue = c(0.01, 0.02, 0.06),
    stringsAsFactors = FALSE
  )
}

test_that(".tkoi_clean_expression() drops blank genes and keeps the first row of each gene", {
  cleaned = .tkoi_clean_expression(clean_expression_input())
  expect_identical(cleaned, clean_expression_expected())
  expect_s3_class(cleaned, "data.frame", exact = TRUE)
  expect_identical(rownames(cleaned), c("1", "2", "3"))
})

test_that(".tkoi_clean_expression() converts factors and integers", {
  expression = clean_expression_input()
  expression$gene_name = factor(expression$gene_name)
  expression$logfc = as.integer(expression$logfc)
  cleaned = .tkoi_clean_expression(expression)
  expect_type(cleaned$gene_name, "character")
  expect_type(cleaned$logfc, "double")
  expect_identical(cleaned, clean_expression_expected())
})

test_that(".tkoi_clean_expression() accepts tibbles and data.tables", {
  expected = clean_expression_expected()
  expect_identical(.tkoi_clean_expression(dplyr::as_tibble(clean_expression_input())), expected)

  skip_if_not_installed("data.table")
  expect_identical(.tkoi_clean_expression(data.table::as.data.table(clean_expression_input())), expected)
})

test_that(".tkoi_clean_expression() handles tables without usable genes", {
  empty = .tkoi_clean_expression(clean_expression_input()[0, ])
  expect_named(empty, c("gene_name", "logfc", "pvalue"))
  expect_identical(nrow(empty), 0L)

  blank = .tkoi_clean_expression(clean_expression_input()[c(3, 4, 8), ])
  expect_named(blank, c("gene_name", "logfc", "pvalue"))
  expect_identical(nrow(blank), 0L)
})

# .tkoi_prepare_network() ----------------------------------------------------------

test_that(".tkoi_prepare_network() rejects unusable networks", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())

  expect_error(.tkoi_prepare_network(as.raw(1:10)), "Encrypted networks are no longer supported")
  expect_error(.tkoi_prepare_network(data.frame(from = 1, to = 2)), "must be an igraph object")
  expect_error(.tkoi_prepare_network(list()), "must be an igraph object")
  expect_error(.tkoi_prepare_network(NULL), "must be an igraph object")
  expect_error(.tkoi_prepare_network(igraph::make_empty_graph(0, directed = FALSE)), "no vertices")

  unnamed = igraph::make_ring(4)
  igraph::V(unnamed)$labels = "['Gene']"
  expect_error(.tkoi_prepare_network(unnamed), "vertex attribute `name`", fixed = TRUE)

  unlabelled = igraph::make_ring(4)
  igraph::V(unlabelled)$name = letters[1:4]
  expect_error(.tkoi_prepare_network(unlabelled), "vertex attribute `labels`", fixed = TRUE)
})

test_that(".tkoi_prepare_network() extracts the matrix and vertex metadata", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = small_network()
  prepared = .tkoi_prepare_network(graph)

  expect_named(prepared, c("csr", "vertex_names", "node_type", "identifier", "valency", "n"))
  expect_equal(prepared$n, 5)
  expect_identical(prepared$vertex_names, paste0("node", 1:5))
  expect_identical(prepared$node_type, c("Gene", "Gene", "Disease", "Pathway", "Gene"))
  expect_identical(prepared$identifier, paste0("ID:", 1:5))
  expect_identical(prepared$valency, c(10, 20, 30, 40, 50))
  expect_identical(prepared$csr, small_network_csr())
})

test_that(".tkoi_prepare_network() uses the weight edge attribute", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = small_network()
  igraph::E(graph)$weight = c(1, 2, 0.5, 3, 4, 2)
  prepared = .tkoi_prepare_network(graph)

  expect_identical(prepared$csr, small_network_csr(c(1, 2, 0.5, 3, 4, 2)))
  expect_false(isTRUE(all.equal(prepared$csr$val, small_network_csr()$val)))

  # PageRank from the prepared matrix follows igraph's weighted PageRank.
  seeds = seed_sets(c(1, 3), c(1, 2), list(c(2, 5)), c(1, 1))
  result = run_ppr(prepared$csr, seeds)
  reference = igraph_ppr(graph, seeds)
  expect_equal(cbind(result$observed, null_matrix(result)), reference, tolerance = 1e-10)
})

test_that(".tkoi_prepare_network() falls back to `label`, vertex names, and igraph::degree()", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = igraph::make_graph(c(1, 2, 2, 3, 3, 3, 1, 3), n = 4, directed = FALSE)
  igraph::V(graph)$name = c("a", "b", "c", "d")
  igraph::V(graph)$label = c("['Gene']", "['Compound']", "Disease", NA)
  prepared = .tkoi_prepare_network(graph)

  expect_identical(prepared$node_type, c("Gene", "Compound", "Disease", "Unknown"))
  expect_identical(prepared$identifier, c("a", "b", "c", "d"))
  expect_equal(prepared$valency, igraph::degree(graph))
  # The loop counts twice in the degree.
  expect_equal(as.numeric(prepared$valency), c(2, 2, 4, 0))

  # `labels` takes precedence over `label`.
  igraph::V(graph)$labels = "['Pathway']"
  expect_identical(.tkoi_prepare_network(graph)$node_type, rep("Pathway", 4))
})

test_that(".tkoi_prepare_network() reuses the cached matrix for the same graph", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = small_network()
  original_csr = small_network_csr()
  extended = igraph::add_edges(graph, c(1, 5))
  expect_false(identical(igraph::graph_id(extended), igraph::graph_id(graph)))
  edges = igraph::as_edgelist(extended, names = FALSE)
  extended_csr = .tkoi_build_csr(edges[, 1], edges[, 2], numeric(0), 5L)
  counter = local_build_counter()

  first = .tkoi_prepare_network(graph)
  expect_identical(counter$calls, 1L)
  expect_identical(first$csr, original_csr)
  expect_identical(.tkoi_cache$network$key, igraph::graph_id(graph))
  expect_identical(.tkoi_cache$network$csr, first$csr)

  # The same graph, or a copy of it, is served from the cache.
  expect_identical(.tkoi_prepare_network(graph), first)
  copy = graph
  expect_identical(.tkoi_prepare_network(copy), first)
  expect_identical(counter$calls, 1L)

  # A structurally different graph is rebuilt and replaces the cache entry.
  rebuilt = .tkoi_prepare_network(extended)
  expect_identical(counter$calls, 2L)
  expect_identical(rebuilt$csr, extended_csr)
  expect_identical(.tkoi_cache$network$key, igraph::graph_id(extended))

  # Only the most recent network is kept.
  expect_identical(.tkoi_prepare_network(graph), first)
  expect_identical(counter$calls, 3L)

  # Clearing the cache empties it.
  .tkoi_clear_cache()
  expect_identical(ls(.tkoi_cache), character(0))
  expect_identical(.tkoi_prepare_network(graph), first)
  expect_identical(counter$calls, 4L)
})

test_that(".tkoi_prepare_network() rebuilds the matrix when edge weights change", {
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = small_network()
  weight = c(1, 2, 0.5, 3, 4, 2)
  unweighted_csr = small_network_csr()
  weighted_csr = small_network_csr(weight)
  nudged_csr = small_network_csr(replace(weight, 2, 2 + 1e-12))
  counter = local_build_counter()

  expect_identical(.tkoi_prepare_network(graph)$csr, unweighted_csr)
  expect_identical(counter$calls, 1L)

  # igraph keeps the graph ID when only attributes change, so the cache must
  # compare the weights too.
  reweighted = graph
  igraph::E(reweighted)$weight = weight
  expect_identical(igraph::graph_id(reweighted), igraph::graph_id(graph))
  expect_identical(.tkoi_prepare_network(reweighted)$csr, weighted_csr)
  expect_identical(counter$calls, 2L)
  expect_identical(.tkoi_prepare_network(reweighted)$csr, weighted_csr)
  expect_identical(counter$calls, 2L)

  # Even a tiny change is noticed.
  igraph::E(reweighted)$weight[2] = 2 + 1e-12
  expect_identical(.tkoi_prepare_network(reweighted)$csr, nudged_csr)
  expect_identical(counter$calls, 3L)

  # Removing the weights rebuilds the unweighted matrix.
  expect_identical(.tkoi_prepare_network(igraph::delete_edge_attr(reweighted, "weight"))$csr, unweighted_csr)
  expect_identical(counter$calls, 4L)
})

test_that(".tkoi_prepare_network() reads vertex metadata afresh when the matrix is cached", {
  # These attributes do not change the graph ID, and they used to be served
  # stale from the cache.
  .tkoi_clear_cache()
  withr::defer(.tkoi_clear_cache())
  graph = small_network()
  counter = local_build_counter()
  first = .tkoi_prepare_network(graph)

  edited = graph
  igraph::V(edited)$labels = c("['Disease']", NA, "", "['Gene']", "['Compound']")
  prepared = .tkoi_prepare_network(edited)
  expect_identical(prepared$node_type, c("Disease", "Unknown", "Unknown", "Gene", "Compound"))

  igraph::V(edited)$name = paste0("renamed", 1:5)
  prepared = .tkoi_prepare_network(edited)
  expect_identical(prepared$vertex_names, paste0("renamed", 1:5))

  igraph::V(edited)$identifier = paste0("NEW:", 1:5)
  prepared = .tkoi_prepare_network(edited)
  expect_identical(prepared$identifier, paste0("NEW:", 1:5))

  igraph::V(edited)$degree = c(5, 4, 3, 2, 1)
  prepared = .tkoi_prepare_network(edited)
  expect_identical(prepared$valency, c(5, 4, 3, 2, 1))

  # Every call above reused the cached matrix.
  expect_identical(igraph::graph_id(edited), igraph::graph_id(graph))
  expect_identical(counter$calls, 1L)
  expect_identical(prepared$csr, first$csr)
  expect_identical(prepared$node_type, c("Disease", "Unknown", "Unknown", "Gene", "Compound"))
  expect_identical(prepared$vertex_names, paste0("renamed", 1:5))
  expect_identical(prepared$identifier, paste0("NEW:", 1:5))

  # Without `identifier` and `degree`, names and igraph::degree() are used.
  bare = edited |>
    igraph::delete_vertex_attr("identifier") |>
    igraph::delete_vertex_attr("degree")
  prepared = .tkoi_prepare_network(bare)
  expect_identical(prepared$identifier, paste0("renamed", 1:5))
  expect_equal(prepared$valency, igraph::degree(bare))

  # The original graph still gives its own metadata.
  expect_identical(.tkoi_prepare_network(graph), first)
})
