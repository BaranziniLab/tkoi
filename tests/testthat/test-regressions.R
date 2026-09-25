# Targeted checks for behaviours that plausible bugs could break without
# failing any other test (found by planting bugs and running the suite).

regression_result = local({
  cache = NULL
  function() {
    if (is.null(cache)) {
      cache <<- toy_run()
    }
    cache
  }
})

# Ensembl IDs of toy-network genes that have exactly one Ensembl ID.
toy_ensembl = function(n) {
  universe = toy_universe()
  unique_ensembl = universe$ensembl[!duplicated(universe$ensembl) &
    !universe$ensembl %in% universe$ensembl[duplicated(universe$ensembl)] &
    nzchar(universe$ensembl)]
  utils::head(unique_ensembl, n)
}

test_that("seed thresholds are inclusive and rows with missing values are not seeds", {
  ensembl = toy_ensembl(6)
  expression = data.frame(
    gene_name = ensembl,
    logfc = c(0.25, -0.25, 1, 0.2499999, 1, NA),
    pvalue = c(0.05, 0.05, 0.0500001, 0.01, NA, 0.01)
  )
  seeds = .tkoi_seed_genes(expression, pvalue_threshold = 0.05, logfc_threshold = 0.25)
  expected = tkoi::genes$id[match(ensembl[1:2], tkoi::genes$ensembl)]
  expect_setequal(seeds$id, expected)
})

test_that("a gene whose first row has a missing p-value is not a seed", {
  ensembl = toy_ensembl(2)
  expression = data.frame(
    gene_name = c(ensembl[1], ensembl[1], ensembl[2]),
    logfc = c(2, 2, 2),
    pvalue = c(NA, 1e-6, 1e-6)
  )
  seeds = .tkoi_seed_genes(expression, pvalue_threshold = 0.05, logfc_threshold = 0.25)
  expect_identical(seeds$id, tkoi::genes$id[match(ensembl[2], tkoi::genes$ensembl)])
})

test_that("p-values keep full relative precision in the far upper tail", {
  tables = regression_result()@network_summary_statistics
  beta = unlist(lapply(tables, function(x) x$beta))
  p_value = unlist(lapply(tables, function(x) x$p_value))
  # Beyond beta of about 38 the p-value is below the smallest double and is 0.
  representable = is.finite(beta) & stats::pnorm(beta, lower.tail = FALSE, log.p = TRUE) > -700
  expect_true(any(p_value[representable] < 1e-16))
  expect_true(all(p_value[representable] > 0))
  log_expected = stats::pnorm(beta[representable], lower.tail = FALSE, log.p = TRUE)
  expect_lt(max(abs(log(p_value[representable]) - log_expected)), 1e-10)
})

test_that("the 1-norm error of a PageRank vector is about 2 * tolerance even with extreme weights", {
  withr::local_seed(11)
  g = igraph::sample_gnm(400, 2400)
  weight = exp(stats::runif(igraph::ecount(g), log(1e-6), log(1e6)))
  edges = igraph::as_edgelist(g, names = FALSE)
  csr = .tkoi_build_csr(edges[, 1], edges[, 2], weight, 400L)
  seeds = c(3L, 50L, 199L)
  seed_weight = c(0.2, 0.5, 0.3)
  tol = 1e-12
  result = .tkoi_ppr_null(
    csr, seeds, seed_weight, matrix(integer(0), 3, 0), seed_weight,
    0.85, tol, 5000L, 1L, 1L, FALSE
  )

  # Dense reference: (I - d A D^-1) z = u, normalized.
  a = matrix(0, 400, 400)
  for (e in seq_len(nrow(edges))) {
    i = edges[e, 1]
    j = edges[e, 2]
    a[i, j] = a[i, j] + weight[e]
    a[j, i] = a[j, i] + weight[e]
  }
  degree = colSums(a)
  transition = sweep(a, 2, ifelse(degree > 0, degree, 1), "/")
  u = numeric(400)
  u[seeds] = seed_weight
  z = solve(diag(400) - 0.85 * transition, u)
  reference = z / sum(z)

  expect_lte(result$rel_residual, tol)
  expect_lt(sum(abs(result$observed - reference)), 10 * tol)
})

test_that("a gene with logfc exactly 0 is drawn in the up-regulated panel", {
  result = regression_result()
  expression = as.data.frame(result@expression_data)
  gene_ids = result@network_summary_statistics$Gene$node_id
  in_table = expression$gene_name %in% tkoi::genes$ensembl[tkoi::genes$id %in% gene_ids]
  target = expression$gene_name[in_table][1]
  expression$logfc[expression$gene_name == target] = 0
  result@expression_data = expression

  plot_data = make_gene_exploration_plot(result)$data
  expect_identical(unique(plot_data$direction[plot_data$gene_name == target]), "Up-regulated")
})

test_that("plot_network with degree_expansion = 1 adds no bridging nodes", {
  result = regression_result()
  net = toy_network()
  expression = .tkoi_clean_expression(result@expression_data)
  significant = expression$gene_name[expression$pvalue <= result@pvalue_threshold &
    abs(expression$logfc) >= result@logfc_threshold]
  sig_ids = intersect(tkoi::genes$id[tkoi::genes$ensembl %in% significant], igraph::V(net)$name)

  # A target adjacent to a significant gene that also shares a neighbour with
  # it, so a two-hop path exists that must not be drawn at one hop.
  candidates = result@network_summary_statistics$BiologicalProcess$node_id
  target = NULL
  for (node in candidates) {
    neighbours = names(igraph::neighbors(net, node, mode = "all"))
    adjacent_genes = intersect(neighbours, sig_ids)
    if (length(adjacent_genes) == 0) next
    shared = intersect(neighbours, names(igraph::neighbors(net, adjacent_genes[1], mode = "all")))
    if (length(setdiff(shared, c(node, adjacent_genes))) > 0) {
      target = node
      break
    }
  }
  skip_if(is.null(target), "no suitable target in the toy network")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  subnet = plot_network(result, target, degree_expansion = 1, subnetwork = net)
  adjacent = intersect(names(igraph::neighbors(net, target, mode = "all")), sig_ids)
  expect_equal(igraph::vcount(subnet), length(unique(c(target, adjacent))))
})

test_that("the non-convergence warning counts every unconverged PageRank vector", {
  expect_warning(
    toy_run(maximum_iteration = 1, n_permutation = 6),
    "^7 PageRank vector\\(s\\) did not reach"
  )
})
