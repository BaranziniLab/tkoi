test_that("analysis retains custom topology and attributes after saving", {
  graph = toy_network()
  graph = igraph::delete_edges(graph, 1)
  igraph::E(graph)$weight = rep(2, igraph::ecount(graph))
  result = withr::with_seed(42, run_tkoi(
    toy_expression(), subnetwork = graph, n_permutation = 2,
    n_cores = 1, verbose = FALSE
  ))
  expect_identical(get_analysis_graph(result), graph)
  path = tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(result, path)
  restored = get_analysis_graph(readRDS(path))
  expect_identical(igraph::as_edgelist(restored), igraph::as_edgelist(graph))
  expect_identical(igraph::vertex_attr(restored), igraph::vertex_attr(graph))
  expect_identical(igraph::edge_attr(restored), igraph::edge_attr(graph))
  igraph::E(graph)$weight = 3
  expect_true(all(igraph::E(get_analysis_graph(result))$weight == 2))
})

test_that("missing or inconsistent analysis graphs are never replaced", {
  expect_error(get_analysis_graph(list()), "must be a tKOIList")
  result = toy_run(n_permutation = 2)
  result@subnetwork = NULL
  expect_error(get_analysis_graph(result), "original subnetwork")
  expect_error(plot_network(result, "missing"), "original subnetwork")
  result@subnetwork = igraph::delete_vertices(toy_network(), 1)
  expect_error(get_analysis_graph(result), "do not match")
})
