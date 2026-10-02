shiny_app_environment = local({
  environment = NULL
  function() {
    for (package in .tkoi_app_dependencies()) testthat::skip_if_not_installed(package)
    if (is.null(environment)) {
      app_root = system.file("shiny", "tkoi", package = "tkoi", mustWork = TRUE)
      environment <<- new.env(parent = globalenv())
      environment$app_root = app_root
      withr::with_dir(app_root, {
        sys.source("global.R", envir = environment)
        sys.source("server.R", envir = environment)
      })
    }
    environment
  }
})

shiny_test_graph = local({
  graph = NULL
  universe = NULL
  function() {
    if (is.null(graph)) {
      genes = as.data.frame(tkoi::genes)
      genes = genes[!is.na(genes$ensembl) & genes$degree >= 10 & genes$degree <= 30, ]
      common_degree = as.numeric(names(sort(table(genes$degree), decreasing = TRUE))[1])
      universe <<- head(genes[genes$degree == common_degree & !duplicated(genes$ensembl), ], 40)
      vertices = unique(unlist(lapply(
        igraph::ego(tkoi::tkoi_net, order = 1, nodes = universe$id, mode = "all"),
        names
      )))
      graph <<- igraph::induced_subgraph(tkoi::tkoi_net, vertices)
    }
    list(graph = graph, genes = universe)
  }
})
shiny_test_expression = function() {
  genes = shiny_test_graph()$genes
  n = nrow(genes)
  data.frame(
    gene_name = genes$ensembl,
    logfc = rep(c(-1.5, 0.25, 0.8, 2), length.out = n),
    pvalue = seq(0.001, 0.04, length.out = n)
  )
}
shiny_test_app = function() {
  app_environment = shiny_app_environment()
  instance = new.env(parent = app_environment)
  instance$example_data = shiny_test_expression()
  sys.source(file.path(app_environment$app_root, "server.R"), envir = instance)
  instance$create_server(shiny_test_graph()$graph)
}
shiny_set_controls = function(session, ...) {
  session$setInputs(
    pvalue_threshold = 0.02,
    logfc_threshold = 0.25,
    indirect_link_threshold = 3,
    topology_similarity = 0.9,
    n_permutation = 6,
    damping_factor = 0.85,
    maximum_iteration = 500,
    n_cores = 2,
    ...
  )
}
