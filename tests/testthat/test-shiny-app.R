test_that("the app loads the installed tkoi namespace and class", {
  app_environment = shiny_app_environment()
  testthat::expect_identical(get("run_tkoi", envir = app_environment, inherits = TRUE), tkoi::run_tkoi)
  testthat::expect_identical(methods::getClass("tKOIList")@package, "tkoi")
  testthat::expect_identical(app_environment$tkoi_version, as.character(utils::packageVersion("tkoi")))
  testthat::expect_gte(utils::packageVersion("tkoi"), package_version("1.1.0"))
})
test_that("supported upload formats preserve the same expression input", {
  app_environment = shiny_app_environment()
  expected = shiny_test_expression()
  csv = withr::local_tempfile(fileext = ".csv")
  tsv = withr::local_tempfile(fileext = ".tsv")
  xlsx = withr::local_tempfile(fileext = ".xlsx")
  write.csv(expected, csv, row.names = FALSE)
  write.table(expected, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  openxlsx::write.xlsx(expected, xlsx)
  for (file in c(csv, tsv, xlsx)) {
    actual = app_environment$read_expression_data(list(name = basename(file), datapath = file))
    testthat::expect_equal(actual, expected)
  }
  testthat::expect_equal(app_environment$read_expression_data(NULL, expected), expected)
  invalid = withr::local_tempfile(fileext = ".csv")
  write.csv(data.frame(gene_name = "ENSG0001"), invalid, row.names = FALSE)
  testthat::expect_error(
    app_environment$read_expression_data(list(name = "invalid.csv", datapath = invalid)),
    "Missing expression columns"
  )
})
test_that("example analysis matches the actual tkoi engine for every result table", {
  # Observe the RNG state at the library boundary. Shiny generates random
  # progress IDs before analysis; the reference must use the engine's seed.
  shiny_app_environment()
  real_run = tkoi::run_tkoi
  boundary = new.env()
  testthat::local_mocked_bindings(
    run_tkoi = function(...) {
      boundary$rng = get(".Random.seed", envir = globalenv())
      real_run(...)
    },
    .package = "tkoi"
  )
  withr::local_preserve_seed()
  set.seed(101)
  shiny::testServer(shiny_test_app(), {
    testthat::expect_null(tkoi_output())
    shiny_set_controls(session, run_tkoi_analysis = 1)
    actual = tkoi_output()
    testthat::expect_s4_class(actual, "tKOIList")
    testthat::expect_equal(actual@expression_data, shiny_test_expression())
    testthat::expect_equal(actual@n_permutation, 6)
    testthat::expect_named(actual@pagerank_data, c("node_id", "pagerank", "null_mean", "null_sd"))
    assign(".Random.seed", boundary$rng, envir = globalenv()) # nolint: object_name_linter.
    reference = real_run(
      expression_data = shiny_test_expression(), subnetwork = shiny_test_graph()$graph,
      pvalue_threshold = 0.02, logfc_threshold = 0.25, indirect_link_threshold = 3,
      topology_similarity = 0.9, n_permutation = 6, damping_factor = 0.85,
      maximum_iteration = 500, n_cores = 2, keep_permutations = FALSE, verbose = FALSE
    )
    testthat::expect_equal(actual@pagerank_data, reference@pagerank_data, tolerance = 0)
    testthat::expect_equal(actual@network_summary_statistics, reference@network_summary_statistics, tolerance = 0)
    session$setInputs(result_category = "BiologicalProcess", lookup_category = "Gene")
    testthat::expect_match(output$tkoi_result_ui$html, "result_category")
    testthat::expect_match(output$tkoi_network_visualization_lookup$html, "lookup_category")
    testthat::expect_true(is.character(output$tkoi_result_table))
    testthat::expect_true(is.character(output$tkoi_result_table_lookup))
    excel = withr::local_tempfile(fileext = ".xlsx")
    export_excel(excel)
    testthat::expect_setequal(openxlsx::getSheetNames(excel), names(actual@network_summary_statistics))
    first = names(actual@network_summary_statistics)[1]
    exported = openxlsx::read.xlsx(excel, sheet = first)
    testthat::expect_equal(nrow(exported), nrow(actual@network_summary_statistics[[first]]))
    significant = withr::local_tempfile(fileext = ".xlsx")
    export_excel(significant, significant_only = TRUE)
    for (name in openxlsx::getSheetNames(significant)) {
      data = openxlsx::read.xlsx(significant, sheet = name)
      if (nrow(data)) testthat::expect_true(all(data$fdr <= 0.05))
    }
  })
})
test_that("the current upload is used immediately, including after another upload", {
  shiny_app_environment()
  withr::local_preserve_seed()
  set.seed(102)
  first = shiny_test_expression()[1:15, ]
  second = shiny_test_expression()[21:35, ]
  rownames(second) = NULL
  one = withr::local_tempfile(fileext = ".csv")
  two = withr::local_tempfile(fileext = ".tsv")
  write.csv(first, one, row.names = FALSE)
  write.table(second, two, sep = "\t", row.names = FALSE, quote = FALSE)
  shiny::testServer(shiny_test_app(), {
    shiny_set_controls(session, upload_data = data.frame(name = "first.csv", datapath = one), run_tkoi_analysis = 1)
    testthat::expect_equal(tkoi_output()@expression_data, first)
    session$setInputs(pvalue_threshold = 0.05,
      upload_data = data.frame(name = "second.tsv", datapath = two), run_tkoi_analysis = 2)
    testthat::expect_equal(tkoi_output()@expression_data, second)
    testthat::expect_equal(processed_data()$gene_name, second$gene_name)
    testthat::expect_equal(tkoi_output()@pvalue_threshold, 0.05)
    session$setInputs(upload_data = NULL, n_permutation = 9, run_tkoi_analysis = 3)
    testthat::expect_equal(tkoi_output()@expression_data, shiny_test_expression())
    testthat::expect_equal(tkoi_output()@n_permutation, 9)
  })
})
test_that("invalid input clears obsolete results and a corrected run recovers", {
  shiny_app_environment()
  withr::local_preserve_seed()
  set.seed(103)
  shiny::testServer(shiny_test_app(), {
    shiny_set_controls(session, run_tkoi_analysis = 1)
    testthat::expect_s4_class(tkoi_output(), "tKOIList")
    session$setInputs(n_permutation = 1, run_tkoi_analysis = 2)
    testthat::expect_null(tkoi_output())
    testthat::expect_null(analysis_parameters())
    session$setInputs(n_permutation = 5, run_tkoi_analysis = 3)
    testthat::expect_s4_class(tkoi_output(), "tKOIList")
    testthat::expect_equal(tkoi_output()@n_permutation, 5)
  })
})

test_that("startup rejects a stale tkoi namespace in an existing R session", {
  app_environment = shiny_app_environment()
  stale_session = new.env(parent = globalenv())
  stale_session$getNamespaceVersion = function(package) package_version("1.0.0")
  testthat::expect_error(
    withr::with_dir(app_environment$app_root, sys.source("global.R", envir = stale_session)),
    "Restart R before launching the tKOI app"
  )
  testthat::expect_identical(
    app_environment$tkoi_version, as.character(getNamespaceVersion("tkoi"))
  )
})
