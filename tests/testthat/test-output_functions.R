###########################
# Author: Wanjun Gu
# Email: wanjun.gu@ucsf.edu
# Date: 2026-04-20
###########################

library(tkoi)

expression_data = data.table::fread(testthat::test_path("ra.csv"))

set.seed(42)
tkoi_result = run_tkoi(
  expression_data  = expression_data,
  subnetwork       = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold  = 0.25,
  topology_similarity   = 0.9,
  indirect_link_threshold = 3,
  n_permutation    = 10,
  damping_factor   = 0.85,
  maximum_iteration = 500
)

test_that("export_gene_exploration_data returns a data frame with expected columns", {
  gene_data = export_gene_exploration_data(tkoi_result)
  expect_true(is.data.frame(gene_data))
  expected_cols = c("gene_name", "gene_symbol", "id", "identifier",
                    "experimental_logfc", "experimental_pvalue",
                    "pagerank", "tkoi_beta", "tkoi_pvalue", "tkoi_fdr")
  expect_equal(names(gene_data), expected_cols)
  expect_true(nrow(gene_data) > 0)
})

test_that("visualize_topn returns a ggplot for BiologicalProcess", {
  plt = visualize_topn(tkoi_result, category = "BiologicalProcess", top_n = 10)
  expect_s3_class(plt, "ggplot")
})

test_that("visualize_topn returns a ggplot for Gene", {
  plt = visualize_topn(tkoi_result, category = "Gene", top_n = 10)
  expect_s3_class(plt, "ggplot")
})

test_that("visualize_topn returns a ggplot for Pathway", {
  plt = visualize_topn(tkoi_result, category = "Pathway", top_n = 10)
  expect_s3_class(plt, "ggplot")
})

test_that("export_network_summary_statistics writes an xlsx file", {
  tmp = tempfile(fileext = ".xlsx")
  export_network_summary_statistics(tkoi_result, filename = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
  file.remove(tmp)
})
