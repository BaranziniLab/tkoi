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

test_that("run_tkoi returns a tKOIList object", {
  expect_s4_class(tkoi_result, "tKOIList")
})

test_that("tKOIList slots are populated", {
  expect_true(nrow(tkoi_result@expression_data) > 0)
  expect_true(nrow(tkoi_result@pagerank_data) > 0)
  expect_true(length(tkoi_result@network_summary_statistics) > 0)
})

test_that("pagerank_data has required columns", {
  cols = names(tkoi_result@pagerank_data)
  expect_true("node_id" %in% cols)
  expect_true("pagerank" %in% cols)
  expect_true(any(grepl("^perm\\.", cols)))
})

test_that("network_summary_statistics contains key node types", {
  node_types = names(tkoi_result@network_summary_statistics)
  expect_true("Gene" %in% node_types)
  expect_true("BiologicalProcess" %in% node_types)
  expect_true("Pathway" %in% node_types)
  expect_true("Disease" %in% node_types)
})

test_that("network enrichment statistics have expected columns", {
  gene_stats = tkoi_result@network_summary_statistics$Gene
  expect_true(all(c("node_id", "pagerank", "beta", "p_value", "fdr") %in% names(gene_stats)))
})

test_that("stored parameters match inputs", {
  expect_equal(tkoi_result@n_permutation, 10)
  expect_equal(tkoi_result@pvalue_threshold, 0.05)
  expect_equal(tkoi_result@logfc_threshold, 0.25)
  expect_equal(tkoi_result@damping_factor, 0.85)
})
