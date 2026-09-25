# Slow tests on the full tkoi_net and on clusterProfiler. They load the
# 0.5 GB network and take a few minutes, so they only run when
# TKOI_EXTENDED_TESTS=true, and never on CRAN (so outside devtools::test(),
# also set NOT_CRAN=true).

skip_unless_extended = function() {
  skip_if_not(identical(Sys.getenv("TKOI_EXTENDED_TESTS"), "true"), "set TKOI_EXTENDED_TESTS=true to run")
  skip_on_cran()
}

extended_node_types = sort(
  c(
    "Anatomy", "BiologicalProcess", "CellType", "CellularComponent", "ClinicalLab", "Complex",
    "Compound", "Disease", "EC", "Gene", "MiRNA", "MolecularFunction", "Pathway", "Protein",
    "ProteinDomain", "ProteinFamily", "PwGroup", "Reaction"
  ),
  method = "radix"
)

extended_statistic_columns = c(
  "node_id", "node_type", "identifier", "pagerank", "beta", "p_value", "fdr",
  "direct_links", "indirect_links", "valency"
)

extended_seed = 20240917

# The first 90 genes of the real differential expression table.
extended_expression = function() {
  utils::read.csv(test_path("ra.csv"), nrows = 90, stringsAsFactors = FALSE)
}

# run_tkoi() and the released algorithm on the full network, each computed
# once and shared by the tests below.
extended_cache = new.env(parent = emptyenv())

extended_run = function() {
  if (is.null(extended_cache$run)) {
    extended_cache$run = withr::with_seed(extended_seed, run_tkoi(
      expression_data = extended_expression(),
      n_permutation = 3,
      n_cores = 2,
      verbose = FALSE
    ))
  }
  extended_cache$run
}

extended_legacy_run = function() {
  if (is.null(extended_cache$legacy)) {
    extended_cache$legacy = withr::with_seed(extended_seed, legacy_run_tkoi(
      expression_data = extended_expression(),
      subnetwork = tkoi::tkoi_net,
      n_permutation = 3
    ))
  }
  extended_cache$legacy
}

test_that("tkoi_net is a plain undirected igraph of the documented size", {
  skip_unless_extended()
  network = tkoi::tkoi_net

  expect_false(is.raw(network))
  expect_true(igraph::is_igraph(network))
  expect_false(igraph::is_directed(network))
  expect_equal(igraph::vcount(network), 939059)
  expect_equal(igraph::ecount(network), 10622200)
  expect_true(all(c("name", "identifier", "source", "labels", "degree") %in% igraph::vertex_attr_names(network)))
  expect_true("edge_type" %in% igraph::edge_attr_names(network))
})

test_that("tkoi_net vertex attributes are consistent", {
  skip_unless_extended()
  network = tkoi::tkoi_net
  vertex_names = igraph::V(network)$name

  expect_type(vertex_names, "character")
  expect_identical(anyDuplicated(vertex_names), 0L)
  expect_type(igraph::V(network)$identifier, "character")
  expect_type(igraph::V(network)$source, "character")

  labels = igraph::V(network)$labels
  expect_true(all(grepl("^\\['[A-Za-z]+'\\]$", labels)))
  expect_identical(sort(unique(gsub("[][']", "", labels)), method = "radix"), extended_node_types)

  expect_equal(as.numeric(igraph::V(network)$degree), as.numeric(igraph::degree(network)))
})

test_that("every tkoi::genes id is a Gene vertex of tkoi_net with the listed degree", {
  skip_unless_extended()
  network = tkoi::tkoi_net
  genes = as.data.frame(tkoi::genes)

  rows = match(genes$id, igraph::V(network)$name)
  expect_false(anyNA(rows))
  expect_true(all(igraph::V(network)$labels[rows] == "['Gene']"))
  expect_equal(as.numeric(igraph::V(network)$degree[rows]), as.numeric(genes$degree))
  expect_identical(igraph::V(network)$identifier[rows], as.character(genes$identifier))
})

test_that("run_tkoi works on the full network and reports all 18 node types", {
  skip_unless_extended()
  result = extended_run()
  network = tkoi::tkoi_net

  expect_s4_class(result, "tKOIList")
  expect_equal(result@n_permutation, 3)

  pagerank_data = result@pagerank_data
  expect_identical(names(pagerank_data), c("node_id", "pagerank", "perm.1", "perm.2", "perm.3"))
  expect_identical(pagerank_data$node_id, igraph::V(network)$name)
  expect_identical(rownames(pagerank_data), pagerank_data$node_id)
  for (column in c("pagerank", "perm.1", "perm.2", "perm.3")) {
    expect_true(all(pagerank_data[[column]] >= 0), info = column)
    expect_equal(sum(pagerank_data[[column]]), 1, tolerance = 1e-8)
  }

  statistics = result@network_summary_statistics
  expect_identical(names(statistics), extended_node_types)
  for (type in names(statistics)) {
    table = statistics[[type]]
    expect_s3_class(table, "tbl_df")
    expect_identical(names(table)[seq_along(extended_statistic_columns)], extended_statistic_columns)
    expect_gt(nrow(table), 0)
    expect_true(all(table$node_type == type), info = type)
    fdr = table$fdr[!is.nan(table$fdr)]
    expect_true(all(fdr >= 0 & fdr <= 1), info = type)
  }
  expect_true(all(statistics$Compound$identifier %in% tkoi::human_metabolites))
})

test_that("run_tkoi seeds on the full network are the mapped, filtered genes", {
  skip_unless_extended()
  result = extended_run()
  expression = extended_expression()

  selected = expression[expression$pvalue <= 0.05 & abs(expression$logfc) >= 0.25, ]
  seeds = unique(tkoi::genes$id[tkoi::genes$ensembl %in% selected$gene_name])
  expect_gt(length(seeds), 0)

  genes = result@network_summary_statistics$Gene
  seed_rows = genes[genes$node_id %in% seeds, ]
  expect_identical(sort(unique(seed_rows$node_id)), sort(seeds))
  # A seed is within one hop of itself.
  expect_true(all(seed_rows$direct_links >= 1))
  expect_true(all(seed_rows$indirect_links >= seed_rows$direct_links))
})

test_that("run_tkoi matches the released algorithm on the full network", {
  skip_unless_extended()
  # Same check as the toy-network parity tests (helper-parity.R): identical
  # PageRank, beta and p_value up to numerical precision, with the deliberate
  # 1.1.0 changes (unannotated nodes kept, repaired Compound family, constant
  # nulls untestable).
  # On tkoi_net, igraph's solver (used by 1.0.0) itself has per-node relative
  # errors up to ~3e-7 on nodes with tiny PageRank (the new solver: ~8e-8), so
  # the relative PageRank bound is 1e-6 here instead of 1e-8.
  expect_matches_legacy(extended_run(), extended_legacy_run(), pagerank_relative = 1e-6)
})

test_that("run_gene_enrichment adds the GO comparison to a toy result", {
  skip_unless_extended()
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  result = toy_run()
  enriched = suppressMessages(run_gene_enrichment(result))

  expect_s4_class(enriched, "tKOIList")
  comparison = enriched@gene_enrichment_comparison
  expect_named(comparison, c("enrichment_result", "comparison_scatter1", "comparison_scatter2"))

  enrichment = comparison$enrichment_result
  expect_s3_class(enrichment, "data.frame")
  expect_identical(
    names(enrichment),
    c(
      "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count",
      "Namespace", "tkoi_node_id", "tkoi_pagerank", "tkoi_beta", "tkoi_p_value", "tkoi_fdr", "definition"
    )
  )
  expect_gt(nrow(enrichment), 0)
  expect_true(all(startsWith(enrichment$ID, "GO:")))
  expect_true(all(enrichment$Namespace %in% c("Biological Process", "Cellular Component", "Molecular Function")))

  # GO terms that are nodes of the toy network carry that node's statistics.
  go_nodes = do.call(rbind, lapply(
    c("BiologicalProcess", "CellularComponent", "MolecularFunction"),
    function(type) {
      table = result@network_summary_statistics[[type]]
      if (is.null(table)) NULL else as.data.frame(table)[, c("identifier", "node_id", "beta")]
    }
  ))
  matched = enrichment[!is.na(enrichment$tkoi_node_id), ]
  rows = match(matched$ID, go_nodes$identifier)
  expect_false(anyNA(rows))
  expect_identical(matched$tkoi_node_id, go_nodes$node_id[rows])
  expect_identical(matched$tkoi_beta, go_nodes$beta[rows])

  for (plot in comparison[c("comparison_scatter1", "comparison_scatter2")]) {
    expect_true(inherits(plot, "ggplot"))
    expect_no_error(suppressWarnings(ggplot2::ggplot_build(plot)))
  }

  # The other slots are left untouched.
  expect_identical(enriched@network_summary_statistics, result@network_summary_statistics)
  expect_identical(enriched@pagerank_data, result@pagerank_data)
})
