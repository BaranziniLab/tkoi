#' tkoi: Transcriptomic Knowledge-Graph Omics Integration
#'
#' Network-aware gene enrichment analysis. Transcriptomic signals are
#' propagated through a human biological knowledge graph with personalized
#' PageRank, and a degree-matched permutation null identifies enriched
#' biological concepts. The main entry point is [run_tkoi()].
#'
#' @section Performance:
#' The PageRank solver, neighborhood counting, and null sampling are
#' implemented in C++ and run on multiple threads. See the `n_cores`
#' argument of [run_tkoi()] and the `tkoi.n_cores` option.
#'
#' @keywords internal
#' @useDynLib tkoi, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

# Column names used inside dplyr verbs.
utils::globalVariables(c(
  "BgRatio", "Count", "Description", "GeneRatio", "ID", "Namespace",
  "beta", "definition", "direction", "ensembl", "fdr", "gene_name", "geneID",
  "id", "identifier", "logfc", "met_threshold", "name", "node_id", "node_type", "p.adjust",
  "p_value", "pagerank", "pvalue", "qvalue", "tkoi_beta"
))
