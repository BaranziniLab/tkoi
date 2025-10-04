setClass(
  Class = "tKOIList",
  slots = list(
    expression_data = "data.frame",
    pagerank_data = "data.frame",
    network_summary_statistics = "list",
    subnetwork = "ANY",
    pvalue_threshold = "numeric",
    logfc_threshold = "numeric",
    topology_similarity = "numeric",
    n_permutation = "numeric",
    damping_factor = "numeric",
    maximum_iteration = "numeric",
    gene_enrichment_comparison = "ANY"
  ),
  prototype = list(
    subnetwork = NULL,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.25,
    topology_similarity = 0.9,
    n_permutation = 100,
    damping_factor = 0.85,
    maximum_iteration = 500,
    gene_enrichment_comparison = list()
  )
)

setMethod(
  f = "show",
  signature = "tKOIList",
  definition = function(object){
    parameter_set = data.frame(parameter = c(object@pvalue_threshold,
                                             object@logfc_threshold,
                                             object@topology_similarity,
                                             object@n_permutation,
                                             object@damping_factor,
                                             object@maximum_iteration))
    rownames(parameter_set) = c("P.Value Threshold",
                                "Log Fold Change Threshold",
                                "Topology Similarity",
                                "N Permutation",
                                "Damping Factor",
                                "Maximum Iteration")

    cat("[tKOIList object]\n")
    cat(glue::glue("{dim(object@expression_data)[1]} genes included in the experiment."))
    cat("\nTotal object size: ")
    print(utils::object.size(object), units = "auto")
    cat("Network enrichment analysis run with these parameters:\n")
    print(knitr::kable(parameter_set))
  }
)
