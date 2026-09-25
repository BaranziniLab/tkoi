#' Run Gene Enrichment and Compare with TKOI Data
#'
#' This function performs gene enrichment analysis using the `clusterProfiler` package and
#' integrates the results
#' with the `tKOIList` object, specifically comparing the enrichment results with the TKOI data.
#' It generates
#' scatter plots visualizing the relationship between TKOI scores and enrichment results.
#'
#' @param tkoi_list An object of class `tKOIList` that contains the input data for analysis. The
#' object should
#' include the following slots:
#'   - `expression_data`: A `data.frame` with columns including `gene_name` and `pvalue`.
#'   - `network_summary_statistics`: A `list` containing network-level summary data for
#'     Biological Process (BP),
#'     Cellular Component (CC), and Molecular Function (MF), each stored as a `data.frame`.
#'   - `pvalue_threshold`: A `numeric` value specifying the significance threshold for filtering
#'     genes.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts genes with p-values at or below the threshold stored in `tkoi_list`
#'    (the log fold change threshold is not applied here).
#' 2. Conducts Gene Ontology (GO) enrichment analysis for Biological Process (BP), Cellular
#'   Component (CC),
#'    and Molecular Function (MF) using `clusterProfiler::enrichGO`.
#' 3. Merges the enrichment results with TKOI network statistics to create a unified dataset.
#' 4. Generates scatter plots to visualize the relationship between TKOI network enrichment
#'   effect size
#' and gene enrichment q-values.
#'
#' @return A `tKOIList` object with the `gene_enrichment_comparison` slot populated. This slot
#' is a `list`
#' containing:
#'   - `enrichment_result`: A `data.frame` with the merged results of the enrichment analysis
#'     and TKOI network data.
#'   - `comparison_scatter1`: A `ggplot` object visualizing a scatter plot without facets.
#'   - `comparison_scatter2`: A `ggplot` object visualizing a scatter plot with facets for each
#'     namespace.
#'
#' @examples
#' \dontrun{
#' tkoi_result = run_tkoi(expression_data)
#' tkoi_result = run_gene_enrichment(tkoi_result)
#'
#' tkoi_result@gene_enrichment_comparison$enrichment_result
#' tkoi_result@gene_enrichment_comparison$comparison_scatter1
#' }
#' @export
run_gene_enrichment = function(tkoi_list) {
  if (!methods::is(tkoi_list, "tKOIList")) {
    stop("`tkoi_list` must be a tKOIList object from run_tkoi().", call. = FALSE)
  }
  # Same gene rows as run_tkoi(): no blank names, first row per gene.
  expression_data = .tkoi_clean_expression(tkoi_list@expression_data)
  enrichment_gene_data = expression_data$gene_name[
    !is.na(expression_data$pvalue) & expression_data$pvalue <= tkoi_list@pvalue_threshold
  ]
  if (length(enrichment_gene_data) == 0) {
    stop("No genes pass `pvalue_threshold`; nothing to enrich.", call. = FALSE)
  }

  namespaces = c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")
  enrichment_result = lapply(names(namespaces), function(ontology) {
    message(glue::glue("Running GO enrichment on: {namespaces[[ontology]]}\n"))
    enrichment = clusterProfiler::enrichGO(
      gene = enrichment_gene_data,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENSEMBL",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
    if (is.null(enrichment)) {
      return(NULL)
    }
    result = enrichment@result
    result$Namespace = rep(namespaces[[ontology]], nrow(result))
    result
  })
  enrichment_result = do.call(rbind.data.frame, enrichment_result)
  if (is.null(enrichment_result) || nrow(enrichment_result) == 0) {
    stop("GO enrichment returned no terms for these genes.", call. = FALSE)
  }

  go_types = c("BiologicalProcess", "CellularComponent", "MolecularFunction")
  tkoi_go_result = do.call(rbind.data.frame, lapply(go_types, function(type) {
    x = tkoi_list@network_summary_statistics[[type]]
    if (is.null(x)) {
      return(NULL)
    }
    x[, c("identifier", "node_id", "pagerank", "beta", "p_value", "fdr", "definition")]
  }))

  enrichment_result = enrichment_result |>
    dplyr::left_join(tkoi_go_result, by = dplyr::join_by("ID" == "identifier")) |>
    dplyr::select(
      ID, Description, GeneRatio, BgRatio, pvalue,
      p.adjust, qvalue, geneID, Count, Namespace,
      node_id, pagerank, beta, p_value, fdr, definition
    ) |>
    dplyr::rename(
      tkoi_node_id = node_id,
      tkoi_pagerank = pagerank,
      tkoi_beta = beta,
      tkoi_p_value = p_value,
      tkoi_fdr = fdr
    ) |>
    dplyr::filter(!is.na(Namespace))

  comparison_plot = function(facet) {
    plt = ggplot2::ggplot(
      enrichment_result,
      ggplot2::aes(x = tkoi_beta, y = -log(qvalue, base = 10))
    ) +
      ggplot2::geom_point(ggplot2::aes(fill = Namespace), color = "black", size = 3, shape = 21) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = -log(0.05, base = 10)), linetype = "dashed") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "black") +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "") +
      ggplot2::xlab("tKOI Network Enrichment") +
      ggplot2::ylab("Gene Enrichment -log10(qvalue)") +
      ggplot2::ylim(c(0, NA)) +
      ggplot2::theme_minimal(base_size = 14)
    if (facet) {
      plt = plt + ggplot2::facet_grid(. ~ Namespace)
    }
    plt
  }

  tkoi_list@gene_enrichment_comparison = list(
    enrichment_result = enrichment_result,
    comparison_scatter1 = comparison_plot(facet = FALSE),
    comparison_scatter2 = comparison_plot(facet = TRUE)
  )
  tkoi_list
}
