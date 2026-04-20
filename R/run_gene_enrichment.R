#' Run Gene Enrichment and Compare with TKOI Data
#'
#' This function performs gene enrichment analysis using the `clusterProfiler` package and integrates the results
#' with the `tKOIList` object, specifically comparing the enrichment results with the TKOI data. It generates
#' scatter plots visualizing the relationship between TKOI scores and enrichment results.
#'
#' @param tkoi_list An object of class `tKOIList` that contains the input data for analysis. The object should
#' include the following slots:
#'   - `expression_data`: A `data.frame` with columns including `gene_name` and `pvalue`.
#'   - `network_summary_statistics`: A `list` containing network-level summary data for Biological Process (BP),
#'     Cellular Component (CC), and Molecular Function (MF), each stored as a `data.frame`.
#'   - `pvalue_threshold`: A `numeric` value specifying the significance threshold for filtering genes.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts genes with p-values below the threshold specified in the `tkoi_list`.
#' 2. Conducts Gene Ontology (GO) enrichment analysis for Biological Process (BP), Cellular Component (CC),
#'    and Molecular Function (MF) using `clusterProfiler::enrichGO`.
#' 3. Computes pairwise term similarities using `enrichplot::pairwise_termsim`.
#' 4. Merges the enrichment results with TKOI network statistics to create a unified dataset.
#' 5. Generates scatter plots to visualize the relationship between TKOI network enrichment effect size
#' and gene enrichment q-values.
#'
#' @return A `tKOIList` object with the `gene_enrichment_comparison` slot populated. This slot is a `list`
#' containing:
#'   - `enrichment_result`: A `data.frame` with the merged results of the enrichment analysis and TKOI network data.
#'   - `comparison_scatter1`: A `ggplot` object visualizing a scatter plot without facets.
#'   - `comparison_scatter2`: A `ggplot` object visualizing a scatter plot with facets for each namespace.
#'
#' @examples
#' \dontrun{
#' # Example usage of the function
#' library(clusterProfiler)
#' library(enrichplot)
#' library(org.Hs.eg.db)
#' library(ggplot2)
#' library(dplyr)
#'
#' # Create a dummy tKOIList object
#' tkoi_list = new("tKOIList",
#'                 expression_data = data.frame(
#'                   gene_name = c("gene1", "gene2", "gene3"),
#'                   pvalue = c(0.01, 0.02, 0.2)
#'                 ),
#'                 network_summary_statistics = list(
#'                   BiologicalProcess = data.frame(
#'                     identifier = c("GO:0008150", "GO:0009987"),
#'                     node_id = c(1, 2),
#'                     pagerank = c(0.2, 0.3),
#'                     beta = c(0.5, 0.4),
#'                     p_value = c(0.01, 0.02),
#'                     fdr = c(0.05, 0.1),
#'                     definition = c("process1", "process2")
#'                   ),
#'                   CellularComponent = data.frame(),
#'                   MolecularFunction = data.frame()
#'                 ))
#'
#' # Run the enrichment function
#' result = run_gene_enrichment(tkoi_list)
#' }
#' @export
run_gene_enrichment = function(tkoi_list){

  enrichment_gene_data = dplyr::filter(tkoi_list@expression_data,
                                       pvalue <= tkoi_list@pvalue_threshold)$gene_name

  message("Running GO enrichment on: Biological Process\n")
  enrichment_bp = clusterProfiler::enrichGO(gene = enrichment_gene_data,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            keyType = "ENSEMBL",
                                            ont = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05) |>
    enrichplot::pairwise_termsim()
  enrichment_bp_result = enrichment_bp@result |>
    dplyr::mutate(Namespace = "Biological Process")

  message("Running GO enrichment on: Cellular Component\n")
  enrichment_cc = clusterProfiler::enrichGO(gene = enrichment_gene_data,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            keyType = "ENSEMBL",
                                            ont = "CC",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05) |>
    enrichplot::pairwise_termsim()
  enrichment_cc_result = enrichment_cc@result |>
    dplyr::mutate(Namespace = "Cellular Component")

  message("Running GO enrichment on: Molecular Function\n")
  enrichment_mf = clusterProfiler::enrichGO(gene = enrichment_gene_data,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            keyType = "ENSEMBL",
                                            ont = "MF",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05) |>
    enrichplot::pairwise_termsim()
  enrichment_mf_result = enrichment_mf@result |>
    dplyr::mutate(Namespace = "Molecular Function")

  enrichment_result = rbind.data.frame(enrichment_bp_result,
                                       enrichment_cc_result,
                                       enrichment_mf_result)

  tkoi_go_result = rbind.data.frame(tkoi_list@network_summary_statistics$BiologicalProcess,
                                    tkoi_list@network_summary_statistics$CellularComponent,
                                    tkoi_list@network_summary_statistics$MolecularFunction)

  enrichment_result = enrichment_result |>
    dplyr::left_join(tkoi_go_result,
                     by = dplyr::join_by("ID" == "identifier")) |>
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue,
                  p.adjust, qvalue, geneID, Count, Namespace,
                  node_id, pagerank, beta, p_value, fdr, definition) |>
    dplyr::rename("tkoi_node_id" = "node_id") |>
    dplyr::rename("tkoi_pagerank" = "pagerank") |>
    dplyr::rename("tkoi_beta" = "beta") |>
    dplyr::rename("tkoi_p_value" = "p_value") |>
    dplyr::rename("tkoi_fdr" = "fdr") |>
    dplyr::filter(!is.na(Namespace))

  plt1 = ggplot2::ggplot(enrichment_result, ggplot2::aes(x = tkoi_beta,
                                                         y = -log(qvalue, base = 10))) +
    ggplot2::geom_point(ggplot2::aes(fill = Namespace), color = "black", size = 3, shape = 21) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log(0.05, base = 10)), linetype = "dashed") +
    ggplot2::geom_smooth(method = "lm", color = "black") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "") +
    ggplot2::xlab("tKOI Network Enrichment") +
    ggplot2::ylab("Gene Enrichment -log10(qvalue)") +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::theme_minimal(base_size = 14)

  plt2 = ggplot2::ggplot(enrichment_result, ggplot2::aes(x = tkoi_beta,
                                                         y = -log(qvalue, base = 10))) +
    ggplot2::geom_point(ggplot2::aes(fill = Namespace), color = "black", size = 3, shape = 21) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log(0.05, base = 10)), linetype = "dashed") +
    ggplot2::geom_smooth(method = "lm", color = "black") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "") +
    ggplot2::xlab("tKOI Network Enrichment") +
    ggplot2::ylab("Gene Enrichment -log10(qvalue)") +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::facet_grid(.~Namespace) +
    ggplot2::theme_minimal(base_size = 14)

  gene_enrichment_comparison = list(
    enrichment_result = enrichment_result,
    comparison_scatter1 = plt1,
    comparison_scatter2 = plt2
  )

  tkoi_list@gene_enrichment_comparison = gene_enrichment_comparison
  return(tkoi_list)

}
