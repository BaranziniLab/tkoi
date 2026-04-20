#' Export Gene Exploration Data
#'
#' This function extracts and exports a combined dataset of experimental gene expression data
#' and tKOI network statistics, including key metrics such as log fold changes, p-values,
#' PageRank scores, and adjusted false discovery rates (FDR).
#'
#' @param tkoi_list An object of class `tKOIList`. This object must contain the following slots:
#'   - `expression_data`: A `data.frame` containing experimental data with columns `gene_name`, `logfc`, and `pvalue`.
#'   - `network_summary_statistics`: A list with a `data.frame` for gene-level network statistics, including `id`, `pagerank`, `beta`, `p_value`, and `fdr`.
#'   - `pvalue_threshold`: A numeric value specifying the p-value threshold.
#'   - `logfc_threshold`: A numeric value specifying the log fold change threshold.
#'
#' @details
#' The function performs the following steps:
#' 1. Merges the `expression_data` from the `tKOIList` object with gene metadata using `inner_join`.
#' 2. Merges the result with gene-level network data using `right_join`.
#' 3. Selects key columns from the merged dataset, including:
#'    - `gene_name`: Gene names from the experimental data.
#'    - `name`: Gene symbols or common names.
#'    - `id` and `identifier`: Unique identifiers for genes in the network.
#'    - `logfc`: Experimental log fold change values.
#'    - `pvalue`: Experimental p-values.
#'    - `pagerank`: tKOI PageRank scores.
#'    - `beta`: tKOI z-scores for network enrichment.
#'    - `p_value`: tKOI unadjusted p-values.
#'    - `fdr`: tKOI false discovery rates (FDR).
#' 4. Renames columns for clarity and standardization.
#'
#' @return A `data.frame` with the following columns:
#'   - `gene_name`: Gene names from the experimental dataset.
#'   - `gene_symbol`: Gene symbols or common names from the metadata.
#'   - `id`: Unique node identifiers in the network.
#'   - `identifier`: External identifiers for the genes (e.g., Ensembl IDs).
#'   - `experimental_logfc`: Log fold change values from the experimental data.
#'   - `experimental_pvalue`: P-values from the experimental data.
#'   - `pagerank`: Personalized PageRank scores from the tKOI network analysis.
#'   - `tkoi_beta`: tKOI z-scores for network enrichment.
#'   - `tkoi_pvalue`: Unadjusted p-values from the tKOI network analysis.
#'   - `tkoi_fdr`: False discovery rates (FDR) from the tKOI network analysis.
#'
#' @examples
#' \dontrun{
#' # Create a dummy tKOIList object
#' tkoi_list = new("tKOIList",
#'                 expression_data = data.frame(
#'                   gene_name = c("gene1", "gene2", "gene3"),
#'                   logfc = c(1.2, -0.8, 0.5),
#'                   pvalue = c(0.01, 0.03, 0.2)
#'                 ),
#'                 network_summary_statistics = list(
#'                   Gene = data.frame(
#'                     id = c(1, 2, 3),
#'                     pagerank = c(0.05, 0.03, 0.02),
#'                     beta = c(2.1, -1.5, 0.9),
#'                     p_value = c(0.01, 0.02, 0.05),
#'                     fdr = c(0.02, 0.03, 0.1)
#'                   )
#'                 ),
#'                 pvalue_threshold = 0.05,
#'                 logfc_threshold = 0.5)
#'
#' # Export the gene exploration data
#' gene_data = export_gene_exploration_data(tkoi_list)
#' print(gene_data)
#' }
#'
#' @seealso \code{\link[dplyr]{dplyr}}, \code{\link[tkoi]{make_gene_exploration_plot}}
#'
#' @export
export_gene_exploration_data = function(tkoi_list){
  expression_data = tkoi_list@expression_data
  network_data = tkoi_list@network_summary_statistics$Gene
  p_value_treshold = tkoi_list@pvalue_threshold
  logfc_threshold = tkoi_list@logfc_threshold

  gene_data = dplyr::inner_join(expression_data, tkoi::genes,
                                by = dplyr::join_by(gene_name == ensembl)) |>
    dplyr::select(-identifier, -name) |>
    dplyr::right_join(network_data, by = dplyr::join_by("id" == "node_id")) |>
    dplyr::select(gene_name, name, id, identifier, logfc, pvalue, pagerank, beta, p_value, fdr)

  names(gene_data) = c("gene_name", "gene_symbol", "id", "identifier",
                       "experimental_logfc", "experimental_pvalue",
                       "pagerank", "tkoi_beta", "tkoi_pvalue", "tkoi_fdr")

  return(gene_data)
}
