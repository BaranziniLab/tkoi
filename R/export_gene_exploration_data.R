#' Export Gene Exploration Data
#'
#' Combines the differential expression table of a [run_tkoi()] result with
#' the tKOI network statistics of every gene node, so that experimental and
#' network evidence can be compared gene by gene.
#'
#' @param tkoi_list A `tKOIList` object returned by [run_tkoi()]. Two slots
#'   are used:
#'   - `expression_data`: the input table, with columns `gene_name` (Ensembl
#'     gene IDs), `logfc`, and `pvalue`.
#'   - `network_summary_statistics`: its `Gene` table, with columns
#'     `node_id`, `identifier`, `name`, `pagerank`, `beta`, `p_value`, and
#'     `fdr`.
#'
#' @details
#' `expression_data` is cleaned the same way as in [run_tkoi()]: rows with a
#' missing or blank `gene_name` are dropped, and only the first row of each
#' gene is used. Genes are matched to network nodes through their Ensembl IDs
#' (see [genes]). Every row of the `Gene` table is kept, so network genes
#' that are not in `expression_data` have `NA` in `gene_name`,
#' `experimental_logfc`, and `experimental_pvalue`.
#'
#' @return A data frame with one row per row of the `Gene` table and columns:
#'   - `gene_name`: Ensembl gene ID from `expression_data`.
#'   - `gene_symbol`: Gene symbol.
#'   - `id`: Node ID of the gene in the network (`node_id` in the `Gene`
#'     table).
#'   - `identifier`: NCBI Entrez Gene ID, as a character string.
#'   - `experimental_logfc`: `logfc` from `expression_data`.
#'   - `experimental_pvalue`: `pvalue` from `expression_data`.
#'   - `pagerank`: Observed personalized PageRank.
#'   - `tkoi_beta`: tKOI network enrichment z-score (`beta`).
#'   - `tkoi_pvalue`: Unadjusted one-sided tKOI p-value.
#'   - `tkoi_fdr`: tKOI false discovery rate (Benjamini-Hochberg, among gene
#'     nodes).
#'
#' @examples
#' \dontrun{
#' expression_data = data.table::fread(
#'   system.file("extdata", "example_data.csv", package = "tkoi")
#' )
#'
#' set.seed(1)
#' tkoi_result = run_tkoi(expression_data = expression_data)
#'
#' gene_data = export_gene_exploration_data(tkoi_result)
#' head(gene_data)
#'
#' # Genes supported by both the experiment and the network
#' subset(gene_data, experimental_pvalue <= 0.05 & tkoi_fdr <= 0.05)
#' }
#'
#' @seealso [make_gene_exploration_plot()] to plot the same data,
#'   [export_network_summary_statistics()] to export every node type.
#'
#' @export
export_gene_exploration_data = function(tkoi_list) {
  # Same gene rows as run_tkoi(): no blank names, first row per gene.
  expression_data = .tkoi_clean_expression(tkoi_list@expression_data)
  network_data = tkoi_list@network_summary_statistics$Gene

  gene_data = dplyr::inner_join(expression_data, tkoi::genes,
    by = dplyr::join_by(gene_name == ensembl)
  ) |>
    dplyr::select(-identifier, -name) |>
    dplyr::right_join(network_data, by = dplyr::join_by("id" == "node_id")) |>
    dplyr::select(gene_name, name, id, identifier, logfc, pvalue, pagerank, beta, p_value, fdr)

  names(gene_data) = c(
    "gene_name", "gene_symbol", "id", "identifier",
    "experimental_logfc", "experimental_pvalue",
    "pagerank", "tkoi_beta", "tkoi_pvalue", "tkoi_fdr"
  )

  return(gene_data)
}
