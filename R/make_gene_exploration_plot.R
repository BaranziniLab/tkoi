#' Create Gene Exploration Plot
#'
#' Plots, for every gene of the experiment that is a gene node of the network,
#' the experimental p-value against the tKOI false discovery rate (FDR), with
#' up- and down-regulated genes in separate panels. Genes supported by both
#' the experiment and the network appear in the top right of each panel.
#'
#' @param tkoi_list A `tKOIList` object returned by [run_tkoi()]. It uses the
#'   `expression_data` slot (columns `gene_name`, `logfc`, and `pvalue`), the
#'   `Gene` table of `network_summary_statistics` (columns `node_id`, `name`,
#'   and `fdr`), and the `pvalue_threshold` and `logfc_threshold` of the run.
#' @param sig_color Fill color of genes whose `abs(logfc)` is at least the
#'   run's `logfc_threshold`. Default `"#F39B7FB2"`.
#' @param non_sig_color Fill color of the other genes. Default `"gray"`.
#'
#' @details
#' `expression_data` is cleaned the same way as in [run_tkoi()]: rows with a
#' missing or blank `gene_name` are dropped, and only the first row of each
#' gene is used. Genes are matched to network nodes through their Ensembl IDs
#' (see [genes]); network genes without expression data are not plotted.
#'
#' In the plot:
#'   - The x-axis is `-log10(pvalue)` (experimental p-value).
#'   - The y-axis is `-log10(fdr)` (tKOI FDR).
#'   - Genes with `logfc >= 0` are in the `Up-regulated` panel, the others in
#'     the `Down-regulated` panel.
#'   - The dashed vertical line marks the run's `pvalue_threshold` and the
#'     dashed horizontal line an FDR of 0.05.
#'
#' @return A `ggplot` object.
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
#' plt = make_gene_exploration_plot(tkoi_result, sig_color = "#F39B7FB2", non_sig_color = "gray")
#' plt
#' }
#'
#' @seealso [export_gene_exploration_data()] for the underlying table.
#'
#' @export
make_gene_exploration_plot = function(tkoi_list,
                                      sig_color = "#F39B7FB2",
                                      non_sig_color = "gray") {
  # Same gene rows as run_tkoi(): no blank names, first row per gene.
  expression_data = .tkoi_clean_expression(tkoi_list@expression_data)
  network_data = tkoi_list@network_summary_statistics$Gene
  p_value_treshold = tkoi_list@pvalue_threshold
  logfc_threshold = tkoi_list@logfc_threshold

  gene_data = dplyr::inner_join(expression_data, tkoi::genes,
    by = dplyr::join_by(gene_name == ensembl)
  ) |>
    dplyr::select(-identifier) |>
    dplyr::right_join(dplyr::select(network_data, -name), by = dplyr::join_by("id" == "node_id")) |>
    dplyr::mutate(direction = ifelse(logfc >= 0, "Up-regulated", "Down-regulated")) |>
    dplyr::filter(!is.na(direction)) |>
    dplyr::mutate(fill = ifelse(abs(logfc) >= logfc_threshold, sig_color, non_sig_color))

  plt = ggplot2::ggplot(gene_data, ggplot2::aes(
    x = -log(pvalue, base = 10),
    y = -log(fdr, base = 10)
  )) +
    ggplot2::geom_point(
      fill = gene_data[["fill"]],
      color = "black", size = 3, shape = 21
    ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -log(0.05, base = 10)), linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = -log(p_value_treshold, base = 10)), linetype = "dashed") +
    ggplot2::xlab("-log(Experimental P.Value)") +
    ggplot2::ylab("-log(tKOI FDR)") +
    ggplot2::facet_grid(. ~ direction) +
    ggplot2::theme_classic(base_size = 14)

  return(plt)
}
