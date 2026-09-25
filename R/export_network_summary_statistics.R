#' Export Network Summary Statistics to Excel
#'
#' Writes the `network_summary_statistics` slot of a [run_tkoi()] result to an
#' Excel workbook, with one sheet per node type (for example `Gene`,
#' `BiologicalProcess`, `Disease`).
#'
#' @param tkoi_result A `tKOIList` object returned by [run_tkoi()].
#' @param filename Path of the Excel file to write. An existing file is
#'   overwritten. Default `"tkoi_result.xlsx"`.
#'
#' @return The path of the written file, invisibly, as returned by
#'   [writexl::write_xlsx()]. The function is called for its side effect of
#'   writing the workbook.
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
#' path = export_network_summary_statistics(
#'   tkoi_result,
#'   filename = file.path(tempdir(), "tkoi_network_statistics.xlsx")
#' )
#' path
#' }
#'
#' @seealso [export_gene_exploration_data()] for a gene-level table.
#'
#' @importFrom writexl write_xlsx
#' @export
export_network_summary_statistics = function(tkoi_result, filename = "tkoi_result.xlsx") {
  network_summary_statistics = tkoi_result@network_summary_statistics
  write_xlsx(network_summary_statistics, path = filename)
}
