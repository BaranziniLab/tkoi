#' Retrieve the Graph Used for a tKOI Analysis
#'
#' Returns the exact igraph object supplied to [run_tkoi()], including its
#' vertex and edge attributes. Save the result with [saveRDS()] to retain
#' this network for later contextualization or agent traversal.
#'
#' @param tkoi_result A `tKOIList` returned by [run_tkoi()].
#'
#' @details
#' Results created before tkoi 1.3.0 may not contain an analysis graph. This
#' function refuses those results instead of substituting the currently
#' installed `tkoi_net`, which may differ from the network actually analyzed.
#' Rerun the analysis with its original graph to create a complete result.
#'
#' @return The stored igraph object.
#' @examples
#' \dontrun{
#' result = run_tkoi(expression_data, subnetwork = tkoi::tkoi_net)
#' saveRDS(result, "analysis.rds")
#' graph = get_analysis_graph(readRDS("analysis.rds"))
#' get_neighboring_nodes(result@pagerank_data$node_id[1], 1, subnetwork = graph)
#' }
#' @export
get_analysis_graph = function(tkoi_result) {
  if (!methods::is(tkoi_result, "tKOIList")) {
    stop("`tkoi_result` must be a tKOIList.", call. = FALSE)
  }
  graph = tkoi_result@subnetwork
  if (!igraph::is_igraph(graph)) {
    stop(
      "This result does not contain its analysis graph. Rerun with tkoi >= 1.3.0 ",
      "and the original subnetwork; do not substitute another network.",
      call. = FALSE
    )
  }
  if (!identical(as.character(igraph::V(graph)$name), as.character(tkoi_result@pagerank_data$node_id))) {
    stop("The stored graph vertices do not match the analysis result.", call. = FALSE)
  }
  graph
}
