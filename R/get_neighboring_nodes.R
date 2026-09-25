#' Retrieve Neighboring Nodes in a Knowledge Graph
#'
#' Returns every node within \code{degree_expansion} hops of a node in the
#' tKOI knowledge graph, including the node itself.
#'
#' @param gene_node_id A node ID (a vertex name of \code{subnetwork}), for
#'   example a gene's \code{id} in \code{\link{genes}}.
#' @param degree_expansion Number of hops to expand. \code{1} returns the node
#'   and its direct neighbors; \code{2} also returns their neighbors.
#' @param subnetwork An igraph object to search. Default \code{tkoi::tkoi_net}.
#'
#' @details
#' The search uses \code{igraph::ego()} and ignores edge direction.
#'
#' @return A character vector of node IDs.
#'
#' @examples
#' \dontrun{
#' gene_id = tkoi::genes$id[tkoi::genes$name == "TP53"]
#'
#' # The gene and its direct neighbors
#' get_neighboring_nodes(gene_id, degree_expansion = 1)
#'
#' # Everything within two hops
#' get_neighboring_nodes(gene_id, degree_expansion = 2)
#' }
#'
#' @seealso \code{\link[igraph]{ego}}, \code{\link{plot_network}}
#'
#' @export
get_neighboring_nodes = function(gene_node_id, degree_expansion, subnetwork = tkoi::tkoi_net) {
  if (!igraph::is_igraph(subnetwork)) {
    stop("`subnetwork` must be an igraph object.", call. = FALSE)
  }
  if (!is.character(gene_node_id) || length(gene_node_id) != 1 || is.na(gene_node_id)) {
    stop("`gene_node_id` must be a single node ID.", call. = FALSE)
  }
  if (!is.numeric(degree_expansion) || length(degree_expansion) != 1 ||
    is.na(degree_expansion) || degree_expansion < 0 ||
    degree_expansion != round(degree_expansion)) {
    stop("`degree_expansion` must be a non-negative whole number.", call. = FALSE)
  }
  if (!gene_node_id %in% igraph::V(subnetwork)$name) {
    stop("`gene_node_id` is not a vertex of `subnetwork`.", call. = FALSE)
  }

  neighboring_nodes = igraph::ego(
    subnetwork,
    order = degree_expansion,
    nodes = gene_node_id,
    mode = "all"
  )[[1]]
  names(neighboring_nodes)
}
