#' Retrieve Neighboring Nodes in a Knowledge Graph
#'
#' This function returns all neighboring nodes for a given gene node ID in the tKOI knowledge graph.
#' The neighborhood is defined by the number of hops (`degree_expansion`) from the input node.
#'
#' @param gene_node_id A character string representing the unique node ID of a gene in the `tkoi::tkoi_net` igraph object.
#' @param degree_expansion An integer indicating the radius (in number of hops) to expand the search for neighbors.
#' For example, `degree_expansion = 1` returns first-degree (direct) neighbors, while `degree_expansion = 2` includes both direct and second-degree (indirect) neighbors.
#' @param subnetwork An `igraph` object (typically a subgraph of `tkoi_net`) in which
#'   the neighborhood search will be performed.
#'
#' @details
#' The function leverages the `igraph::ego()` function to extract the subgraph containing all nodes within a specified path length
#' from the target node. It then returns the node IDs (as character vector) of these neighbors.
#' This is useful in the tKOI pipeline for assessing local network topology around significantly altered genes.
#'
#' @return A character vector of node IDs representing the neighbors of the specified gene node.
#'
#' @examples
#' \dontrun{
#' # Get direct neighbors of a gene node with ID "ENSG00000123456"
#' get_neighboring_nodes("ENSG00000123456", degree_expansion = 1)
#'
#' # Get both direct and indirect (2-hop) neighbors
#' get_neighboring_nodes("ENSG00000123456", degree_expansion = 2)
#' }
#'
#' @seealso \code{\link[igraph]{ego}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @export
get_neighboring_nodes = function(gene_node_id, degree_expansion, subnetwork){
  neighboring_nodes = igraph::ego(subnetwork, order = degree_expansion,
                                  nodes = gene_node_id, mode = "all")[[1]]
  neighboring_nodes = names(neighboring_nodes)
  return(neighboring_nodes)
}
