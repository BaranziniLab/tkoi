#' Plot a Local Network Around a Target Node from tKOI Results
#'
#' This function visualizes a subnetwork of biologically relevant nodes around a specified
#' target node from a tKOI result. It includes all significant genes (based on p-value and logFC thresholds)
#' within a user-defined network neighborhood, as well as all connecting paths between the target and these genes.
#' The resulting network is colored by log fold change (logFC), sized by beta value, and plotted with
#' an igraph layout.
#'
#' @param tkoi_result An object of class `tkoi_result` containing differential expression results,
#'        thresholds, and network information from a tKOI analysis.
#' @param target_node_id A character string specifying the node ID (e.g., GO term) to center the network around.
#' @param degree_expansion Integer specifying the maximum graph distance from the target node to explore.
#'        Defaults to 2.
#' @param network_layout_type A character string indicating the layout algorithm to use for plotting.
#'        Options include: `"kk"` (Kamada-Kawai), `"fr"` (Fruchterman-Reingold),
#'        `"gem"` (Graph Embedding), `"graphopt"` (Graph Optimization),
#'        `"lgl"` (Large Graph Layout), and `"mds"` (Multidimensional Scaling).
#'
#' @details
#' The function first identifies differentially expressed genes that meet user-defined thresholds
#' and are within the specified graph neighborhood of the target node.
#' It then computes all simple paths between the target node and each significant gene,
#' compiles a list of nodes involved, and builds a subgraph.
#'
#' Each node is colored by logFC (gradient), or gray if not a gene. The target node is highlighted in orange.
#' Node sizes are scaled to beta values. The layout algorithm used can be selected with `network_layout_type`.
#'
#' @return A network plot is rendered directly using `plot.igraph()`.
#'
#' @import igraph
#' @importFrom dplyr filter inner_join select mutate left_join coalesce
#' @importFrom purrr map reduce
#' @importFrom grDevices colorRampPalette
#' @importFrom scales rescale
#'
#' @examples
#' \dontrun{
#' plot_network(tkoi_result = my_tkoi_result,
#'              target_node_id = "4:c77f6410-bc08-43ba-a172-0503ab1c93db:1234567",
#'              degree_expansion = 2,
#'              network_layout_type = "kk")
#' }
#'
#' @export
plot_network = function(tkoi_result, target_node_id, degree_expansion = 2,
                        network_layout_type = "kk"){

  sig_genes = tkoi_result@expression_data |>
    dplyr::filter(pvalue <= tkoi_result@pvalue_threshold) |>
    dplyr::filter(abs(logfc) >= tkoi_result@logfc_threshold) |>
    dplyr::inner_join(tkoi::genes, by = dplyr::join_by("gene_name" == "ensembl"))

  neighboring_nodes = tkoi::get_neighboring_nodes(target_node_id,
                                                  degree_expansion = degree_expansion)

  sig_genes = sig_genes |>
    dplyr::filter(id %in% neighboring_nodes)

  concept_cluster = c(sig_genes$id, target_node_id)
  for(gene_id in sig_genes$id){
    concept_cluster = c(concept_cluster,
                        igraph::all_simple_paths(tkoi::tkoi_net, from = target_node_id,
                                                 to = gene_id, cutoff = degree_expansion) |>
                          unlist() |>
                          names()
    )
  }

  concept_cluster = unique(concept_cluster)

  subnet = igraph::induced_subgraph(tkoi::tkoi_net, vids = concept_cluster)
  subnet = igraph::delete_vertices(subnet, which(igraph::degree(subnet) == 0))

  nodes = data.frame(
    target_node_id = igraph::vertex_attr(subnet)$name,
    identifier = igraph::vertex_attr(subnet)$identifier,
    node_type = igraph::vertex_attr(subnet)$labels
  ) |>
    dplyr::left_join(dplyr::select(tkoi::go_annotation, identifier, name), by = "identifier") |>
    dplyr::left_join(dplyr::mutate(dplyr::select(tkoi::genes, identifier, ensembl, name),
                                   identifier = as.character(identifier)),
                     by = "identifier",
                     suffix = c("_go", "_gene")) |>
    dplyr::mutate(label = dplyr::coalesce(name_gene, name_go, identifier)) |>
    dplyr::select(-name_gene, -name_go) |>
    dplyr::mutate(node_type = gsub("^\\['|']$", "", node_type)) |>
    dplyr::left_join(tkoi_result@expression_data, by = dplyr::join_by("ensembl" == "gene_name"))

  beta_df = purrr::map(tkoi_result@network_summary_statistics, function(x){
    dplyr::select(x, node_id, beta)
  }) |>
    purrr::reduce(rbind.data.frame)
  nodes = dplyr::left_join(nodes, beta_df, by = dplyr::join_by("target_node_id" == "node_id"))

  gradient_colors = c("#79ADD6", "white", "#BB2A2D")
  color_vector = rep(NA, nrow(nodes))
  color_vector[nodes$target_node_id == target_node_id] = "orange"
  is_gene = nodes$node_type == "Gene" & nodes$target_node_id != target_node_id
  logfc_vals = nodes$logfc[is_gene]
  scaled_vals = (logfc_vals - min(logfc_vals, na.rm = TRUE)) /
    (max(logfc_vals, na.rm = TRUE) - min(logfc_vals, na.rm = TRUE))

  color_vector[is_gene] = grDevices::colorRampPalette(gradient_colors)(100)[
    ceiling(scaled_vals * 99) + 1
  ]
  color_vector[is.na(color_vector)] = "gray"

  vertex_ids = V(subnet)$name
  beta_vector = nodes$beta[match(vertex_ids, nodes$target_node_id)]
  beta_vector[is.na(beta_vector)] = min(beta_vector, na.rm = TRUE)
  scaled_beta = scales::rescale(beta_vector, to = c(10, 20), na.rm = TRUE)

  V(subnet)$name = nodes$label

  if(network_layout_type == "kk"){
    network_layout = igraph::layout_with_kk(subnet)
  }

  if(network_layout_type == "fr"){
    network_layout = igraph::layout_with_fr(subnet)
  }

  if(network_layout_type == "gem"){
    network_layout = igraph::layout_with_gem(subnet)
  }

  if(network_layout_type == "graphopt"){
    network_layout = igraph::layout_with_graphopt(subnet)
  }

  if(network_layout_type == "lgl"){
    network_layout = igraph::layout_with_lgl(subnet)
  }

  if(network_layout_type == "mds"){
    network_layout = igraph::layout_with_mds(subnet)
  }

  plot(
    subnet,
    layout = network_layout,
    vertex.color = color_vector,
    vertex.size = scaled_beta,
    vertex.label.cex = 0.8,
    vertex.label.color = "black",
    vertex.label.family = "Arial",
    asp = 0.8
  )

}
