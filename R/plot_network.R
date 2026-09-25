#' Plot a Local Network Around a Target Node from tKOI Results
#'
#' Draws the part of the knowledge graph that links a target node (for
#' example an enriched GO term) to the significant genes near it. Genes are
#' colored by log fold change (blue down, white zero, red up), the target is
#' orange, other nodes are gray, and node size follows the tKOI effect size
#' (\code{beta}).
#'
#' @param tkoi_result A \code{tKOIList} returned by \code{\link{run_tkoi}}.
#' @param target_node_id Node ID (vertex name) of the node to center on.
#' @param degree_expansion Maximum number of hops between the target and a
#'   gene. Default \code{2}.
#' @param network_layout_type Layout algorithm: \code{"kk"}
#'   (Kamada-Kawai, the default), \code{"fr"} (Fruchterman-Reingold),
#'   \code{"gem"}, \code{"graphopt"}, \code{"lgl"}, or \code{"mds"}.
#' @param subnetwork The igraph network used for the analysis. Default
#'   \code{tkoi::tkoi_net}.
#'
#' @details
#' Significant genes pass the p-value and log fold change thresholds stored
#' in \code{tkoi_result}. The plot shows every node on a simple path of at
#' most \code{degree_expansion} edges between the target and one of these
#' genes. For \code{degree_expansion <= 2} these nodes are found directly from
#' neighbor sets; longer paths use \code{igraph::all_simple_paths()}, which
#' can be slow around highly connected nodes.
#'
#' @return The plotted igraph subgraph, invisibly.
#'
#' @import igraph
#' @importFrom grDevices colorRampPalette
#' @importFrom scales rescale
#'
#' @examples
#' \dontrun{
#' top_term = tkoi_result@network_summary_statistics$BiologicalProcess$node_id[1]
#' plot_network(tkoi_result, target_node_id = top_term, degree_expansion = 2)
#' }
#'
#' @export
plot_network = function(
  tkoi_result,
  target_node_id,
  degree_expansion = 2,
  network_layout_type = c("kk", "fr", "gem", "graphopt", "lgl", "mds"),
  subnetwork = tkoi::tkoi_net
) {
  network_layout_type = match.arg(network_layout_type)
  if (!igraph::is_igraph(subnetwork)) {
    stop("`subnetwork` must be an igraph object.", call. = FALSE)
  }
  vertex_names = igraph::V(subnetwork)$name
  if (!is.character(target_node_id) || length(target_node_id) != 1 ||
    !target_node_id %in% vertex_names) {
    stop("`target_node_id` must be a single vertex name of `subnetwork`.", call. = FALSE)
  }

  # Same gene rows as run_tkoi(): no blank names, first row per gene.
  expression_data = .tkoi_clean_expression(tkoi_result@expression_data)
  sig_genes = expression_data[
    !is.na(expression_data$pvalue) & !is.na(expression_data$logfc) &
      expression_data$pvalue <= tkoi_result@pvalue_threshold &
      abs(expression_data$logfc) >= tkoi_result@logfc_threshold, ,
    drop = FALSE
  ] |>
    dplyr::inner_join(tkoi::genes, by = dplyr::join_by("gene_name" == "ensembl"))

  neighboring_nodes = get_neighboring_nodes(
    target_node_id,
    degree_expansion = degree_expansion,
    subnetwork = subnetwork
  )
  gene_ids = setdiff(intersect(unique(sig_genes$id), neighboring_nodes), target_node_id)
  if (length(gene_ids) == 0) {
    stop("No significant genes lie within `degree_expansion` hops of the target.", call. = FALSE)
  }

  target = match(target_node_id, vertex_names)
  genes_idx = match(gene_ids, vertex_names)
  if (degree_expansion <= 2) {
    # Simple paths of length <= 2: the direct edge, or target - w - gene.
    target_neighbors = as.integer(igraph::neighbors(subnetwork, target, mode = "all"))
    gene_neighbors = igraph::adjacent_vertices(subnetwork, genes_idx, mode = "all")
    bridges = if (degree_expansion == 2) {
      unlist(lapply(gene_neighbors, function(x) intersect(as.integer(x), target_neighbors)))
    } else {
      integer(0)
    }
    cluster = unique(c(target, genes_idx, bridges))
  } else {
    paths = igraph::all_simple_paths(
      subnetwork,
      from = target,
      to = genes_idx,
      cutoff = degree_expansion,
      mode = "all"
    )
    cluster = unique(c(target, genes_idx, unlist(lapply(paths, as.integer))))
  }

  subnet = igraph::induced_subgraph(subnetwork, vids = cluster)
  subnet = igraph::delete_vertices(subnet, which(igraph::degree(subnet) == 0))
  sub_names = igraph::V(subnet)$name
  sub_identifier = igraph::V(subnet)$identifier
  if (is.null(sub_identifier)) {
    sub_identifier = sub_names
  }
  sub_type = .tkoi_vertex_types(subnet)
  if (is.null(sub_type)) {
    sub_type = rep("Unknown", length(sub_names))
  }

  # One label per vertex: gene symbol, else GO term name, else identifier.
  gene_rows = match(sub_names, tkoi::genes$id)
  go_rows = match(sub_identifier, tkoi::go_annotation$identifier)
  label = dplyr::coalesce(tkoi::genes$name[gene_rows], tkoi::go_annotation$name[go_rows], sub_identifier)

  expression_rows = match(tkoi::genes$ensembl[gene_rows], expression_data$gene_name)
  logfc = expression_data$logfc[expression_rows]

  beta_table = do.call(rbind, lapply(tkoi_result@network_summary_statistics, function(x) {
    data.frame(node_id = x$node_id, beta = x$beta)
  }))
  beta = beta_table$beta[match(sub_names, beta_table$node_id)]

  # Colors: genes on a blue-white-red logFC ramp, target orange, others gray.
  color_vector = rep("gray", length(sub_names))
  is_gene = sub_type == "Gene" & sub_names != target_node_id & !is.na(logfc)
  if (any(is_gene)) {
    # Symmetric around zero, so white always means no change.
    ramp = grDevices::colorRampPalette(c("#79ADD6", "white", "#BB2A2D"))(100)
    max_abs = max(abs(logfc[is_gene]))
    scaled = if (max_abs > 0) (logfc[is_gene] / max_abs + 1) / 2 else rep(0.5, sum(is_gene))
    color_vector[is_gene] = ramp[ceiling(scaled * 99) + 1]
  }
  color_vector[sub_names == target_node_id] = "orange"

  # Sizes: beta rescaled to 10-20. Nodes without a beta get the smallest
  # size; infinite betas are capped at the finite range.
  finite_beta = is.finite(beta)
  if (any(finite_beta)) {
    beta[!is.na(beta) & beta == Inf] = max(beta[finite_beta])
    beta[!is.finite(beta)] = min(beta[finite_beta])
    vertex_size = scales::rescale(beta, to = c(10, 20))
    vertex_size[!is.finite(vertex_size)] = 15
  } else {
    vertex_size = rep(15, length(beta))
  }

  igraph::V(subnet)$name = label
  layout_function = switch(network_layout_type,
    kk = igraph::layout_with_kk,
    fr = igraph::layout_with_fr,
    gem = igraph::layout_with_gem,
    graphopt = igraph::layout_with_graphopt,
    lgl = igraph::layout_with_lgl,
    mds = igraph::layout_with_mds
  )
  network_layout = layout_function(subnet)

  plot(
    subnet,
    layout = network_layout,
    vertex.label = label,
    vertex.color = color_vector,
    vertex.size = vertex_size,
    vertex.label.cex = 0.8,
    vertex.label.color = "black",
    asp = 0.8
  )
  invisible(subnet)
}
