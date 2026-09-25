# Reference implementation: run_tkoi() as released before the native engine
# (tkoi 1.0.0, sequential n_cores = 1 path), kept verbatim apart from:
#   * purrr::map() replaced by the equivalent lapply() (purrr is no longer a dependency);
#   * no decryption step (networks are plain igraph objects now);
#   * `gene_universe` so toy networks can restrict null sampling to the genes
#     they contain (the default, tkoi::genes, is the released behavior);
#   * no progress bar or messages.
# Used only to check that the rewrite reproduces the released statistics.
legacy_run_tkoi = function(
    expression_data,
    subnetwork,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.25,
    indirect_link_threshold = 3,
    topology_similarity = 0.9,
    n_permutation = 100,
    damping_factor = 0.85,
    maximum_iteration = 500,
    gene_universe = tkoi::genes
) {
  gene_data = expression_data |>
    dplyr::distinct(gene_name, .keep_all = TRUE) |>
    dplyr::inner_join(tkoi::genes, by = dplyr::join_by(gene_name == ensembl)) |>
    dplyr::filter(pvalue <= pvalue_threshold) |>
    dplyr::filter(abs(logfc) >= logfc_threshold)

  abs_logfc = abs(gene_data$logfc)
  gene_data[["prob"]] = abs_logfc / sum(abs_logfc)
  gene_data = gene_data |> dplyr::select(id, prob, degree)

  direct_ego = igraph::ego(subnetwork, order = 1, nodes = gene_data$id, mode = "all")
  direct_counts = table(unlist(lapply(direct_ego, names)))
  direct_links_table = data.frame(
    node_id = names(direct_counts),
    direct_links = as.numeric(direct_counts),
    stringsAsFactors = FALSE
  )

  indirect_ego = igraph::ego(subnetwork, order = 2, nodes = gene_data$id, mode = "all")
  indirect_counts = table(unlist(lapply(indirect_ego, names)))
  indirect_links_table = data.frame(
    node_id = names(indirect_counts),
    indirect_links = as.numeric(indirect_counts),
    stringsAsFactors = FALSE
  )

  valency_table = data.frame(
    node_id = igraph::vertex_attr(subnetwork, "name"),
    valency = igraph::vertex_attr(subnetwork, "degree")
  )

  vertex_names = igraph::V(subnetwork)$name
  n_vertices = length(vertex_names)

  prob_vec = numeric(n_vertices)
  gene_idx = match(gene_data$id, vertex_names)
  valid_idx = !is.na(gene_idx)
  prob_vec[gene_idx[valid_idx]] = gene_data$prob[valid_idx]

  observed_pr = igraph::page_rank(
    subnetwork,
    personalized = prob_vec,
    algo = "prpack",
    damping = damping_factor,
    directed = FALSE,
    options = list(maxiter = maximum_iteration, eps = 1e-20)
  )$vector

  all_gene_ids = gene_universe$id
  all_gene_degrees = gene_universe$degree
  lower_bounds = topology_similarity * gene_data$degree
  upper_bounds = (2 - topology_similarity) * gene_data$degree
  n_genes = nrow(gene_data)

  candidate_pools = lapply(seq_len(n_genes), function(i) {
    all_gene_ids[all_gene_degrees >= lower_bounds[i] & all_gene_degrees <= upper_bounds[i]]
  })

  do_perm = function(perm_idx) {
    substitute_genes = character(n_genes)
    for (i in seq_len(n_genes)) {
      pool = setdiff(candidate_pools[[i]], substitute_genes)
      if (length(pool) == 0L) pool = candidate_pools[[i]]
      substitute_genes[i] = sample(pool, 1)
    }

    perm_prob = numeric(n_vertices)
    pidx = match(substitute_genes, vertex_names)
    pv = !is.na(pidx)
    perm_prob[pidx[pv]] = gene_data$prob[pv]

    igraph::page_rank(
      subnetwork,
      personalized = perm_prob,
      algo = "prpack",
      damping = damping_factor,
      directed = FALSE,
      options = list(maxiter = maximum_iteration, eps = 1e-20)
    )$vector
  }

  perm_results = lapply(seq_len(n_permutation), do_perm)

  perm_matrix = matrix(
    unlist(perm_results),
    nrow = n_vertices,
    ncol = n_permutation,
    dimnames = list(vertex_names, paste0("perm.", seq_len(n_permutation)))
  )

  perm_mean = rowMeans(perm_matrix)
  centered = sweep(perm_matrix, 1L, perm_mean)
  perm_sd = sqrt(rowSums(centered^2) / (n_permutation - 1L))

  z_scores = (observed_pr - perm_mean) / perm_sd
  p_values = exp(pnorm(z_scores, lower.tail = FALSE, log.p = TRUE))

  node_stats = data.frame(
    node_id = vertex_names,
    pagerank = observed_pr,
    beta = z_scores,
    p_value = p_values,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  gene_pagerank = cbind(
    data.frame(node_id = vertex_names, pagerank = observed_pr,
               row.names = NULL, stringsAsFactors = FALSE),
    as.data.frame(perm_matrix)
  )

  annotated_node_stats = data.frame(
    node_id = igraph::vertex_attr(subnetwork)$name,
    node_type = igraph::vertex_attr(subnetwork)$label,
    identifier = igraph::vertex_attr(subnetwork)$identifier
  ) |>
    dplyr::inner_join(node_stats, by = "node_id") |>
    dplyr::group_split(node_type) |>
    lapply(function(x) {
      x |>
        dplyr::mutate(fdr = p.adjust(p_value, method = "fdr")) |>
        dplyr::arrange(fdr, dplyr::desc(beta))
    })

  node_types = unlist(lapply(annotated_node_stats, function(x) x$node_type[1]))
  node_types = gsub(pattern = "[[]", replacement = "", node_types)
  node_types = gsub(pattern = "[]]", replacement = "", node_types)
  node_types = gsub(pattern = "[']", replacement = "", node_types)
  names(annotated_node_stats) = node_types

  annotated_node_stats = lapply(annotated_node_stats, function(x) {
    node_type = x$node_type[1]
    node_type = gsub(pattern = "[[]", replacement = "", node_type)
    node_type = gsub(pattern = "[]]", replacement = "", node_type)
    node_type = gsub(pattern = "[']", replacement = "", node_type)
    x$node_type = node_type
    x = dplyr::left_join(x, direct_links_table, by = "node_id") |>
      dplyr::left_join(indirect_links_table, by = "node_id") |>
      dplyr::left_join(valency_table, by = "node_id") |>
      dplyr::mutate(met_treshold = as.numeric(indirect_links >= indirect_link_threshold)) |>
      dplyr::arrange(dplyr::desc(met_treshold), fdr) |>
      dplyr::select(-met_treshold)
    return(x)
  })

  annotate = function(type, annotation) {
    if (is.null(annotated_node_stats[[type]])) return(NULL)
    dplyr::inner_join(annotated_node_stats[[type]], annotation, by = "identifier",
                      relationship = "many-to-many")
  }
  annotated_node_stats$Anatomy = annotate("Anatomy", tkoi::anatomy_annotation)
  annotated_node_stats$BiologicalProcess = annotate("BiologicalProcess", tkoi::go_annotation)
  annotated_node_stats$CellType = annotate("CellType", tkoi::celltype_annotation)
  annotated_node_stats$CellularComponent = annotate("CellularComponent", tkoi::go_annotation)
  annotated_node_stats$ClinicalLab = annotate("ClinicalLab", tkoi::clinicallab_annotation)
  annotated_node_stats$Complex = annotate("Complex", tkoi::complex_annotation)
  annotated_node_stats$Compound = annotate("Compound", tkoi::compound_annotation)
  if (!is.null(annotated_node_stats$Compound)) {
    annotated_node_stats$Compound = dplyr::filter(annotated_node_stats$Compound,
                                                  identifier %in% tkoi::human_metabolites)
  }
  annotated_node_stats$Disease = annotate("Disease", tkoi::disease_annotation)
  annotated_node_stats$EC = annotate("EC", tkoi::ec_annotation)
  annotated_node_stats$Gene = annotate(
    "Gene",
    dplyr::mutate(tkoi::genes, identifier = as.character(identifier)) |>
      dplyr::select(-id, -degree)
  )
  annotated_node_stats$MiRNA = annotate("MiRNA", tkoi::mirna_annotation)
  annotated_node_stats$MolecularFunction = annotate("MolecularFunction", tkoi::go_annotation)
  annotated_node_stats$Pathway = annotate("Pathway", tkoi::pathway_annotation)
  annotated_node_stats$Protein = annotate("Protein", tkoi::protein_annotation)
  annotated_node_stats$ProteinDomain = annotate("ProteinDomain", tkoi::proteindomain_annotation)
  annotated_node_stats$ProteinFamily = annotate("ProteinFamily", tkoi::proteinfamily_annotation)
  annotated_node_stats$PwGroup = annotate("PwGroup", tkoi::pwgroup_annotation)
  annotated_node_stats$Reaction = annotate("Reaction", tkoi::reaction_annotation)

  list(pagerank_data = gene_pagerank, network_summary_statistics = annotated_node_stats)
}
