#' Run tKOI Analysis (Optimized)
#'
#' Vectorized and pre-computation optimized version of \code{\link{run_tkoi}}.
#' Produces identical output to \code{run_tkoi} when \code{n_cores = 1} and
#' the same RNG seed is used. With \code{n_cores > 1}, permutations run in
#' parallel (macOS/Linux only); results are statistically equivalent but
#' numerically different due to independent per-worker RNG streams.
#'
#' Key optimizations over \code{run_tkoi}:
#' \enumerate{
#'   \item \strong{Batch ego calls}: metapath filters use one \code{igraph::ego()} call
#'         for all genes instead of one call per gene.
#'   \item \strong{Vectorized prob vectors}: probability vectors are built with
#'         \code{match()} + direct indexing instead of dplyr join chains.
#'   \item \strong{Pre-computed candidate pools}: degree-based gene pools are computed
#'         once before the permutation loop, eliminating \code{n_perm x n_genes}
#'         \code{dplyr::filter} calls.
#'   \item \strong{Pre-allocated perm matrix}: permutation PageRank results are stored
#'         in a pre-allocated matrix instead of a growing list.
#'   \item \strong{Vectorized statistics}: z-scores and p-values are computed with
#'         \code{rowMeans}/\code{sweep}/\code{rowSums} instead of row-by-row
#'         \code{lapply}.
#'   \item \strong{Optional parallel permutations}: \code{n_cores > 1} runs
#'         permutations in parallel via \code{parallel::mclapply}.
#' }
#'
#' Per-step wall-clock timings are printed on completion and stored as
#' \code{attr(result, "step_timings")}.
#'
#' @param expression_data A data frame with columns \code{gene_name}, \code{logfc}, \code{pvalue}.
#' @param subnetwork An igraph object. Default is \code{tkoi::tkoi_net}.
#' @param pvalue_threshold Numeric p-value filter threshold. Default \code{0.05}.
#' @param logfc_threshold Numeric log fold-change filter threshold. Default \code{0.25}.
#' @param indirect_link_threshold Minimum 2-hop gene connections for prioritization. Default \code{3}.
#' @param topology_similarity Degree similarity for substitute gene sampling. Default \code{0.9}.
#' @param n_permutation Number of permutations. Default \code{100}.
#' @param damping_factor PageRank damping factor. Default \code{0.85}.
#' @param maximum_iteration PageRank max iterations. Default \code{500}.
#' @param n_cores Number of cores for parallel permutations. Default \code{1} (sequential,
#'   identical output to \code{run_tkoi}). Values > 1 require macOS or Linux.
#'
#' @return An S4 object of class \code{tKOIList} (identical structure to \code{run_tkoi}).
#'   Per-step timings (in seconds) are attached as \code{attr(result, "step_timings")}.
#'
#' @seealso \code{\link{run_tkoi}}
#' @export
run_tkoi_fast = function(
    expression_data,
    subnetwork = tkoi::tkoi_net,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.25,
    indirect_link_threshold = 3,
    topology_similarity = 0.9,
    n_permutation = 100,
    damping_factor = 0.85,
    maximum_iteration = 500,
    n_cores = 1
){
  t_wall = proc.time()
  timings = list()

  # ── Step 1: decrypt network + gene filtering ──────────────────────────────
  t0 = proc.time()

  if(is.raw(subnetwork)){
    subnetwork = sodium::data_decrypt(subnetwork, tkoi:::network_attributes)
    subnetwork = memDecompress(subnetwork, type = "gzip")
    subnetwork = unserialize(subnetwork)
  }

  gene_data = expression_data |>
    dplyr::distinct(gene_name, .keep_all = TRUE) |>
    dplyr::inner_join(tkoi::genes, by = dplyr::join_by(gene_name == ensembl)) |>
    dplyr::filter(pvalue <= pvalue_threshold) |>
    dplyr::filter(abs(logfc) >= logfc_threshold)

  abs_logfc = abs(gene_data$logfc)
  gene_data[["prob"]] = abs_logfc / sum(abs_logfc)
  gene_data = gene_data |> dplyr::select(id, prob, degree)

  timings$gene_filtering = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  [%.2fs] Gene filtering (%d genes)", timings$gene_filtering, nrow(gene_data)))

  # ── Step 2: metapath filters — batch ego (one call vs n_genes calls) ──────
  t0 = proc.time()
  message("Generating metapath filters...")

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

  timings$metapath_filters = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  [%.2fs] Metapath filters", timings$metapath_filters))

  # ── Step 3: observed PageRank — vectorized prob vector ────────────────────
  t0 = proc.time()

  vertex_names = igraph::V(subnetwork)$name
  n_vertices = length(vertex_names)

  # build prob vector with match() instead of dplyr join chain
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

  timings$observed_pagerank = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  [%.2fs] Observed PageRank", timings$observed_pagerank))

  # ── Step 4: pre-compute candidate pools (eliminates n_perm×n_genes dplyr calls) ──
  t0 = proc.time()

  all_gene_ids = tkoi::genes$id
  all_gene_degrees = tkoi::genes$degree
  lower_bounds = topology_similarity * gene_data$degree
  upper_bounds = (2 - topology_similarity) * gene_data$degree
  n_genes = nrow(gene_data)

  candidate_pools = lapply(seq_len(n_genes), function(i){
    all_gene_ids[all_gene_degrees >= lower_bounds[i] & all_gene_degrees <= upper_bounds[i]]
  })

  timings$candidate_pool_precompute = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  [%.2fs] Candidate pool pre-computation", timings$candidate_pool_precompute))

  # ── Step 5: permutation test ───────────────────────────────────────────────
  t0 = proc.time()
  message("Running permutation test...")

  # closure captures candidate_pools, gene_data$prob, vertex_names, n_genes, n_vertices
  do_perm = function(perm_idx){
    substitute_genes = character(n_genes)  # pre-allocated; "" not in any pool
    for(i in seq_len(n_genes)){
      pool = setdiff(candidate_pools[[i]], substitute_genes)
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

  if(n_cores > 1){
    perm_results = parallel::mclapply(seq_len(n_permutation), do_perm, mc.cores = n_cores)
  } else {
    pb = txtProgressBar(min = 1, max = n_permutation, style = 3)
    perm_results = lapply(seq_len(n_permutation), function(n){
      setTxtProgressBar(pb, n)
      do_perm(n)
    })
    close(pb)
  }

  # pre-allocated matrix — avoids growing list + bind_cols overhead
  perm_matrix = matrix(
    unlist(perm_results),
    nrow = n_vertices,
    ncol = n_permutation,
    dimnames = list(vertex_names, paste0("perm.", seq_len(n_permutation)))
  )

  timings$permutation = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("\n  [%.2fs] Permutation test (%d perms, %d core(s))",
                  timings$permutation, n_permutation, n_cores))

  # ── Step 6: vectorized statistics (replaces row-by-row lapply) ────────────
  t0 = proc.time()
  message("Calculating network enrichment statistics...")

  perm_mean = rowMeans(perm_matrix)
  centered = sweep(perm_matrix, 1L, perm_mean)        # (n_vertices × n_perm) - mean per row
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

  timings$statistics = as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  [%.2fs] Vectorized statistics", timings$statistics))

  # ── Step 7: node annotation (identical to run_tkoi) ───────────────────────
  t0 = proc.time()
  message("\nAnnotating tKOI analysis results...")

  annotated_node_stats = data.frame(
    node_id = igraph::vertex_attr(subnetwork)$name,
    node_type = igraph::vertex_attr(subnetwork)$label,
    identifier = igraph::vertex_attr(subnetwork)$identifier
  ) |>
    dplyr::inner_join(node_stats, by = "node_id") |>
    dplyr::group_split(node_type) |>
    purrr::map(function(x){
      x |>
        dplyr::mutate(fdr = p.adjust(p_value, method = "fdr")) |>
        dplyr::arrange(fdr, dplyr::desc(beta))
    })

  node_types = unlist(purrr::map(annotated_node_stats, function(x) x$node_type[1]))
  node_types = gsub(pattern = "[[]", replacement = "", node_types)
  node_types = gsub(pattern = "[]]", replacement = "", node_types)
  node_types = gsub(pattern = "[']", replacement = "", node_types)
  names(annotated_node_stats) = node_types

  annotated_node_stats = purrr::map(annotated_node_stats, function(x){
    node_type = x$node_type[1]
    node_type = gsub(pattern = "[[]", replacement = "", node_type)
    node_type = gsub(pattern = "[]]", replacement = "", node_type)
    node_type = gsub(pattern = "[']", replacement = "", node_type)
    x$node_type = node_type
    x = dplyr::left_join(x, direct_links_table, by = "node_id") |>
      dplyr::left_join(indirect_links_table, by = "node_id") |>
      dplyr::left_join(valency_table, by = "node_id") |>
      dplyr::mutate(met_treshold = as.numeric(indirect_links >= indirect_link_threshold)) |>
      dplyr::arrange(desc(met_treshold), fdr) |>
      dplyr::select(-met_treshold)
    return(x)
  })

  annotated_node_stats$Anatomy = dplyr::inner_join(annotated_node_stats$Anatomy, tkoi::anatomy_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$BiologicalProcess = dplyr::inner_join(annotated_node_stats$BiologicalProcess, tkoi::go_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$CellType = dplyr::inner_join(annotated_node_stats$CellType, tkoi::celltype_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$CellularComponent = dplyr::inner_join(annotated_node_stats$CellularComponent, tkoi::go_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$ClinicalLab = dplyr::inner_join(annotated_node_stats$ClinicalLab, tkoi::clinicallab_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Complex = dplyr::inner_join(annotated_node_stats$Complex, tkoi::complex_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Compound = dplyr::inner_join(annotated_node_stats$Compound, tkoi::compound_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Compound = dplyr::filter(annotated_node_stats$Compound, identifier %in% tkoi::human_metabolites)
  annotated_node_stats$Disease = dplyr::inner_join(annotated_node_stats$Disease, tkoi::disease_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$EC = dplyr::inner_join(annotated_node_stats$EC, tkoi::ec_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Gene = dplyr::inner_join(annotated_node_stats$Gene,
                                                dplyr::mutate(tkoi::genes, identifier = as.character(identifier)) |>
                                                  dplyr::select(-id, -degree),
                                                by = "identifier", relationship = "many-to-many")
  annotated_node_stats$MiRNA = dplyr::inner_join(annotated_node_stats$MiRNA, tkoi::mirna_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$MolecularFunction = dplyr::inner_join(annotated_node_stats$MolecularFunction, tkoi::go_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Pathway = dplyr::inner_join(annotated_node_stats$Pathway, tkoi::pathway_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Protein = dplyr::inner_join(annotated_node_stats$Protein, tkoi::protein_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$ProteinDomain = dplyr::inner_join(annotated_node_stats$ProteinDomain, tkoi::proteindomain_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$ProteinFamily = dplyr::inner_join(annotated_node_stats$ProteinFamily, tkoi::proteinfamily_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$PwGroup = dplyr::inner_join(annotated_node_stats$PwGroup, tkoi::pwgroup_annotation, by = "identifier", relationship = "many-to-many")
  annotated_node_stats$Reaction = dplyr::inner_join(annotated_node_stats$Reaction, tkoi::reaction_annotation, by = "identifier", relationship = "many-to-many")

  timings$annotation = as.numeric((proc.time() - t0)["elapsed"])
  timings$total = as.numeric((proc.time() - t_wall)["elapsed"])

  message(sprintf("  [%.2fs] Annotation", timings$annotation))
  message(sprintf("\n── Step timings (run_tkoi_fast) ─────────────────────"))
  for(step in names(timings)){
    message(sprintf("  %-35s %6.2fs", step, timings[[step]]))
  }
  message(sprintf("─────────────────────────────────────────────────────\n"))

  result = new("tKOIList",
               expression_data = expression_data,
               pagerank_data = gene_pagerank,
               network_summary_statistics = annotated_node_stats,
               pvalue_threshold = pvalue_threshold,
               logfc_threshold = logfc_threshold,
               topology_similarity = topology_similarity,
               n_permutation = n_permutation,
               damping_factor = damping_factor,
               maximum_iteration = maximum_iteration
  )

  attr(result, "step_timings") = timings
  return(result)
}
