#' Run tKOI Analysis
#'
#' Performs tKOI (Transcriptomic Knowledge-graph-driven Omics Integration) analysis
#' on a gene expression dataset. The analysis integrates transcriptomic measurements
#' with a biological knowledge graph to identify and interpret biologically enriched
#' subnetworks. The pipeline includes the following steps:
#'
#' \enumerate{
#'   \item \strong{Gene filtering:} Significant genes are selected based on user-defined
#'         thresholds for p-value and log fold change.
#'   \item \strong{Network mapping:} Filtered genes are mapped onto the biological
#'         subnetwork using Ensembl gene identifiers.
#'   \item \strong{Personalized PageRank:} A PageRank propagation is performed using
#'         log fold change-weighted probabilities, quantifying each node's influence.
#'   \item \strong{Permutation testing:} The propagation is repeated \code{n_permutation}
#'         times using substitute genes with similar topological properties to build a
#'         null distribution.
#'   \item \strong{Network enrichment scoring:} Each node's z-score (beta) and p-value
#'         are computed by comparing observed PageRank to the permutation null.
#'   \item \strong{Node annotation:} Nodes are annotated with domain-specific metadata
#'         (GO terms, diseases, cell types, compounds) from curated knowledge resources.
#'   \item \strong{Prioritization:} Results are organized by node type with FDR adjustment
#'         and connectivity statistics for biological prioritization.
#' }
#'
#' @param expression_data A data frame with columns \code{gene_name} (Ensembl IDs),
#'   \code{logfc}, and \code{pvalue}.
#' @param subnetwork An igraph object representing the biological knowledge graph.
#'   Default is \code{tkoi::tkoi_net}.
#' @param pvalue_threshold Numeric threshold for filtering genes by p-value. Default \code{0.05}.
#' @param logfc_threshold Numeric threshold for filtering genes by absolute log fold change.
#'   Default \code{0.25}.
#' @param indirect_link_threshold Minimum number of seed genes that must be 2-hop-connected
#'   to a node for it to be prioritized. Default \code{3}.
#' @param topology_similarity Numeric in (0, 1) controlling how similar in degree a
#'   substitute gene must be to the original during permutation sampling. Default \code{0.9}.
#' @param n_permutation Integer number of permutations for the null distribution. Default \code{100}.
#' @param damping_factor Damping factor for personalized PageRank. Default \code{0.85}.
#' @param maximum_iteration Maximum iterations for PageRank convergence. Default \code{500}.
#' @param n_cores Number of cores for parallel permutation computation. Default \code{1}
#'   (sequential). Values \code{> 1} use \code{parallel::mclapply} and require macOS or Linux;
#'   note that parallel mode produces statistically equivalent but numerically different results
#'   due to independent per-worker RNG streams.
#'
#' @details
#' Probability vectors and metapath filters are built with vectorized operations
#' (\code{match()}, batch \code{igraph::ego()}) rather than per-gene dplyr operations.
#' Candidate gene pools for permutation sampling are pre-computed once before the loop.
#' Enrichment z-scores are computed with \code{rowMeans}/\code{sweep}/\code{rowSums}
#' across the full permutation matrix rather than row-by-row.
#'
#' @return An S4 object of class \code{tKOIList} with slots:
#' \describe{
#'   \item{\code{expression_data}}{The input gene expression data frame.}
#'   \item{\code{pagerank_data}}{A data frame of observed and permutation PageRank scores
#'     for every node in the subnetwork.}
#'   \item{\code{network_summary_statistics}}{A named list of data frames, one per node type,
#'     containing beta, p_value, fdr, and biological annotations.}
#' }
#'
#' @examples
#' \dontrun{
#' expression_data = data.frame(
#'   gene_name = c("ENSG00000000003", "ENSG00000000005"),
#'   logfc = c(1.5, -0.8),
#'   pvalue = c(0.01, 0.02)
#' )
#'
#' result = run_tkoi(
#'   expression_data = expression_data,
#'   subnetwork = tkoi::tkoi_net,
#'   pvalue_threshold = 0.05,
#'   logfc_threshold = 0.25,
#'   n_permutation = 100
#' )
#'
#' result@pagerank_data
#' result@network_summary_statistics
#' }
#'
#' @seealso \code{\link[igraph]{page_rank}}, \code{\link{compute_network_enrichment}},
#'   \code{\link{run_gene_enrichment}}
#'
#' @export
run_tkoi = function(
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
  start_time = Sys.time()

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

  all_gene_ids = tkoi::genes$id
  all_gene_degrees = tkoi::genes$degree
  lower_bounds = topology_similarity * gene_data$degree
  upper_bounds = (2 - topology_similarity) * gene_data$degree
  n_genes = nrow(gene_data)

  candidate_pools = lapply(seq_len(n_genes), function(i){
    all_gene_ids[all_gene_degrees >= lower_bounds[i] & all_gene_degrees <= upper_bounds[i]]
  })

  message("Running permutation test...")

  do_perm = function(perm_idx){
    substitute_genes = character(n_genes)
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

  perm_matrix = matrix(
    unlist(perm_results),
    nrow = n_vertices,
    ncol = n_permutation,
    dimnames = list(vertex_names, paste0("perm.", seq_len(n_permutation)))
  )

  message(" ")
  message("Calculating network enrichment statistics...")

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

  end_time = Sys.time()
  time_diff = end_time - start_time
  message(glue::glue("Analysis finished and took {round(time_diff, 2)} minutes.\n"))

  return(result)
}
