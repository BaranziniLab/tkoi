#' Run tKOI Analysis
#'
#' This function performs a tKOI (Transcriptomic Knowledge-graph-driven Omics Integration) analysis on a gene expression dataset.
#' The analysis integrates transcriptomic measurements with a biological knowledge graph to identify and interpret biologically enriched subnetworks.
#' The pipeline includes the following steps:
#'
#' \enumerate{
#'   \item \strong{Gene filtering:} Significant genes are selected based on user-defined thresholds for p-value and log fold change.
#'   \item \strong{Network mapping:} These genes are mapped onto a predefined biological subnetwork using unique gene identifiers.
#'   \item \strong{Personalized PageRank:} A PageRank-based propagation is performed using log fold change-weighted probabilities, quantifying the influence of each node in the network.
#'   \item \strong{Permutation testing:} To assess statistical significance, the same propagation is repeated multiple times using randomly substituted genes with similar topological properties.
#'   \item \strong{Network enrichment scoring:} Each node's empirical enrichment score (beta) and p-value are computed by comparing observed PageRank to the permutation null distribution.
#'   \item \strong{Node annotation:} Nodes are annotated with domain-specific metadata (e.g., GO terms, diseases, cell types, compounds), drawn from curated knowledge resources linked to the graph.
#'   \item \strong{Prioritization:} Annotated results are organized by node type, with adjusted FDR values and connectivity statistics (e.g., direct/indirect neighbors) used to prioritize biologically relevant findings.
#' }
#'
#' The result is a comprehensive, annotated list of network nodes enriched in the context of the transcriptomic changes observed in the input dataset.

#' @param expression_data A data frame containing the gene expression data with columns for gene identifiers, log fold changes, and p-values.
#' @param subnetwork An igraph object representing the subnetwork to be analyzed. Default is `tkoi::tkoi_net`.
#' @param pvalue_threshold A numeric value specifying the threshold for filtering genes based on p-values. Default is 0.05.
#' @param logfc_threshold A numeric value specifying the threshold for filtering genes based on log fold changes. Default is 0.25.
#' @param indirect_link_threshold A numeric value specifying the threshold for the least amount of genes to be 2-hop-connected to the propagated node. Default is 3.
#' @param topology_similarity A numeric value (between 0 and 1) defining the degree of similarity in topology for selecting substitute genes during permutation. Default is 0.9.
#' @param n_permutation An integer specifying the number of permutations for the network enrichment test. Default is 100.
#' @param damping_factor A numeric value for the damping factor used in the personalized PageRank calculation. Default is 0.85.
#' @param maximum_iteration An integer specifying the maximum number of iterations for the PageRank algorithm. Default is 500.
#'
#' @details
#' The function takes gene expression data and integrates it with a predefined knowledge graph subnetwork. Genes are filtered based on their p-values, and personalized PageRank scores are computed to measure their importance within the network.
#' To assess statistical significance, a permutation test is performed where the genes are replaced by substitutes with similar network topology. The results are annotated with network node information and functional descriptions.
#'
#' The output includes:
#' - **Personalized PageRank scores** for genes.
#' - **Permutation-based network enrichment statistics**.
#' - **Annotated summary statistics** for network nodes.
#'
#' @return An S4 object of class `tKOIList`, containing the following slots:
#' - `expression_data`: The input gene expression data.
#' - `pagerank_data`: A data frame of personalized PageRank scores for nodes in the subnetwork.
#' - `network_summary_statistics`: A list of data frames, one for each node type, containing network enrichment statistics and functional annotations.
#'
#' @examples
#' \dontrun{
#' # Example gene expression data
#' expression_data = data.frame(
#'   gene_name = c("gene1", "gene2", "gene3"),
#'   logfc = c(1.5, -0.8, 2.3),
#'   pvalue = c(0.01, 0.05, 0.02)
#' )
#'
#' # Run tKOI analysis
#' result = run_tkoi(
#'   expression_data = expression_data,
#'   subnetwork = tkoi::tkoi_net,
#'   pvalue_threshold = 0.05,
#'   logfc_threshold = 0.5,
#'   n_permutation = 100
#' )
#'
#' # Access the results
#' result@pagerank_data
#' result@network_summary_statistics
#' }
#'
#' @seealso \code{\link[igraph]{page_rank}}, \code{\link[dplyr]{dplyr}}, \code{\link[tkoi]{compute_stat}}
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
    maximum_iteration = 500
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
  prob = abs_logfc / sum(abs_logfc)
  gene_data[["prob"]] = prob

  gene_data = gene_data |>
    dplyr::select(id, prob, degree)

  message("Generating metapath filters...")
  direct_links_table = sapply(gene_data$id, tkoi::get_neighboring_nodes,
                              degree_expansion = 1, subnetwork = subnetwork) |>
    unlist() |>
    table() |>
    as.data.frame()
  names(direct_links_table) = c("node_id", "direct_links")
  direct_links_table = direct_links_table |>
    dplyr::arrange(desc(direct_links))

  indirect_links_table = sapply(gene_data$id, tkoi::get_neighboring_nodes,
                                degree_expansion = 2, subnetwork = subnetwork) |>
    unlist() |>
    table() |>
    as.data.frame()
  names(indirect_links_table) = c("node_id", "indirect_links")
  indirect_links_table = indirect_links_table |>
    dplyr::arrange(desc(indirect_links))

  valency_table = data.frame(
    node_id = igraph::vertex_attr(subnetwork, "name"),
    valency = igraph::vertex_attr(subnetwork, "degree")
  )

  gene_prob = tkoi::genes |>
    dplyr::select(id) |>
    dplyr::left_join(gene_data, by = "id") |>
    dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))

  prob_vec = data.frame(id = igraph::V(subnetwork)$name) |>
    dplyr::left_join(gene_prob, by = "id") |>
    dplyr::select(id, prob) |>
    dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  prob_vec = prob_vec[["prob"]]

  gene_pagerank = igraph::page_rank(subnetwork, personalized = prob_vec,
                                    algo = "prpack",
                                    damping = damping_factor,
                                    directed = FALSE,
                                    options = list(maxiter = maximum_iteration, eps = 1e-20))$vector
  gene_pagerank = data.frame(node_id = names(gene_pagerank),
                             pagerank = gene_pagerank)
  rownames(gene_pagerank) = NULL


  permute_gene_pagerank = list()
  pb = txtProgressBar(min = 1, max = n_permutation, style = 3)
  message("Running permutation test...")
  for(n in 1:n_permutation){
    setTxtProgressBar(pb, n)
    substitute_gene_list = c()
    for(i in 1:dim(gene_data)[1]){
      gene = gene_data$id[i]
      prob = gene_data$prob[i]
      target_degree = gene_data$degree[i]

      substitute_gene = sample(
        (tkoi::genes |>
           dplyr::filter(!(id %in% substitute_gene_list)) |>
           dplyr::filter(degree <= (2-topology_similarity)*target_degree) |>
           dplyr::filter(degree >= (topology_similarity)*target_degree))$id, 1)
      substitute_gene_list = c(substitute_gene_list, substitute_gene)
    }

    permute_gene_prob = data.frame(
      id = substitute_gene_list,
      prob = gene_data$prob
    )

    permute_prob_vec = data.frame(id = igraph::V(subnetwork)$name) |>
      dplyr::left_join(permute_gene_prob, by = "id") |>
      dplyr::select(id, prob) |>
      dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
    permute_prob_vec = permute_prob_vec[["prob"]]

    permute_gene_pagerank[[n]] = igraph::page_rank(subnetwork, personalized = permute_prob_vec,
                                                   algo = "prpack",
                                                   damping = damping_factor,
                                                   directed = FALSE,
                                                   options = list(maxiter = maximum_iteration, eps = 1e-20))$vector
  }

  names(permute_gene_pagerank) = paste0("perm.", seq(1:n_permutation))
  permute_gene_pagerank = dplyr::bind_cols(permute_gene_pagerank, .name_repair = "unique")

  gene_pagerank = cbind.data.frame(gene_pagerank, permute_gene_pagerank)

  node_stats = gene_pagerank |>
    dplyr::select(node_id, pagerank)
  node_stats$beta = NA
  node_stats$p_value = NA

  pb = txtProgressBar(min = 1, max = nrow(gene_pagerank), style = 3)
  message(" ")
  message("Calculating network enrichment statistics...")
  stats = lapply(1:nrow(gene_pagerank), function(i){
    setTxtProgressBar(pb, i)
    tkoi::compute_network_enrichment(gene_pagerank[i,])
  })
  stats = do.call(rbind, stats)

  node_stats$beta = stats$beta
  node_stats$p_value = stats$p_value

  annotated_node_stats = data.frame(
    node_id = igraph::vertex_attr(subnetwork)$name,
    node_type = igraph::vertex_attr(subnetwork)$label,
    identifier = igraph::vertex_attr(subnetwork)$identifier
  ) |>
    dplyr::inner_join(node_stats, by = "node_id")|>
    dplyr::group_split(node_type) |>
    purrr::map(function(x){
      x = x |>
        dplyr::mutate(fdr = p.adjust(p_value, method = "fdr")) |>
        dplyr::arrange(fdr, dplyr::desc(beta))
      return(x)
    })

  node_types = unlist(purrr::map(annotated_node_stats, function(x){
    x$node_type[1]
  }))
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

  message(" ")
  message("\nAnnotating tKOI analysis results...")
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



