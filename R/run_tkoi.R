#' Run tKOI Analysis
#'
#' Performs tKOI (Transcriptomic Knowledge-graph-driven Omics Integration)
#' analysis on a differential expression table. Transcriptomic signals are
#' propagated through a biological knowledge graph with personalized PageRank,
#' and a degree-matched permutation null identifies enriched concepts. The
#' pipeline has these steps:
#'
#' \enumerate{
#'   \item \strong{Gene filtering:} genes are kept when \code{pvalue <=
#'         pvalue_threshold} and \code{abs(logfc) >= logfc_threshold}.
#'   \item \strong{Network mapping:} kept genes are mapped onto the network
#'         through their Ensembl IDs (see \code{\link{genes}}).
#'   \item \strong{Personalized PageRank:} PageRank is propagated from the
#'         seed genes, each weighted by \code{abs(logfc)}.
#'   \item \strong{Permutation null:} every seed gene is replaced by a random
#'         gene of similar degree, \code{n_permutation} times, and PageRank is
#'         recomputed for each replacement set.
#'   \item \strong{Network enrichment scoring:} each node's z-score
#'         (\code{beta}) and one-sided p-value compare its observed PageRank
#'         with the null.
#'   \item \strong{Node annotation:} nodes are annotated with curated
#'         metadata (GO terms, diseases, cell types, compounds, and more).
#'         Every node is reported; nodes without a curated annotation have
#'         missing annotation columns. Among compounds, only human
#'         metabolites (\code{\link{human_metabolites}}) are reported.
#'   \item \strong{Prioritization:} results are split by node type, adjusted
#'         for multiple testing (Benjamini-Hochberg FDR over the reported
#'         nodes of each type), and ranked.
#' }
#'
#' @param expression_data A data frame with columns \code{gene_name} (Ensembl
#'   gene IDs), \code{logfc}, and \code{pvalue}. Rows with a missing or blank
#'   \code{gene_name} are ignored, and only the first row of each gene is used.
#' @param subnetwork An igraph knowledge graph. Default \code{tkoi::tkoi_net}.
#' @param pvalue_threshold Keep genes with \code{pvalue} at or below this
#'   value. Default \code{0.05}.
#' @param logfc_threshold Keep genes with \code{abs(logfc)} at or above this
#'   value. Default \code{0.25}.
#' @param indirect_link_threshold Nodes within two hops of at least this many
#'   seed genes are ranked first in each result table. This only orders the
#'   tables; no node is excluded. Default \code{3}.
#' @param topology_similarity Number in \code{[0, 1]}. A replacement gene's
#'   degree must lie within \code{[topology_similarity * d, (2 -
#'   topology_similarity) * d]}, where \code{d} is the seed gene's degree in
#'   \code{tkoi::genes} (its degree in \code{tkoi_net}, also for custom
#'   networks).
#'   Default \code{0.9}.
#' @param n_permutation Number of permutations in the null (at least 2).
#'   Default \code{100}.
#' @param damping_factor PageRank damping factor in \code{(0, 1)}. Default
#'   \code{0.85}.
#' @param maximum_iteration Maximum solver iterations per PageRank vector.
#'   The solver needs about 55 iterations on \code{tkoi_net} at the default
#'   \code{tolerance}. Default \code{500}.
#' @param n_cores Number of CPU threads. \code{NULL} (the default) uses every
#'   available core, or \code{getOption("tkoi.n_cores")} when that is set.
#'   Requests above the number of available cores (which respects Slurm and
#'   Linux cgroup CPU limits) are capped. Results do not depend on this value.
#' @param keep_permutations If \code{TRUE} (the default), \code{pagerank_data}
#'   holds every permutation's PageRank vector (\code{perm.1}, \code{perm.2},
#'   ...). If \code{FALSE}, it holds only the null mean and standard deviation
#'   (\code{null_mean}, \code{null_sd}), which saves memory for large
#'   \code{n_permutation}.
#' @param tolerance Relative residual at which the PageRank solver stops. A
#'   vector is accepted when its residual is below \code{tolerance} both in
#'   the symmetric system the solver works on and, in 1-norm, in the original
#'   PageRank system, which bounds the 1-norm error of the PageRank vector by
#'   about \code{2 * tolerance}. The default \code{1e-14} is at least as
#'   accurate as \code{igraph::page_rank()} on \code{tkoi_net}, per node as
#'   well as overall.
#' @param verbose If \code{TRUE} (the default), print progress messages.
#'
#' @details
#' \strong{Speed.} The PageRank solver is written in C++. It solves the
#' symmetric form of the personalized PageRank system with conjugate
#' gradient, up to 16 PageRank vectors per pass over the network, on
#' \code{n_cores} threads. Two-hop neighborhood counts use bitsets, and the
#' permutation null is drawn in C++.
#'
#' \strong{Reproducibility.} The null seed sets are drawn in one sequence from
#' R's random number generator, so \code{set.seed()} fixes the result, and
#' the same seed gives bit-identical statistics for any \code{n_cores}. On
#' \code{tkoi_net}, the same seed also draws the same null seed sets as tkoi
#' 1.0.0 run with \code{n_cores = 1}, and PageRank values, \code{beta}, and
#' \code{p_value} agree with 1.0.0 to numerical precision. Result tables
#' differ from 1.0.0 by design: unannotated nodes are kept, the Compound FDR
#' family is the reported human metabolites, and nodes with a constant null
#' are untestable (see below).
#'
#' \strong{Custom networks.} Seed genes that are not vertices of
#' \code{subnetwork} are dropped before the null is drawn, so the observed run
#' and every null run use the same seeds with the same weights. Replacement
#' genes are drawn only from \code{\link{genes}} that are vertices of
#' \code{subnetwork}. Nodes with a missing type are grouped as
#' \code{"Unknown"}, and nodes of types without a curated annotation are
#' reported without annotation columns.
#'
#' \strong{Memory.} Beyond \code{tkoi_net} itself (about 0.5 GB), the network
#' matrix takes about 0.3 GB and the solver about 0.6 GB. Keeping all
#' permutations takes \code{8 * vcount(subnetwork) * (n_permutation + 1)}
#' bytes (about 0.75 GB for 100 permutations on \code{tkoi_net}). Before any
#' PageRank is computed, the run checks this against the available memory
#' (physical memory, or a Linux container limit) and stops with advice if it
#' would not fit.
#'
#' \strong{Untestable nodes.} When a node's null PageRank is the same in every
#' permutation (standard deviation 0), its z-score is undefined, so its
#' \code{beta}, \code{p_value}, and \code{fdr} are \code{NaN}; such nodes
#' are left out of the FDR adjustment and listed last. This covers nodes
#' that no seed gene or replacement gene reaches, and seed genes in small
#' disconnected components that no replacement gene reaches (which would
#' otherwise get an infinite \code{beta}). \code{direct_links} and
#' \code{indirect_links} are \code{NA} for nodes that have no seed gene within
#' one or two hops. A gene whose first row in \code{expression_data} has a
#' missing \code{pvalue} or \code{logfc} is not used as a seed.
#'
#' @return An S4 object of class \code{tKOIList} with slots:
#' \describe{
#'   \item{\code{expression_data}}{The input data frame.}
#'   \item{\code{pagerank_data}}{A data frame with \code{node_id},
#'     observed \code{pagerank}, and either every permutation
#'     (\code{perm.1}, ...) or \code{null_mean} and \code{null_sd}.}
#'   \item{\code{network_summary_statistics}}{A named list of tibbles, one per
#'     node type, with \code{node_id}, \code{node_type}, \code{identifier},
#'     \code{pagerank}, \code{beta}, \code{p_value}, \code{fdr},
#'     \code{direct_links}, \code{indirect_links}, \code{valency}, and
#'     annotation columns (missing for nodes without a curated annotation).
#'     Every node of the network is listed once, except Compound nodes that
#'     are not human metabolites.}
#' }
#'
#' @examples
#' \dontrun{
#' expression_data = data.table::fread(
#'   system.file("extdata", "example_data.csv", package = "tkoi")
#' )
#'
#' set.seed(1)
#' result = run_tkoi(
#'   expression_data = expression_data,
#'   n_permutation = 100
#' )
#'
#' result@network_summary_statistics$BiologicalProcess
#' }
#'
#' @seealso \code{\link{visualize_topn}}, \code{\link{export_gene_exploration_data}},
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
  n_cores = NULL,
  keep_permutations = TRUE,
  tolerance = 1e-14,
  verbose = TRUE
) {
  start_time = Sys.time()
  .tkoi_check_run_arguments(
    expression_data = expression_data,
    pvalue_threshold = pvalue_threshold,
    logfc_threshold = logfc_threshold,
    indirect_link_threshold = indirect_link_threshold,
    topology_similarity = topology_similarity,
    n_permutation = n_permutation,
    damping_factor = damping_factor,
    maximum_iteration = maximum_iteration,
    keep_permutations = keep_permutations,
    tolerance = tolerance,
    verbose = verbose
  )
  n_threads = .tkoi_resolve_cores(n_cores)
  say = function(...) if (verbose) message(...)

  say("Preparing the network...")
  network = .tkoi_prepare_network(subnetwork)

  # Seed genes -------------------------------------------------------------
  # Genes outside `subnetwork` are dropped here, so the observed run and every
  # null run use the same number of seeds with the same weights.
  gene_data = .tkoi_seed_genes(expression_data, pvalue_threshold, logfc_threshold)
  seed_node = match(gene_data$id, network$vertex_names)
  in_network = !is.na(seed_node)
  if (!any(in_network)) {
    stop("None of the selected genes are in `subnetwork`.", call. = FALSE)
  }
  gene_data = gene_data[in_network, , drop = FALSE]
  seed_node = seed_node[in_network]
  if (sum(gene_data$prob) <= 0) {
    stop("All selected genes in `subnetwork` have `logfc` equal to zero.", call. = FALSE)
  }

  # Null seed sets ----------------------------------------------------------
  universe = tkoi::genes
  universe_node = match(universe$id, network$vertex_names)
  universe_ok = !is.na(universe_node)
  universe_node = universe_node[universe_ok]
  universe_degree = universe$degree[universe_ok]

  # Check memory before any heavy work.
  block_width = .tkoi_plan_memory(
    n_vertices = network$n,
    n_nonzero = length(network$csr$col),
    n_seeds = nrow(gene_data),
    n_permutation = n_permutation,
    keep_permutations = keep_permutations
  )

  say(glue::glue("Sampling {n_permutation} degree-matched null seed sets for {nrow(gene_data)} genes..."))
  pools = .tkoi_candidate_pools(gene_data$degree, universe_degree, topology_similarity)
  null_draws = .tkoi_sample_null(pools$ptr, pools$idx, length(universe_node), as.integer(n_permutation))
  # Map sampled universe rows to vertices, keeping the genes x permutations shape.
  null_draws[] = universe_node[null_draws]

  # PageRank ----------------------------------------------------------------
  say(glue::glue(
    "Running personalized PageRank for the observed data and {n_permutation} ",
    "permutations on {n_threads} thread{if (n_threads == 1) '' else 's'}..."
  ))
  progress = NULL
  if (verbose) {
    bar = utils::txtProgressBar(min = 0, max = n_permutation + 1, style = 3)
    on.exit(close(bar), add = TRUE)
    progress = function(done, total) utils::setTxtProgressBar(bar, done)
  }
  pagerank = .tkoi_ppr_null(
    network$csr,
    as.integer(seed_node),
    gene_data$prob,
    null_draws,
    gene_data$prob,
    damping_factor,
    tolerance,
    as.integer(maximum_iteration),
    n_threads,
    block_width,
    keep_permutations,
    progress
  )
  unconverged = sum(pagerank$rel_residual > tolerance)
  if (unconverged > 0) {
    warning(
      glue::glue(
        "{unconverged} PageRank vector(s) did not reach `tolerance` within ",
        "`maximum_iteration` = {maximum_iteration} iterations."
      ),
      call. = FALSE
    )
  }

  # Statistics ----------------------------------------------------------------
  say("\nCalculating network enrichment statistics...")
  observed_pr = pagerank$observed
  z_scores = (observed_pr - pagerank$null_mean) / pagerank$null_sd
  # A null without variation cannot score a node: its z-score would be 0/0
  # or infinite (e.g. a seed gene in a small component that no replacement
  # gene reaches). Report such nodes as untestable (NaN).
  z_scores[!is.na(pagerank$null_sd) & pagerank$null_sd == 0] = NaN
  p_values = exp(stats::pnorm(z_scores, lower.tail = FALSE, log.p = TRUE))

  reach = .tkoi_seed_reach(network$csr$row_ptr, network$csr$col, as.integer(seed_node), n_threads)
  direct_links = as.numeric(reach$direct)
  direct_links[direct_links == 0] = NA
  indirect_links = as.numeric(reach$indirect)
  indirect_links[indirect_links == 0] = NA

  pagerank_data = .tkoi_pagerank_table(network$vertex_names, observed_pr, pagerank, keep_permutations)

  say("Annotating tKOI analysis results...")
  # Among compounds, only human metabolites are reported, and they form the
  # multiple-testing family of the Compound table.
  reported = !(network$node_type == "Compound" & !network$identifier %in% tkoi::human_metabolites)
  node_statistics = .tkoi_node_statistics(
    network = network,
    reported = reported,
    pagerank = observed_pr,
    beta = z_scores,
    p_value = p_values,
    direct_links = direct_links,
    indirect_links = indirect_links,
    indirect_link_threshold = indirect_link_threshold
  )
  node_statistics = .tkoi_annotate(node_statistics)

  result = methods::new(
    "tKOIList",
    expression_data = expression_data,
    pagerank_data = pagerank_data,
    network_summary_statistics = node_statistics,
    pvalue_threshold = pvalue_threshold,
    logfc_threshold = logfc_threshold,
    topology_similarity = topology_similarity,
    n_permutation = n_permutation,
    damping_factor = damping_factor,
    maximum_iteration = maximum_iteration
  )

  time_diff = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  time_label = if (time_diff < 60) {
    glue::glue("{round(time_diff, 2)} seconds")
  } else {
    glue::glue("{round(time_diff / 60, 2)} minutes")
  }
  say(glue::glue("Analysis finished and took {time_label}.\n"))

  result
}

#' Validate run_tkoi() arguments
#' @noRd
.tkoi_check_run_arguments = function(
  expression_data,
  pvalue_threshold,
  logfc_threshold,
  indirect_link_threshold,
  topology_similarity,
  n_permutation,
  damping_factor,
  maximum_iteration,
  keep_permutations,
  tolerance,
  verbose
) {
  if (!is.data.frame(expression_data)) {
    stop("`expression_data` must be a data frame.", call. = FALSE)
  }
  missing_columns = setdiff(c("gene_name", "logfc", "pvalue"), names(expression_data))
  if (length(missing_columns) > 0) {
    stop(
      "`expression_data` is missing column(s): ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.numeric(expression_data$logfc) || !is.numeric(expression_data$pvalue)) {
    stop("`logfc` and `pvalue` in `expression_data` must be numeric.", call. = FALSE)
  }

  is_number = function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x)
  is_whole = function(x) is_number(x) && x == round(x)
  is_flag = function(x) is.logical(x) && length(x) == 1 && !is.na(x)

  if (!is_number(pvalue_threshold) || pvalue_threshold < 0 || pvalue_threshold > 1) {
    stop("`pvalue_threshold` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is_number(logfc_threshold) || logfc_threshold < 0) {
    stop("`logfc_threshold` must be a non-negative number.", call. = FALSE)
  }
  if (!is_number(indirect_link_threshold)) {
    stop("`indirect_link_threshold` must be a number.", call. = FALSE)
  }
  if (!is_number(topology_similarity) || topology_similarity < 0 || topology_similarity > 1) {
    stop("`topology_similarity` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is_whole(n_permutation) || n_permutation < 2 || n_permutation > .Machine$integer.max) {
    stop("`n_permutation` must be a whole number of at least 2.", call. = FALSE)
  }
  if (!is_number(damping_factor) || damping_factor <= 0 || damping_factor >= 1) {
    stop("`damping_factor` must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is_whole(maximum_iteration) || maximum_iteration < 1 || maximum_iteration > .Machine$integer.max) {
    stop("`maximum_iteration` must be a whole number between 1 and .Machine$integer.max.", call. = FALSE)
  }
  if (!is_number(tolerance) || tolerance <= 0 || tolerance >= 1) {
    stop("`tolerance` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is_flag(keep_permutations)) {
    stop("`keep_permutations` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is_flag(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Select seed genes and their PageRank weights
#'
#' Keeps the first row per gene, maps Ensembl IDs to network genes, applies
#' the thresholds, and weights genes by `abs(logfc)`.
#' @noRd
.tkoi_seed_genes = function(expression_data, pvalue_threshold, logfc_threshold) {
  expression = .tkoi_clean_expression(expression_data)
  gene_data = expression |>
    dplyr::inner_join(tkoi::genes, by = dplyr::join_by(gene_name == ensembl)) |>
    dplyr::filter(pvalue <= pvalue_threshold) |>
    dplyr::filter(abs(logfc) >= logfc_threshold)

  if (nrow(gene_data) == 0) {
    stop(
      "No genes pass the thresholds: check that `gene_name` holds Ensembl gene IDs ",
      "and that some genes have `pvalue <= pvalue_threshold` and ",
      "`abs(logfc) >= logfc_threshold`.",
      call. = FALSE
    )
  }
  abs_logfc = abs(gene_data$logfc)
  if (any(!is.finite(abs_logfc))) {
    stop("Selected genes have infinite `logfc` values; remove or cap them.", call. = FALSE)
  }
  if (sum(abs_logfc) <= 0) {
    stop("All selected genes have `logfc` equal to zero.", call. = FALSE)
  }
  gene_data$prob = abs_logfc / sum(abs_logfc)
  gene_data[, c("id", "prob", "degree")]
}

#' Assemble the pagerank_data table without copying the PageRank vectors
#' @noRd
.tkoi_pagerank_table = function(vertex_names, observed, pagerank, keep_permutations) {
  columns = list(node_id = vertex_names, pagerank = observed)
  if (keep_permutations) {
    null_columns = pagerank$null_columns
    names(null_columns) = paste0("perm.", seq_along(null_columns))
    columns = c(columns, null_columns)
  } else {
    columns$null_mean = pagerank$null_mean
    columns$null_sd = pagerank$null_sd
  }
  structure(columns, class = "data.frame", row.names = vertex_names)
}

#' Per-node-type statistics tables
#'
#' Splits the reported nodes by type, adjusts p-values within each type
#' (Benjamini-Hochberg, over the reported nodes of that type), and orders
#' rows: nodes within two hops of at least `indirect_link_threshold` seeds
#' first, then by FDR, then by decreasing beta.
#' @noRd
.tkoi_node_statistics = function(
  network,
  reported,
  pagerank,
  beta,
  p_value,
  direct_links,
  indirect_links,
  indirect_link_threshold
) {
  node_type = network$node_type
  fdr = rep(NA_real_, length(p_value))
  by_type = split(which(reported), node_type[reported])
  for (rows in by_type) {
    fdr[rows] = stats::p.adjust(p_value[rows], method = "fdr")
  }
  met_threshold = as.numeric(indirect_links >= indirect_link_threshold)

  # Same ordering as the original two dplyr::arrange() passes, done once.
  # Keys are rounded to 10 significant digits so that nodes whose statistics
  # are equal up to floating-point noise keep their network order on every
  # platform.
  keep = which(reported)
  ordering = data.frame(
    row = keep,
    node_type = node_type[keep],
    met_threshold = met_threshold[keep],
    fdr = signif(fdr[keep], 10),
    beta = signif(beta[keep], 10)
  ) |>
    dplyr::arrange(node_type, fdr, dplyr::desc(beta)) |>
    dplyr::arrange(node_type, dplyr::desc(met_threshold), fdr)

  rows_by_type = split(ordering$row, ordering$node_type)
  types = sort(unique(node_type[keep]), method = "radix")

  tables = lapply(types, function(type) {
    rows = rows_by_type[[type]]
    dplyr::as_tibble(data.frame(
      node_id = network$vertex_names[rows],
      node_type = node_type[rows],
      identifier = network$identifier[rows],
      pagerank = pagerank[rows],
      beta = beta[rows],
      p_value = p_value[rows],
      fdr = fdr[rows],
      direct_links = direct_links[rows],
      indirect_links = indirect_links[rows],
      valency = network$valency[rows],
      stringsAsFactors = FALSE
    ))
  })
  names(tables) = types
  tables
}

#' Attach curated annotations to each node-type table
#'
#' Each table is left-joined with its annotation dataset by `identifier`
#' (first row per identifier), so every node is kept and nodes without an
#' annotation get missing annotation columns.
#' @noRd
.tkoi_annotate = function(tables) {
  gene_annotation = function() {
    genes = as.data.frame(tkoi::genes)
    genes$identifier = as.character(genes$identifier)
    genes[, setdiff(names(genes), c("id", "degree")), drop = FALSE]
  }
  annotations = list(
    Anatomy = function() tkoi::anatomy_annotation,
    BiologicalProcess = function() tkoi::go_annotation,
    CellType = function() tkoi::celltype_annotation,
    CellularComponent = function() tkoi::go_annotation,
    ClinicalLab = function() tkoi::clinicallab_annotation,
    Complex = function() tkoi::complex_annotation,
    Compound = function() tkoi::compound_annotation,
    Disease = function() tkoi::disease_annotation,
    EC = function() tkoi::ec_annotation,
    Gene = gene_annotation,
    MiRNA = function() tkoi::mirna_annotation,
    MolecularFunction = function() tkoi::go_annotation,
    Pathway = function() tkoi::pathway_annotation,
    Protein = function() tkoi::protein_annotation,
    ProteinDomain = function() tkoi::proteindomain_annotation,
    ProteinFamily = function() tkoi::proteinfamily_annotation,
    PwGroup = function() tkoi::pwgroup_annotation,
    Reaction = function() tkoi::reaction_annotation
  )

  for (type in intersect(names(annotations), names(tables))) {
    annotation = as.data.frame(annotations[[type]]())
    annotation = annotation[!duplicated(annotation$identifier), , drop = FALSE]
    tables[[type]] = dplyr::left_join(
      tables[[type]],
      annotation,
      by = "identifier",
      relationship = "many-to-one"
    )
  }
  tables
}
