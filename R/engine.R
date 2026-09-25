# Internal helpers for run_tkoi(). The heavy lifting happens in
# src/tkoi_engine.cpp; these functions prepare its inputs and shape its
# outputs.

# Cache of the most recently built network matrix.
.tkoi_cache = new.env(parent = emptyenv())

#' Resolve the number of threads to use
#'
#' `NULL` means "all available cores". Requests are capped at the number of
#' cores the machine reports, at a Linux cgroup (v1 or v2) CPU quota or a
#' Slurm CPU allocation when present, and at 2 when R CMD check limits cores.
#'
#' @param n_cores `NULL` or a positive whole number.
#' @return A positive integer.
#' @noRd
.tkoi_resolve_cores = function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores = getOption("tkoi.n_cores", NULL)
  }
  available = .tkoi_available_cores()

  if (is.null(n_cores)) {
    return(available)
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) ||
    n_cores < 1 || n_cores != round(n_cores)) {
    stop("`n_cores` must be NULL or a single positive whole number.", call. = FALSE)
  }
  as.integer(min(n_cores, available))
}

#' Number of cores this process may use
#' @noRd
.tkoi_available_cores = function() {
  available = suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(available) || available < 1) {
    available = max(1L, .tkoi_hardware_threads())
  }

  # Slurm allocation.
  slurm = suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (!is.na(slurm) && slurm >= 1) {
    available = min(available, floor(slurm))
  }

  # Linux cgroup v2 CPU quota ("max 100000" means unlimited).
  cpu_max = .tkoi_read_first_line("/sys/fs/cgroup/cpu.max")
  if (!is.null(cpu_max)) {
    fields = strsplit(cpu_max, " ", fixed = TRUE)[[1]]
    quota = suppressWarnings(as.numeric(fields[1]))
    period = suppressWarnings(as.numeric(fields[2]))
    if (!is.na(quota) && !is.na(period) && period > 0) {
      available = min(available, max(1, floor(quota / period)))
    }
  }

  # Linux cgroup v1 CPU quota (-1 means unlimited).
  quota_v1 = .tkoi_read_first_line("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
  period_v1 = .tkoi_read_first_line("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
  if (!is.null(quota_v1) && !is.null(period_v1)) {
    quota = suppressWarnings(as.numeric(quota_v1))
    period = suppressWarnings(as.numeric(period_v1))
    if (!is.na(quota) && !is.na(period) && quota > 0 && period > 0) {
      available = min(available, max(1, floor(quota / period)))
    }
  }

  limit_env = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit_env) && !identical(toupper(limit_env), "FALSE")) {
    available = min(available, 2L)
  }
  as.integer(max(1, available))
}

#' Memory this process may use, in bytes (NA if unknown)
#'
#' Physical memory, lowered to a Linux cgroup memory limit when one is set.
#' @noRd
.tkoi_memory_limit = function() {
  total = .tkoi_total_memory()
  for (path in c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    line = .tkoi_read_first_line(path)
    limit = if (is.null(line)) NA else suppressWarnings(as.numeric(line))
    if (!is.na(limit) && limit > 0 && (is.na(total) || limit < total)) {
      total = limit
    }
  }
  total
}

#' First line of a small system file, or NULL
#' @noRd
.tkoi_read_first_line = function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    return(NULL)
  }
  line = tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character(0))
  if (length(line) == 0) NULL else trimws(line)
}

#' Prepare a network for the native engine
#'
#' Builds the normalized sparse matrix and reads the vertex metadata the
#' pipeline needs. The matrix for the most recent graph is cached and reused
#' when the graph structure (`igraph::graph_id()`) and its `weight` edge
#' attribute are unchanged, so repeated runs on `tkoi_net` skip that step.
#' Vertex metadata is read afresh on every call.
#'
#' @param subnetwork An igraph object.
#' @return A list with `csr`, `vertex_names`, `node_type`, `identifier`,
#'   `valency`, and `n`.
#' @noRd
.tkoi_prepare_network = function(subnetwork) {
  if (is.raw(subnetwork)) {
    stop(
      "`subnetwork` is a raw vector. Encrypted networks are no longer supported; ",
      "pass an igraph object such as `tkoi::tkoi_net`.",
      call. = FALSE
    )
  }
  if (!igraph::is_igraph(subnetwork)) {
    stop("`subnetwork` must be an igraph object.", call. = FALSE)
  }

  n = igraph::vcount(subnetwork)
  if (n == 0) {
    stop("`subnetwork` has no vertices.", call. = FALSE)
  }
  vertex_names = igraph::V(subnetwork)$name
  if (is.null(vertex_names)) {
    stop("`subnetwork` must have a vertex attribute `name`.", call. = FALSE)
  }
  node_type = .tkoi_vertex_types(subnetwork)
  if (is.null(node_type)) {
    stop("`subnetwork` must have a vertex attribute `labels` with node types.", call. = FALSE)
  }
  vertex_attrs = igraph::vertex_attr_names(subnetwork)
  identifier = if ("identifier" %in% vertex_attrs) {
    igraph::vertex_attr(subnetwork, "identifier")
  } else {
    vertex_names
  }
  valency = if ("degree" %in% vertex_attrs) {
    igraph::vertex_attr(subnetwork, "degree")
  } else {
    igraph::degree(subnetwork)
  }

  # igraph::page_rank() uses a `weight` edge attribute by default; so do we.
  weight = if ("weight" %in% igraph::edge_attr_names(subnetwork)) {
    as.numeric(igraph::E(subnetwork)$weight)
  } else {
    numeric(0)
  }

  key = igraph::graph_id(subnetwork)
  cached = .tkoi_cache$network
  if (!is.null(cached) && identical(cached$key, key) && identical(cached$weight, weight)) {
    csr = cached$csr
  } else {
    edges = igraph::as_edgelist(subnetwork, names = FALSE)
    csr = .tkoi_build_csr(
      as.integer(edges[, 1]),
      as.integer(edges[, 2]),
      weight,
      as.integer(n)
    )
    rm(edges)
    .tkoi_cache$network = list(key = key, weight = weight, csr = csr)
  }

  list(
    csr = csr,
    vertex_names = vertex_names,
    node_type = node_type,
    identifier = identifier,
    valency = valency,
    n = n
  )
}

#' Clean node types of a network, or NULL if it has none
#'
#' Reads the `labels` vertex attribute (or `label`), strips Neo4j formatting,
#' and names missing types "Unknown".
#' @noRd
.tkoi_vertex_types = function(subnetwork) {
  type_attr = intersect(c("labels", "label"), igraph::vertex_attr_names(subnetwork))
  if (length(type_attr) == 0) {
    return(NULL)
  }
  .tkoi_clean_node_type(igraph::vertex_attr(subnetwork, type_attr[1]))
}

#' Drop the cached network
#' @noRd
.tkoi_clear_cache = function() {
  rm(list = ls(.tkoi_cache), envir = .tkoi_cache)
  invisible(NULL)
}

#' Strip Neo4j label formatting such as `"['Gene']"`; missing types become "Unknown"
#' @noRd
.tkoi_clean_node_type = function(x) {
  x = gsub("[][']", "", as.character(x))
  x[is.na(x) | !nzchar(x)] = "Unknown"
  x
}

#' Clean a differential expression table
#'
#' Keeps `gene_name`, `logfc`, and `pvalue`, drops rows with a missing or
#' blank `gene_name` (which would otherwise match every gene without an
#' Ensembl ID), and keeps the first row of each gene (even when that row has
#' a missing value). Every function that reads `expression_data` uses this,
#' so they all see the same genes.
#' @noRd
.tkoi_clean_expression = function(expression_data) {
  expression = data.frame(
    gene_name = as.character(expression_data$gene_name),
    logfc = as.numeric(expression_data$logfc),
    pvalue = as.numeric(expression_data$pvalue),
    stringsAsFactors = FALSE
  )
  expression = expression[!is.na(expression$gene_name) & nzchar(expression$gene_name), , drop = FALSE]
  expression = expression[!duplicated(expression$gene_name), , drop = FALSE]
  rownames(expression) = NULL
  expression
}

#' Degree-matched candidate pools for null sampling
#'
#' For every seed gene, the pool holds all genes in the universe whose degree
#' lies in `[topology_similarity * degree, (2 - topology_similarity) * degree]`,
#' in universe order. Pools are returned in compressed form for the engine.
#'
#' @param seed_degree Degrees of the seed genes.
#' @param universe_degree Degrees of every gene that may be drawn.
#' @param topology_similarity Number in `[0, 1]`.
#' @return A list with `ptr` (0-based offsets) and `idx` (1-based universe
#'   indices).
#' @noRd
.tkoi_candidate_pools = function(seed_degree, universe_degree, topology_similarity) {
  # Pools depend on the degree only, so build one pool per distinct degree.
  unique_degree = unique(seed_degree)
  pools = lapply(unique_degree, function(degree) {
    lower = topology_similarity * degree
    upper = (2 - topology_similarity) * degree
    which(universe_degree >= lower & universe_degree <= upper)
  })
  pools = pools[match(seed_degree, unique_degree)]
  sizes = lengths(pools)
  if (any(sizes == 0)) {
    stop(
      "No degree-matched candidates exist for some seed genes; ",
      "lower `topology_similarity`.",
      call. = FALSE
    )
  }
  list(
    ptr = as.integer(c(0, cumsum(sizes))),
    idx = as.integer(unlist(pools, use.names = FALSE))
  )
}

#' Plan memory use of the PageRank engine
#'
#' Picks the widest block of simultaneous PageRank vectors that fits in
#' memory, and stops early with advice if the requested outputs cannot fit.
#' Results do not depend on the block width.
#'
#' @param n_vertices,n_nonzero Size of the network and of its matrix.
#' @param n_seeds Number of seed genes.
#' @param n_permutation Number of null seed sets.
#' @param keep_permutations Whether every null vector is kept.
#' @return The block width (16, 8, 4, 2, or 1).
#' @noRd
.tkoi_plan_memory = function(n_vertices, n_nonzero, n_seeds, n_permutation, keep_permutations) {
  bytes = 8
  total = .tkoi_memory_limit()
  graph_bytes = n_nonzero * 12 + n_vertices * 24
  # Null draws are held as integer matrices (sampled indices and vertices).
  seed_bytes = 3 * 4 * n_seeds * n_permutation
  output_bytes = n_vertices * bytes * (if (keep_permutations) n_permutation + 1 else 4)
  work_bytes = function(width) 5 * n_vertices * width * bytes

  width = 16
  if (!is.na(total)) {
    budget = 0.75 * total
    fixed = graph_bytes + seed_bytes + output_bytes
    while (width > 1 && fixed + work_bytes(width) > budget) {
      width = width / 2
    }
    needed = fixed + work_bytes(width)
    if (needed > budget) {
      advice = if (keep_permutations) {
        " Set `keep_permutations = FALSE` or reduce `n_permutation`."
      } else {
        " Reduce `n_permutation` or use a smaller network."
      }
      stop(
        sprintf(
          "run_tkoi() needs about %.1f GB but only %.1f GB of memory is available.%s",
          needed / 1e9, total / 1e9, advice
        ),
        call. = FALSE
      )
    }
  }
  as.integer(width)
}
