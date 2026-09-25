# Tests of run_tkoi() on the toy network (see helper-fixtures.R).

# Shared results --------------------------------------------------------------

# The default toy run and the released algorithm run with the same seed are
# used by several tests, so compute each once per session. Tests that mock
# package functions call toy_result() before mocking, so the cached run is
# never computed under a mock.
toy_result = local({
  cache = NULL
  function() {
    if (is.null(cache)) {
      cache <<- toy_run()
    }
    cache
  }
})

toy_legacy_result = local({
  cache = NULL
  function() {
    if (is.null(cache)) {
      cache <<- withr::with_seed(42, legacy_run_tkoi(
        expression_data = toy_expression(),
        subnetwork = toy_network(),
        n_permutation = 10,
        gene_universe = toy_universe()
      ))
    }
    cache
  }
})

# Helpers ---------------------------------------------------------------------

# Seed genes selected and weighted as documented for run_tkoi(): first row per
# gene, joined to tkoi::genes by Ensembl ID, then filtered.
toy_seeds = function(expression_data = toy_expression(), pvalue_threshold = 0.05, logfc_threshold = 0.25) {
  expression_data |>
    dplyr::distinct(gene_name, .keep_all = TRUE) |>
    dplyr::inner_join(as.data.frame(tkoi::genes), by = dplyr::join_by(gene_name == ensembl)) |>
    dplyr::filter(pvalue <= pvalue_threshold, abs(logfc) >= logfc_threshold) |>
    dplyr::mutate(prob = abs(logfc) / sum(abs(logfc)))
}

# run_tkoi() with quiet, single-threaded defaults on any input.
quiet_run = function(expression_data = toy_expression(), subnetwork = toy_network(), seed = 42, n_permutation = 10,
                     ...) {
  withr::with_seed(seed, run_tkoi(
    expression_data,
    subnetwork,
    n_permutation = n_permutation,
    n_cores = 1,
    verbose = FALSE,
    ...
  ))
}

# Personalized PageRank reset vector of the toy seeds on `network`.
toy_reset_vector = function(network = toy_network()) {
  seeds = toy_seeds()
  reset = numeric(igraph::vcount(network))
  reset[match(seeds$id, igraph::V(network)$name)] = seeds$prob
  reset
}

# igraph's personalized PageRank of the toy seeds on `network`.
igraph_toy_pagerank = function(network) {
  reference = igraph::page_rank(
    network,
    personalized = toy_reset_vector(network),
    algo = "prpack",
    damping = 0.85,
    directed = FALSE
  )$vector
  unname(reference)
}

clean_node_type = function(network = toy_network()) {
  gsub("[][']", "", igraph::V(network)$labels)
}

# TRUE where x and y agree to a relative `tolerance` (missing matches missing).
# The documented row order: nodes with indirect_links >= threshold first, then
# the others, then nodes without any seed within two hops; within each group by
# fdr ascending, then beta descending, then network order.
expect_tkoi_row_order = function(table, threshold, network = toy_network()) {
  met = table$indirect_links >= threshold
  group = ifelse(is.na(met), 3, ifelse(met, 1, 2))
  position = match(table$node_id, igraph::V(network)$name)
  expected = order(group, signif(table$fdr, 10), -signif(table$beta, 10), position)
  expect_identical(expected, seq_len(nrow(table)))
}

# Rows of `table` sorted by node_id, for order-insensitive comparisons.
sort_by_node = function(table) {
  table[order(table$node_id, method = "radix"), ]
}

# Annotation dataset joined to each toy node type.
toy_annotation = function(type) {
  switch(type,
    Anatomy = tkoi::anatomy_annotation,
    BiologicalProcess = tkoi::go_annotation,
    CellType = tkoi::celltype_annotation,
    Compound = tkoi::compound_annotation,
    Disease = tkoi::disease_annotation,
    Gene = {
      genes = as.data.frame(tkoi::genes)
      genes$identifier = as.character(genes$identifier)
      genes[, setdiff(names(genes), c("id", "degree"))]
    },
    MolecularFunction = tkoi::go_annotation,
    Pathway = tkoi::pathway_annotation
  )
}

# Arguments that toy_run(...) passes to .tkoi_plan_memory(). The run stops
# right after planning, before any sampling or PageRank.
toy_memory_request = function(...) {
  request = NULL
  with_mocked_bindings(
    expect_error(toy_run(...), "memory request recorded"),
    .tkoi_plan_memory = function(...) {
      request <<- list(...)
      stop("memory request recorded")
    },
    .package = "tkoi"
  )
  request
}

# Smallest memory limit, in bytes, at which .tkoi_plan_memory() accepts
# `request` (found by bisection, so the test does not depend on the formula).
minimum_memory = function(request) {
  plan = get(".tkoi_plan_memory", envir = asNamespace("tkoi"))
  fits = function(limit) {
    with_mocked_bindings(
      tryCatch(
        {
          do.call(plan, request)
          TRUE
        },
        error = function(e) FALSE
      ),
      .tkoi_memory_limit = function() limit,
      .package = "tkoi"
    )
  }
  low = 0
  high = 1e12
  expect_true(fits(high))
  for (step in seq_len(60)) {
    middle = (low + high) / 2
    if (fits(middle)) {
      high = middle
    } else {
      low = middle
    }
  }
  high
}

# Result object and parameters -----------------------------------------------

test_that("run_tkoi() returns a valid tKOIList that stores its parameters", {
  expression = toy_expression()
  result = withr::with_seed(1, run_tkoi(
    expression_data = expression,
    subnetwork = toy_network(),
    pvalue_threshold = 0.04,
    logfc_threshold = 0.3,
    topology_similarity = 0.85,
    n_permutation = 5,
    damping_factor = 0.8,
    maximum_iteration = 300,
    n_cores = 1,
    verbose = FALSE
  ))

  expect_s4_class(result, "tKOIList")
  expect_true(methods::validObject(result, test = TRUE))
  expect_identical(result@expression_data, expression)
  expect_identical(result@pvalue_threshold, 0.04)
  expect_identical(result@logfc_threshold, 0.3)
  expect_identical(result@topology_similarity, 0.85)
  expect_identical(result@n_permutation, 5)
  expect_identical(result@damping_factor, 0.8)
  expect_identical(result@maximum_iteration, 300)
  expect_null(result@subnetwork)
  expect_identical(result@gene_enrichment_comparison, list())
  expect_identical(names(result@pagerank_data), c("node_id", "pagerank", paste0("perm.", 1:5)))
  expect_output(methods::show(result), "tKOIList")
})

test_that("the default solver tolerance is 1e-14", {
  expect_identical(formals(run_tkoi)$tolerance, 1e-14)
})

test_that("pagerank_data has one row per vertex in network order and columns that sum to one", {
  result = toy_result()
  pagerank_data = result@pagerank_data
  vertex_names = igraph::V(toy_network())$name

  expect_s3_class(pagerank_data, "data.frame")
  expect_identical(names(pagerank_data), c("node_id", "pagerank", paste0("perm.", 1:10)))
  expect_identical(nrow(pagerank_data), length(vertex_names))
  expect_identical(pagerank_data$node_id, vertex_names)
  expect_identical(rownames(pagerank_data), vertex_names)

  vectors = as.matrix(pagerank_data[, -1])
  expect_true(is.double(vectors))
  expect_true(all(is.finite(vectors) & vectors >= 0))
  expect_per_node(colSums(vectors), rep(1, 11), 1e-10, label = "column sums")
})

test_that("keep_permutations = FALSE keeps only the null mean and sd, with identical statistics", {
  full = toy_result()
  compact = toy_run(keep_permutations = FALSE)
  vertex_names = igraph::V(toy_network())$name

  expect_true(methods::validObject(compact, test = TRUE))
  expect_identical(names(compact@pagerank_data), c("node_id", "pagerank", "null_mean", "null_sd"))
  expect_identical(compact@pagerank_data$node_id, vertex_names)
  expect_identical(rownames(compact@pagerank_data), vertex_names)
  expect_identical(compact@pagerank_data$pagerank, full@pagerank_data$pagerank)
  expect_identical(compact@network_summary_statistics, full@network_summary_statistics)

  permutations = unname(as.matrix(full@pagerank_data[, grep("^perm\\.", names(full@pagerank_data))]))
  expect_per_node(compact@pagerank_data$null_mean, rowMeans(permutations), 1e-12, floor = 0, label = "null_mean")
  expect_per_node(compact@pagerank_data$null_sd, apply(permutations, 1, stats::sd), 1e-10, floor = 0, label = "null_sd")
})

test_that("set.seed() makes run_tkoi() reproducible and another seed changes the null", {
  first = toy_result()
  again = toy_run(seed = 42)
  other = toy_run(seed = 7)

  expect_identical(again, first)
  expect_identical(other@pagerank_data$pagerank, first@pagerank_data$pagerank)
  expect_false(identical(other@pagerank_data$perm.1, first@pagerank_data$perm.1))
  expect_false(isTRUE(all.equal(
    other@network_summary_statistics$Gene$beta,
    first@network_summary_statistics$Gene$beta
  )))
})

test_that("results are identical for n_cores = 1 and n_cores = 2", {
  skip_if(.tkoi_resolve_cores(2) < 2, "fewer than two cores are available")
  expect_identical(toy_run(n_cores = 2), toy_result())
})

test_that("tolerance controls how precisely PageRank is solved", {
  default = toy_result()
  loose = toy_run(tolerance = 1e-4)
  difference = max(abs(loose@pagerank_data$pagerank - default@pagerank_data$pagerank))

  expect_gt(difference, 1e-12)
  expect_lt(difference, 1e-3)
})

# Parity with the released algorithm -----------------------------------------

test_that("run_tkoi() matches the released algorithm node by node", {
  expect_matches_legacy(toy_result(), toy_legacy_result())
})

test_that("run_tkoi() matches the released algorithm with non-default parameters", {
  expression = toy_expression()
  network = toy_network()
  result = withr::with_seed(9, run_tkoi(
    expression,
    network,
    pvalue_threshold = 0.02,
    logfc_threshold = 0.5,
    indirect_link_threshold = 1,
    topology_similarity = 0.7,
    n_permutation = 4,
    damping_factor = 0.6,
    n_cores = 1,
    verbose = FALSE
  ))
  legacy = withr::with_seed(9, legacy_run_tkoi(
    expression,
    network,
    pvalue_threshold = 0.02,
    logfc_threshold = 0.5,
    indirect_link_threshold = 1,
    topology_similarity = 0.7,
    n_permutation = 4,
    damping_factor = 0.6,
    gene_universe = toy_universe()
  ))

  expect_matches_legacy(result, legacy, indirect_link_threshold = 1)
  for (table in result@network_summary_statistics) {
    expect_tkoi_row_order(table, threshold = 1)
  }
  # The parameters matter: this is not the default run.
  expect_identical(names(result@pagerank_data), c("node_id", "pagerank", paste0("perm.", 1:4)))
  expect_false(isTRUE(all.equal(result@pagerank_data$pagerank, toy_result()@pagerank_data$pagerank)))
})

# Table contents --------------------------------------------------------------

test_that("each table holds every reported node of one type, with its annotation", {
  result = toy_result()
  network = toy_network()
  tables = result@network_summary_statistics
  node_type = clean_node_type(network)
  vertex_names = igraph::V(network)$name
  identifier = igraph::V(network)$identifier

  expect_identical(names(tables), sort(unique(node_type), method = "radix"))
  for (type in names(tables)) {
    table = tables[[type]]
    annotation = toy_annotation(type)
    annotation = annotation[!duplicated(annotation$identifier), ]
    annotation_columns = setdiff(names(annotation), "identifier")
    expect_s3_class(table, "tbl_df")
    expect_identical(names(table), c(statistic_columns, annotation_columns))
    expect_true(all(table$node_type == type), label = paste(type, "node_type"))

    rows = match(table$node_id, vertex_names)
    expect_false(anyNA(rows))
    expect_false(anyDuplicated(table$node_id) > 0, label = paste(type, "one row per node"))
    expect_identical(node_type[rows], table$node_type)
    expect_identical(table$identifier, identifier[rows])
    expect_identical(table$valency, igraph::V(network)$degree[rows])
    expect_identical(table$pagerank, result@pagerank_data$pagerank[rows])

    # Every node of the type is reported (compounds: only human metabolites).
    expected = which(node_type == type)
    if (type == "Compound") {
      expected = expected[identifier[expected] %in% tkoi::human_metabolites]
    }
    expect_setequal(table$node_id, vertex_names[expected])

    # Annotation columns come from the annotation dataset, or are missing.
    annotation_rows = match(table$identifier, annotation$identifier)
    for (column in annotation_columns) {
      expect_identical(table[[column]], annotation[[column]][annotation_rows], label = paste(type, column))
    }
  }

  expect_true(all(tables$Compound$identifier %in% tkoi::human_metabolites))
  expect_true(all(c("name", "namespace", "definition") %in% names(tables$BiologicalProcess)))
  expect_true(all(c("ensembl", "name") %in% names(tables$Gene)))
  # Genes without a curated annotation are kept, with missing annotation.
  expect_identical(nrow(tables$Gene), sum(node_type == "Gene"))
  expect_true(any(is.na(tables$Gene$name)))
  # Compounds that are not human metabolites are not reported.
  expect_lt(nrow(tables$Compound), sum(node_type == "Compound"))
})

test_that("beta, p_value and fdr follow from pagerank_data", {
  result = toy_result()
  network = toy_network()
  node_type = clean_node_type(network)
  derived = compute_network_enrichment(result@pagerank_data)
  # Constant nulls are untestable.
  constant = constant_null(result)
  expect_true(all(is.nan(derived$beta[constant])))
  # BH adjustment within each node type, over the reported nodes of that type.
  reported = !(node_type == "Compound" & !igraph::V(network)$identifier %in% tkoi::human_metabolites)
  derived$fdr = NA_real_
  for (type in unique(node_type)) {
    rows = which(node_type == type & reported)
    derived$fdr[rows] = stats::p.adjust(derived$p_value[rows], method = "fdr")
  }

  for (type in names(result@network_summary_statistics)) {
    table = result@network_summary_statistics[[type]]
    rows = match(table$node_id, result@pagerank_data$node_id)
    for (column in c("beta", "p_value", "fdr")) {
      expect_per_node(table[[column]], derived[[column]][rows], 1e-8, label = paste(type, column))
    }
  }
})

test_that("a seed gene that no replacement gene reaches is untestable, not infinitely significant", {
  network = toy_network()
  seeds = toy_seeds()
  # An isolated two-vertex component: a seed gene and a new Protein node.
  seed = seeds$id[1]
  seed_row = match(seed, igraph::V(network)$name)
  neighbours = igraph::neighbors(network, seed_row, mode = "all")
  isolated = igraph::delete_edges(network, igraph::incident(network, seed_row, mode = "all"))
  isolated = igraph::add_vertices(
    isolated, 1,
    name = "isolated-protein", identifier = "isolated-protein", labels = "['Protein']", degree = 1
  )
  isolated = igraph::add_edges(isolated, c(seed, "isolated-protein"))
  result = withr::with_seed(3, run_tkoi(toy_expression(), isolated, n_permutation = 5, n_cores = 1, verbose = FALSE))

  null = as.matrix(result@pagerank_data[, grep("^perm\\.", names(result@pagerank_data))])
  rows = match(c(seed, "isolated-protein"), result@pagerank_data$node_id)
  expect_true(all(result@pagerank_data$pagerank[rows] > 0))
  expect_true(all(null[rows, ] == 0))
  statistics = do.call(rbind, lapply(result@network_summary_statistics, function(x) {
    x[x$node_id %in% c(seed, "isolated-protein"), c("node_id", "beta", "p_value", "fdr")]
  }))
  expect_identical(nrow(statistics), 2L)
  expect_true(all(is.nan(statistics$beta)))
  expect_true(all(is.nan(statistics$p_value)))
  expect_true(all(is.na(statistics$fdr)))
  # Untestable nodes are listed after every testable node of their type.
  gene_table = result@network_summary_statistics$Gene
  expect_true(all(is.na(utils::tail(gene_table$beta, sum(is.na(gene_table$beta))))))
  expect_gt(length(neighbours), 0)
})

test_that("rows are ordered by the link threshold, then fdr, then beta", {
  for (table in toy_result()@network_summary_statistics) {
    expect_tkoi_row_order(table, threshold = 3)
  }
})

test_that("direct_links and indirect_links count seeds within one and two hops", {
  network = toy_network()
  seeds = toy_seeds()$id
  seeds = seeds[seeds %in% igraph::V(network)$name]
  count_reach = function(order) {
    table(unlist(lapply(igraph::ego(network, order = order, nodes = seeds, mode = "all"), names)))
  }
  direct = count_reach(1)
  indirect = count_reach(2)

  for (table in toy_result()@network_summary_statistics) {
    expect_identical(table$direct_links, unname(as.numeric(direct[table$node_id])))
    expect_identical(table$indirect_links, unname(as.numeric(indirect[table$node_id])))
  }
  all_links = dplyr::bind_rows(lapply(toy_result()@network_summary_statistics, function(x) x[, statistic_columns]))
  expect_false(any(all_links$direct_links == 0, na.rm = TRUE))
  expect_false(any(all_links$indirect_links == 0, na.rm = TRUE))
  expect_true(anyNA(all_links$indirect_links))
})

test_that("indirect_link_threshold only changes the row order", {
  default = toy_result()
  for (threshold in c(1, 1e6)) {
    result = toy_run(indirect_link_threshold = threshold)
    expect_identical(result@pagerank_data, default@pagerank_data)
    expect_identical(names(result@network_summary_statistics), names(default@network_summary_statistics))
    for (type in names(default@network_summary_statistics)) {
      table = result@network_summary_statistics[[type]]
      expect_identical(sort_by_node(table), sort_by_node(default@network_summary_statistics[[type]]))
      expect_tkoi_row_order(table, threshold = threshold)
    }
    reordered = vapply(names(default@network_summary_statistics), function(type) {
      !identical(result@network_summary_statistics[[type]]$node_id, default@network_summary_statistics[[type]]$node_id)
    }, logical(1))
    expect_true(any(reordered), label = paste("threshold", threshold, "reorders some table"))
  }
})

test_that("nodes that no seed can reach get NaN statistics, sort last, and have no link counts", {
  network = toy_network()
  # Four Pathway vertices in a ring, disconnected from the rest of the network.
  # Their identifiers have exactly one annotation row, so each survives the
  # annotation join once.
  pathway_ids = tkoi::pathway_annotation$identifier
  single_row = pathway_ids[!pathway_ids %in% pathway_ids[duplicated(pathway_ids)]]
  ring_ids = utils::head(setdiff(single_row, igraph::V(network)$identifier), 4)
  ring = paste0("unreachable:", 1:4)
  ringed = network |>
    igraph::add_vertices(4, name = ring, identifier = ring_ids, labels = "['Pathway']", degree = 2) |>
    igraph::add_edges(c(rbind(ring, ring[c(2:4, 1)])))
  expect_length(ring_ids, 4)
  expect_identical(igraph::components(ringed)$no, 2)

  result = quiet_run(subnetwork = ringed)
  default = toy_result()

  # No seed or null seed reaches the ring, so its PageRank is zero in every run.
  ring_pagerank = as.matrix(result@pagerank_data[result@pagerank_data$node_id %in% ring, -1])
  expect_identical(dim(ring_pagerank), c(4L, 11L))
  expect_true(all(ring_pagerank == 0))
  expect_per_node(colSums(as.matrix(result@pagerank_data[, -1])), rep(1, 11), 1e-10, label = "column sums")

  pathway = result@network_summary_statistics$Pathway
  expect_identical(utils::tail(pathway$node_id, 4), ring)
  ring_rows = utils::tail(pathway, 4)
  expect_identical(ring_rows$identifier, ring_ids)
  expect_identical(ring_rows$pagerank, rep(0, 4))
  expect_true(all(is.nan(ring_rows$beta)))
  expect_true(all(is.nan(ring_rows$p_value)))
  expect_true(all(is.nan(ring_rows$fdr)))
  expect_identical(ring_rows$direct_links, rep(NA_real_, 4))
  expect_identical(ring_rows$indirect_links, rep(NA_real_, 4))
  expect_identical(ring_rows$valency, rep(2, 4))
  expect_tkoi_row_order(pathway, threshold = 3, network = ringed)

  # Every other node keeps its statistics: the NaN nodes do not enter the FDR.
  expect_identical(names(result@network_summary_statistics), names(default@network_summary_statistics))
  for (type in names(default@network_summary_statistics)) {
    table = result@network_summary_statistics[[type]]
    table = table[!table$node_id %in% ring, ]
    expected = default@network_summary_statistics[[type]]
    expect_identical(sort_by_node(table)$node_id, sort_by_node(expected)$node_id)
    aligned = table[match(expected$node_id, table$node_id), ]
    expect_false(anyNA(aligned$beta), label = paste(type, "beta"))
    for (column in c("pagerank", "beta", "p_value", "fdr")) {
      expect_per_node(aligned[[column]], expected[[column]], 1e-10, label = paste(type, column))
    }
  }
})

# Node types --------------------------------------------------------------------

test_that("vertex attributes are read afresh: new types stay unannotated and missing types become Unknown", {
  # The toy network is prepared (and its matrix cached) first; the edited copy
  # below keeps the same igraph::graph_id().
  default = toy_result()
  network = toy_network()
  labels = igraph::V(network)$labels
  edited = network
  igraph::V(edited)$labels[labels == "['Anatomy']"] = "['Tissue']"
  igraph::V(edited)$labels[labels == "['CellType']"] = NA
  igraph::V(edited)$degree = igraph::V(network)$degree + 1
  expect_identical(igraph::graph_id(edited), igraph::graph_id(network))

  result = quiet_run(subnetwork = edited)
  tables = result@network_summary_statistics

  expect_identical(result@pagerank_data, default@pagerank_data)
  expected_types = c(setdiff(names(default@network_summary_statistics), c("Anatomy", "CellType")), "Tissue", "Unknown")
  expect_identical(names(tables), sort(expected_types, method = "radix"))

  # Types without a curated annotation keep every node and only the statistics.
  renamed = list(Tissue = "['Anatomy']", Unknown = "['CellType']")
  for (type in names(renamed)) {
    table = tables[[type]]
    expect_s3_class(table, "tbl_df")
    expect_identical(names(table), statistic_columns)
    expect_true(all(table$node_type == type), label = paste(type, "node_type"))
    expect_setequal(table$node_id, igraph::V(network)$name[labels == renamed[[type]]])
    expect_identical(nrow(table), sum(labels == renamed[[type]]))
    expect_tkoi_row_order(table, threshold = 3, network = edited)
  }

  # The renamed nodes keep their statistics (FDR is still adjusted over the same
  # nodes), and every valency comes from the edited degree attribute.
  original = list(Tissue = "Anatomy", Unknown = "CellType")
  for (type in names(original)) {
    expected = default@network_summary_statistics[[original[[type]]]]
    aligned = tables[[type]][match(expected$node_id, tables[[type]]$node_id), ]
    for (column in c("node_id", "identifier", "pagerank", "beta", "p_value", "fdr", "direct_links", "indirect_links")) {
      expect_identical(aligned[[column]], expected[[column]], label = paste(type, column))
    }
  }
  for (table in tables) {
    rows = match(table$node_id, igraph::V(network)$name)
    expect_identical(table$valency, igraph::V(network)$degree[rows] + 1)
  }
  for (type in setdiff(names(tables), c("Tissue", "Unknown"))) {
    unchanged = setdiff(names(tables[[type]]), "valency")
    expect_identical(tables[[type]][unchanged], default@network_summary_statistics[[type]][unchanged])
  }
})

test_that("node types may be stored in a `label` vertex attribute", {
  network = toy_network()
  relabelled = igraph::delete_vertex_attr(network, "labels")
  igraph::V(relabelled)$label = igraph::V(network)$labels

  result = quiet_run(subnetwork = relabelled)
  expect_identical(result@pagerank_data, toy_result()@pagerank_data)
  expect_identical(result@network_summary_statistics, toy_result()@network_summary_statistics)
})

# Input handling --------------------------------------------------------------

test_that("tibble input and a factor gene_name give the same statistics as a data frame", {
  expression = toy_expression()
  reference = quiet_run(expression, n_permutation = 3)
  factor_genes = expression
  factor_genes$gene_name = factor(factor_genes$gene_name)
  inputs = list(tibble = dplyr::as_tibble(expression), factor = factor_genes)

  for (name in names(inputs)) {
    result = quiet_run(inputs[[name]], n_permutation = 3)
    expect_identical(result@expression_data, inputs[[name]], label = name)
    expect_identical(result@pagerank_data, reference@pagerank_data, label = name)
    expect_identical(result@network_summary_statistics, reference@network_summary_statistics, label = name)
  }
})

test_that("data.table input gives the same statistics as a data frame", {
  skip_if_not_installed("data.table")
  expression = toy_expression()
  table = data.table::as.data.table(expression)

  result = quiet_run(table, n_permutation = 3)
  reference = quiet_run(expression, n_permutation = 3)
  # The data.frame slot of tKOIList stores a data.table as a plain data frame
  # (data.table registers that coercion with methods::setAs()).
  expect_identical(result@expression_data, as.data.frame(table))
  expect_identical(result@pagerank_data, reference@pagerank_data)
  expect_identical(result@network_summary_statistics, reference@network_summary_statistics)
})

test_that("blank and missing gene names are ignored", {
  expression = toy_expression()
  padded = rbind(
    data.frame(gene_name = c("", NA), logfc = c(3, -3), pvalue = c(1e-8, 1e-8)),
    expression,
    data.frame(gene_name = c(NA, ""), logfc = c(-2, 2), pvalue = c(1e-6, 1e-6))
  )
  result = quiet_run(padded)
  default = toy_result()

  expect_identical(result@expression_data, padded)
  expect_identical(result@pagerank_data, default@pagerank_data)
  expect_identical(result@network_summary_statistics, default@network_summary_statistics)
})

test_that("the first row of a duplicated gene wins and unknown Ensembl IDs are skipped", {
  expression = toy_expression()
  n = nrow(expression)
  # The last three rows are a duplicate of row 1 and two unknown IDs.
  expect_identical(expression$gene_name[n - 2], expression$gene_name[1])
  cleaned = expression[seq_len(n - 3), ]
  cleaned_result = quiet_run(cleaned)
  expect_identical(cleaned_result@pagerank_data, toy_result()@pagerank_data)
  expect_identical(cleaned_result@network_summary_statistics, toy_result()@network_summary_statistics)

  # With the duplicate listed first, its values are used instead.
  duplicate_first = expression[c(n - 2, 2:(n - 3), 1, n - 1, n), ]
  replaced = cleaned
  replaced[1, ] = expression[n - 2, ]
  duplicate_first_result = quiet_run(duplicate_first, seed = 3, n_permutation = 3)
  replaced_result = quiet_run(replaced, seed = 3, n_permutation = 3)
  expect_identical(duplicate_first_result@pagerank_data, replaced_result@pagerank_data)
  expect_identical(duplicate_first_result@network_summary_statistics, replaced_result@network_summary_statistics)
  expect_false(isTRUE(all.equal(duplicate_first_result@pagerank_data$pagerank, toy_result()@pagerank_data$pagerank)))
})

test_that("run_tkoi() validates its arguments", {
  run_with = function(expression_data = toy_expression(), subnetwork = toy_network(), n_permutation = 2,
                      n_cores = 1, verbose = FALSE, ...) {
    run_tkoi(
      expression_data,
      subnetwork = subnetwork,
      n_permutation = n_permutation,
      n_cores = n_cores,
      verbose = verbose,
      ...
    )
  }
  expression = toy_expression()

  expect_error(run_with(as.matrix(expression)), "must be a data frame")
  expect_error(run_with("ENSG00000112200"), "must be a data frame")
  expect_error(run_with(expression[, c("gene_name", "logfc")]), "missing column\\(s\\): pvalue")
  expect_error(run_with(expression[, c("logfc", "pvalue")]), "missing column\\(s\\): gene_name")
  expect_error(run_with(data.frame(gene = expression$gene_name)), "gene_name, logfc, pvalue")
  expect_error(run_with(transform(expression, logfc = as.character(logfc))), "must be numeric")
  expect_error(run_with(transform(expression, pvalue = as.character(pvalue))), "must be numeric")

  for (value in list(-0.1, 1.5, NA_real_, c(0.01, 0.05), "0.05")) {
    expect_error(run_with(pvalue_threshold = value), "`pvalue_threshold` must be a number between 0 and 1")
  }
  for (value in list(-1, Inf, NA_real_, "1")) {
    expect_error(run_with(logfc_threshold = value), "`logfc_threshold` must be a non-negative number")
  }
  for (value in list("3", NA_real_, c(1, 2))) {
    expect_error(run_with(indirect_link_threshold = value), "`indirect_link_threshold` must be a number")
  }
  for (value in list(-0.1, 1.1, NA_real_)) {
    expect_error(run_with(topology_similarity = value), "`topology_similarity` must be a number between 0 and 1")
  }
  # 3e9 is a whole number but does not fit in an integer.
  for (value in list(1, 0, 2.5, NA_real_, "10", 3e9)) {
    expect_error(run_with(n_permutation = value), "`n_permutation` must be a whole number of at least 2")
  }
  for (value in list(0, 1, 1.2, -0.5, NA_real_)) {
    expect_error(run_with(damping_factor = value), "`damping_factor` must be strictly between 0 and 1")
  }
  for (value in list(0, 1.5, -3, NA_real_, 3e9)) {
    expect_error(
      run_with(maximum_iteration = value),
      "`maximum_iteration` must be a whole number between 1 and .Machine$integer.max",
      fixed = TRUE
    )
  }
  for (value in list(0, 1, -1e-12, NA_real_)) {
    expect_error(run_with(tolerance = value), "`tolerance` must be a number strictly between 0 and 1")
  }
  for (value in list(NA, "TRUE", c(TRUE, FALSE), 1)) {
    expect_error(run_with(keep_permutations = value), "`keep_permutations` must be TRUE or FALSE")
  }
  for (value in list(NA, "yes", 1)) {
    expect_error(run_with(verbose = value), "`verbose` must be TRUE or FALSE")
  }
  for (value in list(0, 1.5, "2", c(1, 2))) {
    expect_error(run_with(n_cores = value), "`n_cores` must be NULL or a single positive whole number")
  }

  expect_error(run_with(logfc_threshold = 100), "No genes pass the thresholds")
  expect_error(run_with(expression[0, ]), "No genes pass the thresholds")
  expect_error(run_with(utils::tail(expression, 2)), "No genes pass the thresholds")
  infinite = expression
  infinite$logfc[2] = Inf
  infinite$pvalue[2] = 1e-5
  expect_error(run_with(infinite), "infinite `logfc`")
  zero = expression
  zero$logfc = 0
  expect_error(run_with(zero, logfc_threshold = 0), "`logfc` equal to zero")

  expect_error(run_with(subnetwork = as.raw(1:16)), "Encrypted networks are no longer supported")
  expect_error(run_with(subnetwork = data.frame(from = "a", to = "b")), "must be an igraph object")
})

test_that("maximum_iteration = 1 warns that PageRank did not converge", {
  expect_warning(
    toy_run(maximum_iteration = 1, n_permutation = 2),
    "PageRank vector\\(s\\) did not reach `tolerance` within `maximum_iteration` = 1 iterations"
  )
})

test_that("verbose = TRUE reports progress and verbose = FALSE is silent", {
  output = utils::capture.output({
    messages = testthat::capture_messages({
      verbose_result = withr::with_seed(42, run_tkoi(
        toy_expression(),
        toy_network(),
        n_permutation = 2,
        n_cores = 1,
        verbose = TRUE
      ))
    })
  })
  messages = paste(messages, collapse = "")
  expect_match(messages, "Preparing the network")
  expect_match(messages, "Sampling 2 degree-matched null seed sets")
  expect_match(messages, "Running personalized PageRank for the observed data and 2 permutations on 1 thread")
  expect_match(messages, "Calculating network enrichment statistics")
  expect_match(messages, "Annotating tKOI analysis results")
  expect_match(messages, "Analysis finished and took")
  expect_match(paste(output, collapse = ""), "100%", fixed = TRUE)

  # verbose = FALSE prints no messages and no progress bar.
  expect_no_message({
    quiet_output = utils::capture.output({
      quiet_result = toy_run(n_permutation = 2)
    })
  })
  expect_identical(quiet_output, character())
  expect_identical(verbose_result, quiet_result)
})

# Memory planning -------------------------------------------------------------

test_that("run_tkoi() checks memory before sampling the null or running PageRank", {
  refuse = function(...) stop("called although the memory check failed")
  local_mocked_bindings(
    # 1 MB: less than the toy network matrix alone.
    .tkoi_memory_limit = function() 1e6,
    .tkoi_sample_null = refuse,
    .tkoi_ppr_null = refuse,
    .package = "tkoi"
  )

  expect_error(
    toy_run(),
    paste0(
      "^run_tkoi\\(\\) needs about [0-9.]+ GB but only 0\\.0 GB of memory is available\\. ",
      "Set `keep_permutations = FALSE` or reduce `n_permutation`\\.$"
    )
  )
  expect_error(
    toy_run(keep_permutations = FALSE),
    paste0(
      "^run_tkoi\\(\\) needs about [0-9.]+ GB but only 0\\.0 GB of memory is available\\. ",
      "Reduce `n_permutation` or use a smaller network\\.$"
    )
  )
})

test_that("keep_permutations = FALSE fits in memory where keeping every permutation does not", {
  # What run_tkoi() asks the planner for: the seeds are the in-network genes.
  request = toy_memory_request(n_permutation = 12)
  expect_identical(request$n_vertices, igraph::vcount(toy_network()))
  expect_identical(request$n_seeds, nrow(toy_seeds()))
  expect_identical(request$n_permutation, 12)
  expect_true(request$keep_permutations)

  keep_limit = minimum_memory(request)
  compact_limit = minimum_memory(utils::modifyList(request, list(keep_permutations = FALSE)))
  expect_lt(compact_limit, keep_limit)

  reference = toy_run(n_permutation = 12, keep_permutations = FALSE)
  original_ppr_null = get(".tkoi_ppr_null", envir = asNamespace("tkoi"))
  pagerank_calls = 0
  limit = (compact_limit + keep_limit) / 2
  local_mocked_bindings(
    .tkoi_memory_limit = function() limit,
    .tkoi_ppr_null = function(...) {
      pagerank_calls <<- pagerank_calls + 1
      original_ppr_null(...)
    },
    .package = "tkoi"
  )

  expect_error(toy_run(n_permutation = 12), "Set `keep_permutations = FALSE` or reduce `n_permutation`")
  expect_identical(pagerank_calls, 0)
  compact = toy_run(n_permutation = 12, keep_permutations = FALSE)
  expect_identical(pagerank_calls, 1)
  expect_identical(compact, reference)
})

test_that("results are identical when PageRank vectors are solved one per block", {
  default = toy_result()
  original_ppr_null = get(".tkoi_ppr_null", envir = asNamespace("tkoi"))
  block_widths = integer()
  local_mocked_bindings(
    .tkoi_plan_memory = function(...) 1L,
    .tkoi_ppr_null = function(csr, observed_node, observed_weight, null_node, null_weight, damping, tol, max_iter,
                              n_threads, max_block, ...) {
      block_widths <<- c(block_widths, max_block)
      original_ppr_null(
        csr, observed_node, observed_weight, null_node, null_weight, damping, tol, max_iter, n_threads, max_block, ...
      )
    },
    .package = "tkoi"
  )

  expect_identical(toy_run(), default)
  expect_identical(block_widths, 1L)
})

# Subnetworks -----------------------------------------------------------------

test_that("a weighted subnetwork gives igraph's weighted personalized PageRank", {
  # A fresh copy of the toy network (new graph id) with random edge weights.
  weighted = igraph::induced_subgraph(toy_network(), igraph::V(toy_network()))
  igraph::E(weighted)$weight = withr::with_seed(1, stats::runif(igraph::ecount(weighted), 0.1, 3))

  result = quiet_run(subnetwork = weighted, n_permutation = 2)
  reference = igraph_toy_pagerank(weighted)

  expect_identical(result@pagerank_data$node_id, igraph::V(toy_network())$name)
  expect_per_node(result@pagerank_data$pagerank, reference, 1e-11, label = "pagerank")
  expect_per_node(result@pagerank_data$pagerank, reference, 1e-8, floor = 0, label = "pagerank (relative)")
  expect_false(isTRUE(all.equal(result@pagerank_data$pagerank, toy_result()@pagerank_data$pagerank)))
})

test_that("adding a weight attribute to a network already analysed changes PageRank", {
  network = toy_network()
  # The first run prepares (and caches) the unweighted network.
  unweighted = quiet_run(subnetwork = network, n_permutation = 2)
  # igraph keeps the graph id when attributes change, so the cached matrix must
  # not be reused for the weighted network.
  igraph::E(network)$weight = withr::with_seed(1, stats::runif(igraph::ecount(network), 0.1, 3))
  expect_identical(igraph::graph_id(network), igraph::graph_id(toy_network()))

  weighted = quiet_run(subnetwork = network, n_permutation = 2)
  reference = igraph_toy_pagerank(network)

  expect_false(isTRUE(all.equal(weighted@pagerank_data$pagerank, unweighted@pagerank_data$pagerank)))
  expect_per_node(weighted@pagerank_data$pagerank, reference, 1e-11, label = "pagerank")
  expect_per_node(weighted@pagerank_data$pagerank, reference, 1e-8, floor = 0, label = "pagerank (relative)")
  # Going back to the unweighted network gives the unweighted result again.
  expect_identical(quiet_run(subnetwork = toy_network(), n_permutation = 2), unweighted)
})

test_that("a subnetwork lacking some node types and genes runs and matches the released algorithm", {
  network = toy_network()
  labels = igraph::V(network)$labels
  seeds = toy_seeds()$id
  genes = igraph::V(network)$name[labels == "['Gene']"]
  non_seed_genes = setdiff(genes, seeds)
  dropped_genes = non_seed_genes[seq(1, length(non_seed_genes), by = 2)]
  keep = labels != "['Anatomy']" & !igraph::V(network)$name %in% dropped_genes
  reduced = igraph::induced_subgraph(network, which(keep))

  result = quiet_run(subnetwork = reduced, seed = 5, n_permutation = 3)
  expect_true(methods::validObject(result, test = TRUE))
  tables = result@network_summary_statistics
  expect_false("Anatomy" %in% names(tables))
  expect_identical(
    names(tables),
    c("BiologicalProcess", "CellType", "Compound", "Disease", "Gene", "MolecularFunction", "Pathway")
  )
  expect_identical(result@pagerank_data$node_id, igraph::V(reduced)$name)
  expect_equal(sum(result@pagerank_data$pagerank), 1, tolerance = 1e-10)

  # Null genes are drawn only from genes in the subnetwork.
  legacy = withr::with_seed(5, legacy_run_tkoi(
    toy_expression(),
    reduced,
    n_permutation = 3,
    gene_universe = toy_universe(reduced)
  ))
  expect_matches_legacy(result, legacy)
})

test_that("seed genes missing from the subnetwork are dropped before the null is drawn", {
  network = toy_network()
  expression = toy_expression()
  seeds = toy_seeds(expression)
  # Remove three seed genes whose Ensembl ID maps to a single tkoi::genes entry.
  genes = as.data.frame(tkoi::genes)
  single = !seeds$gene_name %in% genes$ensembl[duplicated(genes$ensembl)]
  removed = utils::head(seeds$id[single], 3)
  expect_length(removed, 3)
  reduced = igraph::delete_vertices(network, removed)
  n_kept = nrow(seeds) - 3L

  # Record the seeds of every PageRank run, then compute them as usual.
  original_ppr_null = get(".tkoi_ppr_null", envir = asNamespace("tkoi"))
  calls = list()
  local_mocked_bindings(
    .tkoi_ppr_null = function(csr, observed_node, observed_weight, null_node, null_weight, ...) {
      calls[[length(calls) + 1]] <<- list(
        n_observed = length(observed_node),
        n_null = nrow(null_node),
        n_sets = ncol(null_node),
        observed_node = observed_node,
        observed_weight = observed_weight,
        null_node = null_node,
        null_weight = null_weight
      )
      original_ppr_null(csr, observed_node, observed_weight, null_node, null_weight, ...)
    },
    .package = "tkoi"
  )
  result = quiet_run(expression, reduced, seed = 11, n_permutation = 6)

  expect_length(calls, 1)
  call = calls[[1]]
  expect_identical(call$n_observed, n_kept)
  expect_identical(call$n_sets, 6L)
  # Every null seed set has as many seeds as the observed run, with the same weights.
  expect_identical(call$n_null, call$n_observed)
  expect_identical(call$null_weight, call$observed_weight)
  expect_false(anyNA(call$null_node))
  vertex_names = igraph::V(reduced)$name
  expect_setequal(vertex_names[call$observed_node], setdiff(seeds$id, removed))
  expect_true(all(vertex_names[call$null_node] %in% toy_universe(reduced)$id))

  # Same as the released algorithm run on the reduced network with only the
  # genes it contains.
  in_network = genes$ensembl[genes$id %in% vertex_names]
  restricted = expression[expression$gene_name %in% in_network, ]
  legacy = withr::with_seed(11, legacy_run_tkoi(
    restricted,
    reduced,
    n_permutation = 6,
    gene_universe = toy_universe(reduced)
  ))
  expect_matches_legacy(result, legacy)
})

test_that("a subnetwork without any selected gene is an error", {
  network = toy_network()
  no_genes = igraph::induced_subgraph(
    network,
    which(igraph::V(network)$labels %in% c("['BiologicalProcess']", "['Pathway']"))
  )
  expect_error(
    run_tkoi(toy_expression(), no_genes, n_permutation = 2, n_cores = 1, verbose = FALSE),
    "None of the selected genes are in `subnetwork`"
  )
})

test_that("an encrypted (raw vector) subnetwork is rejected", {
  expect_error(
    run_tkoi(toy_expression(), serialize(1:10, NULL), n_permutation = 2, n_cores = 1, verbose = FALSE),
    "Encrypted"
  )
})
