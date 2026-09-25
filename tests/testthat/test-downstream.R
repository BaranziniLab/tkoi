# Tests of the functions that consume a run_tkoi() result: the tKOIList show()
# method, export_gene_exploration_data(), make_gene_exploration_plot(),
# visualize_topn(), export_network_summary_statistics(),
# compute_network_enrichment(), run_gene_enrichment(), get_neighboring_nodes()
# and plot_network(). Everything runs on the toy network.

# Shared values ------------------------------------------------------------------
#
# Computed on first use and cached. A filtered run computes only what its tests
# need, and an error is reported by the tests that need the value instead of
# aborting the whole file.

lazily = function(compute) {
  cache = NULL
  function() {
    if (is.null(cache)) {
      cache <<- compute()
    }
    cache
  }
}

toy_result = lazily(function() toy_run())

gene_table = function() toy_result()@network_summary_statistics$Gene

# The input as the user passed it: one gene is listed twice, with different values.
input_expression = function() toy_result()@expression_data

# The rows run_tkoi() uses: no missing or blank gene names, first row per gene.
first_rows = function(expression = input_expression()) {
  named = !is.na(expression$gene_name) & nzchar(expression$gene_name)
  expression[named & !duplicated(expression$gene_name), ]
}

passes_filters = function(x, result = toy_result()) {
  x$pvalue <= result@pvalue_threshold & abs(x$logfc) >= result@logfc_threshold
}

# Toy-network vertex IDs of a set of Ensembl gene IDs.
toy_gene_ids = function(ensembl) {
  intersect(tkoi::genes$id[tkoi::genes$ensembl %in% ensembl], igraph::V(toy_network())$name)
}

# The seed genes of the run: significant genes (first row per gene) in the network.
seed_ids = lazily(function() {
  rows = first_rows()
  toy_gene_ids(rows$gene_name[passes_filters(rows)])
})

# Hops from each BiologicalProcess node to the nearest seed gene, named by node ID.
hops_to_seed = lazily(function() {
  process_ids = toy_result()@network_summary_statistics$BiologicalProcess$node_id
  distance = igraph::distances(toy_network(), v = seed_ids(), to = process_ids, mode = "all")
  stats::setNames(apply(distance, 2, min), process_ids)
})

processes_at = function(hops) names(which(hops_to_seed() == hops))

# One (node_id, beta, p_value) row per node over every node-type table.
all_node_statistics = function(tkoi_list) {
  statistics = lapply(tkoi_list@network_summary_statistics, function(x) {
    data.frame(node_id = x$node_id, beta = x$beta, p_value = x$p_value)
  })
  statistics = do.call(rbind, unname(statistics))
  statistics[!duplicated(statistics$node_id), ]
}

# Draw base-graphics plots on a null device so nothing is written to disk.
local_null_device = function(env = parent.frame()) {
  grDevices::pdf(NULL)
  device = grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device), envir = env)
}

is_ggplot = function(x) inherits(x, "ggplot")

# plot_network() helpers -------------------------------------------------------------

# Runs `code` and returns its value together with the arguments igraph's plot
# method received (the plotted graph as `x`, then the styling arguments).
capture_plot_arguments = function(code) {
  captured = new.env(parent = emptyenv())
  igraph_namespace = asNamespace("igraph")
  # trace() replaces plot.igraph() in the igraph namespace, but S3 dispatch
  # finds the method in the registry, which keeps its own binding (set lazily,
  # on first use). Register the traced function while `code` runs, then restore.
  registered = utils::getS3method("plot", "igraph")
  suppressMessages(trace(
    "plot.igraph",
    where = igraph_namespace,
    tracer = bquote(assign("arguments", c(list(x = x), list(...)), envir = .(captured))),
    print = FALSE
  ))
  on.exit(
    {
      suppressMessages(untrace("plot.igraph", where = igraph_namespace))
      registerS3method("plot", "igraph", registered, envir = baseenv())
    },
    add = TRUE
  )
  registerS3method("plot", "igraph", get("plot.igraph", envir = igraph_namespace), envir = baseenv())

  value = code
  list(value = value, arguments = captured$arguments)
}

# plot_network() replaces vertex names by display labels, so the tests plot a
# copy of the network that also stores each vertex's node ID.
tagged_network = lazily(function() {
  network = toy_network()
  igraph::set_vertex_attr(network, "node_id", value = igraph::V(network)$name)
})

# Node types of a plotted graph, read from `labels` or `label`.
plotted_types = function(graph) {
  attribute = intersect(c("labels", "label"), igraph::vertex_attr_names(graph))[1]
  gsub("[][']", "", igraph::vertex_attr(graph, attribute))
}

# Labels igraph draws: the vertex.label argument, else a `label` vertex
# attribute, else the vertex names.
drawn_labels = function(arguments) {
  if (!is.null(arguments$vertex.label)) {
    return(arguments$vertex.label)
  }
  if ("label" %in% igraph::vertex_attr_names(arguments$x)) {
    return(igraph::V(arguments$x)$label)
  }
  igraph::V(arguments$x)$name
}

# Log fold change of each node from the rows run_tkoi() uses; NA for nodes that
# are not measured genes.
input_logfc = function(node_ids, expression = input_expression()) {
  rows = first_rows(expression)
  ensembl = tkoi::genes$ensembl[match(node_ids, tkoi::genes$id)]
  rows$logfc[match(ensembl, rows$gene_name)]
}

# Colours of a blue-white-red ramp of 100 steps, centred on a log fold change
# of zero and scaled by the largest absolute value.
ramp_colour = function(logfc) {
  blue_white_red = grDevices::colorRampPalette(c("#79ADD6", "white", "#BB2A2D"))(100)
  blue_white_red[ceiling((logfc / max(abs(logfc)) + 1) / 2 * 99) + 1]
}

# A Disease node next to at least two seed genes. Within two hops, its plot has
# up- and down-regulated seed genes, measured genes that are not significant,
# unmeasured genes and nodes of other types.
colour_target = function() {
  diseases = toy_result()@network_summary_statistics$Disease
  diseases$node_id[!is.na(diseases$direct_links) & diseases$direct_links >= 2][1]
}

# Nodes on a path of at most `order` edges between `target` and a seed gene.
# For order <= 2 these paths are always simple, so this is exactly the node
# set plot_network() draws.
nodes_between = function(target, order) {
  network = toy_network()
  seeds = seed_ids()
  from_target = igraph::distances(network, v = target, mode = "all")[1, ]
  genes = setdiff(seeds[from_target[seeds] <= order], target)
  from_genes = igraph::distances(network, v = genes, mode = "all")
  names(which(from_target + apply(from_genes, 2, min) <= order))
}
identifiers = function(node_ids) igraph::V(toy_network())[node_ids]$identifier

# show() ----------------------------------------------------------------------

test_that("show() prints the header, the input size and the parameter table", {
  result = toy_result()
  expect_output(show(result), "[tKOIList object]", fixed = TRUE)
  expect_output(show(result), glue::glue("{nrow(input_expression())} genes included in the experiment."), fixed = TRUE)
  expect_output(show(result), "Total object size:", fixed = TRUE)
  expect_output(show(result), "Network enrichment analysis run with these parameters:", fixed = TRUE)

  printed = paste(utils::capture.output(show(result)), collapse = "\n")
  parameters = c(
    "P.Value Threshold" = result@pvalue_threshold,
    "Log Fold Change Threshold" = result@logfc_threshold,
    "Topology Similarity" = result@topology_similarity,
    "N Permutation" = result@n_permutation,
    "Damping Factor" = result@damping_factor,
    "Maximum Iteration" = result@maximum_iteration
  )
  for (parameter in names(parameters)) {
    row = regmatches(printed, regexpr(paste0("\\|", parameter, " *\\| *[0-9.]+\\|"), printed))
    expect_length(row, 1)
    value = as.numeric(sub(".*\\| *([0-9.]+)\\|$", "\\1", row))
    expect_equal(value, unname(parameters[parameter]))
  }

  # Printing at the console dispatches to the same method.
  expect_output(print(result), "[tKOIList object]", fixed = TRUE)
})

test_that("show() reports the parameters stored in the object", {
  custom = toy_result()
  custom@n_permutation = 7
  custom@damping_factor = 0.5
  expect_output(show(custom), "N Permutation\\s*\\|\\s*7(\\.0+)?\\|")
  expect_output(show(custom), "Damping Factor\\s*\\|\\s*0\\.50*\\|")

  empty = methods::new("tKOIList")
  expect_output(show(empty), "0 genes included in the experiment.", fixed = TRUE)
  expect_output(show(empty), "Maximum Iteration\\s*\\|\\s*500(\\.0+)?\\|")
})

# export_gene_exploration_data() ------------------------------------------------

test_that("export_gene_exploration_data() returns one row per Gene node with the documented columns", {
  exported = export_gene_exploration_data(toy_result())
  genes = gene_table()

  expect_s3_class(exported, "data.frame")
  expect_named(exported, c(
    "gene_name", "gene_symbol", "id", "identifier",
    "experimental_logfc", "experimental_pvalue",
    "pagerank", "tkoi_beta", "tkoi_pvalue", "tkoi_fdr"
  ))
  # The toy input repeats one gene; it must still give a single row per node.
  expect_equal(nrow(exported), nrow(genes))
  expect_false(anyDuplicated(exported$id) > 0)
  expect_setequal(exported$id, genes$node_id)
})

test_that("export_gene_exploration_data() carries the tKOI statistics of each Gene node", {
  exported = export_gene_exploration_data(toy_result())
  genes = gene_table()
  rows = match(exported$id, genes$node_id)

  expect_false(anyNA(rows))
  expect_identical(exported$identifier, genes$identifier[rows])
  expect_identical(exported$gene_symbol, genes$name[rows])
  expect_identical(exported$pagerank, genes$pagerank[rows])
  expect_identical(exported$tkoi_beta, genes$beta[rows])
  expect_identical(exported$tkoi_pvalue, genes$p_value[rows])
  expect_identical(exported$tkoi_fdr, genes$fdr[rows])
})

test_that("export_gene_exploration_data() reports the input values, and NA for genes not in the input", {
  exported = export_gene_exploration_data(toy_result())
  genes = gene_table()
  rows = first_rows()
  ensembl = genes$ensembl[match(exported$id, genes$node_id)]
  in_input = ensembl %in% rows$gene_name

  # Every gene of the input that is a Gene node appears, and only those have values.
  expect_setequal(ensembl[in_input], intersect(rows$gene_name, genes$ensembl))
  expect_identical(exported$gene_name[in_input], ensembl[in_input])
  expect_true(all(is.na(exported$gene_name[!in_input])))
  expect_true(all(is.na(exported$experimental_logfc[!in_input])))
  expect_true(all(is.na(exported$experimental_pvalue[!in_input])))
  expect_gt(sum(!in_input), 0)

  # Values come from the first input row of each gene.
  source_rows = match(exported$gene_name[in_input], rows$gene_name)
  expect_equal(exported$experimental_logfc[in_input], rows$logfc[source_rows])
  expect_equal(exported$experimental_pvalue[in_input], rows$pvalue[source_rows])
  expect_gt(sum(in_input), 0)

  # IDs that are not genes of tkoi::genes are dropped.
  expect_false(any(c("ENSG_NOT_A_GENE", "ENSG00000000000") %in% exported$gene_name))
})

# make_gene_exploration_plot() ---------------------------------------------------

test_that("make_gene_exploration_plot() builds a ggplot faceted into up- and down-regulated genes", {
  local_null_device()
  plot = make_gene_exploration_plot(toy_result())

  expect_true(is_ggplot(plot))
  built = expect_no_error(ggplot2::ggplot_build(plot))
  expect_setequal(as.character(built$layout$layout$direction), c("Down-regulated", "Up-regulated"))
  expect_gt(nrow(ggplot2::layer_data(plot, 1)), 0)
})

test_that("make_gene_exploration_plot() plots each measured gene at (-log10 p-value, -log10 FDR)", {
  local_null_device()
  plot = make_gene_exploration_plot(toy_result())
  built = ggplot2::ggplot_build(plot)
  points = ggplot2::layer_data(plot, 1)
  rows = first_rows()

  genes = gene_table()
  measured = genes[genes$ensembl %in% rows$gene_name, ]
  input_rows = match(measured$ensembl, rows$gene_name)
  expect_equal(nrow(points), nrow(measured))
  expect_equal(sort(points$x), sort(-log10(rows$pvalue[input_rows])))
  expect_equal(sort(points$y), sort(-log10(measured$fdr)))

  # Genes are split by the sign of their log fold change.
  panels = built$layout$layout
  up_panel = panels$PANEL[panels$direction == "Up-regulated"]
  expect_equal(sum(points$PANEL == up_panel), sum(rows$logfc[input_rows] >= 0))
})

test_that("make_gene_exploration_plot() colors genes by the log fold change threshold", {
  local_null_device()
  result = toy_result()
  plot = make_gene_exploration_plot(result, sig_color = "red", non_sig_color = "blue")
  points = ggplot2::layer_data(plot, 1)

  expect_true(all(points$fill %in% c("red", "blue")))
  expected = ifelse(abs(plot$data$logfc) >= result@logfc_threshold, "red", "blue")
  expect_identical(points$fill, expected)
  expect_true(any(points$fill == "red"))
  expect_true(any(points$fill == "blue"))
})

test_that("make_gene_exploration_plot() draws the FDR and p-value threshold lines", {
  local_null_device()
  result = toy_result()
  geoms = function(plot) unname(vapply(plot$layers, function(layer) class(layer$geom)[1], character(1)))

  plot = make_gene_exploration_plot(result)
  expect_identical(geoms(plot), c("GeomPoint", "GeomHline", "GeomVline"))
  expect_equal(unique(ggplot2::layer_data(plot, 2)$yintercept), -log10(0.05))
  expect_equal(unique(ggplot2::layer_data(plot, 3)$xintercept), -log10(result@pvalue_threshold))

  # The vertical line follows the p-value threshold of the run; the FDR line stays at 0.05.
  strict = result
  strict@pvalue_threshold = 0.001
  plot = make_gene_exploration_plot(strict)
  expect_equal(unique(ggplot2::layer_data(plot, 2)$yintercept), -log10(0.05))
  expect_equal(unique(ggplot2::layer_data(plot, 3)$xintercept), 3)
})

test_that("make_gene_exploration_plot() draws without messages", {
  local_null_device()
  plot = make_gene_exploration_plot(toy_result())
  # The plot has no fill scale, so a fill label would make ggplot2 (>= 4.0)
  # print "Ignoring unknown labels" every time the plot is drawn.
  expect_no_message(ggplot2::ggplot_build(plot))
  expect_no_message(print(plot))
})

# visualize_topn() -----------------------------------------------------------------

test_that("visualize_topn() draws at most top_n bars for each node type", {
  local_null_device()
  result = toy_result()
  for (category in c("BiologicalProcess", "Gene", "Pathway", "Disease", "Compound", "CellType")) {
    plot = visualize_topn(result, category = category, top_n = 10)
    expect_true(is_ggplot(plot))
    expect_no_error(ggplot2::ggplot_build(plot))
    bars = ggplot2::layer_data(plot)
    expect_gt(nrow(bars), 0)
    expect_lte(nrow(bars), 10)
    expect_identical(plot$labels$title, category)
  }
})

test_that("visualize_topn() respects top_n", {
  local_null_device()
  result = toy_result()
  bars = function(...) nrow(ggplot2::layer_data(visualize_topn(result, ...)))

  expect_equal(bars(category = "BiologicalProcess", top_n = 1), 1)
  expect_equal(bars(category = "BiologicalProcess", top_n = 5), 5)
  expect_equal(bars(category = "BiologicalProcess"), 25)
  # A top_n beyond the table size draws every node that has an effect size.
  cell_types = result@network_summary_statistics$CellType
  expect_equal(bars(category = "CellType", top_n = 1e4), sum(!is.na(cell_types$beta)))
})

test_that("visualize_topn() plots the top rows' effect sizes and names", {
  local_null_device()
  result = toy_result()
  processes = result@network_summary_statistics$BiologicalProcess
  processes = processes[!is.na(processes$beta), ]
  top = utils::head(processes, 8)

  plot = visualize_topn(result, category = "BiologicalProcess", top_n = 8, lognorm = FALSE)
  bars = ggplot2::layer_data(plot)
  expect_equal(sort(bars$x), sort(top$beta))
  # BiologicalProcess bars are labelled with the term names, ordered by effect size.
  labels = levels(plot$data$identifier)
  expect_setequal(labels, top$name)
  expect_equal(plot$data$beta, sort(plot$data$beta))
  expect_identical(as.character(plot$labels$x), "Network Enrichment Effect Size")

  genes = gene_table()
  plot = visualize_topn(result, category = "Gene", top_n = 8, lognorm = FALSE)
  expect_setequal(as.character(plot$data$identifier), utils::head(genes$name[!is.na(genes$beta)], 8))
})

test_that("visualize_topn() labels Disease and MolecularFunction bars with names, other types with identifiers", {
  local_null_device()
  result = toy_result()
  top_rows = function(category, n) {
    table = result@network_summary_statistics[[category]]
    utils::head(table[!is.na(table$beta), ], n)
  }

  for (category in c("Disease", "MolecularFunction")) {
    top = top_rows(category, 6)
    plot = visualize_topn(result, category = category, top_n = 6, lognorm = FALSE)
    # Named nodes are labelled by name; unannotated nodes fall back to their identifier.
    expect_setequal(levels(plot$data$identifier), unique(dplyr::coalesce(top$name, top$identifier)))
    named = !is.na(top$name)
    expect_false(any(levels(plot$data$identifier) %in% top$identifier[named]))
  }

  top = top_rows("Pathway", 6)
  plot = visualize_topn(result, category = "Pathway", top_n = 6, lognorm = FALSE)
  expect_setequal(levels(plot$data$identifier), top$identifier)
})

test_that("visualize_topn() applies the ranknorm and lognorm transformations", {
  local_null_device()
  result = toy_result()
  processes = result@network_summary_statistics$BiologicalProcess
  beta = processes$beta[!is.na(processes$beta)]
  top_beta = function(...) {
    plot = visualize_topn(result, category = "BiologicalProcess", top_n = 6, ...)
    list(x = sort(ggplot2::layer_data(plot)$x), label = as.character(plot$labels$x))
  }

  plain = top_beta(ranknorm = FALSE, lognorm = FALSE)
  logged = top_beta(ranknorm = FALSE, lognorm = TRUE)
  ranked = top_beta(ranknorm = TRUE, lognorm = FALSE)
  both = top_beta(ranknorm = TRUE, lognorm = TRUE)

  expect_equal(plain$x, sort(utils::head(beta, 6)))
  expect_equal(logged$x, sort(log2(utils::head(beta, 6))))
  expect_equal(ranked$x, sort(utils::head(RNOmni::RankNorm(beta), 6)))
  # With both options, ranknorm is ignored.
  expect_equal(both$x, logged$x)

  expect_identical(plain$label, "Network Enrichment Effect Size")
  for (transformed in list(logged, ranked, both)) {
    expect_identical(transformed$label, "Transformed Network Enrichment Effect Size")
  }
  for (ranknorm in c(FALSE, TRUE)) {
    for (lognorm in c(FALSE, TRUE)) {
      plot = visualize_topn(result, category = "Disease", top_n = 6, ranknorm = ranknorm, lognorm = lognorm)
      expect_no_error(ggplot2::ggplot_build(plot))
    }
  }
})

test_that("visualize_topn() clamps zero FDRs and p-values to the smallest positive FDR", {
  local_null_device()
  result = toy_result()
  processes = result@network_summary_statistics$BiologicalProcess
  with_beta = processes[!is.na(processes$beta), ]
  smallest = min(with_beta$fdr[with_beta$fdr > 0], na.rm = TRUE)

  zero = result
  processes$fdr[1] = 0
  processes$p_value[1] = 0
  zero@network_summary_statistics$BiologicalProcess = processes
  plot = visualize_topn(zero, category = "BiologicalProcess", top_n = 5, lognorm = FALSE)

  row = plot$data[plot$data$node_id == processes$node_id[1], ]
  expect_equal(nrow(row), 1)
  expect_equal(row$fdr, smallest)
  expect_equal(row$p_value, smallest)
  expect_true(all(plot$data$fdr >= smallest))
  expect_true(all(plot$data$p_value >= smallest))
  # The bar color, -log10(FDR), stays finite.
  expect_true(all(is.finite(-log10(plot$data$fdr))))
  expect_no_error(ggplot2::ggplot_build(plot))

  # Without any positive FDR, the limit is the smallest positive double.
  processes$fdr = 0
  zero@network_summary_statistics$BiologicalProcess = processes
  plot = visualize_topn(zero, category = "BiologicalProcess", top_n = 5, lognorm = FALSE)
  expect_equal(unique(plot$data$fdr), .Machine$double.xmin)
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("visualize_topn() rejects unknown categories and lists the available ones", {
  result = toy_result()
  available = paste(names(result@network_summary_statistics), collapse = ", ")
  message = paste0("`category` must be one of: ", available, ".")

  expect_error(visualize_topn(result, category = "NotANodeType"), message, fixed = TRUE)
  # A node type the toy network does not have.
  expect_error(visualize_topn(result, category = "Protein"), message, fixed = TRUE)
  expect_error(visualize_topn(result, category = c("Gene", "Pathway")), message, fixed = TRUE)
  expect_error(visualize_topn(result, category = 1), message, fixed = TRUE)
})

test_that("visualize_topn() rejects an invalid top_n", {
  result = toy_result()
  for (top_n in list(0, -3, NA, NA_real_, "10", c(5, 10), NULL)) {
    expect_error(
      visualize_topn(result, category = "Pathway", top_n = top_n),
      "`top_n` must be a positive number.",
      fixed = TRUE
    )
  }
})

test_that("visualize_topn() keeps one bar per identifier when a node is repeated", {
  local_null_device()
  result = toy_result()
  duplicated_result = result
  processes = result@network_summary_statistics$BiologicalProcess
  duplicated_result@network_summary_statistics$BiologicalProcess = dplyr::bind_rows(processes[1, ], processes)
  pathways = result@network_summary_statistics$Pathway
  duplicated_result@network_summary_statistics$Pathway = dplyr::bind_rows(pathways[c(1, 1), ], pathways)

  plot = expect_no_error(visualize_topn(duplicated_result, category = "BiologicalProcess", top_n = 5))
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_equal(nrow(ggplot2::layer_data(plot)), 4)
  expect_false(anyDuplicated(as.character(plot$data$identifier)) > 0)

  plot = expect_no_error(visualize_topn(duplicated_result, category = "Pathway", top_n = 5, lognorm = FALSE))
  expect_equal(nrow(ggplot2::layer_data(plot)), 3)
  expect_setequal(as.character(plot$data$identifier), pathways$identifier[1:3])
})

test_that("visualize_topn() errors when no node of the type has an effect size", {
  no_beta = toy_result()
  no_beta@network_summary_statistics$Pathway$beta = NA_real_
  expect_error(visualize_topn(no_beta, category = "Pathway"), "No Pathway nodes have an effect size to plot.")
})

# export_network_summary_statistics() -------------------------------------------------

# Sheet names and cell ranges ("A1:L93") of an .xlsx file, read from its XML parts.
xlsx_sheet_ranges = function(path) {
  directory = withr::local_tempdir()
  utils::unzip(path, exdir = directory)
  read_part = function(part) {
    paste(readLines(file.path(directory, part), warn = FALSE, encoding = "UTF-8"), collapse = "")
  }
  elements = function(xml, tag) regmatches(xml, gregexpr(paste0("<", tag, "\\b[^>]*>"), xml))[[1]]
  attribute = function(element, name) sub(paste0(".*[[:space:]]", name, "=\"([^\"]*)\".*"), "\\1", element)

  sheets = elements(read_part("xl/workbook.xml"), "sheet")
  relations = elements(read_part("xl/_rels/workbook.xml.rels"), "Relationship")
  targets = attribute(relations, "Target")[match(attribute(sheets, "r:id"), attribute(relations, "Id"))]
  parts = ifelse(startsWith(targets, "/"), substring(targets, 2), file.path("xl", targets))
  ranges = vapply(parts, function(part) attribute(elements(read_part(part), "dimension"), "ref"), character(1))
  data.frame(sheet = attribute(sheets, "name"), range = unname(ranges))
}

# Number of rows and columns spanned by a range such as "A1:L93".
range_size = function(range) {
  last_cell = sub(".*:", "", range)
  column_letters = strsplit(sub("[0-9]+$", "", last_cell), "")[[1]]
  columns = sum(match(column_letters, LETTERS) * 26^(rev(seq_along(column_letters)) - 1))
  c(rows = as.integer(sub("^[A-Z]+", "", last_cell)), columns = columns)
}

test_that("export_network_summary_statistics() writes one worksheet per node type", {
  result = toy_result()
  path = withr::local_tempfile(fileext = ".xlsx")
  export_network_summary_statistics(result, filename = path)

  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)

  # An .xlsx file is a zip archive; read the sheet names from its workbook part.
  contents = utils::unzip(path, list = TRUE)$Name
  expect_true("xl/workbook.xml" %in% contents)
  expect_length(grep("^xl/worksheets/sheet[0-9]+\\.xml$", contents), length(result@network_summary_statistics))

  sheets = xlsx_sheet_ranges(path)
  expect_identical(sheets$sheet, names(result@network_summary_statistics))
})

test_that("export_network_summary_statistics() overwrites the file and returns its path invisibly", {
  result = toy_result()
  path = withr::local_tempfile(fileext = ".xlsx")
  writeLines("not a workbook", path)

  returned = expect_invisible(export_network_summary_statistics(result, filename = path))
  expect_identical(normalizePath(returned), normalizePath(path))
  expect_identical(xlsx_sheet_ranges(path)$sheet, names(result@network_summary_statistics))
})

test_that("export_network_summary_statistics() writes every row and column of each table", {
  result = toy_result()
  path = withr::local_tempfile(fileext = ".xlsx")
  export_network_summary_statistics(result, filename = path)

  sheets = xlsx_sheet_ranges(path)
  expect_equal(nrow(sheets), length(result@network_summary_statistics))
  for (i in seq_len(nrow(sheets))) {
    table = result@network_summary_statistics[[sheets$sheet[i]]]
    size = range_size(sheets$range[i])
    # One header row, then one row per node.
    expect_equal(unname(size["rows"]) - 1, nrow(table), label = paste(sheets$sheet[i], "data rows"))
    expect_equal(unname(size["columns"]), ncol(table), label = paste(sheets$sheet[i], "columns"))
  }
})

# compute_network_enrichment() ---------------------------------------------------------

test_that("compute_network_enrichment() computes the z-score and upper-tail p-value", {
  enrichment = compute_network_enrichment(list(pagerank = 0.3, perm.1 = 0.1, perm.2 = 0.2, perm.3 = 0.15))
  expect_s3_class(enrichment, "data.frame")
  expect_named(enrichment, c("beta", "p_value"))
  expect_equal(enrichment$beta, 3)
  expect_equal(enrichment$p_value, stats::pnorm(3, lower.tail = FALSE))

  nodes = data.frame(pagerank = c(0.3, 0.1), perm.1 = c(0.1, 0.1), perm.2 = c(0.2, 0.12))
  enrichment = compute_network_enrichment(nodes)
  expect_equal(nrow(enrichment), 2)
  expect_equal(enrichment$beta, c((0.3 - 0.15) / stats::sd(c(0.1, 0.2)), (0.1 - 0.11) / stats::sd(c(0.1, 0.12))))
  expect_equal(enrichment$p_value, stats::pnorm(enrichment$beta, lower.tail = FALSE))
})

test_that("compute_network_enrichment() reproduces the run's beta and p-value for every node", {
  result = toy_result()
  pagerank_data = result@pagerank_data
  enrichment = compute_network_enrichment(pagerank_data)

  expect_s3_class(enrichment, "data.frame")
  expect_named(enrichment, c("beta", "p_value"))
  expect_equal(nrow(enrichment), nrow(pagerank_data))

  statistics = all_node_statistics(result)
  rows = match(statistics$node_id, pagerank_data$node_id)
  expect_false(anyNA(rows))
  expect_equal(unname(enrichment$beta[rows]), statistics$beta, tolerance = 1e-8)
  expect_equal(unname(enrichment$p_value[rows]), statistics$p_value, tolerance = 1e-8)
})

test_that("compute_network_enrichment() gives the same answer for a list, a one-row data frame and a table", {
  result = toy_result()
  pagerank_data = result@pagerank_data
  table_form = compute_network_enrichment(pagerank_data)
  statistics = all_node_statistics(result)
  top_process = result@network_summary_statistics$BiologicalProcess$node_id[1]
  top_gene = gene_table()$node_id[1]

  for (node_id in c(top_process, top_gene)) {
    row = match(node_id, pagerank_data$node_id)
    expected = statistics[statistics$node_id == node_id, ]

    from_list = compute_network_enrichment(as.list(pagerank_data[row, ]))
    from_row = compute_network_enrichment(pagerank_data[row, , drop = FALSE])
    expect_equal(nrow(from_list), 1)
    expect_equal(nrow(from_row), 1)
    expect_equal(unname(from_list$beta), expected$beta, tolerance = 1e-8)
    expect_equal(unname(from_list$p_value), expected$p_value, tolerance = 1e-8)
    expect_equal(unname(from_row$beta), unname(from_list$beta), tolerance = 1e-12)
    expect_equal(unname(from_row$p_value), unname(from_list$p_value), tolerance = 1e-12)
    expect_equal(unname(table_form$beta[row]), unname(from_list$beta), tolerance = 1e-12)
  }
})

test_that("compute_network_enrichment() treats a list of equal-length vectors as a table of nodes", {
  nodes = list(
    pagerank = c(0.3, 0.1, 0.2),
    perm.1 = c(0.1, 0.1, 0.2),
    perm.2 = c(0.2, 0.12, 0.1),
    perm.3 = c(0.15, 0.05, 0.3)
  )
  from_list = compute_network_enrichment(nodes)
  expect_equal(nrow(from_list), 3)
  expect_equal(from_list, compute_network_enrichment(as.data.frame(nodes)))
  expect_equal(from_list$beta[2], (0.1 - mean(c(0.1, 0.12, 0.05))) / stats::sd(c(0.1, 0.12, 0.05)))

  # The whole pagerank_data table as a list, node_id column included.
  pagerank_data = toy_result()@pagerank_data
  from_list = compute_network_enrichment(as.list(pagerank_data))
  from_table = compute_network_enrichment(pagerank_data)
  expect_equal(nrow(from_list), nrow(pagerank_data))
  expect_equal(unname(from_list$beta), unname(from_table$beta))
  expect_equal(unname(from_list$p_value), unname(from_table$p_value))
})

test_that("compute_network_enrichment() needs pagerank and at least two perm values", {
  message = "`node` needs a `pagerank` value and at least two `perm*` values."
  # Data frames.
  expect_error(compute_network_enrichment(data.frame(pagerank = 0.3, perm.1 = 0.1)), message, fixed = TRUE)
  expect_error(compute_network_enrichment(data.frame(pagerank = 0.3)), message, fixed = TRUE)
  expect_error(compute_network_enrichment(data.frame(perm.1 = 0.1, perm.2 = 0.2)), message, fixed = TRUE)
  # Lists of single values.
  expect_error(compute_network_enrichment(list(perm.1 = 0.1, perm.2 = 0.2, perm.3 = 0.3)), message, fixed = TRUE)
  expect_error(compute_network_enrichment(list(pagerank = 0.3, perm.1 = 0.1)), message, fixed = TRUE)
  expect_error(compute_network_enrichment(list(0.3, 0.1, 0.2)), message, fixed = TRUE)
  # Lists of vectors.
  expect_error(
    compute_network_enrichment(list(perm.1 = c(0.1, 0.2), perm.2 = c(0.2, 0.3))),
    message,
    fixed = TRUE
  )
  expect_error(
    compute_network_enrichment(list(pagerank = c(0.3, 0.1), perm.1 = c(0.1, 0.2))),
    message,
    fixed = TRUE
  )

  # The layout of pagerank_data from run_tkoi(keep_permutations = FALSE).
  summary_only = toy_result()@pagerank_data[, c("node_id", "pagerank")]
  summary_only$null_mean = 0
  summary_only$null_sd = 1
  expect_error(compute_network_enrichment(summary_only), message, fixed = TRUE)
})

test_that("compute_network_enrichment() accepts a numeric matrix with pagerank and perm columns", {
  # The function has a branch for matrices, but the argument check reads
  # names(), which is NULL for a matrix, so a valid matrix is rejected.
  nodes = data.frame(pagerank = c(0.3, 0.1), perm.1 = c(0.1, 0.1), perm.2 = c(0.2, 0.12), perm.3 = c(0.4, 0.3))
  expect_equal(compute_network_enrichment(as.matrix(nodes)), compute_network_enrichment(nodes))
})

# run_gene_enrichment() -----------------------------------------------------------------

# A table in the layout of clusterProfiler's enrichResult@result.
go_enrichment_table = function(ids) {
  n = length(ids)
  data.frame(
    ID = ids,
    Description = paste("term", seq_len(n)),
    GeneRatio = paste0(seq_len(n), "/20"),
    BgRatio = paste0(seq_len(n) * 10, "/1000"),
    RichFactor = seq_len(n) / 10,
    FoldEnrichment = seq_len(n) + 1,
    zScore = seq_len(n) + 0.5,
    pvalue = 10^-(seq_len(n) + 2),
    p.adjust = 10^-(seq_len(n) + 1),
    qvalue = 10^-seq_len(n),
    geneID = paste0("ENSG", seq_len(n)),
    Count = seq_len(n),
    row.names = ids
  )
}

# Replaces clusterProfiler::enrichGO() for the calling test. `results` maps an
# ontology (BP, CC, MF) to the enrichResult to return; other ontologies give
# NULL, as enrichGO() does when nothing maps. Returns the calls it received.
local_mock_enrich_go = function(results, env = parent.frame()) {
  calls = new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    # Argument names follow clusterProfiler::enrichGO().
    enrichGO = function(gene, keyType, ont, ...) { # nolint: object_name_linter.
      assign(ont, list(gene = gene, keyType = keyType), envir = calls)
      results[[ont]]
    },
    .package = "clusterProfiler",
    .env = env
  )
  calls
}

test_that("run_gene_enrichment() rejects input that is not a tKOIList", {
  calls = local_mock_enrich_go(list())
  message = "`tkoi_list` must be a tKOIList object from run_tkoi()."
  for (not_a_result in list(list(), toy_expression(), NULL, "tKOIList")) {
    expect_error(run_gene_enrichment(not_a_result), message, fixed = TRUE)
  }
  expect_length(ls(calls), 0)
})

test_that("run_gene_enrichment() errors when no gene passes pvalue_threshold", {
  calls = local_mock_enrich_go(list())
  message = "No genes pass `pvalue_threshold`; nothing to enrich."

  strict = toy_result()
  strict@pvalue_threshold = 0
  expect_error(run_gene_enrichment(strict), message, fixed = TRUE)

  missing_pvalues = toy_result()
  missing_pvalues@expression_data$pvalue = NA_real_
  expect_error(run_gene_enrichment(missing_pvalues), message, fixed = TRUE)
  expect_length(ls(calls), 0)
})

test_that("run_gene_enrichment() joins the GO enrichment tables with the tKOI statistics", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not(methods::isClass("enrichResult"), "clusterProfiler's enrichResult class is not available")
  result = toy_result()
  processes = result@network_summary_statistics$BiologicalProcess
  functions = result@network_summary_statistics$MolecularFunction
  # Two GO terms of the run and one it does not have.
  process_ids = c(processes$identifier[1:2], "GO:9999999")
  function_ids = functions$identifier[1]
  expect_false("GO:9999999" %in% c(processes$identifier, functions$identifier))

  calls = local_mock_enrich_go(list(
    BP = methods::new("enrichResult", result = go_enrichment_table(process_ids)),
    MF = methods::new("enrichResult", result = go_enrichment_table(function_ids))
  ))
  messages = testthat::capture_messages(run_gene_enrichment(result))
  expect_identical(messages, paste0(
    "Running GO enrichment on: ", c("Biological Process", "Cellular Component", "Molecular Function"), "\n"
  ))
  enriched = suppressMessages(run_gene_enrichment(result))

  # One call per ontology, on Ensembl IDs of the genes that pass the p-value threshold.
  expect_setequal(ls(calls), c("BP", "CC", "MF"))
  rows = first_rows()
  expected_genes = unique(rows$gene_name[rows$pvalue <= result@pvalue_threshold])
  for (ontology in c("BP", "CC", "MF")) {
    expect_setequal(calls[[ontology]]$gene, expected_genes)
    expect_identical(calls[[ontology]]$keyType, "ENSEMBL")
  }

  expect_s4_class(enriched, "tKOIList")
  expect_named(
    enriched@gene_enrichment_comparison,
    c("enrichment_result", "comparison_scatter1", "comparison_scatter2")
  )
  table = enriched@gene_enrichment_comparison$enrichment_result
  expect_named(table, c(
    "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count",
    "Namespace", "tkoi_node_id", "tkoi_pagerank", "tkoi_beta", "tkoi_p_value", "tkoi_fdr", "definition"
  ))

  # Enrichment rows are kept in order, labelled with their namespace (CC gave no result).
  expect_identical(table$ID, c(process_ids, function_ids))
  expect_identical(table$Namespace, c(rep("Biological Process", 3), "Molecular Function"))
  source = go_enrichment_table(process_ids)
  expect_equal(table$qvalue[1:3], source$qvalue)
  expect_identical(table$geneID[1:3], source$geneID)

  # tKOI columns come from the node with the same GO identifier; NA when the run has none.
  tkoi_rows = rbind(
    processes[match(process_ids, processes$identifier), ],
    functions[match(function_ids, functions$identifier), ]
  )
  expect_identical(table$tkoi_node_id, tkoi_rows$node_id)
  expect_identical(table$tkoi_pagerank, tkoi_rows$pagerank)
  expect_identical(table$tkoi_beta, tkoi_rows$beta)
  expect_identical(table$tkoi_p_value, tkoi_rows$p_value)
  expect_identical(table$tkoi_fdr, tkoi_rows$fdr)
  expect_identical(table$definition, tkoi_rows$definition)
  expect_true(all(is.na(table[3, c("tkoi_node_id", "tkoi_pagerank", "tkoi_beta", "tkoi_fdr")])))
  expect_false(anyNA(table$tkoi_node_id[-3]))

  scatter1 = enriched@gene_enrichment_comparison$comparison_scatter1
  scatter2 = enriched@gene_enrichment_comparison$comparison_scatter2
  expect_true(is_ggplot(scatter1))
  expect_true(is_ggplot(scatter2))
  expect_s3_class(scatter1$facet, "FacetNull")
  expect_s3_class(scatter2$facet, "FacetGrid")
  expect_identical(scatter1$data, table)
})

test_that("run_gene_enrichment() errors when GO enrichment finds no terms", {
  skip_if_not_installed("clusterProfiler")
  local_mock_enrich_go(list())
  expect_error(
    suppressMessages(run_gene_enrichment(toy_result())),
    "GO enrichment returned no terms for these genes.",
    fixed = TRUE
  )
})

test_that("run_gene_enrichment() enriches the same gene rows as run_tkoi()", {
  skip_if_not_installed("clusterProfiler")
  # The rows every other function uses: no blank gene names, first row per
  # gene. Here a blank name passes the p-value threshold, and a repeated gene
  # fails it in its first row but passes in its second.
  result = toy_result()
  rows = first_rows()
  failing = rows$gene_name[rows$pvalue > result@pvalue_threshold][1]
  edited = result
  edited@expression_data = rbind(
    input_expression(),
    data.frame(gene_name = c("", failing), logfc = c(2, 2), pvalue = c(1e-8, 1e-8))
  )

  calls = local_mock_enrich_go(list())
  expect_error(suppressMessages(run_gene_enrichment(edited)), "GO enrichment returned no terms")
  expected_genes = rows$gene_name[rows$pvalue <= result@pvalue_threshold]
  expect_setequal(calls$BP$gene, expected_genes)
  expect_false("" %in% calls$BP$gene)
  expect_false(failing %in% calls$BP$gene)
})

# get_neighboring_nodes() --------------------------------------------------------------

test_that("get_neighboring_nodes() matches igraph::ego() on the toy network", {
  result = toy_result()
  network = toy_network()
  loop_vertex = igraph::ends(network, igraph::E(network)[igraph::which_loop(network)], names = TRUE)[1, 1]
  nodes = c(
    seed_ids()[1],
    result@network_summary_statistics$BiologicalProcess$node_id[1],
    result@network_summary_statistics$Pathway$node_id[1],
    loop_vertex
  )
  for (node in nodes) {
    for (order in 0:2) {
      neighbors = get_neighboring_nodes(node, degree_expansion = order, subnetwork = network)
      expected = names(igraph::ego(network, order = order, nodes = node, mode = "all")[[1]])
      expect_type(neighbors, "character")
      expect_identical(neighbors, expected)
      expect_true(node %in% neighbors)
      expect_false(anyDuplicated(neighbors) > 0)
    }
    expect_identical(get_neighboring_nodes(node, 0, network), node)
    direct = names(igraph::neighbors(network, node, mode = "all"))
    expect_setequal(get_neighboring_nodes(node, 1, network), unique(c(node, direct)))
  }
})

test_that("get_neighboring_nodes() ignores edge direction and counts hops", {
  path = igraph::make_graph(c("a", "b", "b", "c", "c", "d", "d", "e"), directed = TRUE)
  expect_setequal(get_neighboring_nodes("c", 0, path), "c")
  expect_setequal(get_neighboring_nodes("c", 1, path), c("b", "c", "d"))
  expect_setequal(get_neighboring_nodes("c", 2, path), c("a", "b", "c", "d", "e"))
  expect_setequal(get_neighboring_nodes("a", 3, path), c("a", "b", "c", "d"))
  expect_setequal(get_neighboring_nodes("a", 10, path), c("a", "b", "c", "d", "e"))
})

test_that("get_neighboring_nodes() validates its inputs", {
  network = toy_network()
  node = seed_ids()[1]
  expect_error(
    get_neighboring_nodes(node, 1, subnetwork = list()),
    "`subnetwork` must be an igraph object.",
    fixed = TRUE
  )
  expect_error(
    get_neighboring_nodes(node, 1, subnetwork = igraph::as_data_frame(network)),
    "`subnetwork` must be an igraph object.",
    fixed = TRUE
  )
  for (bad_node in list(1, c(node, node), NA_character_, character(0), NULL)) {
    expect_error(
      get_neighboring_nodes(bad_node, 1, subnetwork = network),
      "`gene_node_id` must be a single node ID.",
      fixed = TRUE
    )
  }
  for (bad_order in list(-1, 1.5, NA, NA_real_, "1", c(1, 2), numeric(0))) {
    expect_error(
      get_neighboring_nodes(node, bad_order, subnetwork = network),
      "`degree_expansion` must be a non-negative whole number.",
      fixed = TRUE
    )
  }
  expect_error(
    get_neighboring_nodes("not-a-node", 1, subnetwork = network),
    "`gene_node_id` is not a vertex of `subnetwork`.",
    fixed = TRUE
  )
  expect_identical(formals(get_neighboring_nodes)$subnetwork, quote(tkoi::tkoi_net))
})

# plot_network() ---------------------------------------------------------------------------

test_that("plot_network() draws the target and the significant genes around it", {
  local_null_device()
  result = toy_result()
  target = processes_at(1)[1]
  processes = result@network_summary_statistics$BiologicalProcess
  target_row = processes[processes$node_id == target, ]

  graph = expect_invisible(plot_network(result, target, degree_expansion = 2, subnetwork = toy_network()))
  expect_true(igraph::is_igraph(graph))
  expect_true(target_row$identifier %in% igraph::V(graph)$identifier)
  expect_true(target_row$name %in% igraph::V(graph)$name)
  expect_equal(sort(igraph::V(graph)$identifier), sort(identifiers(nodes_between(target, 2))))
  expect_false(any(igraph::degree(graph) == 0))

  # Genes are labelled with their symbols.
  genes_in_plot = igraph::V(graph)$identifier[plotted_types(graph) == "Gene"]
  symbols = tkoi::genes$name[match(intersect(nodes_between(target, 2), seed_ids()), tkoi::genes$id)]
  expect_true(all(symbols %in% igraph::V(graph)$name))
  expect_gt(length(genes_in_plot), 0)

  graph = expect_invisible(plot_network(result, target, degree_expansion = 1, subnetwork = toy_network()))
  expect_equal(sort(igraph::V(graph)$identifier), sort(identifiers(nodes_between(target, 1))))
})

test_that("plot_network() colours the target orange, genes by log fold change and other nodes gray", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  plotted = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  arguments = plotted$arguments
  graph = arguments$x
  expect_true(igraph::is_igraph(graph))
  # The graph returned is the graph drawn.
  expect_identical(igraph::V(graph)$node_id, igraph::V(plotted$value)$node_id)
  expect_equal(dim(arguments$layout), c(igraph::vcount(graph), 2))

  colours = arguments$vertex.color
  node_ids = igraph::V(graph)$node_id
  is_target = node_ids == target
  is_gene = plotted_types(graph) == "Gene" & !is_target
  logfc = input_logfc(node_ids)
  measured = is_gene & !is.na(logfc)
  expect_length(colours, igraph::vcount(graph))
  expect_equal(sum(is_target), 1)

  # The plot holds up- and down-regulated genes, unmeasured genes and other node types.
  expect_gt(sum(measured & logfc > 0), 0)
  expect_gt(sum(measured & logfc < 0), 0)
  expect_gt(sum(is_gene & is.na(logfc)), 0)
  expect_gt(sum(!is_gene & !is_target), 0)

  expect_identical(colours[is_target], "orange")
  expect_true(all(colours[!measured & !is_target] == "gray"))
  expect_identical(colours[measured], ramp_colour(logfc[measured]))

  # Up-regulated genes are on the red side of white, down-regulated ones on the
  # blue side, and the most changed gene gets the end of the ramp.
  rgb = grDevices::col2rgb(colours[measured])
  changed = abs(logfc[measured]) >= result@logfc_threshold
  up = changed & logfc[measured] > 0
  down = changed & logfc[measured] < 0
  expect_true(all(rgb["red", up] > rgb["blue", up]))
  expect_true(all(rgb["blue", down] > rgb["red", down]))
  most_changed = which.max(abs(logfc[measured]))
  expect_identical(
    colours[measured][most_changed],
    if (logfc[measured][most_changed] > 0) "#BB2A2D" else "#79ADD6"
  )
})

test_that("plot_network() sizes vertices by effect size between 10 and 20", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  plotted = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  sizes = plotted$arguments$vertex.size
  node_ids = igraph::V(plotted$arguments$x)$node_id
  statistics = all_node_statistics(result)
  beta = statistics$beta[match(node_ids, statistics$node_id)]
  finite = is.finite(beta)

  expect_length(sizes, length(node_ids))
  expect_true(all(sizes >= 10 & sizes <= 20))
  expect_equal(range(sizes), c(10, 20))
  # Sizes are linear in beta. Every node is reported, so every plotted node has a beta.
  expect_true(all(finite))
  expect_equal(sizes[finite], 10 + 10 * (beta[finite] - min(beta[finite])) / diff(range(beta[finite])))

  # A node without an effect size is drawn at the smallest size.
  unscored = setdiff(node_ids, target)[1]
  missing = result
  missing@network_summary_statistics = lapply(missing@network_summary_statistics, function(table) {
    table$beta[table$node_id == unscored] = NA
    table
  })
  replotted = capture_plot_arguments(
    plot_network(missing, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  replotted_ids = igraph::V(replotted$arguments$x)$node_id
  expect_identical(replotted$arguments$vertex.size[replotted_ids == unscored], 10)
})

test_that("plot_network() caps an infinite effect size at the largest finite one", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  infinite = result
  diseases = infinite@network_summary_statistics$Disease
  diseases$beta[diseases$node_id == target] = Inf
  infinite@network_summary_statistics$Disease = diseases

  plotted = capture_plot_arguments(
    plot_network(infinite, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  sizes = plotted$arguments$vertex.size
  node_ids = igraph::V(plotted$arguments$x)$node_id
  statistics = all_node_statistics(infinite)
  beta = statistics$beta[match(node_ids, statistics$node_id)]
  finite = is.finite(beta)

  expect_true(all(is.finite(sizes)))
  expect_true(all(sizes >= 10 & sizes <= 20))
  expect_equal(sizes[node_ids == target], 20)
  expect_equal(sizes[finite], 10 + 10 * (beta[finite] - min(beta[finite])) / diff(range(beta[finite])))
  expect_equal(sizes[which.max(replace(beta, !finite, -Inf))], 20)
})

test_that("plot_network() colours a repeated gene by its first row", {
  local_null_device()
  result = toy_result()
  network = toy_network()
  # The toy input lists one gene twice: first with a small change that fails
  # the filters (logfc near 0), then with logfc = -5. Draw it as the bridge of
  # a three-node chain: target - repeated gene - seed gene.
  repeated = input_expression()$gene_name[duplicated(input_expression()$gene_name)]
  repeated_id = toy_gene_ids(repeated)
  neighbors = names(igraph::neighbors(network, repeated_id, mode = "all"))
  seed = intersect(neighbors, seed_ids())[1]
  target = intersect(neighbors, result@network_summary_statistics$BiologicalProcess$node_id)[1]
  expect_false(is.na(seed))
  expect_false(is.na(target))
  chain = igraph::induced_subgraph(tagged_network(), c(target, repeated_id, seed))

  plotted = capture_plot_arguments(plot_network(result, target, degree_expansion = 2, subnetwork = chain))
  node_ids = igraph::V(plotted$arguments$x)$node_id
  expect_setequal(node_ids, c(target, repeated_id, seed))
  colours = stats::setNames(plotted$arguments$vertex.color, node_ids)

  logfc = input_logfc(c(repeated_id, seed))
  expect_equal(logfc[1], input_expression()$logfc[match(repeated, input_expression()$gene_name)])
  expected = ramp_colour(logfc)
  expect_identical(unname(colours[c(repeated_id, seed)]), expected)
  # With the second row, the repeated gene would be the darkest blue.
  expect_false(colours[[repeated_id]] == "#79ADD6")
})

test_that("plot_network() works on a network without the identifier attribute", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  default = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  no_identifier = igraph::delete_vertex_attr(tagged_network(), "identifier")
  plotted = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = no_identifier)
  )
  graph = plotted$value
  expect_true(igraph::is_igraph(graph))
  expect_false("identifier" %in% igraph::vertex_attr_names(graph))

  # Same nodes, colours and sizes as with identifiers.
  node_ids = igraph::V(graph)$node_id
  expect_identical(node_ids, igraph::V(default$value)$node_id)
  expect_identical(plotted$arguments$vertex.color, default$arguments$vertex.color)
  expect_identical(plotted$arguments$vertex.size, default$arguments$vertex.size)
  expect_true(any(!plotted$arguments$vertex.color %in% c("gray", "orange")))

  # Genes of tkoi::genes keep their symbols; other nodes, the target included,
  # fall back to their vertex names.
  symbols = tkoi::genes$name[match(node_ids, tkoi::genes$id)]
  expect_gt(sum(!is.na(symbols)), 0)
  expect_identical(igraph::V(graph)$name, dplyr::coalesce(symbols, node_ids))
  expect_identical(igraph::V(graph)$name[node_ids == target], target)
})

test_that("plot_network() reads node types from a label attribute", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  default = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = tagged_network())
  )
  relabelled = tagged_network() |>
    igraph::set_vertex_attr("label", value = igraph::V(tagged_network())$labels) |>
    igraph::delete_vertex_attr("labels")
  plotted = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 2, subnetwork = relabelled)
  )

  expect_identical(igraph::V(plotted$value)$node_id, igraph::V(default$value)$node_id)
  # Genes are still coloured by log fold change.
  expect_identical(plotted$arguments$vertex.color, default$arguments$vertex.color)
  expect_true(any(!plotted$arguments$vertex.color %in% c("gray", "orange")))
  expect_identical(plotted$arguments$vertex.size, default$arguments$vertex.size)
  expect_identical(igraph::V(plotted$value)$name, igraph::V(default$value)$name)
})

test_that("plot_network() draws node names, not a label attribute holding node types", {
  local_null_device()
  result = toy_result()
  target = colour_target()
  default = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 1, subnetwork = tagged_network())
  )
  relabelled = tagged_network() |>
    igraph::set_vertex_attr("label", value = igraph::V(tagged_network())$labels) |>
    igraph::delete_vertex_attr("labels")
  plotted = capture_plot_arguments(
    plot_network(result, target, degree_expansion = 1, subnetwork = relabelled)
  )

  # igraph draws a `label` vertex attribute in place of vertex names, so a
  # network that stores its node types in `label` must not show "['Gene']".
  expect_identical(drawn_labels(default$arguments), igraph::V(default$value)$name)
  expect_identical(drawn_labels(plotted$arguments), drawn_labels(default$arguments))
})

test_that("plot_network() supports every layout", {
  local_null_device()
  result = toy_result()
  target = processes_at(1)[1]
  for (layout in c("kk", "fr", "gem", "graphopt", "lgl", "mds")) {
    graph = expect_invisible(
      plot_network(result, target, network_layout_type = layout, subnetwork = toy_network())
    )
    expect_true(igraph::is_igraph(graph))
    expect_true(identifiers(target) %in% igraph::V(graph)$identifier)
  }
  expect_error(
    plot_network(result, target, network_layout_type = "circle", subnetwork = toy_network()),
    "should be one of"
  )
})

test_that("plot_network() errors clearly when no significant gene is within reach", {
  local_null_device()
  result = toy_result()
  network = toy_network()
  seeds = seed_ids()
  message = "No significant genes lie within `degree_expansion` hops of the target."

  # Two hops from the nearest seed: nothing at one hop, a plot at two.
  target = processes_at(2)[1]
  expect_error(plot_network(result, target, degree_expansion = 1, subnetwork = network), message, fixed = TRUE)
  graph = plot_network(result, target, degree_expansion = 2, subnetwork = network)
  expect_equal(sort(igraph::V(graph)$identifier), sort(identifiers(nodes_between(target, 2))))

  # Three hops from the nearest seed. Pick the node with the smallest
  # neighborhood so the path search of degree_expansion = 3 stays fast.
  far = processes_at(3)
  expect_gt(length(far), 0)
  neighborhood = vapply(far, function(x) sum(igraph::degree(network, igraph::ego(network, 2, x)[[1]])), numeric(1))
  target = far[which.min(neighborhood)]
  expect_error(plot_network(result, target, degree_expansion = 2, subnetwork = network), message, fixed = TRUE)

  graph = expect_invisible(plot_network(result, target, degree_expansion = 3, subnetwork = network))
  reachable_genes = seeds[igraph::distances(network, v = target, to = seeds, mode = "all")[1, ] <= 3]
  expect_true(all(identifiers(c(target, reachable_genes)) %in% igraph::V(graph)$identifier))
  # Every node drawn lies on a walk of at most three edges from the target to a gene.
  from_target = igraph::distances(network, v = target, mode = "all")[1, ]
  from_genes = apply(igraph::distances(network, v = reachable_genes, mode = "all"), 2, min)
  expect_true(all(igraph::V(graph)$identifier %in% identifiers(names(which(from_target + from_genes <= 3)))))
})

test_that("plot_network() follows paths of three edges on a small network", {
  local_null_device()
  result = toy_result()
  network = toy_network()
  seeds = seed_ids()
  target = processes_at(3)[1]
  gene = seeds[igraph::distances(network, v = target, to = seeds, mode = "all")[1, ] == 3][1]
  path = igraph::shortest_paths(network, from = target, to = gene, mode = "all")$vpath[[1]]
  chain = igraph::induced_subgraph(network, path)
  expect_equal(igraph::vcount(chain), 4)

  message = "No significant genes lie within `degree_expansion` hops of the target."
  expect_error(plot_network(result, target, degree_expansion = 1, subnetwork = chain), message, fixed = TRUE)
  expect_error(plot_network(result, target, degree_expansion = 2, subnetwork = chain), message, fixed = TRUE)
  graph = expect_invisible(plot_network(result, target, degree_expansion = 3, subnetwork = chain))
  expect_equal(sort(igraph::V(graph)$identifier), sort(igraph::V(chain)$identifier))
  expect_equal(igraph::ecount(graph), igraph::ecount(chain))
})

test_that("plot_network() validates the target and the network", {
  result = toy_result()
  network = toy_network()
  target = processes_at(1)[1]
  message = "`target_node_id` must be a single vertex name of `subnetwork`."
  for (bad_target in list("not-a-node", c(target, target), 1, NA_character_, character(0))) {
    expect_error(plot_network(result, bad_target, subnetwork = network), message, fixed = TRUE)
  }
  expect_error(
    plot_network(result, target, subnetwork = list()),
    "`subnetwork` must be an igraph object.",
    fixed = TRUE
  )
  expect_error(
    plot_network(result, target, degree_expansion = -1, subnetwork = network),
    "`degree_expansion` must be a non-negative whole number.",
    fixed = TRUE
  )
})

# Repeated and blank genes ----------------------------------------------------------------

test_that("export_gene_exploration_data() and make_gene_exploration_plot() use the first row of a repeated gene", {
  local_null_device()
  result = toy_result()
  # The toy input lists one gene twice. run_tkoi() uses its first row, so the
  # exports should report that row once.
  expression = input_expression()
  repeated = expression$gene_name[duplicated(expression$gene_name)]
  expect_length(repeated, 1)
  first = first_rows()[first_rows()$gene_name == repeated, ]

  exported = export_gene_exploration_data(result)
  repeated_rows = exported[exported$gene_name %in% repeated, ]
  expect_equal(nrow(repeated_rows), 1)
  expect_equal(unique(repeated_rows$experimental_logfc), first$logfc)
  expect_equal(unique(repeated_rows$experimental_pvalue), first$pvalue)

  plot = make_gene_exploration_plot(result)
  genes = gene_table()
  measured = genes[genes$ensembl %in% expression$gene_name, ]
  expect_equal(nrow(ggplot2::layer_data(plot, 1)), nrow(measured))
  expect_equal(sum(plot$data$gene_name == repeated), 1)
  expect_equal(plot$data$logfc[plot$data$gene_name %in% repeated], first$logfc)
})

test_that("plot_network() uses the first row of a repeated gene, like run_tkoi()", {
  local_null_device()
  result = toy_result()
  network = toy_network()
  # The toy input lists one gene twice: first with values that fail the
  # filters, then with values that pass. run_tkoi() uses the first row, so the
  # gene is not a seed and must not be drawn as a significant gene.
  expression = input_expression()
  repeated = expression$gene_name[duplicated(expression$gene_name)]
  expect_length(repeated, 1)
  expect_false(passes_filters(first_rows()[first_rows()$gene_name == repeated, ]))
  repeated_id = toy_gene_ids(repeated)
  expect_false(repeated_id %in% seed_ids())

  processes = result@network_summary_statistics$BiologicalProcess
  next_to_repeated = names(igraph::neighbors(network, repeated_id, mode = "all"))
  target = processes$node_id[processes$node_id %in% next_to_repeated & is.na(processes$direct_links)][1]
  expect_false(is.na(target))

  expect_error(
    plot_network(result, target, degree_expansion = 1, subnetwork = network),
    "No significant genes lie within `degree_expansion` hops of the target.",
    fixed = TRUE
  )
})

test_that("blank gene names in the input are ignored downstream, like in run_tkoi()", {
  local_null_device()
  result = toy_result()
  genes = gene_table()
  # Genes without an Ensembl ID have an empty `ensembl` field in tkoi::genes.
  no_ensembl = as.data.frame(tkoi::genes)
  no_ensembl = no_ensembl[!is.na(no_ensembl$ensembl) & no_ensembl$ensembl == "", ][1, ]
  blank_input = result
  blank_input@expression_data = rbind(input_expression(), data.frame(gene_name = "", logfc = 3, pvalue = 1e-6))

  # export_gene_exploration_data(): the blank row must not be attributed to a gene.
  gene_row = genes[1, ]
  gene_row$node_id = no_ensembl$id
  gene_row$identifier = as.character(no_ensembl$identifier)
  gene_row$ensembl = no_ensembl$ensembl
  gene_row$name = no_ensembl$name
  blank_input@network_summary_statistics$Gene = dplyr::bind_rows(genes, gene_row)
  exported = export_gene_exploration_data(blank_input)
  expect_equal(nrow(exported), nrow(genes) + 1)
  expect_true(is.na(exported$experimental_logfc[exported$id == no_ensembl$id]))

  # make_gene_exploration_plot(): no point for that gene.
  plot = make_gene_exploration_plot(blank_input)
  expect_false(no_ensembl$id %in% plot$data$id)
  expect_equal(nrow(ggplot2::layer_data(plot, 1)), sum(genes$ensembl %in% first_rows()$gene_name))

  # plot_network(): the blank row must not make that gene significant.
  target = result@network_summary_statistics$BiologicalProcess$node_id[1]
  small = igraph::graph_from_data_frame(
    data.frame(from = target, to = no_ensembl$id),
    directed = FALSE,
    vertices = data.frame(
      name = c(target, no_ensembl$id),
      identifier = c(igraph::V(toy_network())[target]$identifier, as.character(no_ensembl$identifier)),
      labels = c("['BiologicalProcess']", "['Gene']")
    )
  )
  expect_error(
    plot_network(blank_input, target, degree_expansion = 1, subnetwork = small),
    "No significant genes lie within `degree_expansion` hops of the target.",
    fixed = TRUE
  )
})
