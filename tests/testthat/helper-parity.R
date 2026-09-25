# Comparison of run_tkoi() results with legacy_run_tkoi() (tkoi 1.0.0), shared
# by test-run_tkoi.R and test-extended.R.

statistic_columns = c(
  "node_id", "node_type", "identifier", "pagerank", "beta", "p_value", "fdr",
  "direct_links", "indirect_links", "valency"
)

near_equal = function(x, y, tolerance = 1e-7) {
  both_missing = is.na(x) & is.na(y)
  close = x == y | abs(x - y) <= tolerance * pmax(abs(x), abs(y))
  both_missing | (!is.na(close) & close)
}

# Node-by-node agreement of `x` with the reference `y` (vectors, matrices, or
# numeric data frames of the same shape): the largest |x - y| / max(|y|, floor)
# must be at most `tolerance`. With floor = 1 the check is relative for values
# above 1 and absolute below; with floor = 0 it is relative everywhere. Missing
# values (NA or NaN) must be missing in both.
expect_per_node = function(x, y, tolerance, floor = 1, label = "values") {
  x = unname(as.matrix(x))
  y = unname(as.matrix(y))
  expect_identical(dim(x), dim(y), label = paste(label, "shape"))
  expect_identical(is.na(x), is.na(y), label = paste(label, "missing values"))
  both = !is.na(x) & !is.na(y)
  difference = abs(x[both] - y[both])
  scaled = ifelse(difference == 0, 0, difference / pmax(abs(y[both]), floor))
  expect_lte(max(scaled, 0), tolerance, label = paste(label, "largest scaled difference"))
}

# Rows of `new_table` are in the order of `legacy_table`, except that rows whose
# fdr and beta agree to about 7 significant digits may be swapped.
expect_order_up_to_ties = function(new_table, legacy_table, threshold, label) {
  # Statistics of the rows in the released order, taken from the new table.
  in_legacy_order = new_table[match(legacy_table$node_id, new_table$node_id), ]
  expect_identical(
    new_table$indirect_links >= threshold,
    in_legacy_order$indirect_links >= threshold,
    label = paste(label, "threshold groups")
  )
  tied = near_equal(new_table$fdr, in_legacy_order$fdr, 1e-7) &
    near_equal(new_table$beta, in_legacy_order$beta, 1e-7)
  expect_true(all(tied), label = paste(label, "rows out of order are near-ties"))
}

# Nodes whose null PageRank is constant (sd 0): tkoi 1.1.0 reports them as
# untestable (NaN), where the released algorithm gave NaN or +/-Inf.
constant_null = function(result) {
  pagerank = result@pagerank_data
  null = as.matrix(pagerank[, grep("^perm\\.", names(pagerank)), drop = FALSE])
  stats::setNames(apply(null, 1, stats::sd) == 0, pagerank$node_id)
}

# A run_tkoi() result agrees with legacy_run_tkoi() (tkoi 1.0.0) node by node,
# apart from the deliberate 1.1.0 changes: every node is kept (1.0.0 dropped
# unannotated nodes), Compound nodes are the repaired human metabolites, and
# nodes with a constant null are NaN. So: identical pagerank_data within
# 1e-11 (absolute) and `pagerank_relative` (relative); every 1.0.0 row is present with the same values (beta and
# p_value within 1e-7; fdr too wherever the FDR family is unchanged); the extra
# rows are unannotated nodes; rows keep the 1.0.0 order up to near-ties.
expect_matches_legacy = function(result, legacy, indirect_link_threshold = 3, pagerank_relative = 1e-8) {
  new_pagerank = result@pagerank_data
  legacy_pagerank = legacy$pagerank_data
  expect_identical(names(new_pagerank), names(legacy_pagerank))
  expect_identical(new_pagerank$node_id, legacy_pagerank$node_id)
  expect_per_node(new_pagerank[, -1], legacy_pagerank[, -1], 1e-11, label = "pagerank_data")
  expect_per_node(
    new_pagerank[, -1], legacy_pagerank[, -1], pagerank_relative,
    floor = 0, label = "pagerank_data (relative)"
  )

  untestable = constant_null(result)
  new_tables = result@network_summary_statistics
  legacy_tables = legacy$network_summary_statistics
  expect_true(all(names(legacy_tables) %in% names(new_tables)))
  for (type in names(legacy_tables)) {
    new_table = new_tables[[type]]
    legacy_table = legacy_tables[[type]]
    expect_identical(names(new_table), names(legacy_table), label = paste(type, "columns"))
    expect_true(all(legacy_table$node_id %in% new_table$node_id), label = paste(type, "keeps 1.0.0 rows"))

    aligned = new_table[match(legacy_table$node_id, new_table$node_id), ]
    constant = untestable[legacy_table$node_id]
    expect_per_node(aligned$pagerank, legacy_table$pagerank, 1e-11, label = paste(type, "pagerank"))
    expect_true(all(is.nan(aligned$beta[constant])), label = paste(type, "constant-null nodes are NaN"))
    expect_true(all(!is.finite(legacy_table$beta[constant])), label = paste(type, "1.0.0 gave these no finite beta"))
    for (column in c("beta", "p_value")) {
      expect_per_node(
        aligned[[column]][!constant], legacy_table[[column]][!constant], 1e-7,
        label = paste(type, column)
      )
    }
    # The FDR family is unchanged unless constant-null nodes left it or the
    # type is Compound (now the repaired human metabolites).
    all_of_type = names(untestable)[result@pagerank_data$node_id %in% new_table$node_id]
    if (type != "Compound" && !any(untestable[all_of_type])) {
      expect_per_node(aligned$fdr, legacy_table$fdr, 1e-7, label = paste(type, "fdr"))
    }
    for (column in setdiff(names(legacy_table), c("pagerank", "beta", "p_value", "fdr"))) {
      expect_identical(aligned[[column]], legacy_table[[column]], label = paste(type, column))
    }

    # Rows 1.0.0 did not report are nodes without a curated annotation (or,
    # for Compound, human metabolites whose identifiers were repaired).
    extra = new_table[!new_table$node_id %in% legacy_table$node_id, ]
    annotation_columns = setdiff(names(new_table), statistic_columns)
    if (type != "Compound" && nrow(extra) > 0 && length(annotation_columns) > 0) {
      expect_true(all(is.na(as.matrix(extra[, annotation_columns]))), label = paste(type, "extra rows are unannotated"))
    }

    if (type != "Compound" && !any(untestable[all_of_type])) {
      shared = new_table[new_table$node_id %in% legacy_table$node_id, ]
      expect_order_up_to_ties(shared, legacy_table, indirect_link_threshold, label = type)
    }
  }
}
