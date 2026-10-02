probability_fixture = function() {
  data.frame(
    node_id = paste0("GO:", seq_len(6)),
    name = c("Ordinary", "One", "Small", "Floored", "Legacy", "Unavailable"),
    beta = seq_len(6),
    p_value = c(0.05, 1, 4.89e-15, .Machine$double.xmin, 0, NA_real_),
    fdr = c(0.1, 1, 1.00e-13, .Machine$double.xmin, 0, NA_real_),
    log_p_value = c(log(0.05), 0, log(4.89e-15), -1000, -1200, NA_real_),
    log_fdr = c(log(0.1), 0, log(1.00e-13), -990, -1190, NA_real_),
    p_value_bounded = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    fdr_bounded = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
  )
}

testthat::test_that("probability displays preserve raw probabilities, logs and flags", {
  app_environment = shiny_app_environment()
  original = probability_fixture()
  before = serialize(original, NULL)
  presented = app_environment$format_result_table(original)
  exported = app_environment$probability_export_table(original)
  numeric_fields = c("p_value", "fdr", "log_p_value", "log_fdr", "p_value_bounded", "fdr_bounded")
  testthat::expect_identical(presented[numeric_fields], original[numeric_fields])
  testthat::expect_identical(exported[names(original)], original)
  testthat::expect_identical(serialize(original, NULL), before)
  testthat::expect_match(exported$p_value_display[1], "e-02$")
  testthat::expect_match(exported$p_value_display[2], "e\\+00$")
  testthat::expect_match(exported$p_value_display[3], "e-15$")
  testthat::expect_true(all(grepl("^<=.*e-308$", exported$p_value_display[4:5])))
  testthat::expect_identical(exported$p_value_display[6], "NA")
  testthat::expect_false(any(grepl("^0[.]0+e", exported$p_value_display)))
  html = app_environment$probability_labels(original, "p_value", "html")
  testthat::expect_match(html[3], "\u00d7 10<sup>\u221215</sup>", fixed = TRUE)
  testthat::expect_match(html[4], "^&le;")
  withr::local_options(scipen = 999)
  testthat::expect_identical(app_environment$probability_labels(original, "p_value"), exported$p_value_display)
})

testthat::test_that("legacy zero probabilities use logs or an explicit unresolved label", {
  app_environment = shiny_app_environment()
  recoverable = data.frame(p_value = 0, log_p_value = log(1e-12))
  testthat::expect_match(app_environment$probability_labels(recoverable, "p_value"), "^1[.]00e-12$")
  unknown = data.frame(p_value = 0)
  testthat::expect_identical(app_environment$probability_labels(unknown, "p_value"),
    "Unresolved zero (log unavailable)")
  testthat::expect_identical(unknown$p_value, 0)
})

testthat::test_that("DT retains numeric probability data and escapes source text", {
  app_environment = shiny_app_environment()
  data = probability_fixture()
  data$name[1] = "<script>untrusted()</script>"
  table = app_environment$probability_datatable(data)
  testthat::expect_identical(table$x$data$p_value, data$p_value)
  testthat::expect_identical(table$x$data$fdr, data$fdr)
  testthat::expect_type(table$x$data$p_value, "double")
  testthat::expect_match(table$x$data$.tkoi_p_value_html[3], "<sup>", fixed = TRUE)
  testthat::expect_match(table$x$data$.tkoi_p_value_html[5], "&le;", fixed = TRUE)
  # The runtime DOM test separately verifies that the supplied script is escaped.
  definitions = table$x$options$columnDefs
  testthat::expect_true(any(vapply(definitions, function(x) identical(x$type, "num"), logical(1))))
})

testthat::test_that("Excel exports remain numeric and retain visible bound labels", {
  app_environment = shiny_app_environment()
  data = probability_fixture()
  before = serialize(data, NULL)
  path = withr::local_tempfile(fileext = ".xlsx")
  workbook = openxlsx::createWorkbook()
  app_environment$write_result_worksheet(workbook, "Probabilities", data)
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  actual = openxlsx::read.xlsx(path, sheet = "Probabilities", skipEmptyRows = FALSE)
  testthat::expect_type(actual$p_value, "double")
  testthat::expect_type(actual$fdr, "double")
  testthat::expect_equal(actual$p_value[1:3], data$p_value[1:3], tolerance = 1e-14)
  testthat::expect_gt(actual$p_value[4], 0)
  testthat::expect_identical(actual$p_value[5], 0)
  testthat::expect_equal(actual$log_p_value[1:5], data$log_p_value[1:5])
  testthat::expect_true(all(grepl("^<=.*e-308$", actual$p_value_display[4:5])))
  testthat::expect_identical(serialize(data, NULL), before)
  extracted = withr::local_tempdir()
  utils::unzip(path, files = "xl/styles.xml", exdir = extracted)
  styles = readLines(file.path(extracted, "xl/styles.xml"), warn = FALSE)
  testthat::expect_match(styles, "0.00E+00", fixed = TRUE)
  testthat::expect_match(styles, "&lt;=", fixed = TRUE)
})

testthat::test_that("scientific example CSV notation round-trips at machine precision", {
  app_environment = shiny_app_environment()
  data = app_environment$example_data
  file = withr::local_tempfile(fileext = ".csv")
  data$pvalue = formatC(data$pvalue, format = "e", digits = 16)
  write.csv(data, file, row.names = FALSE)
  imported = read.csv(file)
  expected = app_environment$example_data$pvalue
  tolerance = 4 * .Machine$double.eps * pmax(abs(expected), .Machine$double.xmin)
  testthat::expect_true(all(abs(imported$pvalue - expected) <= tolerance))
  testthat::expect_true(all(grepl("e[+-][0-9]+$", data$pvalue)))
})


testthat::test_that("empty result categories remain valid empty tables", {
  app_environment = shiny_app_environment()
  data = probability_fixture()[FALSE, ]
  table = app_environment$probability_datatable(data)
  testthat::expect_identical(nrow(table$x$data), 0L)
  testthat::expect_identical(table$x$data$p_value, numeric())
  testthat::expect_identical(table$x$data$.tkoi_p_value_html, character())
})

probability_order_definition = function(table, column) {
  target = match(column, names(table$x$data)) - 1L
  matched = Filter(function(x) identical(x$targets, target) && identical(x$type, "num"),
    table$x$options$columnDefs)
  testthat::expect_length(matched, 1L)
  matched[[1]]
}

testthat::test_that("probability headers sort by unfloored log columns in the displayed table", {
  app_environment = shiny_app_environment()
  data = probability_fixture()
  # Two more rows at the double.xmin bound with distinct logs, stored out of log order.
  extra = data[c(4, 4), ]
  extra$node_id = c("GO:7", "GO:8")
  extra$name = c("Shallow bound", "Deep bound")
  extra$log_p_value = c(-990, -5000)
  extra$log_fdr = c(-985, -4990)
  data = rbind(data, extra)
  # Reverse the input columns so displayed indices differ from input indices.
  data = data[, rev(names(data))]
  before = serialize(data, NULL)
  table = app_environment$probability_datatable(data)
  shown = table$x$data
  displayed = app_environment$format_result_table(data)
  testthat::expect_identical(names(shown)[seq_along(displayed)], names(displayed))
  for (column in c("p_value", "fdr")) {
    log_column = paste0("log_", column)
    definition = probability_order_definition(table, column)
    testthat::expect_identical(definition$orderData, match(log_column, names(displayed)) - 1L)
    testthat::expect_identical(names(shown)[definition$orderData + 1L], log_column)
    testthat::expect_identical(definition$type, "num")
    # Display and filter still read the hidden HTML and plain label columns.
    html_index = match(paste0(".tkoi_", column, "_html"), names(shown)) - 1L
    plain_index = match(paste0(".tkoi_", column, "_plain"), names(shown)) - 1L
    testthat::expect_match(definition$render, sprintf("type==='display')return row[%d]", html_index), fixed = TRUE)
    testthat::expect_match(definition$render, sprintf("type==='filter')return row[%d]", plain_index), fixed = TRUE)
    testthat::expect_identical(shown[[column]], data[[column]])
    testthat::expect_identical(shown[[log_column]], data[[log_column]])
    log_order = order(shown[[definition$orderData + 1L]])
    testthat::expect_identical(log_order, order(data[[log_column]]))
    testthat::expect_false(identical(log_order, order(data[[column]])))
  }
  testthat::expect_identical(shown$name[order(shown$log_p_value)][1:4],
    c("Deep bound", "Legacy", "Floored", "Shallow bound"))
  testthat::expect_identical(serialize(data, NULL), before)
})

testthat::test_that("probability headers fall back to numeric sorting without a usable log", {
  app_environment = shiny_app_environment()
  expression = data.frame(gene_name = c("GENE1", "GENE2", "GENE3"), pvalue = c(0.1, 4.89e-15, 0))
  table = app_environment$probability_datatable(expression)
  definition = probability_order_definition(table, "pvalue")
  testthat::expect_null(definition$orderData)
  testthat::expect_identical(definition$type, "num")
  testthat::expect_identical(table$x$data$pvalue, expression$pvalue)
  expression$log_pvalue = c(log(0.1), log(4.89e-15), -800)
  table = app_environment$probability_datatable(expression)
  definition = probability_order_definition(table, "pvalue")
  testthat::expect_identical(names(table$x$data)[definition$orderData + 1L], "log_pvalue")
  partial = probability_fixture()
  partial$log_p_value[1] = NA_real_
  partial$log_fdr[2] = NaN
  table = app_environment$probability_datatable(partial)
  testthat::expect_null(probability_order_definition(table, "p_value")$orderData)
  testthat::expect_null(probability_order_definition(table, "fdr")$orderData)
  testthat::expect_identical(partial$p_value, probability_fixture()$p_value)
})

testthat::test_that("empty tables keep valid probability ordering definitions", {
  app_environment = shiny_app_environment()
  table = app_environment$probability_datatable(probability_fixture()[FALSE, ])
  for (column in c("p_value", "fdr")) {
    definition = probability_order_definition(table, column)
    testthat::expect_identical(names(table$x$data)[definition$orderData + 1L], paste0("log_", column))
  }
  table = app_environment$probability_datatable(data.frame(gene_name = character(), pvalue = numeric()))
  testthat::expect_identical(nrow(table$x$data), 0L)
  testthat::expect_null(probability_order_definition(table, "pvalue")$orderData)
})
