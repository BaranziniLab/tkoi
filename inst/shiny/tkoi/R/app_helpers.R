# File import and presentation only. All analysis belongs to tkoi.
read_expression_data = function(upload, example = example_data) {
  if (is.null(upload)) {
    data = as.data.frame(example)
  } else {
    extension = tolower(tools::file_ext(upload$name))
    data = switch(extension,
      csv = data.table::fread(upload$datapath, data.table = FALSE),
      tsv = data.table::fread(upload$datapath, data.table = FALSE),
      xlsx = as.data.frame(readxl::read_xlsx(upload$datapath)),
      stop("Upload a .csv, .tsv, or .xlsx file.", call. = FALSE)
    )
  }
  missing = setdiff(c("gene_name", "logfc", "pvalue"), names(data))
  if (length(missing)) {
    stop("Missing expression columns: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(data$logfc) || !is.numeric(data$pvalue)) {
    stop("The logfc and pvalue columns must be numeric.", call. = FALSE)
  }
  data
}
# Keep numeric probabilities untouched. Only the trusted display strings contain HTML.
probability_columns = function(data) intersect(c("p_value", "fdr", "pvalue"), names(data))
probability_log_column = function(column) {
  switch(column, p_value = "log_p_value", fdr = "log_fdr", pvalue = "log_pvalue")
}
# Zero-based index of the unfloored log column that ranks `column`, or NULL when
# the log is absent, non-numeric or missing where the probability is available.
probability_order_index = function(display, column) {
  log_column = probability_log_column(column)
  logs = display[[log_column]]
  if (!is.numeric(logs) || any(is.na(logs) != is.na(display[[column]]))) return(NULL)
  match(log_column, names(display)) - 1L
}
probability_labels = function(data, column, format = "plain") {
  log_column = probability_log_column(column)
  bound_column = paste0(column, "_bounded")
  values = data[[column]]
  logs = if (log_column %in% names(data)) data[[log_column]] else NULL
  bounded = if (bound_column %in% names(data)) data[[bound_column]] else NULL
  tkoi::tkoi_format_probability(values, log_p = logs, bounded = bounded,
    format = format, digits = 3)
}
format_result_table = function(data) {
  data = data[, setdiff(names(data), c("node_type", "direct_links", "valency")), drop = FALSE]
  data = dplyr::select(data,
    dplyr::any_of(c("node_id", "identifier", "name", "beta", "p_value", "fdr")),
    dplyr::everything()
  )
  dplyr::mutate(data, dplyr::across(
    dplyr::any_of(c("pagerank", "beta")),
    function(value) formatC(value, format = "e", digits = 3)
  ))
}
probability_datatable = function(data, page_length = 10) {
  display = format_result_table(data)
  probabilities = probability_columns(display)
  hidden = integer()
  html_hidden = integer()
  definitions = list()
  for (column in probabilities) {
    html_column = paste0(".tkoi_", column, "_html")
    plain_column = paste0(".tkoi_", column, "_plain")
    html_labels = probability_labels(data, column, "html")
    display[[html_column]] = if (length(html_labels)) {
      paste0('<span style="white-space:nowrap">', html_labels, "</span>")
    } else {
      character()
    }
    display[[plain_column]] = probability_labels(data, column, "plain")
    html_index = match(html_column, names(display)) - 1L
    plain_index = match(plain_column, names(display)) - 1L
    hidden = c(hidden, html_index, plain_index)
    html_hidden = c(html_hidden, html_index)
    definition = list(
      targets = match(column, names(display)) - 1L,
      type = "num",
      render = DT::JS(sprintf(
        "function(data,type,row){if(type==='display')return row[%d];if(type==='filter')return row[%d];return data;}",
        html_index, plain_index
      ))
    )
    # Numeric probabilities tie at the double.xmin bound; the unfloored log ranks them.
    # orderData also drives DT's server-side ordering. Without a usable log, sort numerically.
    definition$orderData = probability_order_index(display, column)
    definitions[[length(definitions) + 1L]] = definition
  }
  if (length(hidden)) {
    definitions[[length(definitions) + 1L]] = list(targets = hidden, visible = FALSE)
    definitions[[length(definitions) + 1L]] = list(targets = html_hidden, searchable = FALSE)
  }
  DT::datatable(display,
    options = list(pageLength = page_length, scrollX = TRUE, columnDefs = definitions),
    # Gene/concept names and every source-provided field stay escaped.
    escape = setdiff(seq_along(display), hidden + 1L), rownames = FALSE
  )
}
probability_export_table = function(data) {
  for (column in probability_columns(data)) {
    data[[paste0(column, "_display")]] = probability_labels(data, column, "plain")
  }
  data
}
write_result_worksheet = function(workbook, sheet, data) {
  exported = probability_export_table(data)
  openxlsx::addWorksheet(workbook, sheetName = sheet)
  openxlsx::writeData(workbook, sheet = sheet, x = exported)
  if (!nrow(data)) return(invisible(exported))
  for (column in probability_columns(data)) {
    index = match(column, names(exported))
    openxlsx::addStyle(workbook, sheet, style = openxlsx::createStyle(numFmt = "0.00E+00"),
      rows = seq_len(nrow(data)) + 1L, cols = index, gridExpand = TRUE, stack = TRUE)
    labels = exported[[paste0(column, "_display")]]
    # A numeric cell cannot encode a bound. Its number format displays the bound,
    # while the numeric value, log probability and explicit display column remain available.
    literal = which(!is.na(labels) & (grepl("^(<=|\u2264)", labels) | (!is.na(data[[column]]) & data[[column]] == 0)))
    for (label in unique(labels[literal])) {
      selected = literal[labels[literal] == label]
      style = openxlsx::createStyle(numFmt = paste0('"', gsub('"', '""', label, fixed = TRUE), '"'))
      openxlsx::addStyle(workbook, sheet, style = style, rows = selected + 1L,
        cols = index, gridExpand = TRUE, stack = TRUE)
    }
  }
  invisible(exported)
}
