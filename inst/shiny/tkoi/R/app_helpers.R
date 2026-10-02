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
format_result_table = function(data) {
  data = data[, setdiff(names(data), c("node_type", "direct_links", "valency")), drop = FALSE]
  data = dplyr::select(data,
    dplyr::any_of(c("node_id", "identifier", "name", "beta", "p_value", "fdr")),
    dplyr::everything()
  )
  dplyr::mutate(data, dplyr::across(
    dplyr::any_of(c("pagerank", "beta", "p_value", "fdr")),
    function(value) formatC(value, format = "e", digits = 3)
  ))
}
