#' Format Probabilities in Scientific Notation
#'
#' Format probability values without changing their numeric storage or global R
#' printing options. Stored zeros are recovered from available log probabilities;
#' a zero without its logarithm is labelled as unresolved. Values at the floor carry a
#' less-than-or-equal sign when their bound flag or original log value is given.
#' @param x Numeric probabilities.
#' @param log_p Optional unfloored natural log probabilities, aligned with `x`.
#' @param format Either `"plain"` for e notation or `"html"` for a multiplication
#'   sign and an editable HTML superscript exponent.
#' @param digits Number of significant digits, between 1 and 15.
#' @param bounded Optional logical bound flags aligned with `x`.
#' @return A character vector. Unavailable values are displayed as `"NA"`.
#' @export
tkoi_format_probability = function(x, log_p = NULL, format = c("plain", "html"), digits = 3, bounded = NULL) {
  format = match.arg(format)
  if (!is.numeric(x) || length(digits) != 1L || !is.finite(digits) ||
      digits < 1 || digits > 15 || digits != floor(digits)) {
    stop("Use numeric probabilities and 1 to 15 significant digits.", call. = FALSE)
  }
  if ((!is.null(log_p) && (!is.numeric(log_p) || length(log_p) != length(x))) ||
      (!is.null(bounded) && (!is.logical(bounded) || length(bounded) != length(x)))) {
    stop("Log probabilities and bound flags must align with `x`.", call. = FALSE)
  }
  if (is.null(bounded)) bounded = rep(FALSE, length(x))
  if (!is.null(log_p)) {
    bounded = bounded | (!is.na(log_p) & log_p < log(.Machine$double.xmin))
  }
  bounded[is.na(bounded)] = FALSE
  shown = x
  if (!is.null(log_p)) {
    recover = !is.na(x) & x == 0 & !is.na(log_p) & log_p <= 0
    shown[recover] = .tkoi_bounded_probability(log_p[recover])
  }
  shown[bounded & is.finite(x) & x >= 0 & x <= 1] = .Machine$double.xmin
  unresolved = !is.na(shown) & shown == 0
  valid = is.finite(shown) & shown > 0 & shown <= 1
  out = rep("NA", length(x))
  out[unresolved] = "Unresolved zero (log unavailable)"
  out[valid] = sprintf(sprintf("%%.%de", as.integer(digits) - 1L), shown[valid])
  if (format == "html") {
    parts = strsplit(out[valid], "e", fixed = TRUE)
    out[valid] = vapply(parts, function(z) {
      exponent = as.character(as.integer(z[2L]))
      exponent = sub("^-", "\u2212", exponent)
      paste0(z[1L], " \u00d7 10<sup>", exponent, "</sup>")
    }, character(1))
  }
  prefix = if (format == "html") "&le;" else "<="
  out[valid & bounded] = paste0(prefix, out[valid & bounded])
  names(out) = names(x)
  out
}

#' Prepare a Probability Table for Display
#'
#' Add scientific display columns while retaining all original numeric columns.
#' @param x A data frame with probability columns.
#' @param format Display format passed to [tkoi_format_probability()].
#' @param digits Significant digits.
#' @return A data frame with additional `*_display` character columns.
#' @export
tkoi_probability_table = function(x, format = c("plain", "html"), digits = 3) {
  format = match.arg(format)
  columns = intersect(
    c("p_value", "fdr", "pvalue", "p.adjust", "experimental_pvalue", "tkoi_pvalue", "tkoi_fdr"),
    names(x)
  )
  log_columns = c(p_value = "log_p_value", fdr = "log_fdr", tkoi_pvalue = "tkoi_log_p_value", tkoi_fdr = "tkoi_log_fdr")
  flag_columns = c(
    p_value = "p_value_bounded", fdr = "fdr_bounded",
    tkoi_pvalue = "tkoi_pvalue_bounded", tkoi_fdr = "tkoi_fdr_bounded"
  )
  for (column in columns) {
    log_name = unname(log_columns[column])
    flag_name = unname(flag_columns[column])
    log_p = if (!is.na(log_name) && log_name %in% names(x)) x[[log_name]] else NULL
    bounded = if (!is.na(flag_name) && flag_name %in% names(x)) x[[flag_name]] else NULL
    x[[paste0(column, "_display")]] = tkoi_format_probability(x[[column]], log_p, format, digits, bounded)
  }
  x
}

#' @export
print.tkoi_statistics = function(x, ...) {
  display = tkoi_probability_table(as.data.frame(x), format = "plain")
  shown = grep("_display$", names(display), value = TRUE)
  for (column in shown) display[[sub("_display$", "", column)]] = display[[column]]
  display[shown] = NULL
  print(dplyr::as_tibble(display), ...)
  invisible(x)
}

.tkoi_statistics_table = function(x) {
  class(x) = unique(c("tkoi_statistics", class(x)))
  x
}
