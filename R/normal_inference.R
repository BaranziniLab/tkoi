#' Stable Normal Upper-Tail Inference
#'
#' Compute normal upper-tail probabilities from observed scores and saved null
#' means and standard deviations. No randomizations are performed.
#' @param observed Numeric observed score vector.
#' @param null_mean Numeric null mean vector of the same length.
#' @param null_sd Numeric null standard deviation vector of the same length.
#' @return A data frame with `beta`, positive numeric `p_value`, unfloored
#'   `log_p_value` (natural logarithm), `p_value_bounded`, and `inference_status`.
#'   Numeric probabilities below `.Machine$double.xmin` are represented by that
#'   bound. Zero-spread or invalid nulls are untestable and have missing results.
#' @export
tkoi_normal_inference = function(observed, null_mean, null_sd) {
  values = list(observed, null_mean, null_sd)
  if (!all(vapply(values, is.numeric, logical(1))) ||
      length(unique(lengths(values))) != 1L) {
    stop("Observed scores and null moments must be numeric vectors of equal length.", call. = FALSE)
  }
  n = length(observed)
  status = rep("testable", n)
  status[!is.finite(null_mean) | !is.finite(null_sd) | null_sd < 0] = "invalid_null"
  status[is.finite(null_mean) & is.finite(null_sd) & null_sd == 0] = "constant_null"
  status[!is.finite(observed)] = "nonfinite_observed"
  beta = rep(NaN, n)
  ok = status == "testable"
  beta[ok] = (observed[ok] - null_mean[ok]) / null_sd[ok]
  status[ok & !is.finite(beta)] = "nonfinite_beta"
  beta[status != "testable"] = NaN
  log_p = stats::pnorm(beta, lower.tail = FALSE, log.p = TRUE)
  data.frame(
    beta = beta,
    p_value = .tkoi_bounded_probability(log_p),
    log_p_value = log_p,
    p_value_bounded = log_p < log(.Machine$double.xmin),
    inference_status = status
  )
}

# Keep inference in log space; bound only the numeric representation.
.tkoi_bounded_probability = function(log_p) {
  if (!is.numeric(log_p) || any(log_p > 0, na.rm = TRUE)) {
    stop("Log probabilities must be numeric and no greater than zero.", call. = FALSE)
  }
  bounded = !is.na(log_p) & log_p < log(.Machine$double.xmin)
  probability = exp(pmax(log_p, log(.Machine$double.xmin)))
  probability[bounded] = .Machine$double.xmin
  probability
}

# Match stats::p.adjust(method = "BH"), including its default NA semantics:
# the default n is evaluated after removing unavailable probabilities.
.tkoi_bh_log = function(log_p, n = length(log_p)) {
  if (!is.numeric(log_p) || any(log_p > 0, na.rm = TRUE)) {
    stop("Log probabilities must be numeric and no greater than zero.", call. = FALSE)
  }
  out = log_p
  keep = !is.na(log_p)
  if (!all(keep)) log_p = log_p[keep]
  size = length(log_p)
  if (length(n) != 1L || !is.finite(n) || n < size || n != floor(n)) {
    stop("`n` must be a whole number at least the number of available probabilities.", call. = FALSE)
  }
  if (size == 0L || n <= 1L) return(out)
  ord = order(log_p, decreasing = TRUE)
  adjusted = pmin(0, cummin(log_p[ord] + log(n / (size:1L))))
  out[keep] = pmax(log_p, adjusted[order(ord)])
  out
}
