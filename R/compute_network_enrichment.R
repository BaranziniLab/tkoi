#' Compute Z-Score and P-Value for Node Pagerank
#'
#' This function calculates the z-score and p-value for a node's pagerank
#' by comparing it to a set of permutation values. It returns a data frame
#' containing the z-score (`beta`) and p-value (`p_value`).
#'
#' @param node Either a named list / one-row data frame with `pagerank` and
#'   permutation values whose names start with `"perm"` (e.g. `perm.1`,
#'   `perm.2`), or a data frame with one row per node in that layout, such as
#'   the `pagerank_data` slot of a [run_tkoi()] result with
#'   `keep_permutations = TRUE`. A list of equal-length vectors is treated as
#'   such a data frame. `node` needs a `pagerank` value and at least two
#'   `perm*` values, or an error is raised.
#'
#' @return A data frame with one row per node and columns:
#'   - `beta`: The computed z-score.
#'   - `p_value`: The one-tailed p-value derived from the z-score.
#'
#' @details
#' The function calculates the z-score as:
#' \deqn{z = \frac{\text{pagerank} -
#' \text{mean}(\text{perm_values})}{\text{sd}(\text{perm_values})}}
#' (\code{NaN} when the permutation values have no spread, as in
#' \code{\link{run_tkoi}}), and the p-value is calculated as the survival
#' function of the z-score:
#' \deqn{p = 1 - \Phi(z)}
#' where \eqn{\Phi} is the cumulative distribution function of the standard normal distribution.
#'
#' @examples
#' compute_network_enrichment(list(pagerank = 0.3, perm.1 = 0.1, perm.2 = 0.2, perm.3 = 0.15))
#'
#' nodes = data.frame(pagerank = c(0.3, 0.1), perm.1 = c(0.1, 0.1), perm.2 = c(0.2, 0.12))
#' compute_network_enrichment(nodes)
#'
#' @export
compute_network_enrichment = function(node) {
  # A matrix, or a list of equal-length vectors, is a table of nodes.
  if (is.matrix(node) || (is.list(node) && !is.data.frame(node) && any(lengths(node) > 1))) {
    node = as.data.frame(node)
  }
  if (is.null(names(node)) || !"pagerank" %in% names(node) ||
    length(grep("^perm", names(node))) < 2) {
    stop("`node` needs a `pagerank` value and at least two `perm*` values.", call. = FALSE)
  }

  if (is.data.frame(node)) {
    perm_values = as.matrix(node[, grep("^perm", names(node)), drop = FALSE])
    pagerank = node[["pagerank"]]
    perm_mean = rowMeans(perm_values)
    perm_sd = sqrt(rowSums((perm_values - perm_mean)^2) / (ncol(perm_values) - 1))
  } else {
    pagerank = as.numeric(node[["pagerank"]])
    perm_values = as.numeric(unlist(node[grep("^perm", names(node))]))
    perm_mean = mean(perm_values)
    perm_sd = stats::sd(perm_values)
  }

  z_score = (pagerank - perm_mean) / perm_sd
  # A constant null cannot score a node (as in run_tkoi()).
  z_score[!is.na(perm_sd) & perm_sd == 0] = NaN
  p_value = exp(stats::pnorm(z_score, lower.tail = FALSE, log.p = TRUE))
  data.frame(beta = z_score, p_value = p_value)
}
