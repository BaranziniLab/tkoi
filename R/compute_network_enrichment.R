#' Compute Z-Score and P-Value for Node Pagerank
#'
#' This function calculates the z-score and p-value for a node's pagerank
#' by comparing it to a set of permutation values. It returns a data frame
#' containing the z-score (`beta`) and p-value (`p_value`).
#'
#' @param node A named list or data structure containing:
#'   - `pagerank`: A numeric value representing the pagerank of the node.
#'   - Permutation values: One or more numeric values with names matching the
#'     pattern `"^perm.*"` (e.g., `perm1`, `perm2`), used for calculating the z-score.
#'
#' @return A data frame with the following columns:
#'   - `beta`: The computed z-score.
#'   - `p_value`: The one-tailed p-value derived from the z-score.
#'
#' @details
#' The function calculates the z-score as:
#' \deqn{z = \frac{\text{pagerank} - \text{mean}(\text{perm_values})}{\text{sd}(\text{perm_values})}}
#' and the p-value is calculated as the survival function of the z-score:
#' \deqn{p = 1 - \Phi(z)}
#' where \eqn{\Phi} is the cumulative distribution function of the standard normal distribution.
#'
#' @export
compute_network_enrichment = function(node){
  pagerank = node[["pagerank"]]
  perm_values = as.numeric(node[grep("^perm", names(node))])

  z_score = (pagerank - mean(perm_values)) / sd(perm_values)
  p_value = exp(pnorm(z_score, lower.tail = FALSE, log.p = TRUE))

  results = data.frame(
    beta = z_score,
    p_value = p_value
  )
  return(results)
}
