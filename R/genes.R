#' Gene Metadata
#'
#' A dataset containing information about genes, including identifiers, names, and degree in the network.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{id}{Character. Unique identifier for the gene in the network.}
#'   \item{identifier}{Integer. Internal identifier for the gene.}
#'   \item{ensembl}{Character. Ensembl gene identifier.}
#'   \item{name}{Character. Gene symbol or common name.}
#'   \item{degree}{Integer. Degree of the gene node in the network, representing the number of connections.}
#' }
#'
#' @details
#' This dataset is essential for integrating gene-level information with biological networks and functional studies.
#'
#' @examples
#' data(genes)
#' head(genes)
"genes"
