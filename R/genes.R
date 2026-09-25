#' Gene Metadata
#'
#' The genes of \code{\link{tkoi_net}} that can seed an analysis, with their
#' Ensembl IDs and network degree. \code{\link{run_tkoi}} maps
#' \code{gene_name} to \code{ensembl}, and draws degree-matched replacement
#' genes for the permutation null from this table.
#'
#' @format A data frame (data.table) with 17,569 rows and 5 columns:
#' \describe{
#'   \item{id}{Character. Node ID of the gene in \code{tkoi_net}.}
#'   \item{identifier}{Integer. NCBI Entrez Gene ID.}
#'   \item{ensembl}{Character. Ensembl gene ID. A few genes have no Ensembl
#'     ID (empty string), and a few Ensembl IDs map to two genes.}
#'   \item{name}{Character. Gene symbol.}
#'   \item{degree}{Integer. Degree of the gene node in \code{tkoi_net}.}
#' }
#'
#' @details
#' This dataset is essential for integrating gene-level information with biological networks and
#' functional studies.
#'
#' @examples
#' data(genes)
#' head(genes)
"genes"
