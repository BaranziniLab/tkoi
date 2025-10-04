#' Human miRNA Annotations (miRBase-style)
#'
#' A curated dataset mapping human mature miRNA identifiers to their corresponding miRBase precursor accessions.
#'
#' @format A data frame with two columns and multiple rows (1 per miRNA):
#' \describe{
#'   \item{identifier}{Character string representing the mature miRNA ID, including the species prefix (e.g., `"hsa-miR-6129"`).}
#'   \item{accession}{miRBase accession number for the precursor miRNA hairpin (e.g., `"MI0021274"`).}
#' }
#'
#' @details
#' This dataset enables annotation and cross-referencing of mature miRNA names with miRBase precursor IDs.
#' It is useful for linking results from miRNA expression studies with external reference databases such as miRBase.
#'
#' The `identifier` field typically corresponds to results from expression matrices, differential analysis tools, or network analysis pipelines.
#'
#' @examples
#' data(mirna_annotation)
#' head(mirna_annotation)
#' subset(mirna_annotation, grepl("miR-23", identifier))
#'
#' @source Extracted from a parsed version of miRBase release data.
#'
#' @seealso \code{\link[tkoi]{run_tkoi}}, \code{\link[tkoi]{reaction_annotation}}
#'
#' @keywords dataset miRNA annotation miRBase
"mirna_annotation"
