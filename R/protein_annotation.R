#' Protein Annotations (UniProt-Based)
#'
#' A reference annotation table mapping protein identifiers to gene symbols and functional descriptions.
#' The dataset is derived primarily from UniProt and related sources, and is used within the tKOI framework to annotate nodes of type "Protein" in biological knowledge graphs.
#'
#' @format A data frame with 194076 rows and 3 columns:
#' \describe{
#'   \item{identifier}{A character string representing the UniProt protein accession (e.g., "A0A023HJ61").}
#'   \item{name}{The gene symbol or short name for the protein (e.g., "RAB4A"). This may be `NA` for uncharacterized or unnamed proteins.}
#'   \item{description}{A textual description of the protein’s function or classification, including enzyme commission (EC) numbers, if available (e.g., "ATP synthase subunit a").}
#' }
#'
#' @details
#' This dataset supports the interpretation of graph-based analyses by providing biological context for protein nodes.
#' Descriptions are derived from automated or curated functional annotations such as RuleBase and ECO evidence codes.
#'
#' The `identifier` column is intended to be used as a join key for protein-related nodes in the tKOI network output, while the `description` field adds depth to functional interpretation.
#'
#' @source UniProt \url{https://www.uniprot.org/}
#'
#' @examples
#' data(protein_annotation)
#' subset(protein_annotation, grepl("ATP", description, ignore.case = TRUE))
#'
#' @seealso \code{\link[tkoi]{proteinfamily_annotation}}, \code{\link[tkoi]{run_tkoi}}, \code{\link[tkoi]{complex_annotation}}
#'
#' @keywords dataset protein uniprot annotation
"protein_annotation"


