#' Chemical Compound Annotations
#'
#' A curated annotation table for small molecules and chemical compounds used in pharmacological, biochemical, or research contexts.
#' Each entry includes a standardized identifier (e.g., InChIKey or ChEMBL ID) and the corresponding compound name, when available.
#' This dataset is used within the tKOI framework to annotate nodes representing chemical entities or drugs in biological knowledge graphs.
#'
#' @format A data frame with 554526 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the unique compound ID, either as an InChIKey (e.g., "inchikey:NZLJDTKLZIMONR-UHFFFAOYSA-N") or a ChEMBL compound accession (e.g., "chembl.compound:CHEMBL5219790").}
#'   \item{name}{The human-readable compound name, when available (e.g., "ACRIFLAVINE"). Missing values indicate uncharacterized or unnamed entries.}
#' }
#'
#' @details
#' This dataset enables mapping of compound-level features in biomedical networks, including chemical perturbagens, drug candidates, or environmental exposures.
#' It can be joined to tKOI network results using the `identifier` field to enrich nodes of type "Compound" with interpretable names.
#'
#' @source Data compiled from chemical databases such as ChEMBL and InChI registry.
#'
#' @examples
#' data(compound_annotation)
#' subset(compound_annotation, grepl("CETRIMIDE", name))
#'
#' @seealso \code{\link[tkoi]{complex_annotation}}, \code{\link[tkoi]{clinicallab_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset compound drug chemical_annotation
"compound_annotation"


