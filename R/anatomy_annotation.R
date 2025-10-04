#' Anatomical Entity Annotations from UBERON
#'
#' A curated annotation table of anatomical entities sourced from the UBERON ontology. This dataset provides standard identifiers,
#' human-readable names, and textual definitions for anatomical structures. It is used within the tKOI analysis pipeline to enrich
#' network nodes of type "Anatomy" with biological context and ontology-based metadata.
#'
#' @format A data frame with 15,668 rows and 3 columns:
#' \describe{
#'   \item{identifier}{A character string representing the UBERON ontology identifier (e.g., "UBERON:0000002").}
#'   \item{name}{The human-readable name of the anatomical structure (e.g., "uterine cervix").}
#'   \item{definition}{A character string describing the biological definition or functional role of the anatomical entity.}
#' }
#'
#' @details
#' The dataset supports the annotation and interpretation of graph nodes labeled as anatomical concepts in biological networks.
#' It can be joined to results using the `identifier` column to enhance interpretability of enriched or prioritized nodes.
#'
#' @source UBERON Ontology \url{http://uberon.github.io/}
#'
#' @examples
#' data(anatomy_annotation)
#' head(anatomy_annotation)
#'
#' @seealso \code{\link[tkoi]{run_tkoi}}, \code{\link[tkoi]{celltype_annotation}}, \code{\link[tkoi]{disease_annotation}}
#'
#' @keywords dataset anatomy ontology
"anatomy_annotation"


