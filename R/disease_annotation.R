#' Disease Ontology (DO) Annotations
#'
#' A structured annotation table mapping Disease Ontology (DOID) identifiers to human-readable disease names and definitions.
#' This dataset is used in the tKOI framework to annotate network nodes associated with disease terms, enabling interpretation of enriched or prioritized disease-related features in biological networks.
#'
#' @format A data frame with 10002 rows and 3 columns:
#' \describe{
#'   \item{identifier}{A character string representing the Disease Ontology ID (e.g., "DOID:0040003").}
#'   \item{name}{The name of the disease (e.g., "benzylpenicillin allergy").}
#'   \item{definition}{A detailed textual definition of the disease, including etiology, affected systems, or allergy triggers, often with references to biomedical sources.}
#' }
#'
#' @details
#' Disease Ontology (DO) provides a standardized, ontology-based representation of human diseases. This dataset allows integration of disease concepts into biological graphs and supports downstream analysis of disease associations.
#'
#' The `identifier` column serves as a join key to nodes of type "Disease" in tKOI networks, while the `definition` column provides biological and clinical context.
#'
#' @source Disease Ontology \url{https://disease-ontology.org/}
#'
#' @examples
#' data(disease_annotation)
#' subset(disease_annotation, grepl("allergy", name))
#'
#' @seealso \code{\link[tkoi]{clinicallab_annotation}}, \code{\link[tkoi]{compound_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset disease ontology annotation DOID
"disease_annotation"


