#' Pathway Annotations (WikiPathways)
#'
#' A curated annotation table mapping biological pathway identifiers to pathway names. These annotations are derived primarily from WikiPathways and describe metabolic, signaling, and disease-related pathways relevant to cellular and physiological processes.
#' This dataset is used in the tKOI framework to annotate nodes of type "Pathway" in biological knowledge graphs.
#'
#' @format A data frame with 4831 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the pathway ID, typically from WikiPathways (e.g., "WP2875").}
#'   \item{name}{The human-readable name of the pathway (e.g., "Constitutive androstane receptor pathway").}
#' }
#'
#' @details
#' Pathways describe sets of molecular interactions and regulatory events that carry out defined biological functions.
#' This dataset enables biological interpretation and functional enrichment of graph nodes associated with pathway membership or influence.
#'
#' The `identifier` column is used to map pathway node IDs from graph analyses (e.g., \code{\link[tkoi]{run_tkoi}}), while the `name` column provides context for labeling and visualization.
#'
#' @source WikiPathways \url{https://www.wikipathways.org/}
#'
#' @examples
#' data(pathway_annotation)
#' subset(pathway_annotation, grepl("thyroid", name, ignore.case = TRUE))
#'
#' @seealso \code{\link[tkoi]{reaction_annotation}}, \code{\link[tkoi]{disease_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset pathway annotation wikipathways
"pathway_annotation"


