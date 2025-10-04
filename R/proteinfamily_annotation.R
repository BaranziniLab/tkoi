#' Protein Family Annotations (Pfam Clans)
#'
#' A curated annotation table of protein families, primarily derived from Pfam clan classifications.
#' Each entry maps a protein family identifier to a short descriptive name, supporting the annotation of graph nodes
#' that represent protein families in the tKOI knowledge network.
#'
#' @format A data frame with 659 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the Pfam clan or protein family identifier (e.g., "CL0163").}
#'   \item{name}{The short name of the protein family (e.g., "Calcineurin"). These often reflect structural, evolutionary, or functional groupings of related protein domains.}
#' }
#'
#' @details
#' Protein families group related domains based on common evolutionary origins, structure, or function.
#' This dataset allows high-level functional grouping of nodes in biological networks and supports annotation of gene products in systems biology and bioinformatics analyses.
#'
#' The `identifier` field serves as a join key for nodes of type "ProteinFamily" in the tKOI framework.
#'
#' @source Pfam Clans \url{https://pfam.xfam.org/clan}
#'
#' @examples
#' data(proteinfamily_annotation)
#' subset(proteinfamily_annotation, grepl("Calcineurin", name))
#'
#' @seealso \code{\link[tkoi]{proteindomain_annotation}}, \code{\link[tkoi]{protein_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset protein family clan annotation
"proteinfamily_annotation"


