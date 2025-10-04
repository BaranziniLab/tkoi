#' Pathway Group Annotations (Reactome Interaction Sets)
#'
#' A curated annotation table for grouped pathway interaction events from Reactome.
#' Each entry represents a molecular interaction or signaling complex involved in a pathway context, typically grouped by protein complexes or ligand-receptor interactions.
#' These pathway groupings support detailed annotation of graph nodes in the tKOI network that represent composite signaling units.
#'
#' @format A data frame with 6343 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the Reactome stable identifier (e.g., "reactome:R-HSA-9624419").}
#'   \item{name}{A human-readable description of the interaction group or signaling unit (e.g., "EGFR:EGF-like ligands").}
#' }
#'
#' @details
#' This dataset provides molecular resolution for pathway-associated interactions and groupings (e.g., dimerized receptors, ligand-bound complexes),
#' enhancing the interpretability of nodes labeled as "PwGroup" in tKOI analyses. These groups often represent intermediate or context-specific biochemical states within signaling cascades.
#'
#' The `identifier` field can be used to join with graph nodes, while the `name` field is useful for labeling and summarizing pathway-level results.
#'
#' @source Reactome Pathway Database \url{https://reactome.org}
#'
#' @examples
#' data(pwgroup_annotation)
#' subset(pwgroup_annotation, grepl("EGFR", name))
#'
#' @seealso \code{\link[tkoi]{pathway_annotation}}, \code{\link[tkoi]{reaction_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset pathway interaction reactome annotation
"pwgroup_annotation"

