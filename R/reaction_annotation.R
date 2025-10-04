#' Biochemical Reaction Annotations (KEGG-Based)
#'
#' A structured annotation table of biochemical reactions, primarily derived from KEGG reaction identifiers. Each entry maps a unique reaction ID to a human-readable reaction name or formula, representing enzymatic or spontaneous transformations of biological molecules.
#' This dataset supports annotation of nodes labeled as "Reaction" in the tKOI framework.
#'
#' @format A data frame with 24659 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the KEGG reaction ID (e.g., "R04184").}
#'   \item{name}{The reaction name or stoichiometric equation (e.g., "globoside N-acetylhexosaminohydrolase; ..."). This may be `NA` for reactions without an assigned label or formula.}
#' }
#'
#' @details
#' Reactions describe biochemical events such as substrate-to-product transformations catalyzed by enzymes or occurring spontaneously.
#' This dataset allows enrichment and interpretation of pathway-embedded biochemical steps, enhancing systems-level understanding of metabolic or signaling processes.
#'
#' The `identifier` field enables joining with reaction nodes in tKOI network output, and the `name` column provides interpretability in reporting or visualization.
#'
#' @source KEGG Reaction Database \url{https://www.genome.jp/kegg/reaction/}
#'
#' @examples
#' data(reaction_annotation)
#' subset(reaction_annotation, grepl("glucosyl", name, ignore.case = TRUE))
#'
#' @seealso \code{\link[tkoi]{pwgroup_annotation}}, \code{\link[tkoi]{pathway_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset reaction kegg annotation biochemistry
"reaction_annotation"


