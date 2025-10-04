#' Protein Domain Annotations (Pfam)
#'
#' A curated annotation table of protein domains based on Pfam identifiers. This dataset provides structured mappings from Pfam domain IDs to domain names,
#' supporting the annotation and interpretation of domain-level biological features in knowledge graphs used by the tKOI framework.
#'
#' @format A data frame with 14193 rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the Pfam domain identifier (e.g., "PF19543").}
#'   \item{name}{The name or description of the domain, often including structural or functional attributes (e.g., "Glycoside hydrolase 123, N-terminal domain").}
#' }
#'
#' @details
#' Protein domains represent conserved functional or structural units within proteins. This dataset enables enrichment and network-based inference at the domain level.
#' The `identifier` column can be joined to graph nodes representing domains in tKOI output, while the `name` column adds interpretability for visualization and reporting.
#'
#' Domains include a variety of structural motifs, enzyme catalytic regions, and uncharacterized conserved segments (e.g., DUFs).
#'
#' @source Pfam Database \url{https://pfam.xfam.org/}
#'
#' @examples
#' data(proteindomain_annotation)
#' subset(proteindomain_annotation, grepl("hydrolase", name, ignore.case = TRUE))
#'
#' @seealso \code{\link[tkoi]{protein_annotation}}, \code{\link[tkoi]{proteinfamily_annotation}}, \code{\link[tkoi]{run_tkoi}}
#'
#' @keywords dataset protein domain pfam annotation
"proteindomain_annotation"


