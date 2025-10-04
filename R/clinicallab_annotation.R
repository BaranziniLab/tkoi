#' Clinical Laboratory Measurement Annotations (LOINC)
#'
#' A structured annotation table of clinical laboratory tests and biomarkers, derived from the LOINC (Logical Observation Identifiers Names and Codes) system.
#' This dataset maps LOINC identifiers to human-readable test names and is used within the tKOI pipeline to annotate graph nodes representing clinical laboratory features.
#'
#' @format A data frame with multiple rows and 2 columns:
#' \describe{
#'   \item{identifier}{A character string representing the LOINC code (e.g., "2251-7").}
#'   \item{name}{A human-readable description of the laboratory test or measurement, including analyte, method, specimen type, timing, and scale (e.g., "Estriol:MCnc:Pt:Ser/Plas:Qn").}
#' }
#'
#' @details
#' This dataset supports the integration of clinical laboratory data into biomedical knowledge graphs by providing interpretive labels for LOINC-coded test results.
#' These annotations are essential for associating molecular or phenotypic features with clinical diagnostics in translational research or precision medicine applications.
#'
#' The dataset is intended to be joined to graph nodes with `identifier` values corresponding to LOINC codes.
#'
#' @source LOINC — Logical Observation Identifiers Names and Codes \url{https://loinc.org/}
#'
#' @examples
#' data(clinicallab_annotation)
#' subset(clinicallab_annotation, grepl("Estriol", name))
#'
#' @seealso \code{\link[tkoi]{run_tkoi}}, \code{\link[tkoi]{disease_annotation}}, \code{\link[tkoi]{compound_annotation}}
#'
#' @keywords dataset clinical lab LOINC annotation
"clinicallab_annotation"


