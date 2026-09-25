#' human_metabolites: Parsed Human Metabolite Metadata from HMDB
#'
#' A character vector of compound identifiers for metabolites curated from
#' the Human Metabolome Database (HMDB), written like the \code{identifier}
#' of Compound nodes in \code{\link{tkoi_net}}: InChIKeys (for example
#' \code{"inchikey:BRMWTNUJHUMWMS-LURJTMIESA-N"}), ChEBI IDs (for example
#' \code{"CHEBI:50599"}), and ChEMBL IDs (for example
#' \code{"chembl.compound:CHEMBL4159192"}). \code{\link{run_tkoi}} reports
#' only Compound nodes whose \code{identifier} is in this vector, and these
#' nodes form the multiple-testing family of the Compound table.
#'
#' In tkoi 1.1.0 the identifiers were repaired: ChEBI IDs had a doubled
#' prefix (\code{"CHEBI:CHEBI:50599"}), ChEMBL IDs lacked the \code{"CHEMBL"}
#' prefix, placeholder entries (\code{"CHEBI:NA"}, a bare \code{"inchikey:"})
#' were removed, and duplicates were dropped. As a result, 45 Compound nodes
#' identified by ChEBI ID are now reported in addition to the 24,099
#' InChIKey-identified compounds.
#'
#' ## Use Cases
#' This table is useful for:
#'
#' - Mapping metabolites across chemical databases (InChIKey, ChEBI, ChEMBL)
#' - Filtering knowledge graph nodes to retain only human-relevant compounds
#' - Cross-referencing with metabolomics datasets or annotations
#'
#' ## Construction
#' This object was created by parsing `hmdb_metabolites.xml` using the `xml2` package.
#' Chemical identifiers were extracted using XPath queries for `<inchikey>`, `<chebi_id>`,
#' and `<chembl_id>`, and prefixed with their namespace.
#'
#' @format A character vector of 231,466 unique identifiers (217,895
#'   InChIKeys, 13,562 ChEBI IDs, and 9 ChEMBL IDs).
#'
#' @examples
#' data(human_metabolites)
#' head(human_metabolites)
"human_metabolites"
