#' human_metabolites: Parsed Human Metabolite Metadata from HMDB
#'
#' A data frame containing key chemical identifiers for metabolites curated
#' from the Human Metabolome Database (HMDB). This dataset was constructed
#' by parsing the official HMDB XML export and selecting metabolite-level
#' identifiers relevant for downstream bioinformatics and network biology analyses.
#'
#' ## Description
#' Each row corresponds to a unique metabolite entry from HMDB. The following
#' columns are included:
#'
#' - `name`: The primary name of the metabolite as recorded in HMDB
#' - `inchikey`: The standard InChIKey chemical structure identifier
#' - `chebi_id`: The corresponding ChEBI identifier (when available), prefixed with `"CHEBI:"`
#' - `chembl_id`: The corresponding ChEMBL compound ID (when available)
#'
#' All values are character strings. Missing values are represented as `NA`.
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
#' and `<chembl_id>`. ChEBI identifiers were post-processed to include the standard
#' `"CHEBI:"` prefix when present.
#'
#' @format A data frame with `r nrow(human_metabolites)` rows and 4 columns:
#' \describe{
#'   \item{name}{Metabolite name (character)}
#'   \item{inchikey}{Standard InChIKey identifier (character)}
#'   \item{chebi_id}{ChEBI ID with "CHEBI:" prefix, or NA (character)}
#'   \item{chembl_id}{ChEMBL compound ID, or NA (character)}
#' }
#'
#' @examples
#' data(human_metabolites)
#' head(human_metabolites)
"human_metabolites"
