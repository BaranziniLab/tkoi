# Shape of the shipped datasets and removal of the encryption artifacts.
# None of these tests touch tkoi::tkoi_net (see test-extended.R for that).

annotation_datasets = c(
  "anatomy", "celltype", "clinicallab", "complex", "compound", "disease", "ec", "go", "mirna",
  "pathway", "protein", "proteindomain", "proteinfamily", "pwgroup", "reaction"
)

node_types = c(
  "Anatomy", "BiologicalProcess", "CellType", "CellularComponent", "ClinicalLab", "Complex",
  "Compound", "Disease", "EC", "Gene", "MiRNA", "MolecularFunction", "Pathway", "Protein",
  "ProteinDomain", "ProteinFamily", "PwGroup", "Reaction"
)

dependency_names = function(field) {
  if (is.null(field) || is.na(field)) {
    return(character())
  }
  entries = trimws(strsplit(field, ",", fixed = TRUE)[[1]])
  entries = sub("\\s*\\(.*$", "", entries)
  entries[nzchar(entries)]
}

test_that("genes has the documented columns, types and size", {
  genes = tkoi::genes

  expect_s3_class(genes, "data.frame")
  expect_identical(names(genes), c("id", "identifier", "ensembl", "name", "degree"))
  expect_identical(nrow(genes), 17569L)
  expect_type(genes$id, "character")
  expect_type(genes$identifier, "integer")
  expect_type(genes$ensembl, "character")
  expect_type(genes$name, "character")
  expect_type(genes$degree, "integer")
})

test_that("genes has one row per unique, non-missing node id", {
  genes = tkoi::genes

  expect_false(anyNA(genes$id))
  expect_true(all(nzchar(genes$id)))
  expect_identical(anyDuplicated(genes$id), 0L)
  expect_false(anyNA(genes$identifier))
  expect_false(anyNA(genes$degree))
  expect_true(all(genes$degree >= 1L))
})

test_that("genes Ensembl IDs are Ensembl gene IDs or blank", {
  ensembl = tkoi::genes$ensembl

  expect_false(anyNA(ensembl))
  expect_true(all(grepl("^ENSG[0-9]+$", ensembl) | ensembl == ""))
  # Only "a few" genes lack an Ensembl ID.
  expect_lt(sum(ensembl == ""), 0.01 * length(ensembl))
})

# The identifier forms of human_metabolites, as documented on its help page:
# the same forms as Compound node identifiers.
metabolite_forms = c(
  inchikey = "^inchikey:[A-Z]{14}-[A-Z]{10}-[A-Z]$",
  chebi = "^CHEBI:[0-9]+$",
  chembl = "^chembl[.]compound:CHEMBL[0-9]+$"
)

test_that("human_metabolites is a character vector of the documented size", {
  metabolites = tkoi::human_metabolites

  expect_type(metabolites, "character")
  expect_null(dim(metabolites))
  expect_length(metabolites, 231466L)
  expect_false(anyNA(metabolites))
  expect_false(anyDuplicated(metabolites) > 0)
})

test_that("every human_metabolites entry is written like a Compound node identifier", {
  metabolites = tkoi::human_metabolites
  matches = vapply(metabolite_forms, grepl, logical(length(metabolites)), x = metabolites)

  expect_true(all(rowSums(matches) == 1L))
  expect_identical(unname(colSums(matches)), c(217895, 13562, 9))
  # The 1.0.0 formatting problems are gone.
  expect_false(any(grepl("^CHEBI:CHEBI:|^CHEBI:NA$|^inchikey:$", metabolites)))
  # The examples on the help page.
  examples = c("inchikey:BRMWTNUJHUMWMS-LURJTMIESA-N", "CHEBI:50599", "chembl.compound:CHEMBL4159192")
  expect_true(all(examples %in% metabolites))
})

test_that("human metabolites match Compound identifiers of every form", {
  compound_forms = "^inchikey:[A-Z]{14}-[A-Z]{10}-[A-Z]$|^CHEBI:[0-9]+$|^chembl[.]compound:CHEMBL[0-9]+$"
  metabolites = tkoi::human_metabolites
  compounds = tkoi::compound_annotation$identifier

  # compound_annotation has one row per Compound node, in the documented forms.
  expect_gt(mean(grepl(compound_forms, compounds)), 0.999)
  expect_identical(
    as.vector(table(sub(":.*$", "", compounds))[c("inchikey", "CHEBI", "chembl.compound")]),
    c(530612L, 22652L, 1262L)
  )

  # 24,099 InChIKey and 45 ChEBI compounds are human metabolites.
  matched = compounds[compounds %in% metabolites]
  expect_identical(as.vector(table(sub(":.*$", "", matched))[c("inchikey", "CHEBI")]), c(24099L, 45L))
})

test_that("every annotation dataset is a data frame with a character identifier column", {
  for (dataset in annotation_datasets) {
    annotation = getExportedValue("tkoi", paste0(dataset, "_annotation"))
    expect_s3_class(annotation, "data.frame")
    expect_true("identifier" %in% names(annotation), info = dataset)
    expect_type(annotation$identifier, "character")
    expect_gt(nrow(annotation), 0)
    expect_false(anyNA(annotation$identifier), info = dataset)
  }
})

test_that("the package ships exactly the documented datasets", {
  items = utils::data(package = "tkoi")$results[, "Item"]

  expect_setequal(
    items,
    c(paste0(annotation_datasets, "_annotation"), "genes", "human_metabolites", "tkoi_net")
  )
})

test_that("go_annotation covers the three GO namespaces", {
  go = tkoi::go_annotation

  expect_true(all(c("identifier", "name", "namespace", "definition") %in% names(go)))
  expect_true(all(c("biological_process", "molecular_function", "cellular_component") %in% go$namespace))
  expect_true(all(startsWith(go$identifier, "GO:")))
})

test_that("the encryption artifacts are gone", {
  expect_false(exists("network_attributes", envir = asNamespace("tkoi"), inherits = FALSE))
  expect_false(exists("network_attributes", envir = getNamespaceInfo("tkoi", "lazydata"), inherits = FALSE))
  expect_false("network_attributes" %in% utils::data(package = "tkoi")$results[, "Item"])
})

test_that("DESCRIPTION no longer depends on sodium", {
  description = read.dcf(system.file("DESCRIPTION", package = "tkoi"))
  field = function(name) if (name %in% colnames(description)) description[1, name] else NA_character_

  imports = dependency_names(field("Imports"))
  expect_true(length(imports) > 0)
  expect_true("igraph" %in% imports)
  expect_false("sodium" %in% imports)
  for (name in c("Depends", "LinkingTo", "Suggests")) {
    expect_false("sodium" %in% dependency_names(field(name)), info = name)
  }
  expect_false("sodium" %in% dependency_names(utils::packageDescription("tkoi")$Imports))
})

test_that("the toy fixture is an undirected igraph with the tkoi_net vertex attributes", {
  network = toy_network()

  expect_true(igraph::is_igraph(network))
  expect_false(igraph::is_directed(network))
  expect_equal(igraph::vcount(network), 5260)
  expect_true(all(c("name", "identifier", "labels", "degree") %in% igraph::vertex_attr_names(network)))

  vertex_names = igraph::V(network)$name
  expect_type(vertex_names, "character")
  expect_identical(anyDuplicated(vertex_names), 0L)
  expect_type(igraph::V(network)$identifier, "character")
  expect_false(anyNA(igraph::V(network)$identifier))
  expect_type(igraph::V(network)$degree, "double")
  expect_identical(sum(igraph::which_loop(network)), 4L)
})

test_that("the toy fixture labels use the Neo4j label form of the documented node types", {
  labels = igraph::V(toy_network())$labels

  expect_type(labels, "character")
  expect_true(all(grepl("^\\['[A-Za-z]+'\\]$", labels)))
  types = sort(unique(gsub("[][']", "", labels)), method = "radix")
  expect_identical(
    types,
    c("Anatomy", "BiologicalProcess", "CellType", "Compound", "Disease", "Gene", "MolecularFunction", "Pathway")
  )
  expect_true(all(types %in% node_types))
})

test_that("toy fixture degrees are full-network degrees", {
  network = toy_network()
  degree = igraph::V(network)$degree

  # An induced subgraph can only lose edges.
  expect_true(all(degree == round(degree)))
  expect_true(all(degree >= igraph::degree(network)))

  # Gene vertices carry the degree and Entrez ID listed in tkoi::genes.
  universe = toy_universe(network)
  expect_gt(nrow(universe), 0)
  rows = match(universe$id, igraph::V(network)$name)
  expect_false(anyNA(rows))
  expect_equal(degree[rows], as.numeric(universe$degree))
  expect_identical(igraph::V(network)$identifier[rows], as.character(universe$identifier))
  expect_true(all(igraph::V(network)$labels[rows] == "['Gene']"))
})
