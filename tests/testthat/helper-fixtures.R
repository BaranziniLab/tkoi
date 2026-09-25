# Small inputs shared by the tests.
#
# toy_network() is a 5,260-vertex induced subgraph of tkoi_net (250 random
# genes, their GO / pathway / disease / cell type / anatomy / compound
# neighbors, and the edges among them). It keeps the tkoi_net vertex
# attributes (name, identifier, labels, degree), so every pipeline step runs
# on it in well under a second.

toy_network = local({
  cache = NULL
  function() {
    if (is.null(cache)) {
      cache <<- readRDS(test_path("fixtures", "toy_network.rds"))
    }
    cache
  }
})

# Genes of tkoi::genes that are vertices of the toy network.
toy_universe = function(network = toy_network()) {
  universe = as.data.frame(tkoi::genes)
  universe[universe$id %in% igraph::V(network)$name, ]
}

# A differential expression table on toy-network genes. Some genes fail the
# p-value or log fold change filters, one gene is duplicated, and two IDs do
# not exist in tkoi::genes.
toy_expression = function(n_genes = 80, seed = 1) {
  universe = toy_universe()
  universe = universe[!duplicated(universe$ensembl), ]
  withr::with_seed(seed, {
    rows = sample(nrow(universe), n_genes)
    expression = data.frame(
      gene_name = universe$ensembl[rows],
      logfc = round(stats::rnorm(n_genes, sd = 1.2), 3),
      pvalue = signif(stats::runif(n_genes)^3, 3),
      stringsAsFactors = FALSE
    )
  })
  extra = data.frame(
    gene_name = c(expression$gene_name[1], "ENSG_NOT_A_GENE", "ENSG00000000000"),
    logfc = c(-5, 2, 2),
    pvalue = c(1e-10, 1e-10, 1e-10),
    stringsAsFactors = FALSE
  )
  rbind(expression, extra)
}

# run_tkoi() on the toy network with quiet, fast defaults.
toy_run = function(..., seed = 42, n_permutation = 10, n_cores = 1) {
  withr::with_seed(seed, run_tkoi(
    expression_data = toy_expression(),
    subnetwork = toy_network(),
    n_permutation = n_permutation,
    n_cores = n_cores,
    verbose = FALSE,
    ...
  ))
}
