# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

**tKOI** (Transcriptomic Knowledge-Graph Omics Integration) is an R
package for network-aware gene enrichment analysis. It propagates
transcriptomic signals through a human biological knowledge graph via
personalized PageRank, then uses permutation testing to identify
statistically enriched biological concepts (GO terms, diseases, cell
types, pathways).

## Common Commands

``` r
# Install dependencies
devtools::install_deps()

# Rebuild documentation from Roxygen2 comments
devtools::document()

# Load the package locally for development
devtools::load_all()

# Run tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-<name>.R")

# Build package
devtools::build()

# Check package (CRAN-style)
devtools::check()

# Build pkgdown documentation site
pkgdown::build_site()
```

> Do not edit `README.md` directly — it is generated from `README.Rmd`
> via
> [`devtools::build_readme()`](https://devtools.r-lib.org/reference/build_readme.html).

## Architecture

### S4 Object: `tKOIList`

All analysis state lives in a single S4 object defined in
[R/tKOIList.R](R/tKOIList.R). Key slots: - `expression_data` — input
differential expression table - `pagerank_data` — gene-level PageRank
scores after network propagation - `network_summary_statistics` —
per-node enrichment results, split by node type (GO, Disease, CellType,
Pathway, etc.) - `subnetwork` — `igraph` object of the propagated
subnetwork - `gene_enrichment_comparison` — integrated clusterProfiler
GO results

### Main Pipeline: `run_tkoi()`

[R/run_tkoi.R](R/run_tkoi.R) is the primary entry point. Execution
order: 1. Filter genes by `pvalue_threshold` and `logfc_threshold` 2.
Map filtered genes to the knowledge graph via Ensembl IDs 3. Run
personalized PageRank weighted by log fold change 4. Permutation test:
swap each gene with topologically similar substitutes, recompute
PageRank `n_permutation` times 5. Compute Z-scores and FDR per node 6.
Annotate nodes with biological labels from pre-loaded `.rda` datasets

### Pre-loaded Annotation Data

The `data/` directory holds ~20 `.rda` files (GO, disease, cell type,
pathway, reaction, anatomy, compound, protein annotations). These are
loaded automatically by the package and are the primary reference tables
used to annotate network nodes after PageRank.

### Knowledge Graph

The biological network is encrypted (via `sodium`) and bundled inside
the package. It is decrypted at runtime in
[`run_tkoi()`](reference/run_tkoi.md). The graph is an `igraph` object
representing human gene–concept relationships.

### Visualization & Export

- [R/visualize_topn.R](R/visualize_topn.R) — bar plots of top-ranked
  enriched terms
- [R/make_gene_exploration_plot.R](R/make_gene_exploration_plot.R) —
  per-gene scatter plot
- [R/plot_network.R](R/plot_network.R) — network subgraph visualization
- [R/export_gene_exploration_data.R](R/export_gene_exploration_data.R) /
  [R/export_network_summary_statistics.R](R/export_network_summary_statistics.R)
  — table exports

### Reference Workflow

See [example/example.R](example/example.R) for a complete end-to-end
run, and
[vignettes/getting-started-with-tkoi.Rmd](vignettes/getting-started-with-tkoi.Rmd)
for the narrative tutorial.
