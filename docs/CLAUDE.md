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

## Code Style

- Assignment with `=` (never `<-`), native pipe `|>` (never `%>%`),
  tidyverse layout otherwise.
- Keep
  [`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html)
  clean.

## Common Commands

\
`# Install dependencies`\
`devtools``::`[`install_deps`](https://devtools.r-lib.org/reference/install_deps.html)`(``)`\
\
`# Rebuild documentation from Roxygen2 comments`\
`devtools``::`[`document`](https://devtools.r-lib.org/reference/document.html)`(``)`\
\
`# Load the package locally for development`\
`devtools``::`[`load_all`](https://devtools.r-lib.org/reference/load_all.html)`(``)`\
\
`# Run tests`\
`devtools``::`[`test`](https://devtools.r-lib.org/reference/test.html)`(``)`\
\
`# Run a single test file`\
`testthat``::`[`test_file`](https://testthat.r-lib.org/reference/test_file.html)`(``"tests/testthat/test-<name>.R"``)`\
\
`# Include the slow tests (full network, GO enrichment)`\
[`Sys.setenv`](https://rdrr.io/r/base/Sys.setenv.html)`(``TKOI_EXTENDED_TESTS ``=`` ``"true"``)``; ``devtools``::`[`test`](https://devtools.r-lib.org/reference/test.html)`(``)`\
\
`` # Lint (config in .lintr: `=` assignment, `|>` pipe) ``\
`lintr``::`[`lint_package`](https://lintr.r-lib.org/reference/lint.html)`(``)`\
\
`# Build package`\
`devtools``::`[`build`](https://devtools.r-lib.org/reference/build.html)`(``)`\
\
`# Check package (CRAN-style)`\
`devtools``::`[`check`](https://devtools.r-lib.org/reference/check.html)`(``)`\
\
`# Build pkgdown documentation site`\
`pkgdown``::`[`build_site`](https://pkgdown.r-lib.org/reference/build_site.html)`(``)`

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
Map filtered genes to the knowledge graph via Ensembl IDs 3. Draw
`n_permutation` degree-matched null seed sets (C++, R’s RNG, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) applies) 4. Run
personalized PageRank for the observed seeds and every null set in one
batched, multithreaded C++ solve 5. Compute Z-scores, p-values, and
per-node-type FDR 6. Count seeds within 1 and 2 hops (C++ bitsets),
order tables, and annotate nodes from the `.rda` datasets

### Native Engine

[src/tkoi_engine.cpp](src/tkoi_engine.cpp) (Rcpp) holds the
performance-critical code; [R/engine.R](R/engine.R) prepares its
inputs. - `.tkoi_build_csr()` builds the symmetric normalized matrix;
`.tkoi_prepare_network()` caches it per
[`igraph::graph_id()`](https://r.igraph.org/reference/graph_id.html). -
`.tkoi_ppr_null()` solves personalized PageRank with block conjugate
gradient, 16 vectors per pass, on `n_cores` threads. It matches
`igraph::page_rank(algo = "prpack", directed = FALSE)` conventions
(self-loops count twice, `weight` edge attribute used). - Results are
bit-identical for any thread count (fixed chunking and reduction
order). - `.tkoi_sample_null()` reproduces the legacy
[`setdiff()`](https://rdrr.io/r/base/sets.html) +
[`sample()`](https://rdrr.io/r/base/sample.html) loop draw for draw. -
After editing C++ or `// [[Rcpp::export]]` tags, run
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
then
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html). -
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
compiles with `-O0`; benchmark only an installed build
(`R CMD INSTALL`).

### Pre-loaded Annotation Data

The `data/` directory holds ~20 `.rda` files (GO, disease, cell type,
pathway, reaction, anatomy, compound, protein annotations). These are
loaded automatically by the package and are the primary reference tables
used to annotate network nodes after PageRank.

### Knowledge Graph

`data/tkoi_net.rda` is the plain `igraph` network (939,059 nodes,
10,622,200 undirected edges, xz-compressed, ~55 MB). It lazy-loads as
[`tkoi::tkoi_net`](reference/tkoi_net.md). There is no encryption.

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
