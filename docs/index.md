# tKOI

**Transcriptomic Knowledge-Graph Omics Integration for Human Pathway
Analysis**

The `tkoi` package provides an integrative framework that combines
transcriptomic data with a human-specific biological knowledge graph.
This enables network-aware enrichment, functional interpretation, and
gene prioritization via personalized PageRank and ontology-aware
annotation.

## Web Application

For non-power users, please use the web application of
**[tKOI](https://comphealth.ucsf.edu/app/tkoi)** at `tkoi.org`.

## Contextualization Agent

For subsequent analysis upon getting network enrichment statistics,
please use **[tKOIAgent](https://github.com/BaranziniLab/tKOIAgent)**
for contextualization and network treversal.

## Documentations

Please refer to
**[Documentation](https://broccolito.github.io/software/tkoi/index.html)**
for detailed documentation of this R package.

## Installation

To install the development version from GitHub:

``` r
# Install devtools if necessary
install.packages("devtools")

# Install tkoi
devtools::install_github("Broccolito/tkoi")
```

## Core Features

- Personalized PageRank propagation using transcriptomic weights
- Permutation-based enrichment scoring for network nodes
- Functional annotation using Gene Ontology, Disease Ontology, Cell
  Ontology, Reactome, and more
- Modular and extensible S4 object design (`tKOIList`)
- Export and visualization tools for enriched subnetworks
- Seamless compatibility with `clusterProfiler`, `enrichplot`, and
  `ggplot2` \## Getting Started

### Example Workflow

This section walks you through a complete example using the `tkoi`
package—from reading expression data, running the core network analysis,
to visualizing enrichment results.

### Step 1: Load Example Gene Expression Data

The `tkoi` package includes a small example CSV file containing
simulated gene expression results. We’ll read it using `data.table` for
performance.

``` r
library(tkoi)
library(data.table)

# Get the file path of the example expression data
file_path = system.file("extdata", "example_data.csv", package = "tkoi")

# Read the CSV file
expression_data = fread(file_path)
head(expression_data)
```

The file includes columns:

- `gene_name`: Ensembl gene identifiers
- `logfc`: log2 fold-change values
- `pvalue`: associated p-values for differential expression

### Step 2: Run tKOI Network Enrichment Analysis

`tKOI` integrates transcriptomic changes with a biological knowledge
graph using a personalized PageRank algorithm. It also performs
permutations to assess statistical enrichment.

``` r
tkoi_result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,    # Predefined igraph network included with the package
  pvalue_threshold = 0.05,        # p-value filter for differential expression
  logfc_threshold = 0.25,         # Minimum log fold change
  indirect_link_threshold = 3,    # Required indirect connectivity for downstream inclusion
  topology_similarity = 0.9,      # Similarity for selecting matched genes in permutations
  n_permutation = 100,            # Number of random permutations
  damping_factor = 0.85,          # PageRank damping factor
  maximum_iteration = 500         # Max iterations for convergence
)
```

The result is an S4 object (`tKOIList`) that stores PageRank scores,
permutation statistics, and network annotations.

### Step 3: Perform Gene Ontology (GO) Enrichment

You can extend the analysis by integrating GO term enrichment using
`clusterProfiler`. This allows for side-by-side comparisons of
ontology-based and graph-based enrichment.

``` r
tkoi_result = run_gene_enrichment(tkoi_result)
```

This adds a `gene_enrichment_comparison` slot containing GO enrichment
tables and visual summaries.

### Step 4: Visualize GO vs Graph Enrichment

Two visualizations are automatically generated:

#### Scatter Plot (All Terms)

``` r
tkoi_result@gene_enrichment_comparison$comparison_scatter1
```

#### Scatter Plot (Faceted by GO Namespace)

``` r
tkoi_result@gene_enrichment_comparison$comparison_scatter2
```

These plots compare tKOI network enrichment (`beta`) with gene ontology
q-values.

### Step 5: Visualize Differential Genes in the Network

The
[`make_gene_exploration_plot()`](reference/make_gene_exploration_plot.md)
function highlights upregulated and downregulated genes in a scatter
plot based on both experimental and network evidence.

``` r
plt1 = make_gene_exploration_plot(
  tkoi_list = tkoi_result,
  sig_color = "#F39B7FB2",
  non_sig_color = "gray"
)
plt1
```

### Step 6: Export Gene-Level Prioritization Table

This returns a data frame containing logFC, p-values, PageRank scores,
and FDRs for each gene.

``` r
gene_data = export_gene_exploration_data(tkoi_result)
head(gene_data)
```

### Step 7: Visualize Top N Enriched Nodes

Use [`visualize_topn()`](reference/visualize_topn.md) to highlight the
most significantly enriched genes, pathways, or biological concepts
based on network-level statistics.

``` r
plt2 = visualize_topn(
  tkoi_list = tkoi_result,
  category = "Gene",       # Can also be "Pathway", "BiologicalProcess", etc.
  top_n = 25,
  high_color = "#FF5733",  # Strong enrichment
  low_color = "#154360"    # Moderate enrichment
)
plt2
```

### Step 8: (Optional) Save the Analysis Result

Save your full analysis object for future use:

``` r
save(tkoi_result, file = "tkoi_result.rda")
```

## S4 Object Structure

`tKOIList` is an S4 object returned by
[`run_tkoi()`](reference/run_tkoi.md) with the following slots:

- `expression_data`: Input transcriptomic measurements
- `pagerank_data`: Personalized PageRank vectors
- `network_summary_statistics`: Node-level enrichment results
- `gene_enrichment_comparison`: GO enrichment overlay and plots

## Annotation Resources

Built-in annotation tables support functional interpretation of the
knowledge graph:

- `go_annotation`, `disease_annotation`, `celltype_annotation`,
  `anatomy_annotation`
- `compound_annotation`, `protein_annotation`, `complex_annotation`
- `reaction_annotation`, `pathway_annotation`, `pwgroup_annotation`,
  etc.

Inspect them like so:

``` r
data(go_annotation)
head(go_annotation)
```

## License

MIT + file LICENSE

## Citation

Gu, W., Bellucci, G., Peetoom, B., McDonagh, M., & Baranzini, S. (in
preparation). Integrating Large-Scale Knowledge Graphs to Enhance
Transcriptomics Analysis.

## Contact

**Wanjun Gu** <wanjun.gu@ucsf.edu> ORCID:
[0000-0002-7342-7000](https://orcid.org/0000-0002-7342-7000)
