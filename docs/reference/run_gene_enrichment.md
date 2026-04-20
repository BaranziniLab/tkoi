# Run Gene Enrichment and Compare with TKOI Data

This function performs gene enrichment analysis using the
`clusterProfiler` package and integrates the results with the `tKOIList`
object, specifically comparing the enrichment results with the TKOI
data. It generates scatter plots visualizing the relationship between
TKOI scores and enrichment results.

## Usage

``` r
run_gene_enrichment(tkoi_list)
```

## Arguments

- tkoi_list:

  An object of class `tKOIList` that contains the input data for
  analysis. The object should include the following slots:

  - `expression_data`: A `data.frame` with columns including `gene_name`
    and `pvalue`.

  - `network_summary_statistics`: A `list` containing network-level
    summary data for Biological Process (BP), Cellular Component (CC),
    and Molecular Function (MF), each stored as a `data.frame`.

  - `pvalue_threshold`: A `numeric` value specifying the significance
    threshold for filtering genes.

## Value

A `tKOIList` object with the `gene_enrichment_comparison` slot
populated. This slot is a `list` containing:

- `enrichment_result`: A `data.frame` with the merged results of the
  enrichment analysis and TKOI network data.

- `comparison_scatter1`: A `ggplot` object visualizing a scatter plot
  without facets.

- `comparison_scatter2`: A `ggplot` object visualizing a scatter plot
  with facets for each namespace.

## Details

The function performs the following steps:

1.  Extracts genes with p-values below the threshold specified in the
    `tkoi_list`.

2.  Conducts Gene Ontology (GO) enrichment analysis for Biological
    Process (BP), Cellular Component (CC), and Molecular Function (MF)
    using
    [`clusterProfiler::enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

3.  Computes pairwise term similarities using
    [`enrichplot::pairwise_termsim`](https://rdrr.io/pkg/enrichplot/man/pairwise_termsim.html).

4.  Merges the enrichment results with TKOI network statistics to create
    a unified dataset.

5.  Generates scatter plots to visualize the relationship between TKOI
    network enrichment effect size and gene enrichment q-values.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage of the function
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Create a dummy tKOIList object
tkoi_list = new("tKOIList",
                expression_data = data.frame(
                  gene_name = c("gene1", "gene2", "gene3"),
                  pvalue = c(0.01, 0.02, 0.2)
                ),
                network_summary_statistics = list(
                  BiologicalProcess = data.frame(
                    identifier = c("GO:0008150", "GO:0009987"),
                    node_id = c(1, 2),
                    pagerank = c(0.2, 0.3),
                    beta = c(0.5, 0.4),
                    p_value = c(0.01, 0.02),
                    fdr = c(0.05, 0.1),
                    definition = c("process1", "process2")
                  ),
                  CellularComponent = data.frame(),
                  MolecularFunction = data.frame()
                ))

# Run the enrichment function
result = run_gene_enrichment(tkoi_list)
} # }
```
