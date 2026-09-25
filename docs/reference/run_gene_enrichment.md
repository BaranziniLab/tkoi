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

1.  Extracts genes with p-values at or below the threshold stored in
    `tkoi_list` (the log fold change threshold is not applied here).

2.  Conducts Gene Ontology (GO) enrichment analysis for Biological
    Process (BP), Cellular Component (CC), and Molecular Function (MF)
    using
    [`clusterProfiler::enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

3.  Merges the enrichment results with TKOI network statistics to create
    a unified dataset.

4.  Generates scatter plots to visualize the relationship between TKOI
    network enrichment effect size and gene enrichment q-values.

## Examples

``` r
if (FALSE) { # \dontrun{
tkoi_result = run_tkoi(expression_data)
tkoi_result = run_gene_enrichment(tkoi_result)

tkoi_result@gene_enrichment_comparison$enrichment_result
tkoi_result@gene_enrichment_comparison$comparison_scatter1
} # }
```
