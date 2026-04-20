# Create Gene Exploration Plot

This function generates a scatter plot that visualizes the relationship
between experimental p-values and tKOI-adjusted false discovery rates
(FDR) for genes in the network. The plot highlights significantly
upregulated and downregulated genes based on user-defined thresholds for
log fold change and p-values.

## Usage

``` r
make_gene_exploration_plot(
  tkoi_list,
  sig_color = "#F39B7FB2",
  non_sig_color = "gray"
)
```

## Arguments

- tkoi_list:

  An object of class `tKOIList`. This object must contain the following
  slots:

  - `expression_data`: A `data.frame` containing experimental data with
    columns `gene_name`, `logfc`, and `pvalue`.

  - `network_summary_statistics`: A list with a `data.frame` for
    gene-level statistics, including `node_id` and `fdr`.

  - `pvalue_threshold`: A numeric value specifying the p-value
    threshold.

  - `logfc_threshold`: A numeric value specifying the log fold change
    threshold.

- sig_color:

  A character string specifying the color for significant genes in the
  plot. Default is `"#F39B7FB2"`.

- non_sig_color:

  A character string specifying the color for non-significant genes in
  the plot. Default is `"gray"`.

## Value

A `ggplot` object representing the scatter plot.

## Details

The function performs the following steps:

1.  Merges the `expression_data` from the `tKOIList` object with gene
    metadata using `inner_join`.

2.  Merges the result with gene-level network data using `right_join`.

3.  Filters and annotates the data to classify genes as upregulated or
    downregulated based on the sign of `logfc`.

4.  Colors genes based on their significance determined by
    `logfc_threshold`.

5.  Generates a scatter plot with ggplot2, where:

    - The x-axis represents `-log10(pvalue)` (experimental p-value).

    - The y-axis represents `-log10(fdr)` (tKOI-adjusted FDR).

    - Points are colored and categorized into upregulated and
      downregulated facets.

The plot includes vertical and horizontal dashed lines to indicate the
thresholds for p-value and FDR.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a dummy tKOIList object
tkoi_list = new("tKOIList",
                expression_data = data.frame(
                  gene_name = c("gene1", "gene2", "gene3"),
                  logfc = c(1.2, -0.8, 0.5),
                  pvalue = c(0.01, 0.03, 0.2)
                ),
                network_summary_statistics = list(
                  Gene = data.frame(
                    node_id = c(1, 2, 3),
                    fdr = c(0.02, 0.05, 0.1)
                  )
                ),
                pvalue_threshold = 0.05,
                logfc_threshold = 0.5)

# Generate the plot
plt = make_gene_exploration_plot(tkoi_list)
print(plt)
} # }
```
