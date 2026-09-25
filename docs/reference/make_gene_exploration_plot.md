# Create Gene Exploration Plot

Plots, for every gene of the experiment that is a gene node of the
network, the experimental p-value against the tKOI false discovery rate
(FDR), with up- and down-regulated genes in separate panels. Genes
supported by both the experiment and the network appear in the top right
of each panel.

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

  A `tKOIList` object returned by [`run_tkoi()`](run_tkoi.md). It uses
  the `expression_data` slot (columns `gene_name`, `logfc`, and
  `pvalue`), the `Gene` table of `network_summary_statistics` (columns
  `node_id`, `name`, and `fdr`), and the `pvalue_threshold` and
  `logfc_threshold` of the run.

- sig_color:

  Fill color of genes whose `abs(logfc)` is at least the run's
  `logfc_threshold`. Default `"#F39B7FB2"`.

- non_sig_color:

  Fill color of the other genes. Default `"gray"`.

## Value

A `ggplot` object.

## Details

`expression_data` is cleaned the same way as in
[`run_tkoi()`](run_tkoi.md): rows with a missing or blank `gene_name`
are dropped, and only the first row of each gene is used. Genes are
matched to network nodes through their Ensembl IDs (see
[genes](genes.md)); network genes without expression data are not
plotted.

In the plot:

- The x-axis is `-log10(pvalue)` (experimental p-value).

- The y-axis is `-log10(fdr)` (tKOI FDR).

- Genes with `logfc >= 0` are in the `Up-regulated` panel, the others in
  the `Down-regulated` panel.

- The dashed vertical line marks the run's `pvalue_threshold` and the
  dashed horizontal line an FDR of 0.05.

## See also

[`export_gene_exploration_data()`](export_gene_exploration_data.md) for
the underlying table.

## Examples

``` r
if (FALSE) { # \dontrun{
expression_data = data.table::fread(
  system.file("extdata", "example_data.csv", package = "tkoi")
)

set.seed(1)
tkoi_result = run_tkoi(expression_data = expression_data)

plt = make_gene_exploration_plot(tkoi_result, sig_color = "#F39B7FB2", non_sig_color = "gray")
plt
} # }
```
