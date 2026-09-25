# Export Gene Exploration Data

Combines the differential expression table of a
[`run_tkoi()`](run_tkoi.md) result with the tKOI network statistics of
every gene node, so that experimental and network evidence can be
compared gene by gene.

## Usage

``` r
export_gene_exploration_data(tkoi_list)
```

## Arguments

- tkoi_list:

  A `tKOIList` object returned by [`run_tkoi()`](run_tkoi.md). Two slots
  are used:

  - `expression_data`: the input table, with columns `gene_name`
    (Ensembl gene IDs), `logfc`, and `pvalue`.

  - `network_summary_statistics`: its `Gene` table, with columns
    `node_id`, `identifier`, `name`, `pagerank`, `beta`, `p_value`, and
    `fdr`.

## Value

A data frame with one row per row of the `Gene` table and columns:

- `gene_name`: Ensembl gene ID from `expression_data`.

- `gene_symbol`: Gene symbol.

- `id`: Node ID of the gene in the network (`node_id` in the `Gene`
  table).

- `identifier`: NCBI Entrez Gene ID, as a character string.

- `experimental_logfc`: `logfc` from `expression_data`.

- `experimental_pvalue`: `pvalue` from `expression_data`.

- `pagerank`: Observed personalized PageRank.

- `tkoi_beta`: tKOI network enrichment z-score (`beta`).

- `tkoi_pvalue`: Unadjusted one-sided tKOI p-value.

- `tkoi_fdr`: tKOI false discovery rate (Benjamini-Hochberg, among gene
  nodes).

## Details

`expression_data` is cleaned the same way as in
[`run_tkoi()`](run_tkoi.md): rows with a missing or blank `gene_name`
are dropped, and only the first row of each gene is used. Genes are
matched to network nodes through their Ensembl IDs (see
[genes](genes.md)). Every row of the `Gene` table is kept, so network
genes that are not in `expression_data` have `NA` in `gene_name`,
`experimental_logfc`, and `experimental_pvalue`.

## See also

[`make_gene_exploration_plot()`](make_gene_exploration_plot.md) to plot
the same data,
[`export_network_summary_statistics()`](export_network_summary_statistics.md)
to export every node type.

## Examples

``` r
if (FALSE) { # \dontrun{
expression_data = data.table::fread(
  system.file("extdata", "example_data.csv", package = "tkoi")
)

set.seed(1)
tkoi_result = run_tkoi(expression_data = expression_data)

gene_data = export_gene_exploration_data(tkoi_result)
head(gene_data)

# Genes supported by both the experiment and the network
subset(gene_data, experimental_pvalue <= 0.05 & tkoi_fdr <= 0.05)
} # }
```
