# Export Network Summary Statistics to Excel

Writes the `network_summary_statistics` slot of a
[`run_tkoi()`](run_tkoi.md) result to an Excel workbook, with one sheet
per node type (for example `Gene`, `BiologicalProcess`, `Disease`).

## Usage

``` r
export_network_summary_statistics(tkoi_result, filename = "tkoi_result.xlsx")
```

## Arguments

- tkoi_result:

  A `tKOIList` object returned by [`run_tkoi()`](run_tkoi.md).

- filename:

  Path of the Excel file to write. An existing file is overwritten.
  Default `"tkoi_result.xlsx"`.

## Value

The path of the written file, invisibly, as returned by
[`writexl::write_xlsx()`](https://docs.ropensci.org/writexl//reference/write_xlsx.html).
The function is called for its side effect of writing the workbook.

## See also

[`export_gene_exploration_data()`](export_gene_exploration_data.md) for
a gene-level table.

## Examples

``` r
if (FALSE) { # \dontrun{
expression_data = data.table::fread(
  system.file("extdata", "example_data.csv", package = "tkoi")
)

set.seed(1)
tkoi_result = run_tkoi(expression_data = expression_data)

path = export_network_summary_statistics(
  tkoi_result,
  filename = file.path(tempdir(), "tkoi_network_statistics.xlsx")
)
path
} # }
```
