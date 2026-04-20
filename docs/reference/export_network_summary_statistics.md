# Export Network Summary Statistics from tKOI Result

This function exports the `network_summary_statistics` slot from a
`tkoi_result` object to an Excel file.

## Usage

``` r
export_network_summary_statistics(tkoi_result, filename = "tkoi_result.xlsx")
```

## Arguments

- tkoi_result:

  An object of class `tkoi_result`, typically returned by a tKOI
  analysis pipeline. This object must include a slot named
  `@network_summary_statistics`, which contains a data frame of
  network-level metrics.

- filename:

  A character string specifying the name of the output Excel file.
  Defaults to `"tkoi_result.xlsx"`.

## Value

No return value. This function is called for its side effect of writing
a file to disk.

## Examples

``` r
if (FALSE) { # \dontrun{
export_network_summary_statistics(tkoi_result, filename = "my_output.xlsx")
} # }
```
