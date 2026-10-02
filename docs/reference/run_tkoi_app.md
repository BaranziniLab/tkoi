# Launch the tKOI Shiny App

Opens the Shiny interface bundled with the installed tkoi package. The
app accepts differential expression files, runs
[`run_tkoi()`](https://baranzinilab.github.io/tkoi/reference/run_tkoi.md),
displays results, and exports tables and network plots.

## Usage

``` r
run_tkoi_app(
  host = "127.0.0.1",
  port = NULL,
  launch_browser = interactive(),
  ...
)
```

## Arguments

- host:

  Address to bind to. Defaults to `"127.0.0.1"` for local use. Use
  `"0.0.0.0"` when serving the app on a network interface.

- port:

  Port number, or `NULL` to let Shiny choose a port.

- launch_browser:

  Whether to open the app in a browser. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html). A browser
  function accepted by
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) can
  also be supplied.

- ...:

  Additional arguments passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

The value returned by
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html),
invisibly. The app runs until interrupted, for example by pressing
Escape in RStudio or Ctrl+C in R.

## Details

The optional app packages must be installed first:
`install.packages(c("shiny", "shinythemes", "data.table", "readxl", "openxlsx", "plotly", "DT"))`.

## Examples

``` r
if (FALSE) { # \dontrun{
tkoi::run_tkoi_app()
tkoi::run_tkoi_app(host = "0.0.0.0", port = 3838, launch_browser = FALSE)
} # }
```
