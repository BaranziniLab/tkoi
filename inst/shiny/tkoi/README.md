# tKOI Shiny app

Upload differential expression data, run tKOI, explore results, and export tables
or network plots. The app uses the tkoi R package for the analysis, graph and
annotations.

## Standalone release

Extract `tkoi-shiny-VERSION.zip` or `tkoi-shiny-VERSION.tar.gz` into a directory
of your choice. From the extracted directory, install the dependencies once:

```sh
Rscript install.R
```

The `vendor/` directory contains the exact tkoi source package for this release.
The installer uses that local archive and downloads the other R packages from
CRAN and Bioconductor. It does not need a GitHub checkout or GitHub credentials.
R 4.1 or later, internet access for dependency installation, and a C++ compiler
are required. On Linux, install the system libraries required by the R
dependencies, including curl, OpenSSL and libxml2 development headers. On
macOS, use the Xcode Command Line Tools; on Windows, use the matching Rtools.

Launch locally:

```sh
Rscript run_app.R
```

The default address is `http://127.0.0.1:3838`. To serve on a network interface:

```sh
Rscript run_app.R 3838 0.0.0.0
```

For Shiny Server, place the extracted directory under its configured app
directory, for example `/srv/shiny-server/tkoi`, and make the dependencies
available to the account running Shiny Server. For Posit Connect, publish the
extracted directory as a Shiny application after installing its dependencies.
The app entry point is `app.R`; all runtime app files are in this directory.

## R package

The same app is bundled with the installed tkoi package:

```r
tkoi::run_tkoi_app()
tkoi::run_tkoi_app(host = "0.0.0.0", port = 3838, launch_browser = FALSE)
```

## Development

App development now lives in `inst/shiny/tkoi` in the tkoi repository. From that
repository's root:

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/run_app.R
Rscript scripts/build_shiny_release.R
```

Tests are integrated into tkoi's `tests/testthat` suite. The app sources,
example input, assets and integration tests were imported from tKOIvis commit
`a57652c` on 1 October 2026. The original repository and checkout are retained
while existing processes finish. Its standalone RStudio project and auto-sync
scripts are replaced by the tkoi repository's project and development scripts.

## Input

Upload a `.csv`, `.tsv` or `.xlsx` with the columns `gene_name` (Ensembl gene
identifier), `logfc` (numeric log fold change) and `pvalue` (numeric p-value).
Download Example provides a reference input. Both uploads and the example are
analyzed afresh by the tkoi package.

Results can be exported as an `.rda` containing the package's `tKOIList` object
and engine version, or as an Excel workbook with one sheet per node type.

Source code is MIT licensed; see `LICENSE.md`. The knowledge graph is subject
to the licenses of its data sources. Contact Wanjun Gu or Sergio Baranzini at
UCSF about graph data use.

Probability displays use scientific notation with HTML exponents, retain numeric sorting, and export numeric P/Q values alongside log probabilities, bound flags and plain-text display columns. Stored zeros are recovered from available logs; a zero without its logarithm is labelled unresolved. Requires tkoi 1.3.1 or later.
