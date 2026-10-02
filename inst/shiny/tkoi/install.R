# Install a standalone release using its bundled tkoi source package.
local({
  script = commandArgs(trailingOnly = FALSE)
  script = sub("^--file=", "", script[startsWith(script, "--file=")])
  app_dir = if (length(script)) dirname(normalizePath(script[[1]])) else getwd()
  archives = Sys.glob(file.path(app_dir, "vendor", "tkoi_*.tar.gz"))
  if (length(archives) != 1) {
    stop("Extract the standalone release bundle, then run its install.R script.", call. = FALSE)
  }
  metadata_dir = tempfile("tkoi-description-")
  dir.create(metadata_dir)
  on.exit(unlink(metadata_dir, recursive = TRUE), add = TRUE)
  utils::untar(archives[[1]], files = "tkoi/DESCRIPTION", exdir = metadata_dir)
  metadata = read.dcf(file.path(metadata_dir, "tkoi", "DESCRIPTION"))
  core = trimws(strsplit(gsub("\\([^)]*\\)", "", metadata[1, "Imports"]), ",", fixed = TRUE)[[1]])
  app = c("shiny", "shinythemes", "data.table", "readxl", "openxlsx", "plotly", "DT")
  bioc = c("clusterProfiler", "org.Hs.eg.db")
  cran = setdiff(union(core, app), bioc)
  missing = cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    utils::install.packages(missing, repos = "https://cloud.r-project.org")
  }
  missing_bioc = bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
  status = system2(file.path(R.home("bin"), "R"), c(
    "CMD", "INSTALL", "--no-multiarch", shQuote(archives[[1]])
  ))
  if (status != 0L) {
    stop("Installation of the bundled tkoi package failed.", call. = FALSE)
  }
  message("Installed tkoi ", utils::packageVersion("tkoi"), ". Run Rscript run_app.R to launch the app.")
})
