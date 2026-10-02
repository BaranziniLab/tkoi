# Run from the tkoi repository root. An optional argument selects another
# tkoi checkout.
arguments = commandArgs(trailingOnly = TRUE)
local_tkoi = if (length(arguments)) arguments[[1]] else "."
metadata = read.dcf("DESCRIPTION")
dependencies = trimws(strsplit(
  gsub("\\([^)]*\\)", "", paste(metadata[1, c("Imports", "Suggests")], collapse = ",")), ",", fixed = TRUE
)[[1]])
if (file.exists(file.path(local_tkoi, "DESCRIPTION"))) {
  core_imports = trimws(strsplit(
    gsub("\\([^)]*\\)", "", read.dcf(file.path(local_tkoi, "DESCRIPTION"))[1, "Imports"]),
    ",", fixed = TRUE
  )[[1]])
  dependencies = union(dependencies, core_imports)
}
cran = setdiff(dependencies, c("tkoi", "clusterProfiler", "org.Hs.eg.db"))
missing = cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc = c("clusterProfiler", "org.Hs.eg.db")
missing_bioc = bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
if (file.exists(file.path(local_tkoi, "DESCRIPTION"))) {
  # Update the installed package from source, never source its functions in Shiny.
  status = system2(file.path(R.home("bin"), "R"), c(
    "CMD", "INSTALL", "--no-multiarch", shQuote(normalizePath(local_tkoi))
  ))
  if (status != 0L) {
    stop("Installation of the local tkoi package failed.", call. = FALSE)
  }
} else {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("BaranziniLab/tkoi", upgrade = "never")
}
message("Installed tkoi ", utils::packageVersion("tkoi"), ". Restart R, then run tkoi::run_tkoi_app().")
