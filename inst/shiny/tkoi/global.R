# The installed tkoi library owns computation, the S4 class, and graph data.
if (!requireNamespace("tkoi", quietly = TRUE)) {
  stop(
    "The tKOI Shiny app requires tkoi >= 1.2.0. Run Rscript install.R from the release bundle, then restart R.",
    call. = FALSE
  )
}
if (getNamespaceVersion("tkoi") != utils::packageVersion("tkoi")) {
  stop(
    "The loaded tkoi version differs from the installed library. Restart R before launching the tKOI app.",
    call. = FALSE
  )
}
if (getNamespaceVersion("tkoi") < "1.2.0") {
  stop(
    "The tKOI Shiny app requires tkoi >= 1.2.0. Run Rscript install.R from the release bundle, then restart R.",
    call. = FALSE
  )
}
getFromNamespace(".tkoi_check_app_dependencies", "tkoi")()
app_packages = getFromNamespace(".tkoi_app_dependencies", "tkoi")()
suppressPackageStartupMessages(library(tkoi))
invisible(lapply(setdiff(app_packages, "tkoi"), function(package) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}))
tkoi_version = as.character(getNamespaceVersion("tkoi"))
tkoi_defaults = formals(tkoi::run_tkoi)
source("R/app_helpers.R", local = TRUE)
example_environment = new.env(parent = emptyenv())
load("data/example_data.rda", envir = example_environment)
example_data = as.data.frame(example_environment$example_data)
rm(example_environment, app_packages)
