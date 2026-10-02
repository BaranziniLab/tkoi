#' Launch the tKOI Shiny App
#'
#' Opens the Shiny interface bundled with the installed tkoi package. The app
#' accepts differential expression files, runs [run_tkoi()], displays results,
#' and exports tables and network plots.
#'
#' The optional app packages must be installed first:
#' `install.packages(c("shiny", "shinythemes", "data.table", "readxl",
#' "openxlsx", "plotly", "DT"))`.
#'
#' @param host Address to bind to. Defaults to `"127.0.0.1"` for local use.
#'   Use `"0.0.0.0"` when serving the app on a network interface.
#' @param port Port number, or `NULL` to let Shiny choose a port.
#' @param launch_browser Whether to open the app in a browser. Defaults to
#'   `interactive()`. A browser function accepted by [shiny::runApp()] can also
#'   be supplied.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return The value returned by [shiny::runApp()], invisibly. The app runs until
#'   interrupted, for example by pressing Escape in RStudio or Ctrl+C in R.
#'
#' @examples
#' \dontrun{
#' tkoi::run_tkoi_app()
#' tkoi::run_tkoi_app(host = "0.0.0.0", port = 3838, launch_browser = FALSE)
#' }
#'
#' @export
run_tkoi_app = function(host = "127.0.0.1", port = NULL, launch_browser = interactive(), ...) {
  .tkoi_check_app_dependencies()
  app_dir = system.file("shiny", "tkoi", package = "tkoi", mustWork = TRUE)
  invisible(shiny::runApp(appDir = app_dir, host = host, port = port, launch.browser = launch_browser, ...))
}

.tkoi_app_dependencies = function() {
  c("shiny", "shinythemes", "data.table", "dplyr", "ggplot2", "glue", "readxl", "openxlsx", "plotly", "DT")
}

.tkoi_check_app_dependencies = function() {
  packages = .tkoi_app_dependencies()
  missing = packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    quoted = paste(sprintf('"%s"', missing), collapse = ", ")
    stop(
      "The tKOI Shiny app requires additional packages. Install them with:\n",
      "install.packages(c(", quoted, "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
