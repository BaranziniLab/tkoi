test_that("the launcher uses the installed app from any working directory", {
  skip_if_not_installed("shiny")
  arguments = NULL
  local_mocked_bindings(
    .tkoi_check_app_dependencies = function() invisible(TRUE),
    .package = "tkoi"
  )
  local_mocked_bindings(
    runApp = function(...) {
      arguments <<- list(...)
      invisible("app stopped")
    },
    .package = "shiny"
  )
  unrelated = withr::local_tempdir()
  result = withr::with_dir(unrelated, run_tkoi_app(
    host = "0.0.0.0", port = 4387, launch_browser = FALSE, display.mode = "normal"
  ))
  expect_identical(result, "app stopped")
  expect_identical(arguments$appDir, system.file("shiny", "tkoi", package = "tkoi", mustWork = TRUE))
  expect_true(file.exists(file.path(arguments$appDir, "app.R")))
  expect_identical(arguments$host, "0.0.0.0")
  expect_identical(arguments$port, 4387)
  expect_false(arguments$launch.browser)
  expect_identical(arguments$display.mode, "normal")
  expect_identical(getwd(), normalizePath(test_path()))
})

test_that("missing optional app packages produce an installation command", {
  local_mocked_bindings(
    .tkoi_app_dependencies = function() c("tkoi_test_missing_shiny", "tkoi_test_missing_dt"),
    .package = "tkoi"
  )
  expect_error(
    .tkoi_check_app_dependencies(),
    'install.packages\\(c\\("tkoi_test_missing_shiny", "tkoi_test_missing_dt"\\)\\)'
  )
})

test_that("the installed app includes assets, example input and deployment scripts", {
  app = system.file("shiny", "tkoi", package = "tkoi", mustWork = TRUE)
  files = c("app.R", "global.R", "ui.R", "server.R", "R/app_helpers.R", "data/example_data.rda",
    "www/tkoi_logo.png", "www/tkoi_favicon.ico", "install.R", "run_app.R", "README.md", "LICENSE.md")
  expect_true(all(file.exists(file.path(app, files))))
})
