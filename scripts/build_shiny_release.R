# Build the R package and standalone app archives from the tkoi repository root.
local({
  root = normalizePath(getwd())
  metadata = read.dcf(file.path(root, "DESCRIPTION"))
  stopifnot(metadata[1, "Package"] == "tkoi")
  version = metadata[1, "Version"]
  arguments = commandArgs(trailingOnly = TRUE)
  output = if (length(arguments)) arguments[[1]] else file.path(root, "release")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output = normalizePath(output)
  original = getwd()
  on.exit(setwd(original), add = TRUE)
  setwd(output)
  status = system2(file.path(R.home("bin"), "R"), c(
    "CMD", "build", "--no-manual", shQuote(root)
  ))
  if (status != 0L) stop("R package build failed.", call. = FALSE)
  package = file.path(output, paste0("tkoi_", version, ".tar.gz"))
  stopifnot(file.exists(package))

  staging = tempfile("tkoi-shiny-release-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  name = paste0("tkoi-shiny-", version)
  app = file.path(staging, name)
  dir.create(app)
  source = file.path(root, "inst", "shiny", "tkoi")
  files = list.files(source, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  stopifnot(all(file.copy(files, app, recursive = TRUE)))
  dir.create(file.path(app, "vendor"))
  stopifnot(file.copy(package, file.path(app, "vendor")))
  writeLines(c(
    paste("tkoi Shiny release", version),
    paste("tkoi source package:", basename(package)),
    paste("App source commit:", system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE))
  ), file.path(app, "VERSION.txt"))

  setwd(staging)
  zip = file.path(output, paste0(name, ".zip"))
  archive = file.path(output, paste0(name, ".tar.gz"))
  if (file.exists(zip)) unlink(zip)
  utils::zip(zip, name, flags = "-r9Xq")
  utils::tar(archive, files = name, compression = "gzip", tar = "internal")
  stopifnot(file.exists(zip), file.exists(archive))

  setwd(output)
  assets = basename(c(package, zip, archive))
  sha256 = Sys.which("sha256sum")
  if (nzchar(sha256)) {
    checksums = system2(sha256, shQuote(assets), stdout = TRUE)
  } else {
    shasum = Sys.which("shasum")
    if (!nzchar(shasum)) stop("Install sha256sum or shasum to checksum the release.", call. = FALSE)
    checksums = system2(shasum, c("-a", "256", shQuote(assets)), stdout = TRUE)
  }
  writeLines(checksums, "SHA256SUMS")
  cat("\nRelease assets:\n", paste(file.path(output, c(assets, "SHA256SUMS")), collapse = "\n"), "\n")
})
