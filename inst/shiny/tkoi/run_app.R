arguments = commandArgs(trailingOnly = TRUE)
port = if (length(arguments)) as.integer(arguments[[1]]) else 3838L
host = if (length(arguments) >= 2) arguments[[2]] else "127.0.0.1"
script = commandArgs(trailingOnly = FALSE)
script = sub("^--file=", "", script[startsWith(script, "--file=")])
app_dir = if (length(script)) dirname(normalizePath(script[[1]])) else getwd()
shiny::runApp(appDir = app_dir, host = host, port = port, launch.browser = interactive())
