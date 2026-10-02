arguments = commandArgs(trailingOnly = TRUE)
port = if (length(arguments)) as.integer(arguments[[1]]) else 3838L
tkoi::run_tkoi_app(host = "127.0.0.1", port = port, launch_browser = interactive())
