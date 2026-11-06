find_console_root <- function() {
  root <- normalizePath(Sys.getenv("TRANSCRIPTOME_CODE_DIR", unset = getwd()), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(root, "scripts", "console.R")) && file.exists(file.path(root, "_targets.R"))) return(root)
    parent <- dirname(root)
    if (identical(parent, root)) stop("请先在RStudio打开代码项目，或设置TRANSCRIPTOME_CODE_DIR")
    root <- parent
  }
}

open_console <- function() {
  root <- find_console_root()
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "console.R"), envir = env)
  env$launch_console(code_dir = root)
}

close_console <- function() {
  root <- find_console_root()
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "console.R"), envir = env)
  env$stop_console(code_dir = root)
}
