# 从代码根目录加载入口
launch_console <- function(code_dir = getwd(), data_dir = Sys.getenv("TRANSCRIPTOME_PROJECT_DIR", unset = code_dir),
                           config_dir = Sys.getenv("TRANSCRIPTOME_CONFIG_DIR", unset = "config"), viewer = FALSE, open = TRUE) {
  packages <- c("shiny", "callr", "yaml", "filelock", "digest", "data.table", "tibble", "dplyr", "httpuv")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("请先安装控制台依赖：", paste(missing, collapse = ", "))
  helper <- new.env(parent = globalenv())
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  sys.source(file.path(code_dir, "console", "backend.R"), envir = helper)
  context <- helper$console_context_lst(code_dir, data_dir, config_dir)
  registry <- getOption("transcriptome.console.processes", list())
  old <- registry[[code_dir]]
  if (!is.null(old) && old$process$is_alive()) {
    if (!identical(old$data_dir, context$data_dir) || !identical(old$config_dir, context$config_dir)) stop("请先stop_console()，再更换初始项目配置")
    handle <- old
  } else {
    dir.create(context$state_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- tempfile("app_", tmpdir = context$state_dir, fileext = ".log")
    port <- httpuv::randomPort()
    process <- callr::r_bg(function(context, port) {
      env <- new.env(parent = globalenv())
      sys.source(file.path(context$code_dir, "console", "backend.R"), envir = env)
      sys.source(file.path(context$code_dir, "console", "app.R"), envir = env)
      app <- env$console_app(context$code_dir, context$data_dir, context$config_dir)
      shiny::runApp(app, host = "127.0.0.1", port = port, launch.browser = FALSE)
    }, args = list(context, port), wd = code_dir, libpath = .libPaths(), stdout = log_file, stderr = "2>&1",
      supervise = TRUE, cleanup_tree = TRUE, user_profile = FALSE)
    ready <- FALSE
    deadline <- Sys.time() + 25
    while (process$is_alive() && Sys.time() < deadline) {
      connection <- tryCatch(socketConnection("127.0.0.1", port = port, open = "r+b", blocking = TRUE, timeout = 0.2), error = function(e) NULL)
      if (!is.null(connection)) {
        close(connection)
        ready <- TRUE
        break
      }
      Sys.sleep(0.1)
    }
    if (!ready) {
      process$kill_tree()
      stop("控制台未能启动：\n", helper$console_log_tail(log_file))
    }
    handle <- list(process = process, url = paste0("http://127.0.0.1:", port, "/"), log_file = log_file,
      code_dir = code_dir, data_dir = context$data_dir, config_dir = context$config_dir)
    registry[[code_dir]] <- handle
    options(transcriptome.console.processes = registry)
  }
  url <- handle$url
  in_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  if (open) {
    if (in_rstudio && viewer) rstudioapi::viewer(url, height = "maximize") else {
      if (in_rstudio) url <- rstudioapi::translateLocalUrl(url, absolute = TRUE)
      utils::browseURL(url)
    }
  }
  message("控制台已启动；R Console 可继续使用。关闭页面后可再次 launch_console()。")
  message("关闭控制台及其任务：stop_console()；启动日志：", handle$log_file)
  invisible(handle)
}

stop_console <- function(code_dir = getwd()) {
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  registry <- getOption("transcriptome.console.processes", list())
  handle <- registry[[code_dir]]
  if (is.null(handle)) return(invisible(FALSE))
  if (handle$process$is_alive()) {
    handle$process$kill_tree()
    handle$process$wait(timeout = 2000)
    if (handle$process$is_alive()) stop("控制台进程尚未退出")
  }
  registry[[code_dir]] <- NULL
  options(transcriptome.console.processes = registry)
  invisible(TRUE)
}
