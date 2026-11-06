# 从代码根目录加载入口
launch_console <- function(code_dir = getwd(), data_dir = Sys.getenv("TRANSCRIPTOME_PROJECT_DIR", unset = code_dir),
                           config_dir = Sys.getenv("TRANSCRIPTOME_CONFIG_DIR", unset = "config"), viewer = TRUE, open = TRUE,
                           port = getOption("transcriptome.console.port", NULL), background = TRUE) {
  if (!is.null(port) && (!is.numeric(port) || length(port) != 1L || !is.finite(port) || port != round(port) || port < 1 || port > 65535)) {
    stop("port须为1至65535的整数，或NULL（随机端口）")
  }
  if (!is.logical(background) || length(background) != 1L || is.na(background)) stop("background须为TRUE或FALSE")
  packages <- c("shiny", "callr", "yaml", "filelock", "digest", "data.table", "tibble", "dplyr", "httpuv")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("请先安装控制台依赖：", paste(missing, collapse = ", "))
  helper <- new.env(parent = globalenv())
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  sys.source(file.path(code_dir, "console", "backend.R"), envir = helper)
  context <- helper$console_context_lst(code_dir, data_dir, config_dir)
  in_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  registry <- getOption("transcriptome.console.processes", list())
  old <- registry[[code_dir]]
  if (!is.null(old) && old$process$is_alive()) {
    active_port <- if (!is.null(old$port)) old$port else as.integer(sub("^http://127[.]0[.]0[.]1:([0-9]+)/.*$", "\\1", old$url))
    if (!background) stop("后台控制台正在运行；请先stop_console()，再以前台模式启动")
    if (!is.null(port) && port != active_port) stop("当前控制台使用端口", active_port, "；请先stop_console()，再切换到端口", port)
    if (!identical(old$data_dir, context$data_dir) || !identical(old$config_dir, context$config_dir)) stop("请先stop_console()，再更换初始项目配置")
    handle <- old
    handle$port <- active_port
  } else {
    port <- if (is.null(port)) httpuv::randomPort() else as.integer(port)
    probe <- tryCatch(httpuv::startServer("127.0.0.1", port, list(call = function(req) list(status = 200L, body = "")), quiet = TRUE),
      error = function(e) stop("无法监听127.0.0.1:", port, "；端口可能被占用或无权限。", conditionMessage(e), call. = FALSE))
    httpuv::stopServer(probe)
    if (!background) {
      sys.source(file.path(code_dir, "console", "app.R"), envir = helper)
      browser <- if (!open) FALSE else if (in_rstudio && viewer) rstudioapi::viewer else getOption("shiny.launch.browser", interactive())
      message("以前台模式启动控制台：http://127.0.0.1:", port, "/；使用RStudio Stop按钮结束。")
      return(invisible(shiny::runApp(helper$console_app(code_dir, context$data_dir, context$config_dir), host = "127.0.0.1", port = port,
        launch.browser = browser)))
    }
    dir.create(context$state_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- tempfile("app_", tmpdir = context$state_dir, fileext = ".log")
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
      connection <- tryCatch(suppressWarnings(socketConnection("127.0.0.1", port = port, open = "r+b", blocking = TRUE, timeout = 0.2)),
        error = function(e) NULL)
      if (!is.null(connection)) {
        close(connection)
        ready <- TRUE
        break
      }
      Sys.sleep(0.1)
    }
    if (!ready || !process$is_alive()) {
      process$kill_tree()
      stop("控制台未能启动：\n", helper$console_log_tail(log_file))
    }
    handle <- list(process = process, url = paste0("http://127.0.0.1:", port, "/"), log_file = log_file,
      code_dir = code_dir, data_dir = context$data_dir, config_dir = context$config_dir, port = port)
  }
  url <- handle$url
  handle$proxy_path <- if (in_rstudio) tryCatch(rstudioapi::translateLocalUrl(url, absolute = FALSE), error = function(e) NULL) else NULL
  handle$proxy_url <- if (in_rstudio) tryCatch(rstudioapi::translateLocalUrl(url, absolute = TRUE), error = function(e) NULL) else NULL
  registry[[code_dir]] <- handle
  options(transcriptome.console.processes = registry)
  if (open) {
    if (in_rstudio && viewer) rstudioapi::viewer(url, height = "maximize") else {
      if (in_rstudio) {
        if (is.null(handle$proxy_url)) stop("服务已启动，但代理地址转换失败；可改用viewer=TRUE。日志：", handle$log_file)
        url <- handle$proxy_url
      }
      utils::browseURL(url)
    }
  }
  message("控制台服务已监听：", handle$url, "；后台运行，R Console可继续使用。")
  if (!is.null(handle$proxy_path)) message("RStudio代理相对路径：", handle$proxy_path)
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
