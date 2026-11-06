console_worker <- function(request_file) {
  request <- readRDS(request_file)
  context <- request$context
  setwd(context$code_dir)
  Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = context$data_dir, TRANSCRIPTOME_CONFIG_DIR = request$config_dir)
  data.table::setDTthreads(1L)
  message("开始：", format(Sys.time()), "；操作：", request$action)
  config <- console_read_config_lst(request$config_dir)
  report <- console_validate_lst(context, config)
  saveRDS(report, file.path(dirname(request_file), "validation.rds"))
  if (request$action == "build") {
    if (!requireNamespace("targets", quietly = TRUE)) stop("请先安装targets")
    selected <- request$targets
    message("构建目标：", if (is.null(selected)) "全部" else paste(selected, collapse = ", "))
    args <- list(callr_function = NULL, reporter = "verbose", store = context$store)
    if ("as_job" %in% names(formals(targets::tar_make))) args$as_job <- FALSE
    if (!is.null(selected)) args$names <- substitute(tidyselect::all_of(x), list(x = selected))
    # 当前进程已为干净会话
    do.call(targets::tar_make, args)
  }
  message("完成：", format(Sys.time()))
  invisible(TRUE)
}
