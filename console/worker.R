console_worker <- function(request_file) {
  request <- readRDS(request_file)
  context <- request$context
  setwd(context$code_dir)
  Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = context$data_dir, TRANSCRIPTOME_CONFIG_DIR = request$config_dir)
  execution <- if (is.null(request$execution)) console_execution_lst() else request$execution
  Sys.setenv(TRANSCRIPTOME_CREW_LOG_DIR = file.path(dirname(request_file), "crew"),
    TRANSCRIPTOME_INSPECT_ONLY = if (request$action == "plan") "1" else "0")
  data.table::setDTthreads(execution$threads_per_worker)
  message("开始：", format(Sys.time()), "；操作：", request$action)
  message("执行资源：", execution$workers, " worker × ", execution$threads_per_worker, "线程；memory=", execution$memory,
    "；error=", execution$error, "；reporter=", execution$reporter)
  saveRDS(list(options = execution, available_cores = console_available_cores(), data_table_threads = data.table::getDTthreads(),
    environment = Sys.getenv(names(console_execution_env(execution))), R = R.version.string), file.path(dirname(request_file), "execution.rds"))
  selected <- request$targets
  if (request$action == "plan") {
    manifest <- targets::tar_manifest(callr_function = NULL, drop_missing = FALSE)
    if (length(setdiff(selected, manifest$name))) stop("所选target不在当前工程中：", paste(setdiff(selected, manifest$name), collapse = ", "))
    graph_args <- list(targets_only = TRUE, outdated = execution$graph_outdated, callr_function = NULL, store = context$store, physics = FALSE)
    if (!is.null(selected)) graph_args$names <- substitute(tidyselect::all_of(x), list(x = selected))
    graph <- do.call(targets::tar_visnetwork, graph_args)
    pending <- NULL
    if (execution$graph_outdated) {
      args <- list(callr_function = NULL, store = context$store, reporter = "silent")
      if (!is.null(selected)) args$names <- substitute(tidyselect::all_of(x), list(x = selected))
      pending <- do.call(targets::tar_outdated, args)
    }
    manifest <- as.data.frame(manifest)
    manifest$Selected <- if (is.null(selected)) TRUE else manifest$name %in% selected
    manifest$In_graph <- manifest$name %in% graph$x$nodes$id
    manifest$Outdated <- if (execution$graph_outdated) ifelse(manifest$In_graph, manifest$name %in% pending, NA) else NA
    plan <- list(manifest = manifest, graph = graph, outdated = pending, created = Sys.time(),
      snapshot_dir = request$snapshot_dir, targets = selected)
    saveRDS(plan, file.path(dirname(request_file), "plan.rds"))
    message("依赖预览完成；未执行任何分析target。")
    return(invisible(TRUE))
  }
  config <- console_read_config_lst(request$config_dir)
  report <- console_validate_lst(context, config)
  saveRDS(report, file.path(dirname(request_file), "validation.rds"))
  if (request$action == "build") {
    if (!requireNamespace("targets", quietly = TRUE)) stop("请先安装targets")
    on.exit({
      controller <- targets::tar_option_get("controller")
      if (!is.null(controller)) try(controller$terminate(), silent = TRUE)
    }, add = TRUE)
    manifest <- targets::tar_manifest(callr_function = NULL)
    if (length(setdiff(selected, manifest$name))) stop("所选target不在当前工程中：", paste(setdiff(selected, manifest$name), collapse = ", "))
    message("构建目标：", if (is.null(selected)) "全部" else paste(selected, collapse = ", "))
    args <- list(callr_function = NULL, reporter = execution$reporter, store = context$store, shortcut = execution$shortcut)
    if ("use_crew" %in% names(formals(targets::tar_make))) args$use_crew <- execution$workers > 1L
    if ("as_job" %in% names(formals(targets::tar_make))) args$as_job <- FALSE
    if (!is.null(selected)) args$names <- substitute(tidyselect::all_of(x), list(x = selected))
    # 当前进程已为干净会话
    do.call(targets::tar_make, args)
  }
  message("完成：", format(Sys.time()))
  invisible(TRUE)
}
