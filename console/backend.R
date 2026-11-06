# 控制台与计算隔离
console_config_names <- c("paths", "variables", "comparisons", "analysis", "reference", "enrichment", "figures")

console_modules_lst <- function() {
  list(
    input = list(label = "输入与注释", targets = c("input_files", "annotation_files")),
    qc = list(label = "完整 QC", targets = "qc_files"),
    deg = list(label = "差异结果表", targets = "deg_files"),
    deg_plots = list(label = "火山图 / MA / rank", targets = c("deg_plot_files", "deg_combined_files")),
    sets = list(label = "基因集合与表达热图", targets = c("gene_set_files", "expression_heatmap_files")),
    enrichment = list(label = "ORA / GSEA 结果", targets = c("ora_up_files", "ora_down_files", "gsea_files")),
    enrichment_plots = list(label = "富集单图", targets = c("ora_plot_files", "gsea_plot_files")),
    comparison = list(label = "比较富集与组合", targets = c("comparative_files", "comparative_combined_files", "curve_combined_files")),
    pathview = list(label = "Pathview（须配置开启）", targets = "pathview_files")
  )
}

console_module_targets_vec <- function(modules, all = FALSE) {
  if (isTRUE(all)) return(NULL)
  catalog <- console_modules_lst()
  if (!length(modules) || any(!modules %in% names(catalog))) stop("请选择有效模块")
  unique(unlist(lapply(catalog[modules], `[[`, "targets"), use.names = FALSE))
}

console_absolute_path <- function(path, root = getwd(), must_work = TRUE) {
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path))) stop("路径不能为空")
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

console_context_lst <- function(code_dir, data_dir = code_dir, config_dir = "config") {
  code_dir <- console_absolute_path(code_dir)
  data_dir <- console_absolute_path(data_dir, code_dir)
  config_dir <- console_absolute_path(config_dir, code_dir)
  if (!file.exists(file.path(code_dir, "_targets.R"))) stop("代码目录缺少 _targets.R")
  list(code_dir = code_dir, data_dir = data_dir, config_dir = config_dir, state_dir = file.path(code_dir, ".console"),
       store = file.path(code_dir, "_targets"))
}

console_read_config_lst <- function(config_dir) {
  stats::setNames(lapply(console_config_names, function(name) {
    file <- file.path(config_dir, paste0(name, ".yml"))
    value <- yaml::read_yaml(file, eval.expr = FALSE)
    if (!is.list(value) || is.null(names(value))) stop("配置须为映射：", file)
    value
  }), console_config_names)
}

console_effective_config_lst <- function(context, overrides = list()) {
  config <- console_read_config_lst(context$config_dir)
  for (field in intersect(names(overrides), names(config$paths))) config$paths[[field]] <- overrides[[field]]
  if (!is.null(overrides$fdr)) config$analysis$filter$fdr <- as.numeric(overrides$fdr)
  if (!is.null(overrides$abs_logfc)) config$analysis$filter$abs_logfc <- as.numeric(overrides$abs_logfc)
  if (!is.null(overrides$top_terms)) config$figures$display$top_terms <- as.numeric(overrides$top_terms)
  fdr <- config$analysis$filter$fdr
  logfc <- config$analysis$filter$abs_logfc
  top <- config$figures$display$top_terms
  if (length(fdr) != 1L || !is.finite(fdr) || fdr <= 0 || fdr > 1) stop("FDR须在 (0, 1] 内")
  if (length(logfc) != 1L || !is.finite(logfc) || logfc < 0) stop("|log2FC|阈值须为非负数")
  if (length(top) != 1L || !is.finite(top) || top < 1 || top != round(top)) stop("展示条目数须为正整数")
  config
}

console_snapshot_path <- function(context, config) {
  # 相同配置复用同一路径
  content <- lapply(config, yaml::as.yaml)
  key <- digest::digest(content, algo = "sha256")
  path <- file.path(context$state_dir, "configs", key)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (name in names(content)) {
    file <- file.path(path, paste0(name, ".yml"))
    if (file.exists(file)) {
      if (!identical(readLines(file, warn = FALSE), strsplit(content[[name]], "\n", fixed = TRUE)[[1]])) stop("配置快照被修改：", file)
    } else writeLines(content[[name]], file, sep = "", useBytes = TRUE)
  }
  normalizePath(path)
}

console_validate_lst <- function(context, config, check_packages = TRUE) {
  helpers <- new.env(parent = globalenv())
  for (name in c("io", "metadata", "design")) sys.source(file.path(context$code_dir, "R", paste0(name, ".R")), envir = helpers)
  inputs <- c("count_file", "sample_file", "gtf_file", "gene2go_file")
  paths <- vapply(config$paths[inputs], console_absolute_path, character(1), root = context$data_dir)
  if (length(paths) != length(inputs) || any(dir.exists(paths)) || any(file.access(paths, 4) != 0)) stop("输入文件缺失或不可读")
  message("检查样本元数据与比较方案")
  metadata <- helpers$read_sample_metadata_tbl(paths[["sample_file"]], config$variables)
  header <- names(data.table::fread(paths[["count_file"]], nrows = 0L, data.table = FALSE, check.names = FALSE))
  if (!all(c("Geneid", "Length") %in% header) || anyDuplicated(header)) stop("count表须有Geneid、Length及唯一的样本列")
  count_samples <- setdiff(header, c("Geneid", "Length"))
  if (!setequal(metadata$Sample, count_samples)) stop("样本表与count表头不匹配；count独有：", paste(setdiff(count_samples, metadata$Sample), collapse = ", "),
    "；样本表独有：", paste(setdiff(metadata$Sample, count_samples), collapse = ", "))
  jobs <- helpers$make_comparison_jobs_lst(metadata, config$comparisons)
  comparisons <- helpers$comparison_manifest_tbl(jobs)
  min_n <- config$analysis$model$min_replicates
  if (!is.numeric(min_n) || length(min_n) != 1L || !is.finite(min_n) || min_n < 2) stop("min_replicates须至少为2")
  if (any(comparisons$N_experiment < min_n | comparisons$N_control < min_n)) stop("比较的生物学重复数不足")
  if (!identical(config$analysis$model$method, "edgeR") || !identical(config$analysis$model$strategy, "legacy_pairwise")) stop("当前仅支持edgeR / legacy_pairwise")
  output <- console_absolute_path(config$paths$output_dir, context$data_dir, must_work = FALSE)
  parent <- output
  while (!dir.exists(parent) && dirname(parent) != parent) parent <- dirname(parent)
  if (file.access(parent, 2) != 0) stop("输出目录的父目录不可写：", parent)
  packages <- c("shiny", "callr", "yaml", "filelock", "digest", "data.table", "tibble", "dplyr", "targets", "edgeR", "ggplot2",
    "clusterProfiler", "enrichplot", "GO.db", "KEGGREST", "patchwork", "pheatmap", "rstudioapi")
  installed <- vapply(packages, function(name) nzchar(system.file(package = name)), logical(1))
  package_table <- data.frame(Package = packages, Installed = installed)
  message("输入检查通过：", nrow(metadata), "个样本；", nrow(comparisons), "个比较。count数值由分析模块完整检查。")
  if (check_packages && any(!installed)) message("部分分析依赖未安装，见环境表；未自动安装软件包。")
  list(metadata = as.data.frame(metadata), comparisons = as.data.frame(comparisons), packages = package_table,
       inputs = data.frame(Input = names(paths), Path = unname(paths)), output_dir = output,
       note = "检查读取count表头，不加载完整矩阵；完整数值和注释内容由targets分析时验证。")
}

console_new_controller <- function() {
  controller <- new.env(parent = emptyenv())
  controller$process <- NULL
  controller$lock <- NULL
  controller$run_dir <- NULL
  controller$state <- "IDLE"
  controller$error <- ""
  controller$context <- NULL
  controller$report <- NULL
  controller
}

console_poll_lst <- function(controller) {
  process <- controller$process
  if (!is.null(process) && !process$is_alive() && controller$state == "RUNNING") {
    result <- tryCatch(process$get_result(), error = function(e) e)
    controller$state <- if (inherits(result, "error")) "FAILED" else "SUCCEEDED"
    controller$error <- if (inherits(result, "error")) conditionMessage(result) else ""
    file <- file.path(controller$run_dir, "validation.rds")
    if (file.exists(file)) controller$report <- readRDS(file)
    saveRDS(list(state = controller$state, error = controller$error, finished = Sys.time()), file.path(controller$run_dir, "result.rds"))
    if (!is.null(controller$lock)) filelock::unlock(controller$lock)
    controller$lock <- NULL
  }
  list(state = controller$state, error = controller$error, run_dir = controller$run_dir,
       pid = if (is.null(process)) NA_integer_ else process$get_pid())
}

console_start_run <- function(controller, context, config, action = "check", modules = character(), all = FALSE) {
  console_poll_lst(controller)
  if (controller$state == "RUNNING") stop("本工程已有运行中的任务")
  if (!action %in% c("check", "build")) stop("未知操作")
  selected <- if (action == "build") console_module_targets_vec(modules, all) else character()
  dir.create(context$state_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(context$state_dir, "build.lock"), timeout = 0)
  if (is.null(lock)) stop("另一控制台正在使用本工程；请等待任务完成")
  success <- FALSE
  on.exit(if (!success) filelock::unlock(lock), add = TRUE)
  snapshot <- console_snapshot_path(context, config)
  run_dir <- tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = file.path(context$state_dir, "runs"))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  request <- list(context = context, config_dir = snapshot, action = action, targets = selected, created = Sys.time())
  saveRDS(request, file.path(run_dir, "request.rds"))
  process <- callr::r_bg(function(code_dir, request_file) {
    env <- new.env(parent = globalenv())
    sys.source(file.path(code_dir, "console", "backend.R"), envir = env)
    sys.source(file.path(code_dir, "console", "worker.R"), envir = env)
    env$console_worker(request_file)
  }, args = list(context$code_dir, file.path(run_dir, "request.rds")), wd = context$code_dir, libpath = .libPaths(),
    stdout = file.path(run_dir, "run.log"), stderr = "2>&1", supervise = TRUE, cleanup_tree = TRUE, user_profile = FALSE)
  controller$process <- process
  controller$lock <- lock
  controller$run_dir <- run_dir
  controller$context <- context
  controller$state <- "RUNNING"
  controller$error <- ""
  controller$report <- NULL
  success <- TRUE
  invisible(run_dir)
}

console_stop_run <- function(controller) {
  console_poll_lst(controller)
  if (controller$state != "RUNNING") return(invisible(FALSE))
  controller$process$kill_tree()
  controller$process$wait(timeout = 2000)
  if (controller$process$is_alive()) stop("进程尚未停止，暂不释放工程锁")
  controller$state <- "STOPPED"
  controller$error <- "用户停止；已完成的targets缓存保留。"
  saveRDS(list(state = controller$state, error = controller$error, finished = Sys.time()), file.path(controller$run_dir, "result.rds"))
  filelock::unlock(controller$lock)
  controller$lock <- NULL
  invisible(TRUE)
}

console_log_tail <- function(file, max_bytes = 65536L, max_lines = 200L) {
  if (!file.exists(file)) return("尚无日志")
  con <- file(file, "rb")
  on.exit(close(con))
  size <- file.info(file)$size
  start <- max(0, size - max_bytes)
  seek(con, where = start, origin = "start")
  bytes <- readBin(con, "raw", n = min(size, max_bytes))
  if (start > 0 && any(bytes == as.raw(10))) bytes <- bytes[-seq_len(which(bytes == as.raw(10))[1])]
  text <- iconv(rawToChar(bytes), from = "UTF-8", to = "UTF-8", sub = "")
  paste(tail(strsplit(text, "\n", fixed = TRUE)[[1]], max_lines), collapse = "\n")
}

console_result_files_tbl <- function(output_dir, max_files = 3000L) {
  if (!dir.exists(output_dir)) return(data.frame(Path = character(), Format = character(), Bytes = numeric()))
  root <- normalizePath(output_dir)
  files <- list.files(root, pattern = "\\.(png|pdf|tsv|csv|txt|md)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  real <- normalizePath(files, mustWork = FALSE)
  keep <- startsWith(real, paste0(root, "/")) & !dir.exists(files)
  files <- head(files[keep], max_files)
  data.frame(Path = substring(files, nchar(root) + 2L), Format = tolower(tools::file_ext(files)), Bytes = file.info(files)$size)
}

console_result_path <- function(output_dir, relative_path) {
  root <- normalizePath(output_dir, mustWork = TRUE)
  path <- normalizePath(file.path(root, relative_path), mustWork = TRUE)
  if (!startsWith(path, paste0(root, "/")) || dir.exists(path)) stop("结果文件超出输出目录")
  path
}
