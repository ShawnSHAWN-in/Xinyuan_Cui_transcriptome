# 执行资源不参与统计配置
configure_targets_runtime <- function() {
  read_integer <- function(name, default, minimum = 1L) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
    if (length(value) != 1L || !is.finite(value) || value < minimum || value > .Machine$integer.max || value != round(value)) stop("执行参数无效：", name)
    as.integer(value)
  }
  workers <- read_integer("TRANSCRIPTOME_WORKERS", 1L)
  threads <- read_integer("TRANSCRIPTOME_THREADS_PER_WORKER", 1L)
  gc_interval <- read_integer("TRANSCRIPTOME_GC_INTERVAL", 0L, 0L)
  memory <- Sys.getenv("TRANSCRIPTOME_TARGETS_MEMORY", "transient")
  error <- Sys.getenv("TRANSCRIPTOME_TARGETS_ERROR", "stop")
  if (!memory %in% c("transient", "auto", "persistent")) stop("无效的targets内存策略")
  if (!error %in% c("stop", "continue")) stop("无效的targets错误策略")
  data.table::setDTthreads(threads)
  controller <- NULL
  if (workers > 1L && Sys.getenv("TRANSCRIPTOME_INSPECT_ONLY", "0") != "1") {
    if (!requireNamespace("crew", quietly = TRUE)) stop("并行执行需要安装crew")
    log_dir <- Sys.getenv("TRANSCRIPTOME_CREW_LOG_DIR", "")
    if (!nzchar(log_dir)) log_dir <- NULL
    controller <- crew::crew_controller_local(workers = workers,
      options_local = crew::crew_options_local(log_directory = log_dir, log_join = TRUE))
  }
  targets::tar_option_set(controller = controller, memory = memory, garbage_collection = gc_interval, error = error)
  invisible(list(workers = workers, threads_per_worker = threads, memory = memory, garbage_collection = gc_interval, error = error))
}
