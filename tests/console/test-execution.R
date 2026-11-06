testthat::test_that("执行选项受CPU预算及合法取值约束", {
  value <- console_execution_lst(list(workers = 2, threads_per_worker = 2, reporter = "timestamp", error = "continue"), available = 4)
  testthat::expect_equal(value$workers, 2L)
  testthat::expect_equal(value$threads_per_worker, 2L)
  testthat::expect_identical(value$error, "continue")
  testthat::expect_error(console_execution_lst(list(workers = 3, threads_per_worker = 2), available = 4), "可用CPU")
  testthat::expect_error(console_execution_lst(list(workers = 0)), "整数")
  testthat::expect_error(console_execution_lst(list(workers = 2^32)), "整数")
  testthat::expect_error(console_execution_lst(list(workers = 50000, threads_per_worker = 50000), available = 4), "可用CPU")
  testthat::expect_error(console_execution_lst(list(threads_per_worker = 1.5)), "整数")
  testthat::expect_error(console_execution_lst(list(shortcut = "yes")), "TRUE或FALSE")
  testthat::expect_error(console_execution_lst(list(memory = "unknown")), "取值错误")
  testthat::expect_error(console_execution_lst(list(error = "null")), "取值错误")
  testthat::expect_error(console_execution_lst(list(nonsense = 1)), "未知执行选项")
  testthat::expect_identical(console_execution_env(value)[["OMP_NUM_THREADS"]], "2")
})

testthat::test_that("依赖预览不执行分析且不要求输入文件存在", {
  fixture <- console_fixture()
  unlink(file.path(fixture$root, "counts.tsv"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "plan", all = TRUE)
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_s3_class(controller$plan$graph, "visNetwork")
  testthat::expect_true("qc_files" %in% controller$plan$manifest$name)
  testthat::expect_true(all(controller$plan$manifest$In_graph))
  testthat::expect_true(all(is.na(controller$plan$manifest$Outdated)))
  testthat::expect_false(file.exists(file.path(fixture$root, "output", "done.tsv")))
  console_start_run(controller, fixture$context, fixture$config, action = "plan", all = TRUE, execution = list(graph_outdated = TRUE))
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_true("qc_files" %in% controller$plan$outdated)
  testthat::expect_false(file.exists(file.path(fixture$root, "output", "done.tsv")))
})

testthat::test_that("目标子图保留上游且排除无关目标", {
  fixture <- console_fixture()
  writeLines(c("list(targets::tar_target(upstream, 1), targets::tar_target(unrelated, stop('must not run')),",
    "targets::tar_target(qc_files, upstream + 1))"), file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "plan", target_names = "qc_files")
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  in_graph <- controller$plan$manifest$name[controller$plan$manifest$In_graph]
  testthat::expect_setequal(in_graph, c("qc_files", "upstream"))
  console_start_run(controller, fixture$context, fixture$config, action = "build", target_names = "qc_files")
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_equal(targets::tar_read(qc_files, store = fixture$context$store), 2)
  console_start_run(controller, fixture$context, fixture$config, action = "build", target_names = "missing_target")
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "FAILED")
  testthat::expect_match(state$error, "不在当前工程")
})

testthat::test_that("真实工程能够只解析完整依赖图", {
  fixture <- console_fixture()
  file.copy(file.path(project_root, "_targets.R"), file.path(fixture$root, "_targets.R"), overwrite = TRUE)
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "plan", all = TRUE)
  state <- wait_console(controller, timeout = 45)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_true(all(c("deg_fit", "gsea", "qc_files", "summary_files") %in% controller$plan$manifest$name))
  testthat::expect_false(dir.exists(file.path(fixture$root, "output")))
})

testthat::test_that("线程及运行选项传入实际执行且不破坏缓存", {
  testthat::expect_gte(console_available_cores(), 2L)
  fixture <- console_fixture()
  file <- file.path(fixture$root, "thread_probe.rds")
  writeLines(c("sys.source('R/runtime.R', envir = environment())", "configure_targets_runtime()", "list(targets::tar_target(qc_files, {",
    paste0("saveRDS(list(threads = data.table::getDTthreads(), omp = Sys.getenv('OMP_NUM_THREADS')), ", deparse(file), ")"),
    deparse(file), "}, format = 'file'))"), file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc",
    execution = list(threads_per_worker = 2L, reporter = "timestamp", memory = "auto", garbage_collection = 1L))
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_equal(readRDS(file)$threads, 2L)
  testthat::expect_identical(readRDS(file)$omp, "2")
  options <- readRDS(file.path(controller$run_dir, "request.rds"))$execution
  testthat::expect_identical(options$memory, "auto")
  testthat::expect_identical(options$reporter, "timestamp")
  before <- file.info(file)$mtime
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc", execution = list(workers = 2L))
  state <- wait_console(controller, timeout = 45)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_identical(before, file.info(file)$mtime)
})

testthat::test_that("crew真实并行执行两个独立任务", {
  fixture <- console_fixture()
  file <- file.path(fixture$root, "parallel_probe.rds")
  command <- "{ start <- as.numeric(Sys.time()); Sys.sleep(2); list(pid = Sys.getpid(), start = start, end = as.numeric(Sys.time()), threads = data.table::getDTthreads()) }"
  writeLines(c("sys.source('R/runtime.R', envir = environment())", "configure_targets_runtime()", "list(",
    paste0("targets::tar_target(one, ", command, "),"), paste0("targets::tar_target(two, ", command, "),"),
    paste0("targets::tar_target(qc_files, {saveRDS(list(one, two), ", deparse(file), "); ", deparse(file), "}, format = 'file'))")),
    file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc", execution = list(workers = 2L))
  state <- wait_console(controller, timeout = 60)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  results <- readRDS(file)
  testthat::expect_equal(length(unique(vapply(results, `[[`, integer(1), "pid"))), 2L)
  testthat::expect_true(max(vapply(results, `[[`, numeric(1), "start")) < min(vapply(results, `[[`, numeric(1), "end")))
  testthat::expect_true(all(vapply(results, `[[`, integer(1), "threads") == 1L))
})

testthat::test_that("并行停止后worker退出且工程锁可复用", {
  fixture <- console_fixture()
  pid_file <- file.path(fixture$root, "crew.pid")
  writeLines(c("sys.source('R/runtime.R', envir = environment())", "configure_targets_runtime()", "list(targets::tar_target(qc_files, {",
    paste0("writeLines(as.character(Sys.getpid()), ", deparse(pid_file), ")"), "Sys.sleep(120); 'unused'" , "}))"),
    file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc", execution = list(workers = 2L))
  deadline <- Sys.time() + 30
  while (!file.exists(pid_file) && controller$process$is_alive() && Sys.time() < deadline) Sys.sleep(0.1)
  testthat::expect_true(file.exists(pid_file))
  pid <- as.integer(readLines(pid_file))
  console_stop_run(controller)
  deadline <- Sys.time() + 5
  is_alive <- function() tryCatch(ps::ps_is_running(ps::ps_handle(pid)) && ps::ps_status(ps::ps_handle(pid)) != "zombie", error = function(e) FALSE)
  while (is_alive() && Sys.time() < deadline) Sys.sleep(0.1)
  testthat::expect_false(is_alive())
  testthat::expect_identical(console_poll_lst(controller)$state, "STOPPED")
  console_start_run(controller, fixture$context, fixture$config, action = "check")
  testthat::expect_identical(wait_console(controller)$state, "SUCCEEDED")
})

testthat::test_that("继续策略运行独立任务并保留失败状态", {
  fixture <- console_fixture()
  file <- file.path(fixture$root, "survived.txt")
  writeLines(c("sys.source('R/runtime.R', envir = environment())", "configure_targets_runtime()",
    "list(targets::tar_target(bad, stop('expected_failure')),", paste0("targets::tar_target(qc_files, {writeLines('ok', ", deparse(file),
      "); ", deparse(file), "}, format = 'file'))")), file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", all = TRUE, execution = list(error = "continue"))
  state <- wait_console(controller)
  testthat::expect_true(file.exists(file))
  testthat::expect_identical(state$state, "FAILED")
  testthat::expect_match(state$error, "bad")
  console_start_run(controller, fixture$context, fixture$config, action = "build", target_names = "qc_files", execution = list(error = "continue"))
  testthat::expect_identical(wait_console(controller)$state, "SUCCEEDED")
})

testthat::test_that("运行页面提供执行设置和依赖图入口", {
  fixture <- console_fixture()
  app <- console_app(fixture$root)
  html <- app$httpHandler(list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = ""))$content
  testthat::expect_match(html, "threads_per_worker")
  testthat::expect_match(html, "refresh_graph")
  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(data_dir = fixture$root, config_dir = file.path(fixture$root, "config"),
      count_file = "counts.tsv", sample_file = "samples.tsv", gtf_file = "reference.gtf", gene2go_file = "gene2go.tsv",
      output_dir = "output", fdr = 0.05, abs_logfc = 1, top_terms = 12, graph_all = TRUE, refresh_graph = 1)
    deadline <- Sys.time() + 30
    while (app_controller$state == "RUNNING" && Sys.time() < deadline) {
      Sys.sleep(0.1)
      session$elapse(1100)
    }
    testthat::expect_identical(app_controller$state, "SUCCEEDED")
    testthat::expect_match(output$manifest, "qc_files")
    testthat::expect_match(output$plan_note, "图中target数")
  })
})
