testthat::test_that("配置快照稳定且不修改原配置", {
  fixture <- console_fixture()
  source_file <- file.path(fixture$root, "config", "analysis.yml")
  before <- readLines(source_file)
  snapshot <- console_snapshot_path(fixture$context, fixture$config)
  testthat::expect_identical(snapshot, console_snapshot_path(fixture$context, fixture$config))
  modified <- console_effective_config_lst(fixture$context, list(fdr = 0.01, top_terms = 4))
  modified_snapshot <- console_snapshot_path(fixture$context, modified)
  testthat::expect_false(identical(snapshot, modified_snapshot))
  active <- console_stage_config_path(fixture$context, snapshot)
  unchanged <- file.info(file.path(active, "reference.yml"))$mtime
  testthat::expect_identical(active, console_stage_config_path(fixture$context, modified_snapshot))
  testthat::expect_identical(unchanged, file.info(file.path(active, "reference.yml"))$mtime)
  testthat::expect_equal(yaml::read_yaml(file.path(active, "analysis.yml"))$filter$fdr, 0.01)
  testthat::expect_identical(before, readLines(source_file))
  testthat::expect_error(console_effective_config_lst(fixture$context, list(fdr = 2)), "FDR")
  testthat::expect_error(console_effective_config_lst(fixture$context, list(abs_logfc = -1)), "log2FC")
  testthat::expect_error(console_effective_config_lst(fixture$context, list(top_terms = 1.5)), "正整数")
  testthat::expect_error(console_module_targets_vec("missing"), "模块")
  testthat::expect_setequal(console_module_targets_vec("deg_plots"), c("deg_plot_files", "deg_combined_files"))
  source <- paste(readLines(file.path(project_root, "_targets.R")), collapse = "\n")
  target_names <- regmatches(source, gregexpr("tar_target\\([A-Za-z_]+", source))[[1]]
  target_names <- sub("tar_target\\(", "", target_names)
  testthat::expect_true(all(console_module_targets_vec(names(console_modules_lst())) %in% target_names))
})

testthat::test_that("输入校验复用实际分组和比较逻辑", {
  fixture <- console_fixture()
  report <- console_validate_lst(fixture$context, fixture$config)
  testthat::expect_equal(nrow(report$metadata), 36)
  testthat::expect_equal(nrow(report$comparisons), 9)
  testthat::expect_true(all(report$comparisons$N_experiment == 3))
  file <- file.path(fixture$root, "samples.tsv")
  samples <- data.table::fread(file)
  samples$Sample[1] <- "unknown_1"
  data.table::fwrite(samples, file, sep = "\t")
  testthat::expect_error(console_validate_lst(fixture$context, fixture$config), "不匹配")
})

testthat::test_that("日志与结果预览有边界", {
  root <- tempfile()
  dir.create(root)
  log <- file.path(root, "run.log")
  writeLines(rep("中文日志 line", 30000), log)
  text <- console_log_tail(log)
  testthat::expect_lte(length(strsplit(text, "\n")[[1]]), 200)
  testthat::expect_lte(nchar(text, type = "bytes"), 65536)
  outside <- tempfile(fileext = ".tsv")
  writeLines("outside", outside)
  file.symlink(outside, file.path(root, "outside.tsv"))
  writeLines("inside", file.path(root, "inside.tsv"))
  testthat::expect_identical(console_result_files_tbl(root)$Path, "inside.tsv")
  testthat::expect_error(console_result_path(root, "outside.tsv"), "超出")
  testthat::expect_error(console_result_path(root, paste0("../", basename(outside))), "超出")
})

testthat::test_that("后台检查、互斥、失败后可恢复", {
  fixture <- console_fixture()
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config)
  testthat::expect_error(console_start_run(controller, fixture$context, fixture$config), "运行中")
  locked <- callr::r(function(file) is.null(filelock::lock(file, timeout = 0)),
    args = list(file.path(fixture$context$state_dir, "build.lock")))
  testthat::expect_true(locked)
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_equal(nrow(controller$report$metadata), 36)
  config <- fixture$config
  config$paths$count_file <- "missing.tsv"
  console_start_run(controller, fixture$context, config)
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "FAILED")
  console_start_run(controller, fixture$context, fixture$config)
  testthat::expect_identical(wait_console(controller)$state, "SUCCEEDED")
})

testthat::test_that("真实targets构建与缓存复用", {
  fixture <- console_fixture()
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc")
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  file <- file.path(fixture$root, "output", "done.tsv")
  testthat::expect_true(file.exists(file))
  before <- file.info(file)$mtime
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc")
  state <- wait_console(controller)
  testthat::expect_identical(state$state, "SUCCEEDED", info = state$error)
  testthat::expect_identical(before, file.info(file)$mtime)
  testthat::expect_match(console_log_tail(file.path(controller$run_dir, "run.log")), "skip|skipped|up.to.date", perl = TRUE)
})

testthat::test_that("停止会终止计算进程和后代", {
  fixture <- console_fixture()
  child_pid <- file.path(fixture$root, "child.pid")
  lines <- c("list(targets::tar_target(qc_files, {",
    "child <- callr::r_bg(function() Sys.sleep(120), supervise = TRUE)",
    paste0("writeLines(as.character(child$get_pid()), ", deparse(child_pid), ")"),
    "Sys.sleep(120)", "'unused'", "}))")
  writeLines(lines, file.path(fixture$root, "_targets.R"))
  controller <- console_new_controller()
  on.exit(console_stop_run(controller), add = TRUE)
  console_start_run(controller, fixture$context, fixture$config, action = "build", modules = "qc")
  deadline <- Sys.time() + 20
  while (!file.exists(child_pid) && Sys.time() < deadline && controller$process$is_alive()) Sys.sleep(0.05)
  testthat::expect_true(file.exists(child_pid))
  pid <- as.integer(readLines(child_pid))
  console_stop_run(controller)
  testthat::expect_identical(console_poll_lst(controller)$state, "STOPPED")
  testthat::expect_false(controller$process$is_alive())
  Sys.sleep(0.2)
  alive <- tryCatch(ps::ps_is_running(ps::ps_handle(pid)) && ps::ps_status(ps::ps_handle(pid)) != "zombie", error = function(e) FALSE)
  testthat::expect_false(alive)
  console_start_run(controller, fixture$context, fixture$config)
  testthat::expect_identical(wait_console(controller)$state, "SUCCEEDED")
})

testthat::test_that("Shiny页面可渲染且样本表按检查结果展示", {
  fixture <- console_fixture()
  app <- console_app(fixture$root)
  html <- as.character(app$httpHandler(list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = ""))$content)
  testthat::expect_match(html, "转录组分析控制台")
  testthat::expect_match(html, "refresh_results")
  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(data_dir = fixture$root, config_dir = file.path(fixture$root, "config"),
      count_file = "counts.tsv", sample_file = "samples.tsv", gtf_file = "reference.gtf", gene2go_file = "gene2go.tsv",
      output_dir = "output", fdr = 0.05, abs_logfc = 1, top_terms = 12, check = 1)
    deadline <- Sys.time() + 20
    while (app_controller$state == "RUNNING" && Sys.time() < deadline) {
      Sys.sleep(0.1)
      session$elapse(1100)
    }
    testthat::expect_identical(app_controller$state, "SUCCEEDED")
    testthat::expect_match(output$metadata, "Cat_06_1")
    testthat::expect_match(output$comparisons, "Cat_06_vs_NC_06")
    dir.create(file.path(fixture$root, "output"))
    writeLines(c("Gene\tValue", "gene1\t10"), file.path(fixture$root, "output", "preview.tsv"))
    session$setInputs(refresh_results = 1, format = "tsv", search = "preview")
    session$setInputs(result_file = "preview.tsv")
    testthat::expect_match(output$table, "gene1")
    testthat::expect_match(output$result_path, "preview.tsv")
  })
})

testthat::test_that("界面在独立进程启动并可复用关闭", {
  fixture <- console_fixture()
  handle <- launch_console(code_dir = fixture$root, open = FALSE)
  on.exit(stop_console(fixture$root), add = TRUE)
  testthat::expect_true(handle$process$is_alive())
  response <- readLines(url(handle$url), warn = FALSE)
  testthat::expect_true(any(grepl("转录组分析控制台", response, fixed = TRUE)))
  other <- launch_console(code_dir = fixture$root, open = FALSE)
  testthat::expect_identical(handle$process$get_pid(), other$process$get_pid())
  stop_console(fixture$root)
  testthat::expect_false(handle$process$is_alive())
})
