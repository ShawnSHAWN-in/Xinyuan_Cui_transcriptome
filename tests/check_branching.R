# 使用完整冒烟测试的目录
root <- Sys.getenv("TARGETS_TEST_ROOT")
if (!nzchar(root) || !file.exists(file.path(root, "first_run_metadata.tsv"))) stop("请将TARGETS_TEST_ROOT设为完整smoke_targets.R输出的TEST_ROOT")
root <- normalizePath(root)
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = root, TRANSCRIPTOME_CONFIG_DIR = "config")
setwd(root)
statistical_rows <- function(meta, prefix) {
  meta[grepl(paste0("^(", paste(prefix, collapse = "|"), ")_[0-9a-f]+$"), meta$name), , drop = FALSE]
}
assert_unchanged <- function(before, after, prefix) {
  rows <- statistical_rows(before, prefix)
  stopifnot(nrow(rows) > 0L)
  idx <- match(rows$name, after$name)
  stopifnot(!anyNA(idx), identical(rows$time, after$time[idx]))
}
# 固定主题后删除单文件，确认仅重建对应绘图分支
before_delete <- targets::tar_meta(fields = c(name, time))
files <- targets::tar_read(deg_plot_files)
removed <- files[grepl("\\.png$", files)][1]
unlink(removed)
targets::tar_make(names = deg_plot_files, callr_function = NULL, reporter = "balanced")
progress <- targets::tar_progress()
stopifnot(file.exists(removed), sum(grepl("^deg_plot_files_[0-9a-f]+$", progress$name) & progress$progress == "completed") == 1L)
assert_unchanged(before_delete, targets::tar_meta(fields = c(name, time)), c("deg_fit", "ora_up", "ora_down", "gsea"))
before <- targets::tar_meta(fields = c(name, time))
analysis <- yaml::read_yaml("config/analysis.yml")
original_analysis <- analysis
analysis$filter$abs_logfc <- 3
yaml::write_yaml(analysis, "config/analysis.yml")
targets::tar_make(names = ora_up, callr_function = NULL, reporter = "balanced")
after <- targets::tar_meta(fields = c(name, time))
assert_unchanged(before, after, c("deg_fit", "gsea"))
outdated <- targets::tar_outdated(callr_function = NULL)
stopifnot(!any(grepl("^gsea(_[0-9a-f]+)?$", outdated)))
yaml::write_yaml(original_analysis, "config/analysis.yml")
comparison <- yaml::read_yaml("config/comparisons.yml")
comparison$schemes[[2]]$enabled <- TRUE
yaml::write_yaml(comparison, "config/comparisons.yml")
targets::tar_make(names = deg_files, callr_function = NULL, reporter = "balanced")
after <- targets::tar_meta(fields = c(name, time))
assert_unchanged(before, after, "deg_fit")
stopifnot(length(targets::tar_read(comparison_jobs)) == 17L, length(targets::tar_read(deg_fit)) == 17L)
writeLines("TARGETS_FILTER_AND_NEW_BRANCH_PASS", "BRANCHING_PASS.txt")
cat("TARGETS_FILTER_AND_NEW_BRANCH_PASS\n")
