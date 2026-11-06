if (!all(file.exists(file.path("tests", "fixtures", c("legacy_edgeR_expected.tsv", "legacy_enrichment_expected.tsv"))))) {
  stop("请先在代码根目录执行Rscript tests/generate_baselines.R，生成本地合成回归基线")
}
project <- normalizePath(getwd())
source("tests/setup_fixture.R")
root <- Sys.getenv("TARGETS_TEST_ROOT", unset = tempfile("targets_smoke_"))
root <- setup_targets_fixture(project, root, quick = "--quick" %in% commandArgs(TRUE))
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = root, TRANSCRIPTOME_CONFIG_DIR = "config")
setwd(root)
targets::tar_make(callr_function = NULL, reporter = "balanced")
meta <- targets::tar_meta(fields = c(name, time, data, error))
write.table(meta, "first_run_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(length(targets::tar_read(comparison_jobs)) == 9L, length(targets::tar_read(deg_fit)) == 9L,
          length(targets::tar_read(ora_up)) == 18L, length(targets::tar_read(ora_down)) == 18L, length(targets::tar_read(gsea)) == 18L)
expected <- data.table::fread(file.path(project, "tests", "fixtures", "legacy_enrichment_expected.tsv"), data.table = FALSE)
results <- c(targets::tar_read(ora_up), targets::tar_read(ora_down), targets::tar_read(gsea))
for (result in results) {
  kind <- if (result$method == "GSEA") "GSEA" else paste0("ORA_", result$direction)
  ref <- expected[expected$Comparison == result$job$comparison & expected$Database == result$database & expected$Kind == kind, , drop = FALSE]
  tbl <- result$table
  stopifnot(setequal(tbl$ID, ref$ID))
  idx <- match(tbl$ID, ref$ID)
  for (column in intersect(c("pvalue", "p.adjust", "NES"), names(tbl))) {
    stopifnot(max(abs(tbl[[column]] - ref[[column]][idx]), na.rm = TRUE) < 1e-7)
  }
}
writeLines("ALL_ENRICHMENT_BASELINE_PASS", "ENRICHMENT_BASELINE_PASS.txt")
stopifnot(!length(targets::tar_outdated(callr_function = NULL)))
# 原样再跑不重算
before <- targets::tar_meta(fields = c(name, time))
targets::tar_make(callr_function = NULL, reporter = "balanced")
after <- targets::tar_meta(fields = c(name, time))
stopifnot(all(targets::tar_progress()$progress == "skipped"))
# 主题只影响图像
figures <- yaml::read_yaml("config/figures.yml")
figures$style$title_size <- figures$style$title_size + 1
yaml::write_yaml(figures, "config/figures.yml")
outdated <- targets::tar_outdated(callr_function = NULL)
stopifnot(!any(grepl("^(deg_fit|ora_up|ora_down|gsea)(_[0-9a-f]+)?$", outdated)))
targets::tar_make(names = qc_files, callr_function = NULL, reporter = "balanced")
updated <- targets::tar_meta(fields = c(name, time))
checks <- grep("^(deg_fit|ora_up|ora_down|gsea)_[0-9a-f]+$", before$name)
idx <- match(before$name[checks], updated$name)
stopifnot(identical(before$time[checks], updated$time[idx]))
# 删除图像后只重建文件任务
figures$style$title_size <- figures$style$title_size - 1
yaml::write_yaml(figures, "config/figures.yml")
targets::tar_make(names = deg_plot_files, callr_function = NULL, reporter = "balanced")
files <- targets::tar_read(deg_plot_files)
removed_file <- files[grepl("\\.png$", files)][1]
unlink(removed_file)
stat_before <- targets::tar_meta(fields = c(name, time))
targets::tar_make(names = deg_plot_files, callr_function = NULL, reporter = "balanced")
stat_after <- targets::tar_meta(fields = c(name, time))
stopifnot(file.exists(removed_file))
progress <- targets::tar_progress()
stopifnot(sum(grepl("^deg_plot_files_[0-9a-f]+$", progress$name) & progress$progress == "completed") == 1L)
checks <- grep("^(deg_fit|ora_up|ora_down|gsea)_[0-9a-f]+$", stat_before$name)
stopifnot(identical(stat_before$time[checks], stat_after$time[match(stat_before$name[checks], stat_after$name)]))
writeLines("TARGETS_SMOKE_INCREMENTAL_PASS", "PASS.txt")
cat("TARGETS_SMOKE_INCREMENTAL_PASS\nTEST_ROOT=", root, "\n", sep = "")
