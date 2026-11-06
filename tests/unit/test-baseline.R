testthat::test_that("九个edgeR比较复现工程化前的数值", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
  counts <- read_count_data_lst(fixture_path("counts.tsv"))
  expected <- data.table::fread(fixture_path("legacy_edgeR_expected.tsv"), data.table = FALSE)
  for (job in jobs) {
    fit <- fit_comparison_lst(counts, job, config("analysis")$model, legacy_files)
    ref <- expected[expected$Comparison == job$comparison, ]
    testthat::expect_identical(fit$table$Geneid, ref$Geneid)
    for (field in c("logFC", "p.value", "padj", "logCPM", "baseMean")) {
      testthat::expect_equal(fit$table[[field]], ref[[field]], tolerance = 1e-8, info = paste(job$comparison, field))
    }
  }
})
