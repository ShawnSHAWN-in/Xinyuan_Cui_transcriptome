testthat::test_that("legacy样本表生成原九比较", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
  testthat::expect_equal(nrow(metadata), 36L)
  testthat::expect_length(jobs, 9L)
  expected <- as.vector(outer(c("Cat", "Leopard", "Tiger", "NC"), c("06", "12", "24"), paste, sep = "_"))
  testthat::expect_identical(levels(metadata$Condition_Time), expected)
  testthat::expect_identical(levels(metadata$Time_f), c("06", "12", "24"))
  testthat::expect_setequal(metadata$Time_h, c(6, 12, 24))
  testthat::expect_identical(jobs[[1]]$comparison, "Cat_06_vs_NC_06")
  testthat::expect_identical(jobs[[9]]$comparison, "Tiger_24_vs_NC_24")
  testthat::expect_equal(comparison_manifest_tbl(jobs)$N_control, rep(3L, 9))
})

testthat::test_that("任意分类变量与时间可以组合", {
  metadata <- expand.grid(Sex = c("Female", "Male"), Visit = c("T0", "T1", "T2"), Rep = 1:2, stringsAsFactors = FALSE)
  metadata$Sample <- paste0("sample", seq_len(nrow(metadata)))
  file <- tempfile(fileext = ".tsv")
  write_table_path(metadata, file)
  vars <- list(sample_column = "Sample", legacy_parse = FALSE, variables = list(Sex = list(type = "factor", levels = c("Female", "Male")),
    Visit = list(type = "factor", levels = c("T0", "T1", "T2")), Rep = list(type = "integer")), derived = list(),
    combinations = list(Sex_Visit = list(variables = c("Sex", "Visit"), order = c("Visit", "Sex"), separator = ":")))
  tbl <- read_sample_metadata_tbl(file, vars)
  comparisons <- list(schemes = list(list(id = "sex", enabled = TRUE, type = "pairwise", group = "Sex_Visit", variable = "Sex",
    reference = "Female", treatments = "Male", strata = "Visit")))
  jobs <- make_comparison_jobs_lst(tbl, comparisons)
  testthat::expect_length(jobs, 3L)
  testthat::expect_equal(jobs[[1]]$comparison, "Male:T0_vs_Female:T0")
  testthat::expect_identical(levels(tbl$Sex_Visit), c("Female:T0", "Male:T0", "Female:T1", "Male:T1", "Female:T2", "Male:T2"))
})

testthat::test_that("不同时间比较及显式比较受配置驱动", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  cmp <- config("comparisons")
  cmp$schemes[[1]]$enabled <- FALSE
  cmp$schemes[[2]]$enabled <- TRUE
  jobs <- make_comparison_jobs_lst(metadata, cmp)
  testthat::expect_length(jobs, 8L)
  testthat::expect_equal(jobs[[1]]$comparison, "Cat_12_vs_Cat_06")
  cmp$schemes <- list(list(id = "manual", enabled = TRUE, type = "explicit", group = "Condition_Time", strata = character(),
                          pairs = list(list(numerator = "Tiger_24", denominator = "Cat_06"))))
  testthat::expect_equal(make_comparison_jobs_lst(metadata, cmp)[[1]]$comparison, "Tiger_24_vs_Cat_06")
})

testthat::test_that("缺失水平与不可识别的设计及时失败", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  cmp <- config("comparisons")
  testthat::expect_error(make_comparison_jobs_lst(metadata[as.character(metadata$Condition) != "NC", ], cmp), "缺少")
  cmp$schemes <- c(cmp$schemes[1], cmp$schemes[1])
  testthat::expect_error(make_comparison_jobs_lst(metadata, cmp), "重复")
  testthat::expect_error(cast_variable("UNKNOWN", list(type = "factor", levels = c("A", "B")), "Factor"), "未知水平")
  testthat::expect_error(cast_variable("1.5", list(type = "integer"), "Replicate"), "整数")
})

testthat::test_that("不同分析方案的图像遵循统一变量顺序", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  cmp <- config("comparisons")
  cmp$schemes[[1]]$enabled <- FALSE
  cmp$schemes[[2]]$enabled <- TRUE
  jobs <- make_comparison_jobs_lst(metadata, cmp)
  ordered <- order_comparison_jobs_lst(jobs, c("Time_f", "Condition"))
  testthat::expect_identical(unname(vapply(ordered, `[[`, character(1), "comparison")),
    c("Cat_12_vs_Cat_06", "Leopard_12_vs_Leopard_06", "Tiger_12_vs_Tiger_06", "NC_12_vs_NC_06",
      "Cat_24_vs_Cat_06", "Leopard_24_vs_Leopard_06", "Tiger_24_vs_Tiger_06", "NC_24_vs_NC_06"))
})
