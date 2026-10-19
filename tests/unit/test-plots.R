testthat::test_that("Up与Down共用坐标并有固定偏移", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
  tbl <- data.frame(Scheme = "baseline", Database = "GO", Method = "ORA", Job = jobs[[1]]$id, Comparison = jobs[[1]]$comparison,
    ID = "GO:0006260", Description = "DNA replication", Direction = c("Up", "Down"), p.adjust = c(0.001, 0.002), GeneRatio = c("2/10", "3/12"))
  selection <- list(scheme = "baseline", database = "GO", ids = "GO:0006260", labels = c("GO:0006260" = "DNA replication"))
  p <- plot_comparative_plot(tbl, selection, jobs, "ORA", config("analysis")$filter, config("figures")$display, config("figures")$style)
  testthat::expect_identical(levels(p$data$Direction), c("Up", "Down"))
  testthat::expect_equal(p$data$X, c(0.78, 1.22))
  testthat::expect_equal(p$data$Magnitude, c(0.2, 0.25))
  testthat::expect_equal(sort(unique(p$data$ID)), factor("GO:0006260", levels = "GO:0006260"))
  testthat::expect_silent(ggplot2::ggplot_build(p))
})

testthat::test_that("条目筛选包含方向及方法的代表项", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
  tbl <- data.frame(Scheme = "baseline", Database = "GO", Method = "ORA", Job = jobs[[1]]$id,
    ID = c("UP_TERM", "DOWN_TERM"), Description = c("up", "down"), Direction = c("Up", "Down"), p.adjust = c(1e-10, 0.04))
  display <- config("figures")$display
  display$compare_top_each <- 1L
  selected <- select_display_terms_lst(tbl, config("analysis")$filter, display, jobs)
  testthat::expect_setequal(selected$baseline__GO$ids, c("UP_TERM", "DOWN_TERM"))
  testthat::expect_length(selected$baseline__KEGG$ids, 0L)
})

testthat::test_that("导出限制像素并关闭异常设备", {
  root <- tempfile("plot_")
  settings <- list(dpi = 300, max_pixels = 200000)
  p <- empty_result_plot("Smoke plot")
  paths <- save_single_plot_files(p, root, settings, width = 24, height = 22)
  image <- png::readPNG(paths[2], native = TRUE)
  testthat::expect_lte(prod(dim(image)), settings$max_pixels)
  device <- grDevices::dev.cur()
  testthat::expect_error(save_panel_files(function(i) stop("test failure"), tempfile(), settings = settings), "test failure")
  testthat::expect_identical(grDevices::dev.cur(), device)
})

testthat::test_that("数据库不可用不标为不显著", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
  catalog <- collect_enrichment_catalog_tbl(list(), list(), list())
  backgrounds <- list(GO = list(status = "SUCCESS"), KEGG = list(status = "UNAVAILABLE"))
  selected <- select_display_terms_lst(catalog, config("analysis")$filter, config("figures")$display, jobs, backgrounds)
  testthat::expect_equal(selected$baseline__KEGG$status, "UNAVAILABLE")
  p <- plot_comparative_plot(catalog, selected$baseline__KEGG, jobs, "ORA", config("analysis")$filter,
                             config("figures")$display, config("figures")$style)
  rendered <- ggplot2::ggplot_build(p)$data[[1]]
  testthat::expect_true(any(grepl("UNAVAILABLE", rendered$label)))
})
