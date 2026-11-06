testthat::test_that("差异图复用旧版Top8标签和配色", {
  colors <- c(Up = "#C94C4C", Down = "#377EB8", Not = "#C9C9C9")
  testthat::expect_identical(legacy_deg_colors(), colors)
  context <- new_legacy_context(legacy_files, list(direction_colors = colors))
  thresholds <- config("analysis")$filter
  tbl <- data.frame(Geneid = paste0("G", 1:20), SYMBOL = paste0("Symbol", 1:20), logFC = seq(-4, 4, length.out = 20),
    padj = seq(0.0001, 0.06, length.out = 20), baseMean = 2^(seq_len(20) / 2), logCPM = seq_len(20) / 2)
  result <- list(job = list(scheme = "baseline", comparison = "A_vs_B"), table = tbl)
  for (type in c("volcano", "ma", "rank")) {
    actual <- plot_deg_plot(result, thresholds, type, context)
    expected <- context$run_deg_visualization(list(deg_tbl = tbl, contrast_names = "A_vs_B"), type = type, label_top_n = 8,
      out_dir = NULL, pvalue_cutoff = thresholds$fdr, logFC_cutoff = thresholds$abs_logfc)
    testthat::expect_identical(actual$data, expected$data, info = type)
    testthat::expect_length(actual$layers, length(expected$layers))
    for (i in seq_along(actual$layers)) {
      testthat::expect_identical(actual$layers[[i]]$data, expected$layers[[i]]$data)
      testthat::expect_identical(actual$layers[[i]]$aes_params, expected$layers[[i]]$aes_params)
    }
    actual_colors <- ggplot2::ggplot_build(actual)$data[[1]]$colour
    expected_colors <- ggplot2::ggplot_build(expected)$data[[1]]$colour
    testthat::expect_identical(actual_colors, expected_colors, info = type)
  }
})

testthat::test_that("QC恢复filterByExpr分组过滤与未缩放PCA", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  counts <- read_count_data_lst(fixture_path("counts.tsv"))
  settings <- read_qc_settings_lst(config("figures"))
  qc <- calc_qc_lst(counts, metadata, settings)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts$counts))
  group <- metadata[[settings$group_by]][match(colnames(counts$counts), metadata$Sample)]
  keep <- edgeR::filterByExpr(dge, group = group)
  matrix <- t(edgeR::cpm(dge, log = TRUE, prior.count = 2)[keep, , drop = FALSE])
  matrix <- matrix[, apply(matrix, 2, stats::sd) > 0, drop = FALSE]
  expected <- stats::prcomp(matrix, center = TRUE, scale. = FALSE)
  testthat::expect_identical(qc$keep_pca, keep)
  testthat::expect_equal(unname(as.matrix(qc$scores[, colnames(expected$x)])), unname(expected$x))
  testthat::expect_equal(qc$correlation, cor(qc$log_cpm[keep, , drop = FALSE]))
  testthat::expect_equal(qc$library$NormalizationFactor, dge$samples$norm.factors)
  metadata <- metadata[rev(seq_len(nrow(metadata))), ]
  shuffled <- calc_qc_lst(counts, metadata, settings)
  testthat::expect_identical(shuffled$keep_pca, qc$keep_pca)
  testthat::expect_equal(shuffled$variance, qc$variance)
})

testthat::test_that("QC映射支持分类变量替换", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  metadata$Treatment_custom <- metadata$Condition
  metadata$Time_custom <- metadata$Time_f
  settings <- list(group_by = "Condition_Time", color_by = "Treatment_custom", shape_by = "Time_custom", strata = "Time_custom")
  style <- config("figures")$style
  style$variable_colors <- list(Treatment_custom = unlist(style$condition_colors))
  actual <- qc_group_colors_vec(metadata, settings, style)
  testthat::expect_identical(actual, qc_group_colors_vec(metadata, read_qc_settings_lst(config("figures")), style))
  settings$shape_by <- "Time_h"
  testthat::expect_error(qc_group_colors_vec(metadata, settings, style), "分类变量")
})

testthat::test_that("四模板与主PCA在列顺序改变后仍共用坐标", {
  metadata <- read_sample_metadata_tbl(fixture_path("samples.tsv"), config("variables"))
  counts <- read_count_data_lst(fixture_path("counts.tsv"))
  counts$counts <- counts$counts[, rev(seq_len(ncol(counts$counts))), drop = FALSE]
  settings <- read_qc_settings_lst(config("figures"))
  qc <- calc_qc_lst(counts, metadata, settings)
  context <- new_legacy_context(legacy_files)
  ordered <- metadata[match(colnames(qc$counts), metadata$Sample), , drop = FALSE]
  design <- new_qc_design(qc, ordered, settings, context)
  testthat::expect_true(inherits(design, "Object_experiment_design"))
  pca <- context$run_pca_style_visualization(design$filter_genes(rownames(qc$counts)[qc$keep_pca]), assay_name = "log_cpm",
    styles = "minimal_points", log2_transform = FALSE, scale. = FALSE, out_dir = NULL)
  expected <- qc$scores[match(pca$pca_df$Sample, qc$scores$Sample), ]
  testthat::expect_equal(pca$pca_df$PC1, expected$PC1)
  testthat::expect_equal(pca$pca_df$PC2, expected$PC2)
})

testthat::test_that("缺失SYMBOL时保留旧版Geneid标签回退", {
  context <- new_legacy_context(legacy_files, list(direction_colors = legacy_deg_colors()))
  result <- list(job = list(scheme = "baseline", comparison = "A_vs_B"),
    table = data.frame(Geneid = c("G1", "G2", "G3"), SYMBOL = c(NA, "", "Symbol3"), logFC = c(3, -2, 1), padj = c(.001, .002, .003)))
  p <- plot_deg_plot(result, config("analysis")$filter, "volcano", context)
  testthat::expect_identical(p$data$SYMBOL, c("G1", "G2", "Symbol3"))
})

testthat::test_that("导出失败保留已有完整图像", {
  prefix <- tempfile("atomic_plot_")
  settings <- list(dpi = 48, max_pixels = 100000)
  files <- save_single_plot_files(empty_result_plot("Existing complete figure"), prefix, settings, 4, 3)
  before <- tools::md5sum(files)
  testthat::expect_error(save_panel_files(function(i) stop("drawing failed"), prefix, settings = settings), "drawing failed")
  testthat::expect_identical(tools::md5sum(files), before)
  broken <- tempfile()
  writeBin(as.raw(rep(0, 32)), broken)
  testthat::expect_error(commit_plot_path(broken, files[2]), "未完整关闭")
  testthat::expect_identical(tools::md5sum(files), before)
})
