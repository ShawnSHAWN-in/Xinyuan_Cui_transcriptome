if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
inputs <- read_module("01_inputs")
check_module_inputs(inputs)
out_dir <- module_dir("03_sample_qc")
combined_dir <- module_dir("09_combined_figures")
design <- inputs$design
counts <- design$get_assay("counts")
sample_info <- restore_factors_tbl(inputs$sample_info)
dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts))
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
keep_pca <- edgeR::filterByExpr(dge, group = sample_info$Group)
if (sum(keep_pca) < 3L) stop("表达过滤后不足3个基因，不能进行样本关系分析")
design$set_assay("log_cpm", log_cpm)
run_project_visualization(design, type = "normalized_matrix", assay_name = "tpm", out_dir = out_dir)
save_table(tibble::rownames_to_column(as.data.frame(log_cpm), "Geneid"), file.path(out_dir, "00_normalized_logCPM.tsv"))
save_table(tibble::rownames_to_column(as.data.frame(edgeR::cpm(dge)), "Geneid"), file.path(out_dir, "00_normalized_CPM.tsv"))
summary <- data.frame(Samples = colnames(counts), CountSum = colSums(counts), DetectedGenes = colSums(counts > 0),
                      ZeroProportion = colMeans(counts == 0), NormalizationFactor = dge$samples$norm.factors) |>
  dplyr::left_join(sample_info, by = "Samples") |> restore_factors_tbl()
summary$Samples <- factor(summary$Samples, levels = sample_info$Samples)
save_table(summary, file.path(out_dir, "01_sample_summary.tsv"))
save_table(data.frame(Geneid = rownames(counts), RetainedForPCA = keep_pca), file.path(out_dir, "02_pca_gene_filter.tsv"))
plots <- list()
for (metric in c("CountSum", "DetectedGenes", "ZeroProportion")) {
  plots[[metric]] <- ggplot2::ggplot(summary, ggplot2::aes(Samples, .data[[metric]], fill = Treatment)) + ggplot2::geom_col() +
    ggplot2::facet_grid(~Time, scales = "free_x", space = "free_x") + ggplot2::scale_fill_manual(values = cfg$treatment_colors) +
    ggplot2::labs(x = NULL, y = metric) + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
  save_plot_pair(plots[[metric]], file.path(out_dir, paste0("QC_", metric)), 14, 5)
}
combine_plots(plots, file.path(combined_dir, "01_sample_overview"), ncol = 1, width = 15, height = 13)
group_colors <- stats::setNames(vapply(cfg$group_levels, function(group) {
  treatment <- sub("_[^_]+$", "", group)
  time_index <- match(sub("^.*_", "", group), cfg$time_levels)
  rgb <- grDevices::colorRamp(c("#FFFFFF", cfg$treatment_colors[[treatment]]))(c(0.55, 0.78, 1)[time_index])
  grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
}, character(1)), cfg$group_levels)

# 复用箱线图模块
boxplots <- list()
for (assay in c("counts", "log_cpm")) {
  p <- run_project_visualization(design, type = "boxplot", assay_name = assay, out_dir = NULL, log2_transform = assay == "counts")
  p$data$Sample <- factor(p$data$Sample, levels = sample_info$Samples)
  p$data$Group <- factor(p$data$Group, levels = cfg$group_levels)
  p <- p + ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  boxplots[[assay]] <- p
  save_plot_pair(p, file.path(out_dir, paste0("expression_boxplot_", assay)), 16, 6)
}
combine_plots(boxplots, file.path(combined_dir, "02_expression_distributions"), 1, 16, 11)

# 四种PCA展示共用坐标
pca_design <- design$filter_genes(rownames(counts)[keep_pca])
pca <- run_pca_style_visualization(pca_design, assay_name = "log_cpm", styles = "all", log2_transform = FALSE,
                                   scale. = cfg$pca_scale, out_dir = NULL)
pca_plots <- lapply(pca$style_results, function(item) {
  p <- item$plot
  if (is.null(p)) return(NULL)
  p$data$Group <- factor(p$data$Group, levels = cfg$group_levels)
  p <- p + ggplot2::scale_color_manual(values = group_colors, drop = FALSE) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  has_fill <- !is.null(p$mapping$fill) || any(vapply(p$layers, function(layer) !is.null(layer$mapping$fill), logical(1)))
  if (has_fill) p <- p + ggplot2::scale_fill_manual(values = group_colors, drop = FALSE)
  p
})
for (name in names(pca_plots)) {
  if (is.null(pca_plots[[name]])) record_status("02_qc", paste0("PCA_", name), "FAILED", pca$style_results[[name]]$error)
  save_plot_pair(pca_plots[[name]], file.path(out_dir, paste0("PCA_", name)), 9, 7)
}
combine_plots(pca_plots, file.path(combined_dir, "03_PCA_four_styles"), 2, 20, 15)
save_table(pca$style_summary, file.path(out_dir, "PCA_style_status.tsv"))
coordinates <- dplyr::left_join(pca$pca_df[, c("Sample", "PC1", "PC2")], sample_info, by = c("Sample" = "Samples")) |> restore_factors_tbl()
axes <- list(x = paste0("PC1 (", pca$variance_pct[1], "%)"), y = paste0("PC2 (", pca$variance_pct[2], "%)"))
main_pca <- ggplot2::ggplot(coordinates, ggplot2::aes(PC1, PC2, color = Treatment, shape = Time)) + ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = cfg$treatment_colors, drop = FALSE) + ggplot2::labs(x = axes$x, y = axes$y, title = "Global PCA: logCPM") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot_pair(main_pca, file.path(out_dir, "PCA_treatment_time"), 9, 7)
save_plot_pair(main_pca + ggplot2::facet_wrap(~Time, nrow = 1), file.path(out_dir, "PCA_global_coordinates_by_time"), 16, 6)
save_table(coordinates, file.path(out_dir, "03_pca_coordinates.tsv"))
save_table(data.frame(PC = seq_along(pca$variance_pct), ExplainedPercent = pca$variance_pct), file.path(out_dir, "04_pca_variance.tsv"))

# 局部PCA另行计算
local_pca <- lapply(cfg$time_levels, function(time) {
  samples <- sample_info$Samples[sample_info$Time == time]
  matrix <- t(log_cpm[keep_pca, samples, drop = FALSE])
  matrix <- matrix[, apply(matrix, 2, stats::sd) > 0, drop = FALSE]
  fit <- stats::prcomp(matrix, scale. = cfg$pca_scale)
  variance <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  data <- data.frame(Samples = rownames(fit$x), PC1 = fit$x[, 1], PC2 = fit$x[, 2]) |> dplyr::left_join(sample_info, by = "Samples")
  ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, color = Treatment)) + ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = cfg$treatment_colors) + ggplot2::labs(title = paste("Independent PCA, time", time),
      x = paste0("PC1 (", round(variance[1], 1), "%)"), y = paste0("PC2 (", round(variance[2], 1), "%)")) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
})
names(local_pca) <- cfg$time_levels
combine_plots(local_pca, file.path(combined_dir, "04_independent_PCA_by_time"), 3, 21, 6)

# 相关性与聚类
correlation <- cor(log_cpm[keep_pca, , drop = FALSE])
save_table(tibble::rownames_to_column(as.data.frame(correlation), "Samples"), file.path(out_dir, "05_sample_correlation.tsv"))
annotation <- as.data.frame(sample_info[, c("Treatment", "Time")])
rownames(annotation) <- sample_info$Samples
cor_plot <- pheatmap::pheatmap(correlation, annotation_col = annotation, annotation_row = annotation,
                               annotation_colors = list(Treatment = cfg$treatment_colors), cluster_rows = FALSE, cluster_cols = FALSE,
                               display_numbers = FALSE, color = grDevices::colorRampPalette(c("#F7FBFF", "#08306B"))(100), silent = TRUE)$gtable
save_plot_pair(cor_plot, file.path(out_dir, "sample_correlation_fixed_order"), 14, 13)
cor_clustered <- run_project_visualization(pca_design, type = "heatmap", assay_name = "log_cpm", out_dir = NULL,
                                          display_numbers = FALSE)
save_plot_pair(cor_clustered, file.path(out_dir, "sample_correlation_clustered"), 14, 13)
dendrogram <- attempt_plot("02_qc", "sample_dendrogram", {
  run_project_visualization(pca_design, type = "dendro", assay_name = "log_cpm", out_dir = NULL) +
    ggplot2::scale_color_manual(values = group_colors) + ggplot2::theme(panel.grid = ggplot2::element_blank())
})
save_plot_pair(dendrogram, file.path(out_dir, "sample_dendrogram"), 12, 12)
combine_plots(list(PCA = main_pca, Correlation = cor_plot, Dendrogram = dendrogram),
               file.path(combined_dir, "05_sample_relationships"), 3, 24, 10)

# 3D使用得分坐标
if (requireNamespace("scatterplot3d", quietly = TRUE) && ncol(pca$pca_res$x) >= 3) {
  draw_3d <- function() {
    scatterplot3d::scatterplot3d(pca$pca_res$x[, 1:3], color = cfg$treatment_colors[as.character(sample_info$Treatment)],
                                 pch = c(16, 17, 15)[as.integer(sample_info$Time)],
                                 xlab = paste0("PC1 (", pca$variance_pct[1], "%)"), ylab = paste0("PC2 (", pca$variance_pct[2], "%)"),
                                 zlab = paste0("PC3 (", pca$variance_pct[3], "%)"))
    graphics::legend("topright", legend = names(cfg$treatment_colors), col = cfg$treatment_colors, pch = 16, cex = 0.8)
  }
  grDevices::pdf(file.path(out_dir, "PCA_3D_scores.pdf"), width = 9, height = 8)
  draw_3d(); grDevices::dev.off()
  grDevices::png(file.path(out_dir, "PCA_3D_scores.png"), width = 9, height = 8, units = "in", res = cfg$dpi)
  draw_3d(); grDevices::dev.off()
} else record_status("02_qc", "PCA_3D", "SKIPPED", "需要scatterplot3d且至少3个主成分")
coverage_plot <- ggplot2::ggplot(inputs$coverage, ggplot2::aes(Source, Proportion)) + ggplot2::geom_col(fill = "#3B75AF", na.rm = TRUE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format()) + ggplot2::labs(y = "Proportion of count-matrix genes", x = NULL) +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot_pair(coverage_plot, file.path(out_dir, "annotation_coverage"), 7, 5)
qc <- list(log_cpm = log_cpm, keep_pca = keep_pca, summary = summary, pca = pca, correlation = correlation,
            plots = list(main_pca = main_pca, correlation = cor_plot, dendrogram = dendrogram))
save_module(qc, "02_qc")
record_status("02_qc", "sample_quality", "SUCCESS", paste(sum(keep_pca), "genes used for PCA"))
