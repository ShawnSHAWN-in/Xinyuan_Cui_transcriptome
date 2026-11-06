read_qc_settings_lst <- function(config) {
  defaults <- list(group_by = "Condition_Time", color_by = "Condition", shape_by = "Time_f", strata = "Time_f", pca_scale = FALSE)
  if (is.null(config$qc)) defaults else utils::modifyList(defaults, config$qc)
}

qc_group_colors_vec <- function(metadata, settings, style) {
  required <- c(settings$group_by, settings$color_by, settings$shape_by, as.character(unlist(settings$strata)))
  if (!all(required %in% names(metadata))) stop("QC映射含未知变量：", paste(setdiff(required, names(metadata)), collapse = ", "))
  if (any(vapply(metadata[required], is.numeric, logical(1)))) stop("QC分组、颜色、形状和分层请使用分类变量；数值时间请配套声明Time_f")
  levels_of <- function(x) if (is.factor(x)) levels(droplevels(x)) else unique(as.character(x))
  groups <- levels_of(metadata[[settings$group_by]])
  palette <- variable_palette(style, settings$group_by, metadata[[settings$group_by]])
  if (!is.null(palette)) return(palette[groups])
  color_levels <- levels_of(metadata[[settings$color_by]])
  palette <- variable_palette(style, settings$color_by, metadata[[settings$color_by]])
  if (is.null(palette)) palette <- stats::setNames(scales::hue_pal()(length(color_levels)), color_levels)
  shapes <- levels_of(metadata[[settings$shape_by]])
  shade <- if (length(shapes) == 3L) c(0.55, 0.78, 1) else if (length(shapes) == 1L) 1 else seq(0.55, 1, length.out = length(shapes))
  fallback <- stats::setNames(scales::hue_pal()(length(groups)), groups)
  stats::setNames(vapply(groups, function(group) {
    block <- metadata[as.character(metadata[[settings$group_by]]) == group, , drop = FALSE]
    color <- unique(as.character(block[[settings$color_by]]))
    shape <- unique(as.character(block[[settings$shape_by]]))
    if (length(color) != 1L) return(fallback[[group]])
    if (length(shape) != 1L) return(palette[[color]])
    rgb <- grDevices::colorRamp(c("#FFFFFF", palette[[color]]))(shade[match(shape, shapes)])
    grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  }, character(1)), groups)
}

new_qc_design <- function(qc, metadata, settings, context) {
  # 仅补充旧绘图接口所需分组
  expression_class <- R6::R6Class("qc_expression", inherit = context$object_experiment_design,
    public = list(initialize = context$read_counts$public_methods$initialize))
  samples <- metadata$Sample
  design <- expression_class$new(count_matrix = qc$counts[, samples, drop = FALSE], gene_length = qc$lengths)
  design$GroupSamples <- as.matrix(data.frame(Samples = samples, Group = as.character(metadata[[settings$group_by]])))
  for (assay in c("log_cpm", "cpm", "tpm", "fpkm")) design$set_assay(assay, qc[[assay]][, samples, drop = FALSE])
  design
}

qc_coverage_tbl <- function(qc, annotation, backgrounds) {
  anno <- annotation$annotation
  usable <- anno[!duplicated(anno$Geneid) & !duplicated(anno$Geneid, fromLast = TRUE), , drop = FALSE]
  duplicate_ids <- unique(usable$ENTREZID[!is.na(usable$ENTREZID) &
    (duplicated(usable$ENTREZID) | duplicated(usable$ENTREZID, fromLast = TRUE))])
  usable <- usable[!usable$ENTREZID %in% duplicate_ids, , drop = FALSE]
  usable <- usable[usable$Geneid %in% rownames(qc$counts), , drop = FALSE]
  kegg_ok <- backgrounds$KEGG$status == "SUCCESS"
  result <- tibble::tibble(Source = c("NCBI_GeneID_any", "NCBI_GeneID_unambiguous", "GO_BP_unambiguous", "KEGG_unambiguous"),
    Annotated = c(sum(rownames(qc$counts) %in% anno$Geneid[!is.na(anno$ENTREZID)]), sum(!is.na(usable$ENTREZID)),
      sum(usable$ENTREZID %in% backgrounds$GO$background$TERM2GENE$GENE),
      if (kegg_ok) sum(usable$ENTREZID %in% backgrounds$KEGG$background$TERM2GENE$GENE) else NA_integer_),
    Status = c("SUCCESS", "SUCCESS", backgrounds$GO$status, backgrounds$KEGG$status), Total = nrow(qc$counts))
  result$Proportion <- result$Annotated / result$Total
  result
}

save_base_plot_files <- function(draw, prefix, settings, width = 9, height = 8) {
  paths <- paste0(prefix, c(".pdf", ".png"))
  dpi <- min(settings$dpi, floor(sqrt(settings$max_pixels / (width * height))))
  if (!is.finite(dpi) || dpi < 1) stop("3D画布像素预算无效")
  for (i in seq_along(paths)) {
    temporary <- tempfile(".plot_", tmpdir = dirname(prefix), fileext = ".part")
    on.exit(unlink(temporary), add = TRUE)
    if (i == 1L) grDevices::pdf(temporary, width = width, height = height) else {
      grDevices::png(temporary, width = width, height = height, units = "in", res = dpi)
    }
    device <- grDevices::dev.cur()
    tryCatch(draw(), finally = grDevices::dev.off(device))
    commit_plot_path(temporary, paths[i])
  }
  paths
}

export_qc_files <- function(qc, metadata, display, style, settings, output_dir, engine_files, qc_settings, annotation, backgrounds) {
  check_packages(c("ggrepel", "ggpubr", "reshape2", "pheatmap", "ggdendro"))
  dir <- file.path(output_dir, "03_sample_qc")
  combined_dir <- file.path(output_dir, "09_combined_figures")
  metadata <- sort_metadata_tbl(metadata, display$order_by)
  groups <- qc_group_colors_vec(metadata, qc_settings, style)
  context <- new_legacy_context(engine_files)
  design_metadata <- metadata[match(colnames(qc$counts), metadata$Sample), , drop = FALSE]
  design <- new_qc_design(qc, design_metadata, qc_settings, context)
  pca_design <- design$filter_genes(rownames(qc$counts)[qc$keep_pca])
  qc_display <- list(color_by = qc_settings$color_by, shape_by = qc_settings$shape_by, facet_by = character())
  color_levels <- if (is.factor(metadata[[qc_settings$color_by]])) levels(droplevels(metadata[[qc_settings$color_by]])) else {
    unique(as.character(metadata[[qc_settings$color_by]]))
  }
  palette <- variable_palette(style, qc_settings$color_by, metadata[[qc_settings$color_by]])
  if (is.null(palette)) palette <- stats::setNames(scales::hue_pal()(length(color_levels)), color_levels)
  statuses <- list()
  paths <- character()
  save_plot <- function(plot, name, width = 9, height = 7, combined = FALSE) {
    save_single_plot_files(plot, file.path(if (combined) combined_dir else dir, name), settings, width, height)
  }
  save_matrix <- function(matrix, name) write_table_path(tibble::rownames_to_column(as.data.frame(matrix), "Geneid"), file.path(dir, name))
  summary <- dplyr::left_join(qc$library, metadata, by = "Sample")
  summary <- summary[match(metadata$Sample, summary$Sample), , drop = FALSE]
  summary$Sample <- factor(summary$Sample, levels = metadata$Sample)
  summary$CountSum <- summary$LibrarySize
  summary$DetectedGenes <- summary$Detected
  paths <- c(paths, write_table_path(summary, file.path(dir, "01_sample_summary.tsv")),
    write_table_path(summary, file.path(dir, "library_statistics.tsv")),
    write_table_path(data.frame(Geneid = rownames(qc$counts), RetainedForPCA = qc$keep_pca), file.path(dir, "02_pca_gene_filter.tsv")))
  for (assay in c("counts", "cpm", "log_cpm", "tpm", "fpkm")) paths <- c(paths, save_matrix(qc[[assay]], paste0(assay, ".tsv")))

  # 文库与检出情况
  overview <- lapply(c("CountSum", "DetectedGenes", "ZeroProportion", "NormalizationFactor"), function(metric) {
    p <- ggplot2::ggplot(summary, ggplot2::aes(Sample, .data[[metric]], fill = .data[[qc_settings$color_by]])) +
      ggplot2::geom_col() + ggplot2::scale_fill_manual(values = palette, drop = FALSE) + ggplot2::labs(x = NULL, y = metric) +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    strata <- as.character(unlist(qc_settings$strata))
    if (length(strata)) p <- p + ggplot2::facet_grid(cols = ggplot2::vars(!!!rlang::syms(strata)), scales = "free_x", space = "free_x")
    p
  })
  names(overview) <- c("CountSum", "DetectedGenes", "ZeroProportion", "NormalizationFactor")
  for (name in names(overview)) paths <- c(paths, save_plot(overview[[name]], paste0("QC_", name), 14, 5))
  paths <- c(paths, save_plot(overview$CountSum, "02_library_size", 14, 5),
    save_panel_files(function(i) overview[[i]], file.path(combined_dir, "01_sample_overview"), 4, 1, 15, 17, settings))
  rm(overview)

  # 复用旧版表达箱线图
  assays <- c("counts", "log_cpm", "tpm", "fpkm")
  draw_boxplot <- function(i) {
    assay <- assays[i]
    p <- context$run_project_visualization(design, type = "boxplot", assay_name = assay, out_dir = NULL, log2_transform = assay != "log_cpm")
    p$data$Sample <- factor(p$data$Sample, levels = metadata$Sample)
    p$data$Group <- factor(p$data$Group, levels = names(groups))
    p + ggplot2::scale_fill_manual(values = groups, drop = FALSE) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  for (i in seq_along(assays)) paths <- c(paths, save_plot(draw_boxplot(i), paste0("expression_boxplot_", assays[i]), 16, 6))
  paths <- c(paths, save_plot(draw_boxplot(2), "03_expression_boxplot", 16, 6),
    save_panel_files(draw_boxplot, file.path(combined_dir, "02_expression_distributions"), 4, 1, 16, 22, settings))

  # 四种PCA共用过滤和坐标
  pca <- context$run_pca_style_visualization(pca_design, assay_name = "log_cpm", styles = "all", log2_transform = FALSE,
    scale. = isTRUE(qc_settings$pca_scale), out_dir = NULL)
  pca_plots <- lapply(pca$style_results, function(item) {
    if (is.null(item$plot)) return(empty_result_plot("PCA style unavailable", item$error))
    p <- item$plot
    p$data$Group <- factor(p$data$Group, levels = names(groups))
    p <- p + ggplot2::scale_color_manual(values = groups, drop = FALSE) + ggplot2::theme(panel.grid = ggplot2::element_blank())
    has_fill <- !is.null(p$mapping$fill) || any(vapply(p$layers, function(layer) !is.null(layer$mapping$fill), logical(1)))
    if (has_fill) p <- p + ggplot2::scale_fill_manual(values = groups, drop = FALSE)
    p
  })
  for (name in names(pca_plots)) paths <- c(paths, save_plot(pca_plots[[name]], paste0("PCA_", name)))
  paths <- c(paths, save_panel_files(function(i) pca_plots[[i]], file.path(combined_dir, "03_PCA_four_styles"), 4, 2, 20, 15, settings),
    write_table_path(pca$style_summary, file.path(dir, "PCA_style_status.tsv")))
  statuses[["PCA_styles"]] <- data.frame(Analysis = paste0("PCA_", pca$style_summary$style_name),
    Status = toupper(pca$style_summary$status), Detail = pca$style_summary$error)
  rm(pca_plots, pca)
  main_pca <- plot_pca_plot(qc, metadata, qc_display, style)
  paths <- c(paths, save_plot(main_pca, "PCA_treatment_time"), save_plot(main_pca, "01_PCA"),
    write_table_path(qc$scores, file.path(dir, "03_pca_coordinates.tsv")), write_table_path(qc$scores, file.path(dir, "PCA_scores.tsv")),
    write_table_path(data.frame(PC = seq_along(qc$variance), ExplainedPercent = qc$variance), file.path(dir, "04_pca_variance.tsv")))
  strata <- as.character(unlist(qc_settings$strata))
  if (length(strata)) paths <- c(paths, save_plot(main_pca + ggplot2::facet_wrap(stats::reformulate(strata), nrow = 1),
    "PCA_global_coordinates_by_time", 16, 6))
  variance <- data.frame(PC = seq_along(qc$variance), ExplainedPercent = qc$variance, CumulativePercent = cumsum(qc$variance))
  variance_plot <- ggplot2::ggplot(variance, ggplot2::aes(PC, ExplainedPercent)) + ggplot2::geom_col(fill = "#3B75AF") +
    ggplot2::geom_line(ggplot2::aes(y = CumulativePercent), color = "#C94C4C") + ggplot2::geom_point(ggplot2::aes(y = CumulativePercent)) +
    ggplot2::labs(title = "PCA explained variance", y = "Explained / cumulative variance (%)") +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  paths <- c(paths, save_plot(variance_plot, "PCA_explained_variance", 9, 6))

  # 分层PCA独立重新计算
  if (length(strata)) {
    strata_table <- unique(as.data.frame(metadata[strata]))
    local_plots <- vector("list", nrow(strata_table))
    for (i in seq_len(nrow(strata_table))) {
      use <- rep(TRUE, nrow(metadata))
      for (variable in strata) use <- use & as.character(metadata[[variable]]) == as.character(strata_table[[variable]][i])
      block <- metadata[use, , drop = FALSE]
      values <- vapply(strata, function(variable) as.character(strata_table[[variable]][i]), character(1))
      label <- paste(paste(strata, values, sep = "="), collapse = ", ")
      name <- paste0("PCA_independent_", safe_id(label))
      matrix <- t(qc$log_cpm[qc$keep_pca, block$Sample, drop = FALSE])
      variable <- if (nrow(matrix) >= 3L) apply(matrix, 2, stats::sd) > 0 else rep(FALSE, ncol(matrix))
      if (sum(variable) < 2L) {
        local_plots[[i]] <- empty_result_plot(paste("Independent PCA", label), "Insufficient samples or variable genes")
        statuses[[name]] <- data.frame(Analysis = name, Status = "SKIPPED", Detail = "至少3个样本和2个可变基因")
      } else {
        fit <- stats::prcomp(matrix[, variable, drop = FALSE], center = TRUE, scale. = isTRUE(qc_settings$pca_scale))
        explained <- 100 * fit$sdev^2 / sum(fit$sdev^2)
        coordinates <- tibble::rownames_to_column(as.data.frame(fit$x), "Sample") |>
          dplyr::left_join(block, by = "Sample")
        local_plots[[i]] <- ggplot2::ggplot(coordinates, ggplot2::aes(PC1, PC2, color = .data[[qc_settings$color_by]])) +
          ggplot2::geom_point(size = 3) + ggplot2::scale_color_manual(values = palette, drop = FALSE) +
          ggplot2::labs(title = paste("Independent PCA,", label), x = sprintf("PC1 (%.1f%%)", explained[1]),
            y = sprintf("PC2 (%.1f%%)", explained[2])) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
        paths <- c(paths, write_table_path(coordinates, file.path(dir, paste0(name, "_coordinates.tsv"))),
          write_table_path(data.frame(PC = seq_along(explained), ExplainedPercent = explained), file.path(dir, paste0(name, "_variance.tsv"))))
        statuses[[name]] <- data.frame(Analysis = name, Status = "SUCCESS", Detail = paste(nrow(block), "samples"))
      }
      paths <- c(paths, save_plot(local_plots[[i]], name))
    }
    pages <- split(seq_along(local_plots), ceiling(seq_along(local_plots) / settings$combine_per_page))
    for (page in seq_along(pages)) {
      indices <- pages[[page]]
      paths <- c(paths, save_panel_files(function(i) local_plots[[indices[i]]],
        file.path(combined_dir, paste0("04_independent_PCA_by_time_", page)), length(indices), min(3, length(indices)),
        21, ceiling(length(indices) / 3) * 6, settings))
    }
    rm(local_plots)
  }

  # 恢复相关性和聚类
  correlation <- qc$correlation[metadata$Sample, metadata$Sample]
  annotation_columns <- unique(c(qc_settings$color_by, qc_settings$shape_by))
  sample_anno <- as.data.frame(metadata[annotation_columns])
  rownames(sample_anno) <- metadata$Sample
  annotation_colors <- stats::setNames(list(palette), qc_settings$color_by)
  cor_plot <- pheatmap::pheatmap(correlation, annotation_col = sample_anno, annotation_row = sample_anno,
    annotation_colors = annotation_colors, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = FALSE,
    color = grDevices::colorRampPalette(c("#F7FBFF", "#08306B"))(100), silent = TRUE)$gtable
  paths <- c(paths, save_plot(cor_plot, "sample_correlation_fixed_order", 14, 13), save_plot(cor_plot, "04_correlation", 14, 13),
    write_table_path(tibble::rownames_to_column(as.data.frame(correlation), "Samples"), file.path(dir, "05_sample_correlation.tsv")))
  cor_clustered <- context$run_project_visualization(pca_design, type = "heatmap", assay_name = "log_cpm", out_dir = NULL,
    display_numbers = FALSE, annotation_colors = list(Group = groups))
  paths <- c(paths, save_plot(cor_clustered, "sample_correlation_clustered", 14, 13))
  dendrogram <- context$run_project_visualization(pca_design, type = "dendro", assay_name = "log_cpm", out_dir = NULL) +
    ggplot2::scale_color_manual(values = groups) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  paths <- c(paths, save_plot(dendrogram, "sample_dendrogram", 12, 12), save_plot(dendrogram, "05_clustering", 12, 12),
    save_panel_files(function(i) list(main_pca, cor_plot, dendrogram)[[i]], file.path(combined_dir, "05_sample_relationships"), 3, 3, 24, 10, settings))

  # 3D使用得分坐标
  if (requireNamespace("scatterplot3d", quietly = TRUE) && all(c("PC1", "PC2", "PC3") %in% names(qc$scores))) {
    data <- qc$scores[match(metadata$Sample, qc$scores$Sample), , drop = FALSE]
    shape_levels <- if (is.factor(metadata[[qc_settings$shape_by]])) levels(droplevels(metadata[[qc_settings$shape_by]])) else {
      unique(as.character(metadata[[qc_settings$shape_by]]))
    }
    shapes <- stats::setNames(rep(c(16, 17, 15, 18, 0:14, 19:25), length.out = length(shape_levels)), shape_levels)
    draw_3d <- function() {
      scatterplot3d::scatterplot3d(data[, c("PC1", "PC2", "PC3")], color = palette[as.character(data[[qc_settings$color_by]])],
        pch = shapes[as.character(data[[qc_settings$shape_by]])], xlab = sprintf("PC1 (%.2f%%)", qc$variance[1]),
        ylab = sprintf("PC2 (%.2f%%)", qc$variance[2]), zlab = sprintf("PC3 (%.2f%%)", qc$variance[3]))
      graphics::legend("topright", legend = names(palette), col = palette, pch = 16, cex = 0.8)
      graphics::legend("bottomleft", legend = names(shapes), pch = shapes, cex = 0.8)
    }
    paths <- c(paths, save_base_plot_files(draw_3d, file.path(dir, "PCA_3D_scores"), settings))
    statuses[["PCA_3D"]] <- data.frame(Analysis = "PCA_3D_scores", Status = "SUCCESS", Detail = "Principal component scores")
  } else statuses[["PCA_3D"]] <- data.frame(Analysis = "PCA_3D_scores", Status = "SKIPPED", Detail = "需要scatterplot3d且至少3个主成分")
  coverage <- qc_coverage_tbl(qc, annotation, backgrounds)
  coverage_plot <- ggplot2::ggplot(coverage, ggplot2::aes(Source, Proportion)) + ggplot2::geom_col(fill = "#3B75AF", na.rm = TRUE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) + ggplot2::labs(y = "Proportion of count-matrix genes", x = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  paths <- c(paths, save_plot(coverage_plot, "annotation_coverage", 9, 5),
    write_table_path(coverage, file.path(dir, "06_annotation_coverage.tsv")),
    write_table_path(dplyr::bind_rows(statuses), file.path(dir, "QC_status.tsv")))
  plots <- paths[tools::file_ext(paths) %in% c("pdf", "png")]
  sizes <- file.info(plots)$size
  if (anyNA(sizes) || any(sizes == 0)) stop("QC图像缺失或为空：", paste(basename(plots[is.na(sizes) | sizes == 0]), collapse = ", "))
  c(paths, write_table_path(data.frame(File = plots, Format = tools::file_ext(plots)), file.path(dir, "QC_figure_manifest.tsv")))
}
