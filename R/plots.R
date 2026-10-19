plot_pca_plot <- function(qc, metadata, display, style) {
  color_by <- display$color_by
  shape_by <- display$shape_by
  if (!all(c(color_by, shape_by) %in% names(metadata))) stop("PCA映射变量不存在")
  use_shape <- length(shape_by) == 1L && nzchar(shape_by)
  if (use_shape && is.numeric(metadata[[shape_by]])) stop("形状变量必须为分类变量，请使用Time_f")
  df <- dplyr::left_join(qc$scores[, c("Sample", "PC1", "PC2")], metadata, by = "Sample")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::labs(title = "PCA", x = sprintf("PC1 (%.1f%%)", qc$variance[1]), y = sprintf("PC2 (%.1f%%)", qc$variance[2])) + theme_pca(style)
  palette <- variable_palette(style, color_by, df[[color_by]])
  if (!is.null(palette)) p <- p + ggplot2::scale_color_manual(values = palette, drop = FALSE)
  if (use_shape) {
    shape_levels <- if (is.factor(df[[shape_by]])) levels(df[[shape_by]]) else unique(as.character(df[[shape_by]]))
    n_shapes <- length(shape_levels)
    if (n_shapes > 25L) stop("形状水平超过25，请调整PCA映射")
    p <- p + ggplot2::aes(shape = .data[[shape_by]]) + ggplot2::scale_shape_manual(values = stats::setNames(seq_len(n_shapes) - 1L, shape_levels), drop = FALSE)
  }
  facets <- as.character(unlist(display$facet_by))
  if (length(facets)) {
    if (!all(facets %in% names(df))) stop("PCA分面变量不存在")
    p <- p + ggplot2::facet_wrap(stats::reformulate(facets))
  }
  p
}

plot_deg_plot <- function(result, thresholds, style, type) {
  tbl <- result$table
  tbl$Direction <- ifelse(is.finite(tbl$padj) & tbl$padj < thresholds$fdr & abs(tbl$logFC) >= thresholds$abs_logfc,
                          ifelse(tbl$logFC > 0, "Up", "Down"), "Not")
  tbl$Significance <- -log10(pmax(tbl$padj, .Machine$double.xmin))
  if (type == "volcano") {
    p <- ggplot2::ggplot(tbl, ggplot2::aes(logFC, Significance, color = Direction)) +
      ggplot2::geom_point(size = 1.2, alpha = 0.7) + ggplot2::geom_vline(xintercept = c(-1, 1) * thresholds$abs_logfc, linetype = 4) +
      ggplot2::geom_hline(yintercept = -log10(thresholds$fdr), linetype = 4) + ggplot2::labs(x = "log2 fold change", y = "−log10(FDR)")
  } else if (type == "ma") {
    p <- ggplot2::ggplot(tbl, ggplot2::aes(logCPM, logFC, color = Direction)) + ggplot2::geom_point(size = 1.2, alpha = 0.7) +
      ggplot2::labs(x = "Average logCPM", y = "log2 fold change")
  } else {
    tbl <- tbl[order(-tbl$logFC), , drop = FALSE]
    tbl$Rank <- seq_len(nrow(tbl))
    p <- ggplot2::ggplot(tbl, ggplot2::aes(Rank, logFC, color = Direction)) + ggplot2::geom_point(size = 1.2, alpha = 0.7)
  }
  p + ggplot2::scale_color_manual(values = unlist(style$direction_colors), drop = FALSE) +
    ggplot2::labs(title = result$job$comparison) + theme_project(style)
}

ratio_vec <- function(values) {
  vapply(as.character(values), function(x) {
    parts <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
    if (length(parts) == 2L && parts[2] > 0) parts[1] / parts[2] else NA_real_
  }, numeric(1))
}

plot_comparative_plot <- function(catalog, selection, jobs, method, thresholds, display, style, layout = "triangles") {
  title <- paste(if (selection$database == "GO") "GO-BP" else "KEGG", method, selection$scheme)
  tbl <- catalog[catalog$Scheme == selection$scheme & catalog$Database == selection$database & catalog$Method == method &
                   catalog$ID %in% selection$ids & is.finite(catalog$p.adjust) & catalog$p.adjust < thresholds$fdr, , drop = FALSE]
  if (!nrow(tbl)) {
    status <- selection$status
    detail <- if (!is.null(status) && status != "SUCCESS") paste("Database", status) else "No significant result"
    return(empty_result_plot(title, detail))
  }
  manifest <- comparison_manifest_tbl(jobs)
  manifest <- manifest[manifest$Scheme == selection$scheme, , drop = FALSE]
  variables <- intersect(as.character(unlist(display$compare_by)), names(manifest))
  if (length(variables)) manifest <- sort_metadata_tbl(manifest, variables)
  ids <- selection$ids
  tbl$ID <- factor(tbl$ID, levels = rev(ids))
  tbl$Direction <- factor(tbl$Direction, levels = c("Up", "Down"))
  tbl$Column <- match(tbl$Job, manifest$Job)
  tbl$X <- tbl$Column + if (method == "ORA" && layout == "triangles") ifelse(tbl$Direction == "Up", -0.22, 0.22) else 0
  tbl$Significance <- -log10(pmax(tbl$p.adjust, .Machine$double.xmin))
  tbl$Magnitude <- if (method == "ORA") ratio_vec(tbl$GeneRatio) else abs(tbl$NES)
  labels <- paste0(unname(selection$labels[ids]), " [", ids, "]")
  names(labels) <- ids
  p <- ggplot2::ggplot(tbl, ggplot2::aes(X, ID))
  if (length(variables)) {
    blocks <- as.character(manifest[[variables[1]]])
    cuts <- which(blocks[-1] != blocks[-length(blocks)]) + 0.5
    if (length(cuts)) p <- p + ggplot2::geom_vline(xintercept = cuts, linetype = 2, color = "grey75", linewidth = 0.3)
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(shape = Direction, size = Magnitude, fill = Significance), color = "white", stroke = 0.2) +
    ggplot2::scale_shape_manual(values = c(Up = 24, Down = 25), drop = FALSE) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey30", colour = "grey30", size = 4)),
                    size = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "grey40", colour = "grey40"))) +
    ggplot2::scale_fill_gradientn(colours = unlist(style$significance_colors), name = "−log10(FDR)") +
    ggplot2::scale_size_continuous(range = c(2, 7), name = if (method == "ORA") "GeneRatio" else "|NES|") +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(manifest)), labels = sub("_vs_", "\nvs ", manifest$Comparison, fixed = TRUE), limits = c(0.5,
      nrow(manifest) + 0.5)) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(labels[x], 55), drop = FALSE) +
    ggplot2::labs(title = title, caption = if (method == "ORA") "▲ Up-DEG set; ▼ Down-DEG set" else "▲ Positive NES; ▼ Negative NES") +
      theme_enrichment(style) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = if (nrow(manifest) > 6) 8 else style$base_size))
  if (layout == "facets") p <- p + ggplot2::facet_grid(cols = ggplot2::vars(Direction), drop = FALSE)
  p
}

plot_enrichment_pair_plot <- function(up, down, thresholds, display, style, type = "dot") {
  tbl <- dplyr::bind_rows(collect_enrichment_tbl(up), collect_enrichment_tbl(down))
  title <- paste(up$database, "ORA", up$job$comparison)
  if (!nrow(tbl)) return(empty_result_plot(title, up$status))
  tbl <- tbl[is.finite(tbl$p.adjust) & tbl$p.adjust < thresholds$fdr, , drop = FALSE]
  if (!nrow(tbl)) return(empty_result_plot(title))
  ids <- unique(unlist(lapply(split(tbl, tbl$Direction), function(x) utils::head(x$ID[order(x$p.adjust)], display$top_terms))))
  tbl <- tbl[tbl$ID %in% ids, , drop = FALSE]
  tbl$Term <- factor(tbl$ID, levels = rev(ids))
  labels <- stats::setNames(tbl$Description, tbl$ID)
  tbl$Significance <- -log10(pmax(tbl$p.adjust, .Machine$double.xmin))
  tbl$SignedSignificance <- ifelse(tbl$Direction == "Up", 1, -1) * tbl$Significance
  tbl$GeneRatioValue <- ratio_vec(tbl$GeneRatio)
  if (type == "bar") {
    p <- ggplot2::ggplot(tbl, ggplot2::aes(SignedSignificance, Term, fill = Direction)) + ggplot2::geom_col(position = "identity") +
      ggplot2::scale_fill_manual(values = unlist(style$direction_colors)) +
      ggplot2::labs(x = "Direction × −log10(FDR)", caption = "Sign indicates the DEG set; length indicates enrichment significance")
  } else {
    p <- ggplot2::ggplot(tbl, ggplot2::aes(GeneRatioValue, Term, fill = Significance, shape = Direction)) +
      ggplot2::geom_point(ggplot2::aes(size = Count), color = "white") + ggplot2::facet_grid(cols = ggplot2::vars(Direction), drop = FALSE) +
      ggplot2::scale_shape_manual(values = c(Up = 24, Down = 25), drop = FALSE) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = "grey30", colour = "grey30", size = 4)),
                      size = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "grey40", colour = "grey40"))) +
      ggplot2::scale_fill_gradientn(colours = unlist(style$significance_colors), name = "−log10(FDR)") + ggplot2::labs(x = "GeneRatio")
  }
  p + ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(labels[x], 50), drop = FALSE) + ggplot2::labs(title = title, y = NULL) +
    theme_project(style)
}
