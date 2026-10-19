as_combination_plot <- function(plot) {
  if (inherits(plot, c("ggplot", "patchwork"))) return(plot)
  if (inherits(plot, "pheatmap")) plot <- plot$gtable
  if (inherits(plot, c("grob", "gTree", "gtable"))) return(patchwork::wrap_elements(full = plot))
  if (inherits(plot, "gglist")) {
    heights <- attr(plot, "params")$heights
    if (is.null(heights)) heights <- if (length(plot) == 3) c(1, 0.2, 0.6) else rep(1, length(plot))
    return(patchwork::wrap_plots(plot, ncol = 1, heights = heights))
  }
  stop("无法组合的图对象：", paste(class(plot), collapse = ", "))
}

empty_panel_plot <- function(title, detail = "No significant result") {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = detail, size = 4) +
    ggplot2::labs(title = title) + ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

save_plot_pair <- function(plot, prefix, width = NULL, height = NULL, n_terms = NULL) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  if (is.null(width)) width <- if (is.null(n_terms)) 10 else 13
  if (is.null(height)) height <- if (is.null(n_terms)) 7 else min(32, max(5, n_terms * 0.26 + 2))
  plot <- as_combination_plot(plot)
  ggplot2::ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, units = "in", limitsize = FALSE,
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggplot2::ggsave(paste0(prefix, ".png"), plot, width = width, height = height, units = "in", dpi = cfg$dpi, limitsize = FALSE)
  invisible(prefix)
}

save_plot <- function(plot, file_prefix, width, height) save_plot_pair(plot, file_prefix, width, height)

attempt_plot <- function(module, item, expression) {
  tryCatch({
    p <- force(expression)
    if (is.null(p)) {
      record_status(module, item, "NO_RESULT", "没有可绘制的结果")
    } else record_status(module, item, "SUCCESS")
    p
  }, error = function(e) {
    record_status(module, item, "FAILED", conditionMessage(e))
    warning(item, "：", conditionMessage(e), call. = FALSE)
    NULL
  })
}

combine_plots <- function(plots, prefix, ncol = 3, width = NULL, height = NULL) {
  if (!length(plots)) return(NULL)
  labels <- names(plots)
  panels <- lapply(seq_along(plots), function(i) {
    if (is.null(plots[[i]])) empty_panel_plot(labels[i], "See module_status.tsv") else as_combination_plot(plots[[i]])
  })
  p <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
  if (is.null(width)) width <- min(24, ncol * 7)
  if (is.null(height)) height <- min(36, ceiling(length(panels) / ncol) * 6)
  save_plot_pair(p, prefix, width, height)
  p
}

plot_matrix_heatmap <- function(matrix, title, scale_rows = FALSE, sample_info = NULL) {
  if (!nrow(matrix) || !ncol(matrix)) return(empty_panel_plot(title, "No eligible genes"))
  matrix <- as.matrix(matrix)
  if (scale_rows) {
    variable <- apply(matrix, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
    matrix <- matrix[!is.na(variable) & variable, , drop = FALSE]
    if (!nrow(matrix)) return(empty_panel_plot(title, "No variable genes"))
    matrix <- t(scale(t(matrix)))
  }
  finite <- matrix[is.finite(matrix)]
  limit <- if (length(finite)) max(abs(finite)) else 1
  if (limit == 0) limit <- 1
  annotation <- NULL
  colors <- NULL
  if (!is.null(sample_info)) {
    annotation <- as.data.frame(sample_info[match(colnames(matrix), sample_info$Samples), c("Treatment", "Time")])
    rownames(annotation) <- colnames(matrix)
    colors <- list(Treatment = cfg$treatment_colors)
  }
  complete_rows <- apply(matrix, 1, function(x) all(is.finite(x)))
  pheatmap::pheatmap(matrix, color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
                     breaks = seq(-limit, limit, length.out = 101), cluster_rows = nrow(matrix) > 1 && all(complete_rows),
                     cluster_cols = FALSE, annotation_col = annotation, annotation_colors = colors,
                     show_rownames = nrow(matrix) <= 70, border_color = NA, na_col = "#BDBDBD", main = title,
                     fontsize_col = 8, silent = TRUE)$gtable
}

build_effect_heatmap_mat <- function(deg_results, public_data = NULL, contrast_numbers = NULL, go_ids = NULL) {
  comparisons <- names(deg_results)
  if (!is.null(contrast_numbers)) comparisons <- comparisons[contrast_numbers]
  significant <- unique(unlist(lapply(deg_results[comparisons], function(x) x$significant$Geneid), use.names = FALSE))
  if (!is.null(go_ids) && !is.null(public_data)) {
    bg <- public_data$get_go_background("BP")$TERM2GENE
    ids <- unique(bg$GENE[bg$TERM %in% go_ids])
    anno <- public_data$get_gene_annotation()
    significant <- intersect(significant, anno$Geneid[anno$ENTREZID %in% ids])
  }
  matrix <- matrix(NA_real_, nrow = length(significant), ncol = length(comparisons), dimnames = list(significant, comparisons))
  for (comparison in comparisons) {
    tbl <- deg_results[[comparison]]$all
    matrix[, comparison] <- tbl$logFC[match(significant, tbl$Geneid)]
  }
  matrix
}

# 三集合Venn保留精确计数
venn_plot <- function(data, ...) {
  sets <- data[lengths(data) > 0]
  n <- length(sets)
  if (!n) return(empty_panel_plot("Venn", "No significant genes"))
  if (n == 1) return(empty_panel_plot(names(sets)[1], paste("Genes:", length(unique(sets[[1]])))))
  if (n > 3) stop("Venn仅用于2或3集合，请使用UpSet")
  centers <- if (n == 3) data.frame(x = c(0, -0.43, 0.43), y = c(0.45, -0.3, -0.3)) else data.frame(x = c(-0.45, 0.45), y = c(0, 0))
  circle <- dplyr::bind_rows(lapply(seq_len(n), function(i) {
    theta <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = centers$x[i] + 0.72 * cos(theta), y = centers$y[i] + 0.72 * sin(theta), Set = names(sets)[i])
  }))
  all <- unique(unlist(sets, use.names = FALSE))
  membership <- do.call(cbind, lapply(sets, function(s) all %in% s))
  code <- apply(membership, 1, function(x) paste(as.integer(x), collapse = ""))
  positions <- if (n == 3) {
    data.frame(code = c("100", "010", "001", "110", "101", "011", "111"), x = c(0, -0.85, 0.85, -0.34, 0.34, 0, 0),
               y = c(0.95, -0.45, -0.45, 0.2, 0.2, -0.57, -0.04))
  } else data.frame(code = c("10", "01", "11"), x = c(-0.75, 0.75, 0), y = 0)
  positions$Count <- vapply(positions$code, function(x) sum(code == x), integer(1))
  ggplot2::ggplot(circle, ggplot2::aes(x, y, group = Set)) + ggplot2::geom_polygon(ggplot2::aes(fill = Set), alpha = 0.25) +
    ggplot2::geom_path(ggplot2::aes(color = Set), linewidth = 0.6) +
    ggplot2::geom_text(data = positions, ggplot2::aes(x, y, label = Count), inherit.aes = FALSE, size = 4) +
    ggplot2::coord_equal() + ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom")
}

gene_set_intersections_tbl <- function(sets) {
  ids <- unique(unlist(sets, use.names = FALSE))
  if (!length(ids)) return(data.frame())
  membership <- do.call(cbind, lapply(sets, function(s) ids %in% s))
  code <- apply(membership, 1, function(x) paste(as.integer(x), collapse = ""))
  tbl <- as.data.frame(table(code), stringsAsFactors = FALSE)
  names(tbl) <- c("Code", "Count")
  tbl[order(-tbl$Count, tbl$Code), , drop = FALSE]
}

# UpSet辅助接口
upsetr_plot <- function(data, ...) {
  sets <- lapply(data, unique)
  intersections <- gene_set_intersections_tbl(sets)
  if (!nrow(intersections)) return(empty_panel_plot("UpSet", "No significant genes"))
  intersections <- utils::head(intersections, 25)
  intersections$Intersection <- factor(seq_len(nrow(intersections)), levels = seq_len(nrow(intersections)))
  dots <- dplyr::bind_rows(lapply(seq_len(nrow(intersections)), function(i) {
    data.frame(Intersection = intersections$Intersection[i], Set = factor(names(sets), levels = rev(names(sets))),
               Member = strsplit(intersections$Code[i], "", fixed = TRUE)[[1]] == "1")
  }))
  bars <- ggplot2::ggplot(intersections, ggplot2::aes(Intersection, Count)) + ggplot2::geom_col(fill = "#315C80") +
    ggplot2::theme_classic() + ggplot2::labs(x = NULL, y = "Intersection size") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  matrix <- ggplot2::ggplot(dots, ggplot2::aes(Intersection, Set)) + ggplot2::geom_point(color = "grey85", size = 2) +
    ggplot2::geom_point(data = dots[dots$Member, ], color = "#222222", size = 2) +
    ggplot2::theme_classic() + ggplot2::labs(x = "Top intersections", y = NULL)
  patchwork::wrap_plots(bars, matrix, ncol = 1, heights = c(2, max(1, length(sets) / 5)))
}

plot_gene_sets_flower_plot <- function(sets, title) {
  n <- length(sets)
  if (!n) return(empty_panel_plot(title))
  angles <- seq(0, 2 * pi, length.out = n + 1)[seq_len(n)]
  petals <- dplyr::bind_rows(lapply(seq_len(n), function(i) {
    theta <- seq(0, 2 * pi, length.out = 100)
    xx <- 0.75 + 0.65 * cos(theta)
    yy <- 0.25 * sin(theta)
    data.frame(x = xx * cos(angles[i]) - yy * sin(angles[i]), y = xx * sin(angles[i]) + yy * cos(angles[i]), Set = names(sets)[i])
  }))
  unique_counts <- vapply(seq_len(n), function(i) length(setdiff(sets[[i]], unlist(sets[-i], use.names = FALSE))), integer(1))
  common <- if (n) length(Reduce(intersect, sets)) else 0L
  labels <- data.frame(x = 1.2 * cos(angles), y = 1.2 * sin(angles), Label = paste0(names(sets), "\nUnique=", unique_counts))
  ggplot2::ggplot(petals, ggplot2::aes(x, y, group = Set, fill = Set)) + ggplot2::geom_polygon(alpha = 0.7, color = "white") +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = Label), inherit.aes = FALSE, size = 2.4) +
    ggplot2::annotate("label", x = 0, y = 0, label = paste0("Common\n", common), size = 3) +
    ggplot2::coord_equal(xlim = c(-1.65, 1.65), ylim = c(-1.65, 1.65), clip = "off") + ggplot2::theme_void() +
    ggplot2::labs(title = title, caption = "Petals: genes unique to each set; center: intersection of all sets") + ggplot2::guides(fill = "none")
}

collect_comparative_tbl <- function(enrichment, db, method, contrast_range = NULL) {
  comparisons <- names(enrichment)
  if (!is.null(contrast_range)) comparisons <- comparisons[contrast_range]
  rows <- lapply(comparisons, function(comparison) {
    item <- enrichment[[comparison]][[db]]
    if (is.null(item)) return(NULL)
    if (method == "ORA") {
      dplyr::bind_rows(lapply(c("Up", "Down"), function(direction) {
        tbl <- result_table(item$ORA[[direction]], TRUE)
        if (nrow(tbl)) { tbl$Direction <- direction; tbl$Comparison <- comparison }
        tbl
      }))
    } else {
      tbl <- result_table(item$GSEA, TRUE)
      if (nrow(tbl)) { tbl$Direction <- ifelse(tbl$NES > 0, "Up", "Down"); tbl$Comparison <- comparison }
      tbl
    }
  })
  dplyr::bind_rows(rows)
}

ratio_numeric <- function(x) {
  vapply(as.character(x), function(value) {
    parts <- as.numeric(strsplit(value, "/", fixed = TRUE)[[1]])
    if (length(parts) == 2 && parts[2] > 0) parts[1] / parts[2] else NA_real_
  }, numeric(1))
}

plot_comparative_enrichment_plot <- function(tbl, db, method, selected_ids = NULL) {
  title <- paste(if (db == "GO") "GO-BP" else db, method, "comparative enrichment")
  if (!nrow(tbl)) return(empty_panel_plot(title, "No significant enrichment; see status tables"))
  if (is.null(selected_ids)) selected_ids <- select_representative_ids(tbl)
  tbl <- tbl[tbl$ID %in% selected_ids, , drop = FALSE]
  descriptions <- tbl$Description[match(selected_ids, tbl$ID)]
  supplied_labels <- attr(selected_ids, "labels")
  if (!is.null(supplied_labels)) descriptions <- unname(supplied_labels[selected_ids])
  descriptions[is.na(descriptions)] <- selected_ids[is.na(descriptions)]
  labels <- stats::setNames(paste0(descriptions, " [", selected_ids, "]"), selected_ids)
  tbl$ID <- factor(tbl$ID, levels = rev(selected_ids))
  tbl$Comparison <- factor(tbl$Comparison, levels = cfg$contrast_levels)
  tbl$Direction <- factor(tbl$Direction, levels = c("Up", "Down"))
  tbl$Significance <- -log10(pmax(tbl$p.adjust, .Machine$double.xmin))
  tbl$Magnitude <- if (method == "ORA") ratio_numeric(tbl$GeneRatio) else abs(tbl$NES)
  position <- if (method == "ORA") ggplot2::position_dodge(width = 0.4) else "identity"
  ggplot2::ggplot(tbl, ggplot2::aes(Comparison, ID)) +
    ggplot2::geom_vline(xintercept = c(3.5, 6.5), color = "grey85", linetype = 2, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(shape = Direction, size = Magnitude, fill = Significance, group = Direction),
                       position = position, color = "#444444", stroke = 0.25) +
    ggplot2::scale_shape_manual(values = c(Up = 24, Down = 25), drop = FALSE) +
    ggplot2::scale_fill_gradient(low = "#FEE8C8", high = "#B30000", name = "−log10(FDR)") +
    ggplot2::scale_size_continuous(name = if (method == "ORA") "GeneRatio" else "|NES|", range = c(2, 7)) +
    ggplot2::scale_y_discrete(labels = function(ids) stringr::str_wrap(labels[ids], width = 55), drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE, labels = function(x) sub("_vs_", "\nvs ", x, fixed = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL, title = title) + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
}

label_cnet_plot <- function(p, public_data) {
  annotation <- public_data$get_gene_annotation()
  replace_labels <- function(data) {
    if (!is.data.frame(data)) return(data)
    for (field in intersect(c("name", "label"), names(data))) {
      values <- as.character(data[[field]])
      idx <- match(values, annotation$ENTREZID)
      values[!is.na(idx)] <- annotation$SYMBOL[idx[!is.na(idx)]]
      data[[field]] <- values
    }
    data
  }
  p$data <- replace_labels(p$data)
  for (i in seq_along(p$layers)) p$layers[[i]]$data <- replace_labels(p$layers[[i]]$data)
  p
}

plot_ora_cnet <- function(enrich_result, public_data, showCategory = 5, logFC = NULL, title = NULL, ...) {
  res <- subset_enrichment(enrich_result)
  if (is.null(res)) return(NULL)
  set.seed(cfg$seed)
  p <- enrichplot::cnetplot(res, showCategory = showCategory, foldChange = logFC, colorEdge = TRUE)
  p <- label_cnet_plot(p, public_data)
  p + ggplot2::labs(title = title)
}

plot_pathview_result <- function(res, ranks, out_dir, limit) {
  if (!cfg$run_pathview || !requireNamespace("pathview", quietly = TRUE)) return(NULL)
  tbl <- result_table(res, TRUE)
  if (!nrow(tbl)) return(NULL)
  tbl <- tbl[order(tbl$p.adjust), , drop = FALSE]
  ids <- utils::head(tbl$ID, cfg$pathview_top)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (id in ids) {
    pathway <- sub("^[[:alpha:]]+", "", id)
    pathview::pathview(gene.data = ranks, pathway.id = pathway, species = cfg$kegg_code, gene.idtype = "entrez",
                       kegg.dir = out_dir, limit = list(gene = limit), same.layer = FALSE,
                       low = list(gene = "#2166AC"), mid = list(gene = "#F7F7F7"), high = list(gene = "#B2182B"))
  }
  files <- list.files(out_dir, pattern = "pathview.*\\.png$", full.names = TRUE)
  if (length(files)) files else NULL
}
