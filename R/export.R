save_panel_files <- function(draw_panel, prefix, n_panels = 1L, ncol = 1L, width = 10, height = 7, settings) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  dpi <- min(settings$dpi, floor(sqrt(settings$max_pixels / (width * height))))
  if (!is.finite(dpi) || dpi < 1 || n_panels < 1) stop("导出尺寸或像素预算无效")
  paths <- paste0(prefix, c(".pdf", ".png"))
  for (i in seq_along(paths)) {
    if (i == 1L) {
      device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
      device(paths[i], width = width, height = height)
    } else {
      args <- list(filename = paths[i], width = width, height = height, units = "in", res = dpi)
      if (capabilities("cairo")) args$type <- "cairo"
      do.call(grDevices::png, args)
    }
    device_id <- grDevices::dev.cur()
    tryCatch({
      grDevices::dev.control(displaylist = "inhibit")
      grid::grid.newpage()
      grid::grid.display.list(FALSE)
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(ceiling(n_panels / ncol), ncol)))
      for (j in seq_len(n_panels)) {
        grid::pushViewport(grid::viewport(layout.pos.row = (j - 1L) %/% ncol + 1L, layout.pos.col = (j - 1L) %% ncol + 1L))
        tryCatch({
          panel <- draw_panel(j)
          if (!is.null(panel)) grid::grid.draw(as_panel_grob(panel))
          rm(panel)
        }, finally = grid::popViewport())
        invisible(gc())
      }
      grid::popViewport()
    }, finally = grDevices::dev.off(device_id))
  }
  paths
}

save_single_plot_files <- function(plot, prefix, settings, width = 10, height = 7) {
  save_panel_files(function(i) plot, prefix, width = width, height = height, settings = settings)
}

export_metadata_files <- function(metadata, jobs, output_dir) {
  c(write_table_path(metadata, file.path(output_dir, "01_inputs", "sample_metadata.tsv")),
    write_table_path(comparison_manifest_tbl(jobs), file.path(output_dir, "01_inputs", "comparison_manifest.tsv")))
}

export_annotation_files <- function(go, output_dir) {
  dir <- file.path(output_dir, "02_annotation", "GO_BP")
  c(write_table_path(go$annotation, file.path(dir, "gene_annotation.tsv")), write_table_path(go$direct_go, file.path(dir, "gene2go_direct.tsv")),
    write_table_path(go$background$TERM2GENE, file.path(dir, "TERM2GENE.tsv")), write_table_path(go$background$TERM2NAME, file.path(dir, "TERM2NAME.tsv")))
}

export_deg_files <- function(result, selected, output_dir) {
  dir <- file.path(output_dir, "04_differential_expression", safe_id(result$job$id))
  c(write_object_path(result, file.path(dir, "all_tested.rds")), write_table_path(result$table, file.path(dir, "all_tested.tsv")),
    write_table_path(selected$table, file.path(dir, "significant.tsv")),
    write_table_path(selected$table[selected$table$Result_expression == "Up", ], file.path(dir, "Up.tsv")),
    write_table_path(selected$table[selected$table$Result_expression == "Down", ], file.path(dir, "Down.tsv")))
}

export_enrichment_files <- function(result, output_dir) {
  dir <- file.path(output_dir, "07_enrichment", result$database, safe_id(result$job$id), paste(result$method, result$direction, sep = "_"))
  status <- data.frame(Job = result$job$id, Scheme = result$job$scheme, Comparison = result$job$comparison, Database = result$database,
                       Method = result$method, Direction = result$direction, Status = result$status)
  paths <- c(write_object_path(result, file.path(dir, "result.rds")), write_table_path(result$table, file.path(dir, "all_terms.tsv")),
             write_table_path(status, file.path(dir, "status.tsv")))
  if (result$method == "ORA") paths <- c(paths, write_table_path(data.frame(ENTREZID = result$universe), file.path(dir, "tested_background.tsv")),
    write_table_path(data.frame(ENTREZID = result$foreground), file.path(dir, "foreground.tsv"))) else {
    paths <- c(paths, write_table_path(data.frame(ENTREZID = names(result$ranking), RankMetric = unname(result$ranking)), file.path(dir, "ranking.tsv")))
  }
  paths
}

read_enrichment_record_lst <- function(files) {
  row <- data.table::fread(files[basename(files) == "status.tsv"], data.table = FALSE)
  as.list(c(as.list(row[1, ]), list(object_file = files[basename(files) == "result.rds"])))
}

export_qc_files <- function(qc, metadata, display, style, settings, output_dir) {
  dir <- file.path(output_dir, "03_sample_qc")
  ordered <- sort_metadata_tbl(metadata, display$order_by)
  paths <- write_table_path(dplyr::left_join(qc$library, metadata, by = "Sample"), file.path(dir, "library_statistics.tsv"))
  paths <- c(paths, write_table_path(qc$scores, file.path(dir, "PCA_scores.tsv")))
  for (name in c("cpm", "tpm", "log_cpm")) paths <- c(paths, write_table_path(tibble::rownames_to_column(as.data.frame(qc[[name]]), "Geneid"),
                                                                                              file.path(dir, paste0(name, ".tsv"))))
  paths <- c(paths, save_single_plot_files(plot_pca_plot(qc, metadata, display, style), file.path(dir, "01_PCA"), settings, 8, 7))
  library <- dplyr::left_join(qc$library, metadata, by = "Sample")
  library$Sample <- factor(library$Sample, levels = ordered$Sample)
  p <- ggplot2::ggplot(library, ggplot2::aes(Sample, LibrarySize, fill = .data[[display$color_by]])) + ggplot2::geom_col() +
    ggplot2::labs(title = "Library size", x = NULL) + theme_enrichment(style)
  palette <- variable_palette(style, display$color_by, metadata[[display$color_by]])
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_manual(values = palette, drop = FALSE)
  paths <- c(paths, save_single_plot_files(p, file.path(dir, "02_library_size"), settings, max(9, nrow(metadata) * 0.28), 6))
  long <- as.data.frame(as.table(qc$log_cpm), stringsAsFactors = FALSE)
  names(long) <- c("Geneid", "Sample", "Expression")
  long <- dplyr::left_join(long, metadata, by = "Sample")
  long$Sample <- factor(long$Sample, levels = ordered$Sample)
  p <- ggplot2::ggplot(long, ggplot2::aes(Sample, Expression, fill = .data[[display$color_by]])) +
    ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::labs(title = "logCPM distribution", y = "logCPM", x = NULL) + theme_enrichment(style)
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_manual(values = palette, drop = FALSE)
  paths <- c(paths, save_single_plot_files(p, file.path(dir, "03_expression_boxplot"), settings, max(9, nrow(metadata) * 0.28), 6))
  correlation <- qc$correlation[ordered$Sample, ordered$Sample]
  paths <- c(paths, save_single_plot_files(pheatmap::pheatmap(correlation, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
    main = "Sample correlation", silent = TRUE)$gtable, file.path(dir, "04_correlation"), settings, 12, 11))
  tree <- stats::hclust(stats::dist(t(qc$log_cpm)))
  dendro <- ggdendro::dendro_data(tree)
  p <- ggplot2::ggplot(ggdendro::segment(dendro)) + ggplot2::geom_segment(ggplot2::aes(x, y, xend = xend, yend = yend)) +
    ggplot2::scale_x_continuous(breaks = dendro$labels$x, labels = dendro$labels$label) +
    ggplot2::labs(title = "Sample clustering", x = NULL, y = "Distance") + theme_enrichment(style)
  paths <- c(paths, save_single_plot_files(p, file.path(dir, "05_clustering"), settings, max(10, nrow(metadata) * 0.28), 7))
  paths
}

export_deg_plot_files <- function(result, thresholds, style, settings, output_dir) {
  dir <- file.path(output_dir, "04_differential_expression", safe_id(result$job$id))
  unlist(lapply(c("volcano", "rank", "ma"), function(type) save_single_plot_files(plot_deg_plot(result, thresholds, style, type),
    file.path(dir, paste0("plot_", type)), settings, 8, 7)), use.names = FALSE)
}

export_deg_combined_files <- function(results, thresholds, style, settings, output_dir, display = NULL) {
  if (!is.null(display)) results <- order_comparison_results_lst(results, display$compare_by)
  paths <- character()
  schemes <- unique(vapply(results, function(x) x$job$scheme, character(1)))
  for (scheme in schemes) {
    items <- Filter(function(x) x$job$scheme == scheme, results)
    pages <- split(seq_along(items), ceiling(seq_along(items) / settings$combine_per_page))
    for (type in c("volcano", "rank", "ma")) {
      for (page in seq_along(pages)) {
        indices <- pages[[page]]
        paths <- c(paths, save_panel_files(function(i) plot_deg_plot(items[[indices[i]]], thresholds, style, type),
          file.path(output_dir, "09_combined_figures", paste0("DEG_", safe_id(scheme), "_", type, "_", page)), length(indices),
          settings$combine_columns, 24, ceiling(length(indices) / settings$combine_columns) * 7, settings))
      }
    }
  }
  paths
}

export_ora_plot_files <- function(up, down, annotation, thresholds, display, style, settings, engine_files, plot_helper, output_dir) {
  prefix <- file.path(output_dir, "07_enrichment", up$database, safe_id(up$job$id), "ORA_Up_Down")
  paths <- unlist(lapply(c("dot", "bar"), function(type) save_single_plot_files(plot_enrichment_pair_plot(up, down, thresholds, display, style, type),
    paste0(prefix, "_", type), settings, 15, min(22, max(7, display$top_terms * 0.55 + 3)))), use.names = FALSE)
  context <- new_legacy_context(engine_files, list(fdr_cutoff = thresholds$fdr), plot_helper)
  for (type in c("network", "upset")) {
    if (!isTRUE(display[[paste0("enrichment_", if (type == "network") "networks" else "upset")]])) next
    paths <- c(paths, save_panel_files(function(i) {
      item <- list(up, down)[[i]]
      res <- significant_result(item$result, thresholds$fdr)
      if (is.null(res)) return(empty_result_plot(paste(item$job$comparison, item$direction), item$status))
      p <- if (type == "network") context$plot_ora_cnet(res, legacy_public_data(annotation), showCategory = 5, title = item$direction) else {
        context$plot_ora_upset(res, legacy_public_data(annotation), n = min(8, nrow(res@result)), title = item$direction, out_dir = NULL)
      }
      apply_legacy_theme(p, style)
    }, paste0(prefix, "_", type), n_panels = 2L, ncol = 2L, width = 18, height = 8, settings = settings))
  }
  paths
}

export_gsea_plot_files <- function(result, annotation, thresholds, display, style, settings, engine_files, plot_helper, output_dir) {
  dir <- file.path(output_dir, "07_enrichment", result$database, safe_id(result$job$id), "GSEA_plots")
  context <- new_legacy_context(engine_files, list(fdr_cutoff = thresholds$fdr, top_terms = display$top_terms,
    direction_colors = unlist(style$direction_colors)), plot_helper)
  res <- significant_result(result$result, thresholds$fdr)
  p <- if (is.null(res)) empty_result_plot(result$job$comparison, result$status) else context$plot_gsea_double_bar(res, contrast_name = result$job$comparison)
  paths <- save_single_plot_files(apply_legacy_theme(p, style), file.path(dir, "01_Up_Down_bar"), settings, 13, 8)
  if (isTRUE(display$enrichment_networks)) {
    p <- if (is.null(res)) empty_result_plot(result$job$comparison, result$status) else {
      context$plot_gsea_cnet(res, legacy_public_data(annotation), contrast_name = result$job$comparison, showCategory = 5, out_dir = NULL)
    }
    paths <- c(paths, save_single_plot_files(apply_legacy_theme(p, style), file.path(dir, "02_network"), settings, 12, 9))
  }
  if (!is.null(res)) {
    ids <- utils::head(res@result$ID[order(res@result$p.adjust)], display$curve_top_each)
    draw_curve <- function(i) {
      curve <- context$plot_gsea_gseaplot2(significant_result(res, thresholds$fdr, ids[i]), pvalue_cutoff = thresholds$fdr, out_dir = NULL)
      if (!length(curve)) stop("GSEA曲线未生成：", ids[i])
      apply_legacy_theme(curve[[1]], style)
    }
    for (i in seq_along(ids)) paths <- c(paths, save_panel_files(function(j) draw_curve(i), file.path(dir, paste0("curve_", safe_id(ids[i]))),
      width = 9, height = 7, settings = settings))
    if (length(ids)) paths <- c(paths, save_panel_files(draw_curve, file.path(dir, "03_representative_curves"), length(ids), min(3, length(ids)),
      width = 21, height = ceiling(length(ids) / 3) * 7, settings = settings))
  }
  paths
}

export_comparative_files <- function(catalog, selection, jobs, thresholds, display, style, settings, output_dir) {
  paths <- character()
  dir <- file.path(output_dir, "08_comparative_enrichment", safe_id(selection$scheme), selection$database)
  paths <- c(paths, write_table_path(data.frame(ID = selection$ids, Description = unname(selection$labels[selection$ids])), file.path(dir,
    "displayed_terms.tsv")))
  for (method in c("ORA", "GSEA")) {
    for (layout in c("triangles", "facets")) {
      p <- plot_comparative_plot(catalog, selection, jobs, method, thresholds, display, style, layout)
      width <- if (layout == "facets") 24 else 16
      paths <- c(paths, save_single_plot_files(p, file.path(dir, paste(method, layout, sep = "_")), settings, width,
                                              min(26, max(6, length(selection$ids) * 0.4 + 3))))
    }
  }
  block <- catalog[catalog$Scheme == selection$scheme & catalog$Database == selection$database & catalog$Method == "GSEA", , drop = FALSE]
  if (length(selection$ids) && nrow(block)) {
    manifest <- comparison_manifest_tbl(jobs)
    comparisons <- manifest$Job[manifest$Scheme == selection$scheme]
    grid <- expand.grid(ID = selection$ids, Job = comparisons, stringsAsFactors = FALSE)
    grid <- dplyr::left_join(grid, block[, c("ID", "Job", "NES", "p.adjust")], by = c("ID", "Job"))
    grid$ID <- factor(grid$ID, levels = rev(selection$ids))
    grid$Job <- factor(grid$Job, levels = comparisons)
    grid$Significant <- ifelse(!is.na(grid$p.adjust) & grid$p.adjust < thresholds$fdr, "*", "")
    p <- ggplot2::ggplot(grid, ggplot2::aes(Job, ID, fill = NES)) + ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = Significant), size = 4) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey75") +
      ggplot2::scale_y_discrete(labels = selection$labels) + ggplot2::labs(title = "GSEA NES; * significant FDR", caption = "Grey: NES unavailable") +
      theme_enrichment(style)
    paths <- c(paths, save_single_plot_files(p, file.path(dir, "GSEA_NES_heatmap"), settings, 16, min(24, max(6, length(selection$ids) * 0.4 + 3))))
    paths <- c(paths, write_table_path(grid, file.path(dir, "GSEA_NES_matrix_long.tsv")))
  }
  paths
}

export_comparative_combined_files <- function(catalog, selections, jobs, thresholds, display, style, settings, output_dir) {
  paths <- character()
  for (scheme in unique(vapply(jobs, `[[`, character(1), "scheme"))) {
    selected <- Filter(function(x) x$scheme == scheme, selections)
    panels <- expand.grid(Database = c("GO", "KEGG"), Method = c("ORA", "GSEA"), stringsAsFactors = FALSE)
    paths <- c(paths, save_panel_files(function(i) {
      selection <- selected[[paste(scheme, panels$Database[i], sep = "__")]]
      if (is.null(selection)) return(empty_result_plot(paste(panels$Database[i], panels$Method[i])))
      plot_comparative_plot(catalog, selection, jobs, panels$Method[i], thresholds, display, style)
    }, file.path(output_dir, "09_combined_figures", paste0(safe_id(scheme), "_four_comparative_panels")),
      n_panels = 4, ncol = 2, width = 26, height = 26, settings = settings))
  }
  paths
}
