export_gene_set_files <- function(selected, engine_files, plot_helper, style, settings, output_dir, display = NULL) {
  if (!is.null(display)) selected <- order_comparison_results_lst(selected, display$compare_by)
  context <- new_legacy_context(engine_files, plot_helper = plot_helper)
  paths <- character()
  for (scheme in unique(vapply(selected, function(x) x$job$scheme, character(1)))) {
    items <- Filter(function(x) x$job$scheme == scheme, selected)
    dir <- file.path(output_dir, "05_gene_sets", safe_id(scheme))
    sets <- stats::setNames(lapply(items, function(x) unique(x$table$Geneid)), vapply(items, function(x) x$job$comparison, character(1)))
    genes <- unique(unlist(sets, use.names = FALSE))
    membership <- data.frame(Geneid = genes)
    for (name in names(sets)) membership[[name]] <- genes %in% sets[[name]]
    paths <- c(paths, write_table_path(membership, file.path(dir, "membership.tsv")))
    similarity <- matrix(0, length(sets), length(sets), dimnames = list(names(sets), names(sets)))
    for (i in seq_along(sets)) for (j in seq_along(sets)) {
      union <- length(union(sets[[i]], sets[[j]]))
      similarity[i, j] <- if (union) length(intersect(sets[[i]], sets[[j]])) / union else NA_real_
    }
    paths <- c(paths, write_table_path(tibble::rownames_to_column(as.data.frame(similarity), "Comparison"), file.path(dir, "Jaccard.tsv")))
    p <- if (any(is.finite(similarity))) pheatmap::pheatmap(similarity, cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey75",
      border_color = NA, main = "DEG Jaccard similarity", silent = TRUE)$gtable else empty_result_plot("DEG Jaccard similarity", "No significant genes")
    paths <- c(paths, save_single_plot_files(p, file.path(dir, "01_Jaccard"), settings, 12, 10))
    for (type in c("upset", "flower")) {
      paths <- c(paths, save_panel_files(function(i) {
        direction <- c("Up", "Down")[i]
        dir_sets <- stats::setNames(lapply(items, function(x) unique(x$table$Geneid[x$table$Result_expression == direction])), names(sets))
        plot <- if (type == "upset") context$upsetr_plot(dir_sets) else context$plot_gene_sets_flower_plot(dir_sets, direction)
        apply_legacy_theme(plot, style)
      }, file.path(dir, paste0("02_Up_Down_", type)), 2, 2, 20, 9, settings))
    }
    keys <- vapply(items, function(x) {
      fields <- x$job$strata
      if (!length(fields)) "all" else paste(unlist(x$job$labels[fields]), collapse = "_")
    }, character(1))
    for (key in unique(keys)) {
      block <- items[keys == key]
      if (length(block) < 2L || length(block) > 3L) next
      paths <- c(paths, save_panel_files(function(i) {
        direction <- c("Up", "Down")[i]
        dir_sets <- stats::setNames(lapply(block, function(x) unique(x$table$Geneid[x$table$Result_expression == direction])),
                                   vapply(block, function(x) x$job$comparison, character(1)))
        context$venn_plot(dir_sets) + ggplot2::labs(title = paste(key, direction))
      }, file.path(dir, paste0("03_Venn_", safe_id(key))), 2, 2, 18, 8, settings))
    }
  }
  paths
}

make_heatmap_grob <- function(matrix, metadata, style, title, scale_rows = TRUE, annotation_columns = c("Condition", "Time_f")) {
  if (!nrow(matrix)) return(as_panel_grob(empty_result_plot(title, "No eligible genes")))
  if (scale_rows) {
    sd <- apply(matrix, 1, stats::sd, na.rm = TRUE)
    matrix <- matrix[is.finite(sd) & sd > 0, , drop = FALSE]
    if (!nrow(matrix)) return(as_panel_grob(empty_result_plot(title, "No variable genes")))
    matrix <- t(scale(t(matrix)))
  }
  columns <- intersect(as.character(unlist(annotation_columns)), names(metadata))
  annotation <- as.data.frame(metadata[, columns, drop = FALSE])
  rownames(annotation) <- metadata$Sample
  colors <- if ("Condition" %in% columns && all(unique(as.character(metadata$Condition)) %in% names(style$condition_colors))) {
    list(Condition = unlist(style$condition_colors))
  } else NA
  pheatmap::pheatmap(matrix, cluster_cols = FALSE, cluster_rows = nrow(matrix) > 1 && all(is.finite(matrix)), annotation_col = annotation,
    annotation_colors = colors, show_rownames = nrow(matrix) <= 70, border_color = NA, na_col = "grey75", main = title,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100), silent = TRUE)$gtable
}

export_expression_heatmap_files <- function(qc, selected, metadata, annotation_go, display, style, settings, output_dir) {
  selected <- order_comparison_results_lst(selected, display$compare_by)
  ordered <- sort_metadata_tbl(metadata, display$order_by)
  dir <- file.path(output_dir, "06_expression_heatmaps")
  all_sig <- dplyr::bind_rows(lapply(selected, function(x) x$table))
  all_genes <- unique(all_sig$Geneid)
  top <- unique(unlist(lapply(selected, function(x) unlist(lapply(c("Up", "Down"), function(direction) {
    tbl <- x$table[x$table$Result_expression == direction, , drop = FALSE]
    utils::head(tbl$Geneid[order(tbl$padj)], display$top_genes_each)
  }))), use.names = FALSE))
  max_genes <- if (is.null(display$heatmap_max_genes)) 500L else display$heatmap_max_genes
  if (length(top) > max_genes) {
    variability <- apply(qc$log_cpm[top, , drop = FALSE], 1, stats::mad)
    top <- names(sort(variability, decreasing = TRUE))[seq_len(max_genes)]
  }
  paths <- c(write_table_path(data.frame(Geneid = all_genes), file.path(dir, "all_significant_genes.tsv")),
             write_table_path(data.frame(Geneid = top), file.path(dir, "displayed_genes.tsv")))
  for (assay in c("log_cpm", "cpm", "tpm")) {
    matrix <- qc[[assay]][top, ordered$Sample, drop = FALSE]
    paths <- c(paths, save_single_plot_files(make_heatmap_grob(matrix, ordered, style, paste(assay, "row Z-score"),
      annotation_columns = display$heatmap_annotations),
      file.path(dir, paste0("01_", assay, "_row_Z")), settings, 13, min(22, max(7, length(top) * 0.12 + 3))))
  }
  effects <- matrix(NA_real_, length(all_genes), length(selected), dimnames = list(all_genes, vapply(selected, function(x) x$job$id, character(1))))
  for (i in seq_along(selected)) {
    tbl <- selected[[i]]$table
    effects[, i] <- tbl$logFC[match(all_genes, tbl$Geneid)]
  }
  paths <- c(paths, write_table_path(tibble::rownames_to_column(as.data.frame(effects), "Geneid"), file.path(dir, "all_significant_logFC.tsv")))
  display_effects <- effects[top, , drop = FALSE]
  if (nrow(display_effects)) {
    p <- pheatmap::pheatmap(display_effects, cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey75", border_color = NA,
      show_rownames = nrow(display_effects) <= 70, main = "Significant logFC; grey: not significant here", silent = TRUE)$gtable
  } else p <- empty_result_plot("Significant logFC", "No significant genes")
  paths <- c(paths, save_single_plot_files(p, file.path(dir, "02_significant_logFC"), settings, 14, min(22, max(7, length(top) * 0.12 + 3))))
  for (id in as.character(unlist(display$theme_go_ids))) {
    bg <- annotation_go$background$TERM2GENE
    genes <- annotation_go$annotation$Geneid[annotation_go$annotation$ENTREZID %in% bg$GENE[bg$TERM == id]]
    genes <- intersect(top, genes)
    p <- make_heatmap_grob(qc$log_cpm[genes, ordered$Sample, drop = FALSE], ordered, style, paste("GO-BP", id),
      annotation_columns = display$heatmap_annotations)
    paths <- c(paths, save_single_plot_files(p, file.path(dir, paste0("03_", safe_id(id))), settings, 13, 9))
  }
  paths
}

export_curve_combined_files <- function(records, jobs, thresholds, display, style, settings, engine_files, plot_helper, output_dir) {
  context <- new_legacy_context(engine_files, list(fdr_cutoff = thresholds$fdr), plot_helper)
  paths <- character()
  for (scheme in unique(vapply(jobs, `[[`, character(1), "scheme"))) {
    scheme_jobs <- order_comparison_jobs_lst(Filter(function(x) x$scheme == scheme, jobs), display$compare_by)
    for (db in c("GO", "KEGG")) {
      block <- Filter(function(x) x$Scheme == scheme && x$Database == db, records)
      curve_ids <- unlist(lapply(block, function(record) {
        value <- readRDS(record$object_file)
        tbl <- value$table
        tbl <- tbl[is.finite(tbl$p.adjust) & tbl$p.adjust < thresholds$fdr, , drop = FALSE]
        utils::head(tbl$ID[order(tbl$p.adjust)], display$curve_top_each)
      }), use.names = FALSE)
      ids <- utils::head(names(sort(table(curve_ids), decreasing = TRUE)), 3)
      for (id in ids) {
        pages <- split(seq_along(scheme_jobs), ceiling(seq_along(scheme_jobs) / settings$combine_per_page))
        for (page in seq_along(pages)) {
          page_jobs <- scheme_jobs[pages[[page]]]
          paths <- c(paths, save_panel_files(function(i) {
            job <- page_jobs[[i]]
            matches <- Filter(function(x) x$Job == job$id, block)
            if (!length(matches)) return(empty_result_plot(job$comparison, "Database unavailable"))
            value <- readRDS(matches[[1]]$object_file)
            res <- significant_result(value$result, thresholds$fdr, id)
            if (is.null(res)) return(empty_result_plot(job$comparison, "Not significant or not tested"))
            plots <- context$plot_gsea_gseaplot2(res, pvalue_cutoff = thresholds$fdr, out_dir = NULL)
            if (!length(plots)) stop("曲线未生成：", job$id, " ", id)
            p <- apply_legacy_theme(plots[[1]], style)
            grob <- as_panel_grob(p)
            grid::grid.text(job$comparison, y = 0.99, just = "top", gp = grid::gpar(fontsize = style$title_size))
            grid::pushViewport(grid::viewport(y = 0.47, height = 0.92))
            grid::grid.draw(grob, recording = FALSE)
            grid::popViewport()
            NULL
          }, file.path(output_dir, "09_combined_figures", paste(safe_id(scheme), db, safe_id(id), "curves", page, sep = "_")),
          length(page_jobs), settings$combine_columns, 24, ceiling(length(page_jobs) / settings$combine_columns) * 7, settings))
        }
      }
    }
  }
  if (!length(paths)) paths <- write_table_path(data.frame(Status = "NO_SIGNIFICANT_CURVES"), file.path(output_dir, "09_combined_figures", "curve_status.tsv"))
  paths
}

export_pathview_files <- function(records, thresholds, display, reference, settings, output_dir) {
  root <- file.path(output_dir, "07_enrichment", "Pathview")
  if (!isTRUE(display$pathview)) return(write_table_path(data.frame(Status = "DISABLED"), file.path(root, "status.tsv")))
  check_packages(c("pathview", "png"))
  files <- character()
  status <- list()
  for (record in Filter(function(x) x$Database == "KEGG" && x$Status == "SUCCESS", records)) {
    value <- readRDS(record$object_file)
    tbl <- value$table
    tbl <- tbl[is.finite(tbl$p.adjust) & tbl$p.adjust < thresholds$fdr, , drop = FALSE]
    ids <- utils::head(tbl$ID[order(tbl$p.adjust)], display$pathview_top)
    dir <- file.path(root, safe_id(record$Job))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir <- normalizePath(dir)
    for (id in ids) {
      old <- getwd()
      tryCatch({
        setwd(dir)
        pathview::pathview(gene.data = value$ranking, pathway.id = sub("^[[:alpha:]]+", "", id), species = reference$kegg_code,
          gene.idtype = "entrez", kegg.dir = dir, same.layer = FALSE, limit = list(gene = max(abs(value$ranking), na.rm = TRUE)),
          low = list(gene = "#2166AC"), mid = list(gene = "#F7F7F7"), high = list(gene = "#B2182B"))
      }, finally = setwd(old))
      status[[length(status) + 1L]] <- data.frame(Job = record$Job, ID = id, Status = "SUCCESS")
    }
    files <- c(files, list.files(dir, pattern = "pathview.*\\.png$", full.names = TRUE))
  }
  paths <- files
  for (name in unique(basename(files))) {
    group <- files[basename(files) == name]
    pages <- split(seq_along(group), ceiling(seq_along(group) / settings$combine_per_page))
    for (page in seq_along(pages)) {
      page_files <- group[pages[[page]]]
      paths <- c(paths, save_panel_files(function(i) {
        image <- png::readPNG(page_files[i], native = TRUE)
        grid::grid.text(basename(dirname(page_files[i])), y = 0.98, just = "top")
        panel_ratio <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE) /
          grid::convertHeight(grid::unit(0.9, "npc"), "in", valueOnly = TRUE)
        ratio <- ncol(image) / nrow(image)
        grid::grid.raster(image, y = 0.46, width = min(1, ratio / panel_ratio), height = 0.9 * min(1, panel_ratio / ratio), interpolate = FALSE)
        NULL
      }, file.path(output_dir, "09_combined_figures", paste0("Pathview_RASTER_", safe_id(name), "_", page)), length(page_files),
        settings$combine_columns, 24, ceiling(length(page_files) / settings$combine_columns) * 7, settings))
    }
  }
  rows <- if (length(status)) dplyr::bind_rows(status) else data.frame(Status = "NO_SIGNIFICANT_PATHWAYS")
  c(paths, write_table_path(rows, file.path(root, "status.tsv")))
}
