if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
inputs <- read_module("01_inputs")
deg <- read_module("03_deg")
enrichment <- read_module("06_enrichment")
check_module_inputs(inputs)
check_analysis_settings(deg)
check_enrichment_settings(enrichment, inputs)
out_dir <- module_dir("08_comparative_enrichment")
project <- list(design = inputs$design, public_data = inputs$public_data, enrichment = enrichment$results, deg = deg$comparisons,
                 sample_info = inputs$sample_info, normalized_expression = NULL)
trans <- tranexp_object$new(project)
plots <- list()
tables <- list()
for (db in c("GO", "KEGG")) {
  ora <- collect_comparative_tbl(enrichment$results, db, "ORA")
  gsea <- collect_comparative_tbl(enrichment$results, db, "GSEA")
  all <- dplyr::bind_rows(ora, gsea)
  if (!nrow(all)) {
    unavailable <- is.null(inputs$backgrounds[[db]])
    record_status("08_comparative", db, if (unavailable) "SKIPPED" else "NO_SIGNIFICANT_RESULT",
                   if (unavailable) "数据库背景不可用" else "无显著富集，保留空结果说明")
    selected_ids <- character()
  } else selected_ids <- select_representative_ids(all)
  if (length(selected_ids)) attr(selected_ids, "labels") <- stats::setNames(all$Description[match(selected_ids, all$ID)], selected_ids)
  save_table(data.frame(ID = unname(selected_ids)), file.path(out_dir, paste0(db, "_displayed_terms.tsv")))
  save_table(ora, file.path(out_dir, paste0(db, "_ORA_all_comparisons.tsv")))
  save_table(gsea, file.path(out_dir, paste0(db, "_GSEA_all_comparisons.tsv")))
  tables[[db]] <- list(ORA = ora, GSEA = gsea, selected_ids = selected_ids)
  functions <- if (db == "GO") c("compareCluster_GO_enrich_plot", "compareCluster_GSEA_GO_enrich_plot") else
    c("compareCluster_KEGG_enrich_plot", "compareCluster_GSEA_KEGG_enrich_plot")
  pair <- list()
  for (i in seq_along(functions)) {
    method <- c("ORA", "GSEA")[i]
    p <- trans[[functions[i]]](contrast_range = seq_along(cfg$contrast_levels), selected_ids = selected_ids)
    pair[[method]] <- p
    plots[[paste(db, method, sep = "_")]] <- p
    save_plot_pair(p, file.path(out_dir, paste0(db, "_", method, "_comparative_triangles")), width = 16,
                   height = min(24, max(6, length(selected_ids) * 0.38 + 3)))
  }
  combine_plots(pair, file.path(module_dir("09_combined_figures"), paste0("12_", db, "_ORA_GSEA_comparative")), 2, 24,
                 min(26, max(8, length(selected_ids) * 0.45 + 3)))

  # NES与显著性分别编码
  if (length(selected_ids)) {
    nes <- matrix(NA_real_, length(selected_ids), length(cfg$contrast_levels), dimnames = list(selected_ids, cfg$contrast_levels))
    fdr <- nes
    for (comparison in cfg$contrast_levels) {
      item <- enrichment$results[[comparison]][[db]]
      if (is.null(item)) next
      tbl <- result_table(item$GSEA)
      if (!nrow(tbl)) next
      idx <- match(selected_ids, tbl$ID)
      nes[, comparison] <- tbl$NES[idx]
      fdr[, comparison] <- tbl$p.adjust[idx]
    }
    save_table(tibble::rownames_to_column(as.data.frame(nes), "ID"), file.path(out_dir, paste0(db, "_NES_matrix.tsv")))
    save_table(tibble::rownames_to_column(as.data.frame(fdr), "ID"), file.path(out_dir, paste0(db, "_GSEA_FDR_matrix.tsv")))
    display <- expand.grid(ID = selected_ids, Comparison = cfg$contrast_levels, stringsAsFactors = FALSE)
    display$NES <- as.vector(nes)
    display$FDR <- as.vector(fdr)
    display$ID <- factor(display$ID, levels = rev(selected_ids))
    display$Comparison <- factor(display$Comparison, levels = cfg$contrast_levels)
    display$Significant <- ifelse(!is.na(display$FDR) & display$FDR < cfg$fdr_cutoff, "*", "")
    labels <- attr(selected_ids, "labels")
    p <- ggplot2::ggplot(display, ggplot2::aes(Comparison, ID, fill = NES)) + ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = Significant), size = 5, na.rm = TRUE) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0, na.value = "#BDBDBD") +
      ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(paste0(labels[x], " [", x, "]"), 55)) +
      ggplot2::labs(x = NULL, y = NULL, title = paste(db, "GSEA: NES; * FDR <", cfg$fdr_cutoff),
                    caption = "Grey: NES not estimated or absent from returned results") +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    plots[[paste0(db, "_NES")]] <- p
    save_plot_pair(p, file.path(out_dir, paste0(db, "_GSEA_NES_heatmap")), width = 16, n_terms = length(selected_ids))
  }
}
combine_plots(plots[c("GO_ORA", "GO_GSEA", "KEGG_ORA", "KEGG_GSEA")],
               file.path(module_dir("09_combined_figures"), "13_four_comparative_enrichment_panels"), 2, 24, 32)
if (length(cfg$theme_go_ids)) {
  for (id in cfg$theme_go_ids) {
    thematic <- list(ORA = trans$compareCluster_GO_enrich_plot(GO_term = id), GSEA = trans$compareCluster_GSEA_GO_enrich_plot(GO_term = id))
    combine_plots(thematic, file.path(out_dir, paste0(safe_name(id), "_ORA_GSEA")), 2, 24, 14)
  }
}
save_module(list(plots = plots, tables = tables), "08_comparative")
record_status("08_comparative", "four_methods", "SUCCESS", "GO-BP/KEGG × ORA/GSEA；未重新计算差异或富集")
