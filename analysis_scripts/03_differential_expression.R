if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
inputs <- read_module("01_inputs")
check_module_inputs(inputs)
out_dir <- module_dir("04_differential_expression")
comparisons <- inputs$contrasts$Comparison
results <- stats::setNames(vector("list", length(comparisons)), comparisons)
plots <- list(volcano = list(), rank = list(), ma = list())
summary <- list()
for (i in seq_along(comparisons)) {
  comparison <- comparisons[i]
  message("差异分析：", comparison)
  dir <- file.path(out_dir, comparison)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  res <- run_calculate_differential_gene(inputs$design, contrast_index = i, method = cfg$method, out_dir = NULL)
  sig <- run_calculate_significant_DEG(res, padj_threshold = cfg$fdr_cutoff, abs_logFC_threshold = cfg$abs_logfc_cutoff, out_dir = NULL)
  tbl <- attach_annotation_tbl(res$deg_tbl, inputs$annotation) |> restore_factors_tbl()
  significant <- if (is.null(sig)) tbl[0, , drop = FALSE] else tbl[tbl$Geneid %in% sig$DEG_tbl$Geneid, , drop = FALSE]
  results[[comparison]] <- list(all = tbl, significant = significant, raw = res)
  save_table(tbl, file.path(dir, "01_all_tested_genes.tsv"))
  save_table(significant, file.path(dir, "02_significant_genes.tsv"))
  for (direction in c("Up", "Down")) save_table(significant[significant$Result_expression == direction, , drop = FALSE],
                                                file.path(dir, paste0("03_", direction, "_genes.tsv")))
  openxlsx::write.xlsx(list(Up = significant[significant$Result_expression == "Up", , drop = FALSE],
                           Down = significant[significant$Result_expression == "Down", , drop = FALSE]),
                       file.path(dir, "04_significant_genes.xlsx"), overwrite = TRUE)
  saveRDS(results[[comparison]], file.path(dir, "05_differential_expression.rds"))
  summary[[i]] <- data.frame(Comparison = comparison, Tested = nrow(tbl), Up = sum(significant$Result_expression == "Up"),
                             Down = sum(significant$Result_expression == "Down"), Method = cfg$method,
                             FDR = cfg$fdr_cutoff, AbsLogFC = cfg$abs_logfc_cutoff)
  for (type in names(plots)) {
    plot_input <- res
    plot_input$deg_tbl$padj <- pmax(plot_input$deg_tbl$padj, .Machine$double.xmin)
    p <- attempt_plot("03_deg", paste(comparison, type), {
      run_deg_visualization(plot_input, public_data = inputs$public_data, type = type, label_top_n = 8, out_dir = NULL,
                            pvalue_cutoff = cfg$fdr_cutoff, logFC_cutoff = cfg$abs_logfc_cutoff)
    })
    if (!is.null(p)) p <- p + ggplot2::labs(title = comparison) + ggplot2::theme(panel.grid = ggplot2::element_blank())
    plots[[type]][comparison] <- list(p)
    save_plot_pair(p, file.path(dir, paste0("plot_", type)), 8, 7)
  }
  record_status("03_deg", comparison, "SUCCESS", paste(nrow(tbl), "tested;", nrow(significant), "significant"))
}
summary <- dplyr::bind_rows(summary)
save_table(summary, file.path(out_dir, "00_deg_summary.tsv"))
all_effects <- unlist(lapply(results, function(x) x$all$logFC))
all_fdr <- unlist(lapply(results, function(x) x$all$padj))
effect_limit <- max(abs(all_effects[is.finite(all_effects)])) * 1.05
significance_limit <- max(-log10(pmax(all_fdr[is.finite(all_fdr)], .Machine$double.xmin))) * 1.05
for (type in names(plots)) {
  for (comparison in names(plots[[type]])) {
    p <- plots[[type]][[comparison]]
    if (is.null(p)) next
    if (type == "volcano") p <- p + ggplot2::coord_cartesian(xlim = c(-effect_limit, effect_limit), ylim = c(0, significance_limit))
    if (type %in% c("rank", "ma")) p <- p + ggplot2::coord_cartesian(ylim = c(-effect_limit, effect_limit))
    plots[[type]][[comparison]] <- p
    save_plot_pair(p, file.path(out_dir, comparison, paste0("plot_", type)), 8, 7)
  }
  combine_plots(plots[[type]], file.path(module_dir("09_combined_figures"), paste0("06_DEG_", type, "_nine_contrasts")),
                 ncol = 3, width = 23, height = 22)
}
deg <- list(comparisons = results, summary = summary, plots = plots, settings = cfg[c("method", "fdr_cutoff", "abs_logfc_cutoff")])
save_module(deg, "03_deg")
