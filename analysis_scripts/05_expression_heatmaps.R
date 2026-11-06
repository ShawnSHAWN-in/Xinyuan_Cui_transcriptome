if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
inputs <- read_module("01_inputs")
qc <- read_module("02_qc")
deg <- read_module("03_deg")
check_module_inputs(inputs)
check_analysis_settings(deg)
out_dir <- module_dir("06_expression_heatmaps")
effects <- build_effect_heatmap_mat(deg$comparisons)
available <- intersect(rownames(effects), rownames(qc$log_cpm))
expression <- qc$log_cpm[available, , drop = FALSE]
representatives <- unique(unlist(lapply(deg$comparisons, function(item) {
  tbl <- item$significant
  tbl <- tbl[order(tbl$padj, -abs(tbl$logFC), tbl$Geneid), , drop = FALSE]
  unlist(lapply(c("Up", "Down"), function(direction) utils::head(tbl$Geneid[tbl$Result_expression == direction], cfg$top_genes_each)), use.names = FALSE)
}), use.names = FALSE))
save_table(tibble::rownames_to_column(as.data.frame(effects), "Geneid"), file.path(out_dir, "01_gene_by_comparison_logFC.tsv"))
save_table(tibble::rownames_to_column(as.data.frame(expression), "Geneid"), file.path(out_dir, "02_DEG_sample_logCPM.tsv"))
save_table(data.frame(Geneid = representatives), file.path(out_dir, "03_representative_genes.tsv"))
plots <- list()
plots$all_expression <- plot_matrix_heatmap(expression, "DEG union: sample expression (row Z)", TRUE, inputs$sample_info)
plots$all_effects <- plot_matrix_heatmap(effects, "DEG union: log2FC across comparisons", FALSE)
rep_expression <- expression[intersect(representatives, rownames(expression)), , drop = FALSE]
rep_effects <- effects[intersect(representatives, rownames(effects)), , drop = FALSE]
plots$representative_expression <- plot_matrix_heatmap(rep_expression, "Representative DEG: sample expression (row Z)", TRUE, inputs$sample_info)
plots$representative_effects <- plot_matrix_heatmap(rep_effects, "Representative DEG: log2FC", FALSE)
for (name in names(plots)) save_plot_pair(plots[[name]], file.path(out_dir, name), width = if (grepl("expression", name)) 16 else 13,
                                        height = 15)
combine_plots(plots[c("representative_expression", "representative_effects")],
               file.path(module_dir("09_combined_figures"), "11_representative_expression_and_effects"), 2, 24, 16)

# BP主题由用户配置
themes <- list()
if (length(cfg$theme_go_ids)) {
  project <- list(design = inputs$design, public_data = inputs$public_data, enrichment = NULL, deg = deg$comparisons,
                   sample_info = inputs$sample_info, normalized_expression = qc$log_cpm)
  trans <- tranexp_object$new(project)
  for (id in cfg$theme_go_ids) {
    name <- safe_name(id)
    sample_plot <- trans$pheatmap_go_sample2sample(GO_term = id)
    comparison_plot <- trans$pheatmap_go_group2group(GO_term = id)
    themes[[id]] <- list(sample = sample_plot, effect = comparison_plot)
    save_table(tibble::rownames_to_column(as.data.frame(sample_plot$matrix), "Geneid"), file.path(out_dir, paste0(name, "_sample_logCPM.tsv")))
    save_table(tibble::rownames_to_column(as.data.frame(comparison_plot$matrix), "Geneid"), file.path(out_dir, paste0(name, "_comparison_logFC.tsv")))
    combine_plots(list(Expression = sample_plot$plot, Effect = comparison_plot$plot), file.path(out_dir, paste0(name, "_paired_heatmaps")), 2, 24, 14)
  }
} else record_status("05_heatmaps", "BP_themes", "NOT_CONFIGURED", "在00_config.R设置theme_go_ids后启用；不预设专题基因")
save_module(list(plots = plots, effects = effects, representatives = representatives, themes = themes), "05_heatmaps")
record_status("05_heatmaps", "expression_and_effects", "SUCCESS", paste(nrow(effects), "genes"))
