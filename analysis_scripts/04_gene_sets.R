if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
inputs <- read_module("01_inputs")
deg <- read_module("03_deg")
check_module_inputs(inputs)
check_analysis_settings(deg)
out_dir <- module_dir("05_gene_sets")
combined_dir <- module_dir("09_combined_figures")
legacy_deg <- lapply(deg$comparisons, function(item) list(DEG_tbl = item$significant))
plots <- list()
plots$counts_and_proportions <- if (!any(deg$summary$Up + deg$summary$Down)) empty_panel_plot("DEG counts", "No significant genes") else
  attempt_plot("04_sets", "counts_and_proportions", {
  run_DEG_visualization(legacy_deg, public_data = inputs$public_data, type = "bar", out_dir = NULL)
})
save_plot_pair(plots$counts_and_proportions, file.path(out_dir, "01_DEG_counts_and_proportions"), 18, 7)
get_sets <- function(comparisons, direction) stats::setNames(lapply(comparisons, function(name) {
  tbl <- deg$comparisons[[name]]$significant
  unique(tbl$Geneid[tbl$Result_expression == direction])
}), comparisons)
venn_by_time <- list()
venn_by_treatment <- list()
for (time in cfg$time_levels) {
  names <- cfg$contrast_levels[endsWith(cfg$contrast_levels, paste0("NC_", time))]
  for (direction in c("Up", "Down")) {
    item <- paste(time, direction, sep = "_")
    p <- venn_plot(get_sets(names, direction)) + ggplot2::labs(title = paste("Time", time, direction))
    venn_by_time[[item]] <- p
    save_plot_pair(p, file.path(out_dir, paste0("02_venn_time_", item)), 8, 6)
  }
}
for (treatment in setdiff(cfg$treatment_levels, "NC")) {
  names <- cfg$contrast_levels[startsWith(cfg$contrast_levels, paste0(treatment, "_"))]
  for (direction in c("Up", "Down")) {
    item <- paste(treatment, direction, sep = "_")
    p <- venn_plot(get_sets(names, direction)) + ggplot2::labs(title = paste(treatment, direction))
    venn_by_treatment[[item]] <- p
    save_plot_pair(p, file.path(out_dir, paste0("03_venn_treatment_", item)), 8, 6)
  }
}
combine_plots(venn_by_time, file.path(combined_dir, "07_Venn_by_time_and_direction"), 2, 18, 18)
combine_plots(venn_by_treatment, file.path(combined_dir, "08_Venn_by_treatment_and_direction"), 2, 18, 18)
upset <- flower <- list()
for (direction in c("Up", "Down")) {
  sets <- get_sets(cfg$contrast_levels, direction)
  table <- gene_set_intersections_tbl(sets)
  save_table(table, file.path(out_dir, paste0("04_intersections_", direction, ".tsv")))
  membership <- data.frame(Geneid = unique(unlist(sets, use.names = FALSE)))
  for (name in names(sets)) membership[[name]] <- membership$Geneid %in% sets[[name]]
  save_table(membership, file.path(out_dir, paste0("05_gene_membership_", direction, ".tsv")))
  upset[[direction]] <- upsetr_plot(sets) + patchwork::plot_annotation(title = paste(direction, "DEG intersections"))
  flower[[direction]] <- plot_gene_sets_flower_plot(sets, paste(direction, "DEG common/unique genes"))
  save_plot_pair(upset[[direction]], file.path(out_dir, paste0("06_upset_", direction)), 15, 10)
  save_plot_pair(flower[[direction]], file.path(out_dir, paste0("07_flower_", direction)), 12, 11)
}
combine_plots(upset, file.path(combined_dir, "09_UpSet_up_down"), 1, 16, 20)
combine_plots(flower, file.path(combined_dir, "10_flower_up_down"), 2, 23, 11)
plots$venn_by_time <- venn_by_time
plots$venn_by_treatment <- venn_by_treatment
plots$upset <- upset
plots$flower <- flower
save_module(plots, "04_gene_sets")
record_status("04_gene_sets", "DEG_sets", "SUCCESS")
