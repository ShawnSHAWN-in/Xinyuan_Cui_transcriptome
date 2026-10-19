# 不依赖富集算法的绘图契约
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- normalizePath(file.path(dirname(sub("^--file=", "", arg[1])), ".."))
root <- Sys.getenv("CORE_TEST_ROOT")
if (!nzchar(root)) root <- tempfile("compare_plot_contract_", tmpdir = tempdir())
dir.create(root, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = root)
source(file.path(script_dir, "00_config.R"))
cfg$dpi <- 60L
mock <- stats::setNames(lapply(cfg$contrast_levels, function(comparison) {
  make_item <- function(ids) list(
    ORA = list(Up = data.frame(ID = ids, Description = paste("Synthetic term", 1:3), GeneRatio = c("4/10", "5/10", "2/10"),
                                Count = c(4, 5, 2), p.adjust = c(0.001, 0.01, 0.02)),
               Down = data.frame(ID = ids, Description = paste("Synthetic term", 1:3), GeneRatio = c("2/10", "3/10", "4/10"),
                                  Count = c(2, 3, 4), p.adjust = c(0.02, 0.01, 0.001))),
    GSEA = data.frame(ID = ids, Description = paste("Synthetic term", 1:3), NES = c(1.8, -1.9, 1.5), p.adjust = c(0.001, 0.01, 0.02)))
  list(GO = make_item(c("GO:0006260", "GO:0006412", "GO:0006096")), KEGG = make_item(c("fca00010", "fca00020", "fca00030")))
}), cfg$contrast_levels)
trans <- tranexp_object$new(list(design = NULL, public_data = NULL, enrichment = mock, deg = NULL,
                                 sample_info = NULL, normalized_expression = NULL))
plots <- list(GO_ORA = trans$compareCluster_GO_enrich_plot(), GO_GSEA = trans$compareCluster_GSEA_GO_enrich_plot(),
               KEGG_ORA = trans$compareCluster_KEGG_enrich_plot(), KEGG_GSEA = trans$compareCluster_GSEA_KEGG_enrich_plot())
for (name in names(plots)) {
  stopifnot(identical(levels(plots[[name]]$data$Comparison), cfg$contrast_levels))
  stopifnot(setequal(unique(as.character(plots[[name]]$data$Direction)), c("Up", "Down")))
  save_plot_pair(plots[[name]], file.path(root, "compare_contract", name), 13, 6)
}
stopifnot(nrow(plots$GO_ORA$data) == 54L, nrow(plots$GO_GSEA$data) == 27L)
stopifnot(all(plots$GO_ORA$data$Magnitude %in% c(0.2, 0.3, 0.4, 0.5)))
stopifnot(inherits(venn_plot(list(A = "x", B = "x", C = "y")), "ggplot"))
stopifnot(nrow(gene_set_intersections_tbl(list(A = "x", B = "x"))) == 1L)
combine_plots(plots, file.path(root, "compare_contract", "four_comparative_panels"), 2, 24, 14)
capture.output(warnings(), file = file.path(root, "compare_contract_warnings.txt"))
cat("COMPARATIVE_PLOT_CONTRACT_PASS\n")
