# 独立验证富集与比较图
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- normalizePath(file.path(dirname(sub("^--file=", "", arg[1])), ".."))
root <- Sys.getenv("CORE_TEST_ROOT")
if (!nzchar(root)) stop("请将CORE_TEST_ROOT设置为smoke_core.R输出的TEST_ROOT")
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = root)
source(file.path(script_dir, "00_config.R"))
cfg$dpi <- 60L
cfg$run_pathview <- FALSE
cfg$curve_top_each <- 1L
cfg$top_terms <- 4L
cfg$compare_max_terms <- 10L
cfg$max_gs_size <- 350L
require_enrichment_packages()
inputs <- read_module("01_inputs")
go_ids <- c("GO:0006260", "GO:0006281", "GO:0006412", "GO:0006096", "GO:0005975", "GO:0006810")
sets <- lapply(seq_along(go_ids), function(i) as.character(9000000 + 30 + (i - 1) * 15 + seq_len(15)))
assoc <- dplyr::bind_rows(lapply(seq_along(sets), function(i) data.frame(TERM = go_ids[i], GENE = sets[[i]])))
inputs$backgrounds$GO <- list(TERM2GENE = assoc, TERM2NAME = data.frame(TERM = go_ids, NAME = paste("Synthetic BP test set", 1:6)))
kegg_ids <- paste0("fca", c("00010", "00020", "00030", "00040", "00051", "00052"))
assoc$TERM <- kegg_ids[match(assoc$TERM, go_ids)]
inputs$backgrounds$KEGG <- list(TERM2GENE = assoc, TERM2NAME = data.frame(TERM = kegg_ids, NAME = paste("Synthetic KEGG test set", 1:6)))
save_module(inputs, "01_inputs")
for (file in c("06_enrichment.R", "07_enrichment_plots.R", "08_comparative_enrichment.R", "09_combine_export.R")) {
  message("TEST MODULE: ", file)
  source(file.path(script_dir, file))
}
enrichment <- read_module("06_enrichment")
comparative <- read_module("08_comparative")
stopifnot(length(enrichment$results) == 9L)
for (db in c("GO", "KEGG")) {
  stopifnot(nrow(comparative$tables[[db]]$ORA) > 0, nrow(comparative$tables[[db]]$GSEA) > 0)
  p <- comparative$plots[[paste0(db, "_ORA")]]
  stopifnot(identical(levels(p$data$Comparison), cfg$contrast_levels))
  stopifnot(all(c("Up", "Down") %in% as.character(p$data$Direction)))
}
statuses <- data.table::fread(file.path(cfg$output_dir, "10_objects_and_summary", "module_status.tsv"))
if (any(statuses$Status == "FAILED")) stop("存在失败的图：", paste(statuses$Item[statuses$Status == "FAILED"], collapse = ", "))
writeLines("ENRICHMENT_SMOKE_PASS", file.path(cfg$output_dir, "10_objects_and_summary", "enrichment_test_pass.txt"))
cat("ENRICHMENT_SMOKE_PASS\n")
