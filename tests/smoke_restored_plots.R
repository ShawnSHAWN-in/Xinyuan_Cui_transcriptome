project <- normalizePath(getwd())
setwd(project)
for (file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(file)
if (!all(file.exists(file.path("tests/fixtures", c("counts.tsv", "samples.tsv"))))) {
  stop("请先执行Rscript tests/generate_baselines.R --inputs-only")
}
Sys.setenv(ENGINEERING_PROJECT_ROOT = project)
source("tests/unit/helper-project.R")
metadata <- read_sample_metadata_tbl("tests/fixtures/samples.tsv", config("variables"))
counts <- read_count_data_lst("tests/fixtures/counts.tsv")
figures <- config("figures")
qc_settings <- read_qc_settings_lst(figures)
qc <- calc_qc_lst(counts, metadata, qc_settings)
expected_dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts$counts))
expected_keep <- edgeR::filterByExpr(expected_dge, group = metadata$Condition_Time[match(colnames(counts$counts), metadata$Sample)])
stopifnot(identical(qc$keep_pca, expected_keep), all.equal(qc$correlation, cor(qc$log_cpm[expected_keep, , drop = FALSE])))
jobs <- make_comparison_jobs_lst(metadata, config("comparisons"))
annotation <- tibble::tibble(Geneid = rownames(counts$counts), SYMBOL = paste0("Symbol", seq_len(nrow(counts$counts))),
  ENTREZID = as.character(9000000 + seq_len(nrow(counts$counts))), GENENAME = "Synthetic gene")
results <- lapply(jobs, function(job) attach_gene_annotation_tbl(fit_comparison_lst(counts, job, config("analysis")$model, legacy_files), annotation))
limits <- calculate_deg_plot_limits_lst(results)
context <- new_legacy_context(legacy_files, list(direction_colors = legacy_deg_colors()))
for (type in c("volcano", "rank", "ma")) {
  actual <- plot_deg_plot(results[[1]], config("analysis")$filter, type, context, limits)
  tbl <- as.data.frame(results[[1]]$table)
  tbl$padj <- pmax(tbl$padj, .Machine$double.xmin)
  old <- context$run_deg_visualization(list(deg_tbl = tbl, contrast_names = results[[1]]$job$comparison), type = type,
    label_top_n = 8, out_dir = NULL, pvalue_cutoff = config("analysis")$filter$fdr, logFC_cutoff = config("analysis")$filter$abs_logfc)
  stopifnot(identical(actual$data, old$data), identical(length(actual$layers), length(old$layers)))
  for (i in seq_along(actual$layers)) stopifnot(identical(actual$layers[[i]]$data, old$layers[[i]]$data),
    identical(actual$layers[[i]]$aes_params, old$layers[[i]]$aes_params))
  message("LEGACY_LABEL_COLOR_PASS ", type)
}
root <- Sys.getenv("PLOT_TEST_ROOT", unset = tempfile("restored_plots_"))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
root <- normalizePath(root)
out <- file.path(root, "rendered")
export <- figures$export
export$dpi <- 48
export$max_pixels <- 800000
deg_files <- unlist(lapply(results, function(result) export_deg_plot_files(result, config("analysis")$filter, legacy_files, export, out, limits)))
deg_files <- c(deg_files, export_deg_combined_files(results, config("analysis")$filter, legacy_files, export, out, figures$display, limits))
stopifnot(all(file.exists(deg_files)), all(file.info(deg_files)$size > 0))

go <- data.table::fread("tests/fixtures/gene2go.tsv", data.table = FALSE)
backgrounds <- list(GO = list(status = "SUCCESS", background = list(TERM2GENE = data.frame(TERM = go$GO_ID, GENE = as.character(go$GeneID)))),
  KEGG = list(status = "SUCCESS", background = list(TERM2GENE = data.table::fread("tests/fixtures/KEGG_TERM2GENE.tsv", data.table = FALSE))))
files <- export_qc_files(qc, metadata, figures$display, figures$style, export, out, legacy_files, qc_settings,
  list(annotation = annotation), backgrounds)
stopifnot(all(file.exists(files)))
status <- data.table::fread(file.path(out, "03_sample_qc", "QC_status.tsv"))
stopifnot(all(status$Status == "SUCCESS"))
all_plots <- c(files, deg_files)
pngs <- all_plots[tools::file_ext(all_plots) == "png"]
for (file in pngs) stopifnot(prod(dim(png::readPNG(file, native = TRUE))) <= export$max_pixels)
writeLines(c("RESTORED_DEG_QC_PASS", paste("QC PDF/PNG files:", sum(tools::file_ext(files) %in% c("pdf", "png")))), file.path(root, "PASS.txt"))
cat("RESTORED_DEG_QC_PASS\nQC_PDF_PNG_FILES=", sum(tools::file_ext(files) %in% c("pdf", "png")), "\n", sep = "")
cat("TEST_ROOT=", root, "\n", sep = "")
