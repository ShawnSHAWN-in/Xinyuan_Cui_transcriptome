# 验证逐面板导出及旧对象兼容
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "analysis_scripts/tests/combine_export_memory.R"
script_dir <- normalizePath(file.path(dirname(script_file), ".."))
root <- Sys.getenv("COMBINE_TEST_ROOT", unset = tempfile("combine_memory_"))
dir.create(root)
fixture <- file.path(root, "input.tsv")
writeLines("fixture", fixture)
cfg <- list(output_dir = file.path(root, "03_analysis"), dpi = 300L, combine_max_pixels = as.numeric(Sys.getenv("COMBINE_TEST_PIXELS", unset = "200000")),
            method = "edgeR", fdr_cutoff = 0.05, abs_logfc_cutoff = 1, curve_top_each = 1L,
            reference_species = "Felis catus", tax_id = 9685L, kegg_code = "fca", go_ontology = "BP", go_ancestors = TRUE,
            min_gs_size = 10L, max_gs_size = 500L, seed = 1L, group_levels = "fixture", contrast_levels = paste0("Comparison", 1:9))
source(file.path(script_dir, "lib", "project_helpers.R"))
source(file.path(script_dir, "lib", "plot_helpers.R"))
inputs <- list(signature = input_signature(fixture), backgrounds = list(GO = list()),
               settings = cfg[c("tax_id", "kegg_code", "go_ontology", "go_ancestors", "group_levels", "contrast_levels")],
               sample_info = data.frame(Samples = 1:36), contrasts = data.frame(Comparison = cfg$contrast_levels))
deg <- list(settings = cfg[c("method", "fdr_cutoff", "abs_logfc_cutoff")])
save_module(inputs, "01_inputs")
save_module(deg, "03_deg")
methods::setClass("MemoryTestGSEA", slots = c(result = "data.frame"))
gsea <- methods::new("MemoryTestGSEA", result = data.frame(ID = "GO:0006260", Description = "Synthetic BP", p.adjust = 0.001, NES = 2))
enrichment <- list(results = stats::setNames(lapply(1:9, function(i) list(GO = list(GSEA = gsea), KEGG = NULL)), cfg$contrast_levels),
                   settings = cfg[c("min_gs_size", "max_gs_size", "seed", "go_ancestors", "fdr_cutoff", "method")],
                   signature = list(inputs = inputs$signature, background = digest::digest(inputs$backgrounds),
                                    deg = unname(tools::md5sum(file.path(cfg$output_dir, "10_objects_and_summary", "03_deg.rds")))))
save_module(enrichment, "06_enrichment")
# 旧绘图文件不可读取
writeLines("DO NOT READ: legacy plot RDS", file.path(cfg$output_dir, "10_objects_and_summary", "07_enrichment_plots.rds"))
for (i in 1:12) {
  dir <- file.path(cfg$output_dir, "07_enrichment", "KEGG", cfg$contrast_levels[(i - 1) %% 9 + 1], "plots", "pathview", paste0("Method", i))
  dir.create(dir, recursive = TRUE)
  large_png <- Sys.getenv("COMBINE_TEST_PATHVIEW_PNG")
  if (nzchar(large_png)) file.copy(large_png, file.path(dir, "fca00010.pathview.png")) else {
    png::writePNG(array(rep(c(1, 0, 0, 1), each = 80 * 120), dim = c(80, 120, 4)), file.path(dir, "fca00010.pathview.png"))
  }
}
# 不包含在pathview中的原图不参与组合
png::writePNG(matrix(1, 10, 10), file.path(cfg$output_dir, "07_enrichment", "KEGG", "fca00020.pathview.png"))
curve_calls <- character()
plot_gsea_gseaplot2 <- function(res, pvalue_cutoff, out_dir) {
  stopifnot(identical(pvalue_cutoff, cfg$fdr_cutoff), is.null(out_dir))
  curve_calls <<- c(curve_calls, res@result$ID)
  p <- ggplot2::ggplot(data.frame(x = 1:20000, y = sin((1:20000) / 3000)), ggplot2::aes(x, y)) + ggplot2::geom_line() + ggplot2::theme_classic()
  list("GO:0006260" = structure(list(p, p + ggplot2::labs(y = "Hits"), p + ggplot2::labs(y = "Rank")), class = "gglist"))
}
# 验证重新打开时从01/03/06读取
rm(inputs, deg, enrichment)
existing_names <- ls()
source(file.path(script_dir, "09_combine_export.R"))
stopifnot(identical(existing_names, setdiff(ls(), "existing_names")), length(curve_calls) == 18L, grDevices::dev.cur() == 1L)
out <- file.path(cfg$output_dir, "09_combined_figures")
png_files <- list.files(out, pattern = "\\.png$", full.names = TRUE)
pdf_files <- list.files(out, pattern = "\\.pdf$", full.names = TRUE)
stopifnot(length(png_files) == 3L, length(pdf_files) == 3L, any(grepl("_page02", png_files)))
for (file in png_files) {
  image <- png::readPNG(file, native = TRUE)
  stopifnot(prod(dim(image)) <= cfg$combine_max_pixels)
}
status <- data.table::fread(file.path(cfg$output_dir, "10_objects_and_summary", "module_status.tsv"))
stopifnot(!any(status$Status == "FAILED"), file.exists(file.path(cfg$output_dir, "10_objects_and_summary", "09_memory_usage.tsv")))
file.copy(pdf_files[grepl("14_GO_", pdf_files)], file.path(root, "curve_success.pdf"))
file.copy(file.path(cfg$output_dir, "10_objects_and_summary", "09_memory_usage.tsv"), file.path(root, "success_memory_usage.tsv"))
# 注入绘图失败，设备仍应关闭
plot_gsea_gseaplot2 <- function(...) stop("INJECTED_RENDER_FAILURE")
suppressWarnings(source(file.path(script_dir, "09_combine_export.R")))
status <- data.table::fread(file.path(cfg$output_dir, "10_objects_and_summary", "module_status.tsv"))
stopifnot(any(status$Status == "FAILED" & grepl("INJECTED_RENDER_FAILURE", status$Detail)), grDevices::dev.cur() == 1L)
cat("COMBINE_EXPORT_MEMORY_PASS\nTEST_ROOT=", root, "\n", sep = "")
