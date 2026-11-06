# 从项目顶层顺序执行
script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else "analysis_scripts/run_all.R"
}
code_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
source(file.path(code_dir, "00_config.R"))
modules <- c("01_input_annotation.R", "02_sample_qc.R", "03_differential_expression.R", "04_gene_sets.R", "05_expression_heatmaps.R",
             "06_enrichment.R", "07_enrichment_plots.R", "08_comparative_enrichment.R", "09_combine_export.R")
for (module in modules) {
  message("执行模块：", module)
  source(file.path(code_dir, module))
}
