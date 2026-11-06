# 项目路径
.config_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.config_file)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .config_file <- if (length(.script_arg)) sub("^--file=", "", .script_arg[1]) else "analysis_scripts/00_config.R"
}
script_dir  <- dirname(normalizePath(.config_file, mustWork = TRUE))
project_dir <- Sys.getenv("TRANSCRIPTOME_PROJECT_DIR", unset = dirname(script_dir))
project_dir <- normalizePath(project_dir, mustWork = TRUE)
cfg <- list(
  project_dir       = project_dir,
  script_dir        = script_dir,
  count_file        = file.path(project_dir, "02.Quant_定量结果_Calicivirus/result_count.with_Length.tsv"),
  sample_file       = file.path(project_dir, "02.Quant_定量结果_Calicivirus/01_GroupSamples.tsv"),
  contrast_file     = file.path(project_dir, "02.Quant_定量结果_Calicivirus/03_GroupContrast.txt"),
  gtf_file          = file.path(project_dir, "00.database/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.gtf.gz"),
  gene2go_file      = file.path(project_dir, "gene2go.gz"),
  output_dir        = file.path(project_dir, "03_analysis"),
  reference_species = "Felis catus",
  sample_species    = "Panthera leo",
  tax_id            = 9685L,
  kegg_code         = "fca",
  treatment_levels  = c("Cat", "Leopard", "Tiger", "NC"),
  time_levels       = c("06", "12", "24"),
  treatment_colors  = c(Cat = "#3B75AF", Leopard = "#E69F00", Tiger = "#009E73", NC = "#777777"),
  direction_colors  = c(Up = "#C94C4C", Down = "#377EB8", Not = "#C9C9C9"),
  method            = "edgeR",
  fdr_cutoff        = 0.05,
  abs_logfc_cutoff  = 1,
  go_ontology       = "BP",
  go_ancestors      = TRUE,
  min_gs_size       = 10L,
  max_gs_size       = 500L,
  top_terms         = 12L,
  top_genes_each    = 15L,
  compare_top_each  = 5L,
  compare_max_terms = 30L,
  curve_top_each    = 3L,
  pca_scale         = FALSE,
  dpi               = 300L,
  threads           = 1L,
  seed              = 20261019L,
  run_kegg          = TRUE,
  require_kegg      = FALSE,
  run_pathview      = TRUE,
  pathview_top      = 2L,
  theme_go_ids      = character(),
  optional_packages = c("UpSetR", "ggdendro", "ggpubr", "scatterplot3d", "pathview")
)
cfg$group_levels <- as.vector(outer(cfg$treatment_levels, cfg$time_levels, paste, sep = "_"))
cfg$contrast_levels <- unlist(lapply(cfg$time_levels, function(time) {
  paste0(setdiff(cfg$treatment_levels, "NC"), "_", time, "_vs_NC_", time)
}), use.names = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "lib", "project_helpers.R"))

# 统一载入基座
core_packages <- c("R6", "tidyverse", "data.table", "SummarizedExperiment", "S4Vectors", "rtracklayer", "AnnotationDbi", "GO.db",
                   "edgeR", "ggplot2", "patchwork", "pheatmap", "openxlsx", "ggrepel", "R.utils", "digest", "BiocParallel", "ggpubr")
missing_core <- core_packages[!vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_core)) stop("缺少基础包：", paste(missing_core, collapse = ", "), "；参见 README.md 的依赖说明")
if (utils::packageVersion("dplyr") < "1.1.1") stop("需要dplyr >= 1.1.1")
source(file.path(script_dir, "lib", "Rscript_R6_transcriptome_class.R"))
source(file.path(script_dir, "lib", "compatibility.R"))
source(file.path(script_dir, "lib", "plot_helpers.R"))
set.seed(cfg$seed)
data.table::setDTthreads(cfg$threads)
message("配置完成：", cfg$project_dir)
