new_legacy_context <- function(files, settings = list(), plot_helper = NULL) {
  check_packages(c("R6", "tidyverse", "SummarizedExperiment", "S4Vectors", "BiocManager"))
  defaults <- list(output_dir = tempfile("legacy_cache_"), tax_id = 9685L, kegg_code = "fca", go_ancestors = TRUE,
                   fdr_cutoff = 0.05, abs_logfc_cutoff = 1, min_gs_size = 10L, max_gs_size = 500L, seed = 20261019L,
                   top_terms = 12L, direction_colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "#BDBDBD"))
  env <- new.env(parent = globalenv())
  env$cfg <- utils::modifyList(defaults, settings)
  env$`%>%` <- magrittr::`%>%`
  by_name <- stats::setNames(files, basename(files))
  sys.source(by_name[["project_helpers.R"]], env)
  sys.source(by_name[["Rscript_R6_transcriptome_base.R"]], env)
  expressions <- parse(by_name[["Rscript_R6_transcriptome_class.R"]], keep.source = FALSE)
  for (expr in expressions) {
    if (is.call(expr) && as.character(expr[[1]])[1] %in% c("<-", "=")) eval(expr, env)
  }
  sys.source(by_name[["compatibility.R"]], env)
  if (!is.null(plot_helper)) sys.source(plot_helper, env)
  env
}

legacy_public_data <- function(annotation) {
  list(get_gene_annotation = function() annotation)
}
