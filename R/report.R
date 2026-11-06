export_run_summary_files <- function(files, jobs, metadata, reference, model, thresholds, enrichment, figures, output_dir) {
  files <- unique(as.character(unlist(files, use.names = FALSE)))
  if (any(!file.exists(files))) stop("汇总发现输出文件缺失")
  dir <- file.path(output_dir, "10_objects_and_summary")
  inventory <- tibble::tibble(Path = files, Format = tools::file_ext(files), Bytes = file.info(files)$size)
  status_files <- files[basename(files) == "status.tsv"]
  statuses <- dplyr::bind_rows(lapply(status_files, function(file) {
    tbl <- data.table::fread(file, data.table = FALSE)
    tbl$File <- file
    tbl
  }))
  paths <- c(write_table_path(inventory, file.path(dir, "output_inventory.tsv")),
             write_table_path(statuses, file.path(dir, "analysis_status.tsv")))
  config <- list(reference = reference, model = model, thresholds = thresholds, enrichment = enrichment, figures = figures)
  paths <- c(paths, write_object_path(config, file.path(dir, "effective_configuration.rds")))
  session_file <- file.path(dir, "sessionInfo.txt")
  capture.output(sessionInfo(), file = session_file)
  summary_file <- file.path(dir, "run_summary.md")
  writeLines(c("# targets分析汇总", "", paste0("样本：", nrow(metadata)), paste0("比较任务：", length(jobs)),
    paste0("统计方法：", model$method, " / ", model$strategy), paste0("GO：", reference$go_ontology),
    paste0("FDR < ", thresholds$fdr, "；|log2FC| >= ", thresholds$abs_logfc),
    paste0("图像文件：", sum(inventory$Format %in% c("pdf", "png"))), "",
    "ORA方向为上/下调DEG集合；GSEA方向为NES正负。", "UNAVAILABLE、DISABLED和空结果见analysis_status.tsv。",
    "PNG可能受像素预算限制而降低dpi；PDF保留画布尺寸，Pathview包含栅格元素。"), summary_file)
  c(paths, session_file, summary_file)
}
