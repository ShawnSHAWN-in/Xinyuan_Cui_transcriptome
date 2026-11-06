if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
require_enrichment_packages()
inputs <- read_module("01_inputs")
deg <- read_module("03_deg")
check_module_inputs(inputs)
check_analysis_settings(deg)
out_dir <- module_dir("07_enrichment")
settings <- cfg[c("min_gs_size", "max_gs_size", "seed", "go_ancestors", "fdr_cutoff", "method")]
signature <- list(inputs = inputs$signature, deg = unname(tools::md5sum(file.path(cfg$output_dir, "10_objects_and_summary", "03_deg.rds"))),
                   settings = settings, background = digest::digest(inputs$backgrounds))
results <- stats::setNames(vector("list", length(cfg$contrast_levels)), cfg$contrast_levels)
status_rows <- list()
for (comparison in cfg$contrast_levels) {
  message("富集分析：", comparison)
  checkpoint <- file.path(out_dir, paste0(comparison, "_checkpoint.rds"))
  cached <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
  if (!is.null(cached) && identical(cached$signature, signature)) {
    results[[comparison]] <- cached$result
    message("复用已完成富集：", comparison)
    next
  }
  all <- deg$comparisons[[comparison]]$all
  significant <- deg$comparisons[[comparison]]$significant
  universe <- unique(stats::na.omit(as.character(all$ENTREZID)))
  result <- list(GO = NULL, KEGG = NULL)
  for (db in c("GO", "KEGG")) {
    background <- inputs$backgrounds[[db]]
    if (is.null(background)) {
      record_status("06_enrichment", paste(comparison, db), "SKIPPED", "背景不可用；请查看01_annotation状态")
      next
    }
    dir <- file.path(out_dir, if (db == "GO") "GO_BP" else "KEGG", comparison)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    item <- list(ORA = list(Up = NULL, Down = NULL), GSEA = NULL, foreground = list(), universe = character(), ranking = NULL)
    item$universe <- intersect(universe, as.character(background$TERM2GENE$GENE))
    save_table(data.frame(ENTREZID = universe, HasAnnotation = universe %in% item$universe), file.path(dir, "00_tested_background.tsv"))
    for (direction in c("Up", "Down")) {
      genes <- significant[significant$Result_expression == direction & !is.na(significant$ENTREZID), , drop = FALSE]
      item$foreground[[direction]] <- intersect(unique(as.character(genes$ENTREZID)), item$universe)
      if (nrow(genes) && length(item$foreground[[direction]])) {
        ora <- run_enrichment_analysis(genes, inputs$public_data, analysis_type = "ORA", direction = direction, db = db,
                                        background = background, universe = universe)
        item$ORA[direction] <- list(ora$ORA)
      }
      save_table(data.frame(ENTREZID = item$foreground[[direction]]), file.path(dir, paste0("ORA_", direction, "_foreground.tsv")))
      save_table(result_table(item$ORA[[direction]]), file.path(dir, paste0("ORA_", direction, "_all.tsv")))
      save_table(result_table(item$ORA[[direction]], TRUE), file.path(dir, paste0("ORA_", direction, "_significant.tsv")))
    }
    gsea <- run_enrichment_analysis(all, inputs$public_data, analysis_type = "GSEA", db = db, background = background)
    item$GSEA <- gsea$GSEA
    item$ranking <- gsea$ranking
    save_table(data.frame(ENTREZID = as.character(names(item$ranking)), RankMetric = unname(item$ranking)), file.path(dir, "GSEA_ranking.tsv"))
    save_table(result_table(item$GSEA), file.path(dir, "GSEA_all.tsv"))
    save_table(result_table(item$GSEA, TRUE), file.path(dir, "GSEA_significant.tsv"))
    saveRDS(item, file.path(dir, "enrichment_objects.rds"))
    result[db] <- list(item)
    for (direction in c("Up", "Down")) {
      count <- nrow(result_table(item$ORA[[direction]], TRUE))
      record_status("06_enrichment", paste(comparison, db, "ORA", direction), if (count) "SUCCESS" else "NO_SIGNIFICANT_RESULT", paste(count, "terms"))
    }
    count <- nrow(result_table(item$GSEA, TRUE))
    record_status("06_enrichment", paste(comparison, db, "GSEA"), if (count) "SUCCESS" else "NO_SIGNIFICANT_RESULT", paste(count, "terms"))
  }
  results[[comparison]] <- result
  saveRDS(list(signature = signature, result = result), checkpoint)
}
for (comparison in names(results)) {
  for (db in c("GO", "KEGG")) {
    item <- results[[comparison]][[db]]
    if (is.null(item)) next
    status_rows[[paste(comparison, db)]] <- data.frame(Comparison = comparison, Database = if (db == "GO") "GO_BP" else "KEGG",
      Background = length(item$universe), ORA_Up = nrow(result_table(item$ORA$Up, TRUE)), ORA_Down = nrow(result_table(item$ORA$Down, TRUE)),
      GSEA = nrow(result_table(item$GSEA, TRUE)), RankingGenes = length(item$ranking))
  }
}
save_table(dplyr::bind_rows(status_rows), file.path(out_dir, "00_enrichment_summary.tsv"))
save_module(list(results = results, settings = settings, signature = signature), "06_enrichment")
