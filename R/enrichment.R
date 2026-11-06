enrichment_result_tbl <- function(result) {
  if (is.null(result)) return(tibble::tibble(ID = character(), Description = character(), p.adjust = numeric()))
  tibble::as_tibble(result@result)
}

run_ora_lst <- function(selected, background_job, settings, direction) {
  check_packages("clusterProfiler")
  if (!direction %in% c("Up", "Down")) stop("ORA方向错误")
  background <- background_job$background
  universe <- intersect(selected$universe, as.character(background$TERM2GENE$GENE))
  genes <- unique(as.character(selected$table$ENTREZID[selected$table$Result_expression == direction & !is.na(selected$table$ENTREZID)]))
  genes <- intersect(genes, universe)
  result <- NULL
  status <- background_job$status
  if (status == "SUCCESS") {
    status <- "NO_FOREGROUND"
    if (length(genes)) {
      result <- clusterProfiler::enricher(gene = genes, universe = universe, TERM2GENE = background$TERM2GENE,
        TERM2NAME = background$TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
        minGSSize = settings$min_gs_size, maxGSSize = settings$max_gs_size)
      status <- if (is.null(result) || !nrow(enrichment_result_tbl(result))) "NO_TESTABLE_TERMS" else "SUCCESS"
    }
  }
  list(job = selected$job, database = background_job$database, method = "ORA", direction = direction,
       result = result, table = enrichment_result_tbl(result), foreground = genes, universe = universe, status = status)
}

run_gsea_lst <- function(ranked, background_job, settings) {
  check_packages(c("clusterProfiler", "BiocParallel"))
  result <- NULL
  status <- background_job$status
  bg <- background_job$background
  if (status == "SUCCESS") {
    status <- "NO_TESTABLE_TERMS"
    if (length(ranked$ranking) >= settings$min_gs_size && any(bg$TERM2GENE$GENE %in% names(ranked$ranking))) {
      old_kind <- RNGkind()
      on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
      RNGkind(kind = if (is.null(settings$rng_kind)) "Mersenne-Twister" else settings$rng_kind, normal.kind = "Inversion", sample.kind = "Rejection")
      set.seed(settings$seed)
      result <- clusterProfiler::GSEA(geneList = ranked$ranking, TERM2GENE = bg$TERM2GENE, TERM2NAME = bg$TERM2NAME,
        pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = settings$min_gs_size, maxGSSize = settings$max_gs_size,
        eps = 0, seed = TRUE, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
      if (!is.null(result) && nrow(enrichment_result_tbl(result))) status <- "SUCCESS"
    }
  }
  list(job = ranked$job, database = background_job$database, method = "GSEA", direction = "NES",
       result = result, table = enrichment_result_tbl(result), ranking = ranked$ranking, status = status)
}

collect_enrichment_tbl <- function(result) {
  tbl <- result$table
  if (!nrow(tbl)) return(tibble::tibble())
  tbl$Direction <- if (result$method == "GSEA") ifelse(tbl$NES > 0, "Up", "Down") else result$direction
  tbl$Job <- result$job$id
  tbl$Comparison <- result$job$comparison
  tbl$Scheme <- result$job$scheme
  tbl$Database <- result$database
  tbl$Method <- result$method
  tbl
}

collect_enrichment_catalog_tbl <- function(up, down, gsea) {
  tbl <- dplyr::bind_rows(c(up, down, gsea))
  if (!nrow(tbl)) return(tibble::tibble(Scheme = character(), Database = character(), Method = character(), Job = character(),
    ID = character(), Description = character(), Direction = character(), p.adjust = numeric()))
  tbl
}

select_display_terms_lst <- function(tbl, thresholds, display, jobs, backgrounds = NULL) {
  if (!nrow(tbl)) tbl <- tibble::tibble(Scheme = character(), Database = character(), Method = character(), Job = character(),
                                        ID = character(), Description = character(), Direction = character(), p.adjust = numeric())
  tbl <- tbl[is.finite(tbl$p.adjust) & tbl$p.adjust < thresholds$fdr, , drop = FALSE]
  selected <- list()
  for (scheme in unique(vapply(jobs, `[[`, character(1), "scheme"))) {
    for (db in c("GO", "KEGG")) {
      block <- tbl[tbl$Scheme == scheme & tbl$Database == db, , drop = FALSE]
      if (!nrow(block)) {
        selected[[paste(scheme, db, sep = "__")]] <- list(scheme = scheme, database = db, ids = character(), labels = character(),
          status = if (is.null(backgrounds)) "SUCCESS" else backgrounds[[db]]$status)
        next
      }
      groups <- split(block, interaction(block$Job, block$Method, block$Direction, drop = TRUE))
      candidates <- unique(unlist(lapply(groups, function(x) utils::head(x$ID[order(x$p.adjust, x$ID)], display$compare_top_each))))
      summary <- block[block$ID %in% candidates, , drop = FALSE] |>
        dplyr::group_by(ID) |> dplyr::summarise(Comparisons = dplyr::n_distinct(Job), FDR = min(p.adjust), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(Comparisons), FDR, ID)
      ids <- utils::head(summary$ID, display$compare_max_terms)
      selected[[paste(scheme, db, sep = "__")]] <- list(scheme = scheme, database = db, ids = ids,
        labels = stats::setNames(block$Description[match(ids, block$ID)], ids), status = "SUCCESS")
    }
  }
  selected
}

significant_result <- function(result, fdr, ids = NULL) {
  if (is.null(result)) return(NULL)
  copy <- result
  tbl <- result@result
  keep <- is.finite(tbl$p.adjust) & tbl$p.adjust < fdr
  if (!is.null(ids)) keep <- keep & tbl$ID %in% ids
  copy@result <- tbl[keep, , drop = FALSE]
  if (!nrow(copy@result)) NULL else copy
}
