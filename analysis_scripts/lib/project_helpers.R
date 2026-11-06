module_dir <- function(name) {
  path <- file.path(cfg$output_dir, name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

save_table <- function(tbl, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!ncol(tbl)) tbl <- data.frame(Status = character())
  tbl <- as.data.frame(tbl)
  for (name in names(tbl)) if (is.character(tbl[[name]])) tbl[[name]] <- gsub("[\r\n\t]", " ", tbl[[name]])
  data.table::fwrite(tbl, path, sep = "\t", na = "NA", quote = FALSE)
  invisible(path)
}

save_module <- function(object, name) {
  saveRDS(object, file.path(module_dir("10_objects_and_summary"), paste0(name, ".rds")))
  message("保存模块：", name)
  invisible(object)
}

read_module <- function(name) {
  path <- file.path(cfg$output_dir, "10_objects_and_summary", paste0(name, ".rds"))
  if (!file.exists(path)) stop("缺少前置模块：", name, "；请按 README 顺序执行")
  readRDS(path)
}

record_status <- function(module, item, status, detail = "", path = "") {
  file <- file.path(module_dir("10_objects_and_summary"), "module_status.tsv")
  row <- data.frame(Module = module, Item = item, Status = status, Detail = detail, Path = path, Time = as.character(Sys.time()))
  old <- if (file.exists(file)) data.table::fread(file, quote = "", colClasses = "character", data.table = FALSE) else row[0, ]
  old <- old[!(old$Module == module & old$Item == item), , drop = FALSE]
  save_table(rbind(old, row), file)
  invisible(row)
}

restore_factors_tbl <- function(tbl) {
  if ("Group" %in% names(tbl)) tbl$Group <- factor(as.character(tbl$Group), levels = cfg$group_levels)
  if ("Treatment" %in% names(tbl)) tbl$Treatment <- factor(as.character(tbl$Treatment), levels = cfg$treatment_levels)
  if ("Time" %in% names(tbl)) tbl$Time <- factor(as.character(tbl$Time), levels = cfg$time_levels)
  if ("group" %in% names(tbl)) tbl$group <- factor(as.character(tbl$group), levels = cfg$contrast_levels)
  tbl
}

safe_name <- function(x) gsub("[^[:alnum:]_.-]+", "_", x)

result_table <- function(res, significant = FALSE) {
  if (is.null(res)) return(data.frame(ID = character(), Description = character(), pvalue = numeric(), p.adjust = numeric(), qvalue = numeric()))
  tbl <- if (is.data.frame(res)) as.data.frame(res) else as.data.frame(res@result)
  if (significant && nrow(tbl)) tbl <- tbl[is.finite(tbl$p.adjust) & tbl$p.adjust < cfg$fdr_cutoff, , drop = FALSE]
  tbl
}

subset_enrichment <- function(res, ids = NULL, significant = TRUE) {
  if (is.null(res)) return(NULL)
  tbl <- result_table(res, significant)
  if (!is.null(ids)) tbl <- tbl[tbl$ID %in% ids, , drop = FALSE]
  if (!nrow(tbl)) return(NULL)
  copy <- res
  copy@result <- tbl
  copy
}

require_enrichment_packages <- function() {
  packages <- c("clusterProfiler", "enrichplot", "ggupset")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("富集模块缺少包：", paste(missing, collapse = ", "), "；不会自动安装")
  invisible(TRUE)
}

attach_annotation_tbl <- function(tbl, annotation) {
  anno <- annotation[, c("Geneid", "SYMBOL", "ENTREZID", "GENENAME"), drop = FALSE]
  counts <- table(anno$Geneid)
  ambiguous <- names(counts)[counts > 1]
  if (length(ambiguous)) anno <- anno[!anno$Geneid %in% ambiguous, , drop = FALSE]
  duplicated_ids <- unique(anno$ENTREZID[!is.na(anno$ENTREZID) & (duplicated(anno$ENTREZID) | duplicated(anno$ENTREZID, fromLast = TRUE))])
  if (length(duplicated_ids)) anno <- anno[!anno$ENTREZID %in% duplicated_ids, , drop = FALSE]
  tbl <- tbl[, setdiff(names(tbl), c("SYMBOL", "ENTREZID", "GENENAME")), drop = FALSE]
  dplyr::left_join(tbl, anno, by = "Geneid")
}

input_signature <- function(paths) {
  info <- file.info(paths)
  if (any(is.na(info$size))) stop("输入文件不存在：", paste(paths[is.na(info$size)], collapse = ", "))
  data.frame(Path = normalizePath(paths), Size = info$size, Modified = as.character(info$mtime))
}

check_module_inputs <- function(inputs) {
  current <- input_signature(inputs$signature$Path)
  if (!identical(current, inputs$signature)) stop("输入文件已改变，请从01模块重新执行，避免使用旧结果")
  if (!identical(inputs$settings, cfg[c("tax_id", "kegg_code", "go_ontology", "go_ancestors", "group_levels", "contrast_levels")])) {
    stop("输入或注释配置已改变，请重新执行01模块")
  }
  invisible(TRUE)
}

check_analysis_settings <- function(deg) {
  expected <- cfg[c("method", "fdr_cutoff", "abs_logfc_cutoff")]
  if (!identical(deg$settings, expected)) stop("差异分析参数已改变，请重新执行03模块及后续模块")
  invisible(TRUE)
}

check_enrichment_settings <- function(enrichment, inputs) {
  fields <- c("min_gs_size", "max_gs_size", "seed", "go_ancestors", "fdr_cutoff", "method")
  if (!identical(enrichment$settings, cfg[fields]) || !identical(enrichment$signature$inputs, inputs$signature) ||
      !identical(enrichment$signature$background, digest::digest(inputs$backgrounds))) {
    stop("富集参数或背景已改变，请重新执行06模块，不能将旧结果按新配置出图")
  }
  invisible(TRUE)
}

select_representative_ids <- function(results, max_each = cfg$compare_top_each, max_total = cfg$compare_max_terms) {
  if (!nrow(results)) return(character())
  split_tbl <- split(results, results$Comparison)
  ids <- unique(unlist(lapply(split_tbl, function(tbl) {
    tbl <- tbl[order(tbl$p.adjust, tbl$ID), , drop = FALSE]
    utils::head(tbl$ID, max_each)
  }), use.names = FALSE))
  summary <- results[results$ID %in% ids, , drop = FALSE] |>
    dplyr::group_by(ID) |>
    dplyr::summarise(n_comparisons = dplyr::n_distinct(Comparison), min_fdr = min(p.adjust), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_comparisons), min_fdr, ID)
  utils::head(summary$ID, max_total)
}
