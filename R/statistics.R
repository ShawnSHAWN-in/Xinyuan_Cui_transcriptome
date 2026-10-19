fit_comparison_lst <- function(count_data, job, model, engine_files) {
  if (length(setdiff(names(model), c("method", "strategy", "min_replicates")))) stop("模型含未支持的参数；协变量和交互检验需显式新增模型")
  if (!identical(model$method, "edgeR") || !identical(model$strategy, "legacy_pairwise")) stop("首版回归模型为edgeR legacy_pairwise")
  numbers <- table(job$samples$Group)
  if (any(numbers < model$min_replicates)) stop("生物学重复不足：", job$id)
  if (any(!job$samples$Sample %in% colnames(count_data$counts))) stop("比较样本不在count矩阵：", job$id)
  context <- new_legacy_context(engine_files)
  work <- tempfile("legacy_design_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  sample_file <- write_table_path(data.frame(Samples = job$samples$Sample, Group = job$samples$Group), file.path(work, "samples.tsv"))
  contrast_file <- file.path(work, "contrasts.tsv")
  data.table::fwrite(data.frame(job$numerator, job$denominator), contrast_file, sep = "\t", col.names = FALSE)
  design <- context$object_experiment_design$new(count_matrix = count_data$counts[, job$samples$Sample, drop = FALSE],
                                                gene_length = count_data$lengths, file_GroupSamples = sample_file, file_GroupContrast = contrast_file)
  result <- context$run_calculate_differential_gene(design, contrast_index = 1L, method = "edgeR", out_dir = NULL)
  message("edgeR完成：", job$id)
  list(job = job, table = tibble::as_tibble(result$deg_tbl), method = result$method)
}

attach_gene_annotation_tbl <- function(fit, annotation) {
  cols <- c("Geneid", "SYMBOL", "ENTREZID", "GENENAME")
  anno <- annotation[, cols, drop = FALSE]
  ambiguous <- names(which(table(anno$Geneid) > 1L))
  anno <- anno[!anno$Geneid %in% ambiguous, , drop = FALSE]
  duplicated_ids <- unique(anno$ENTREZID[!is.na(anno$ENTREZID) & (duplicated(anno$ENTREZID) | duplicated(anno$ENTREZID, fromLast = TRUE))])
  anno <- anno[!anno$ENTREZID %in% duplicated_ids, , drop = FALSE]
  tbl <- dplyr::left_join(fit$table, anno, by = "Geneid")
  if (nrow(tbl) && all(is.na(tbl$ENTREZID))) stop("受检Geneid无法获得唯一NCBI GeneID，请核对GTF和定量ID")
  list(job = fit$job, table = tbl, method = fit$method)
}

select_deg_lst <- function(annotated, thresholds) {
  if (thresholds$fdr <= 0 || thresholds$fdr > 1 || thresholds$abs_logfc < 0) stop("DEG阈值无效")
  tbl <- annotated$table
  keep <- is.finite(tbl$padj) & tbl$padj < thresholds$fdr & is.finite(tbl$logFC) & abs(tbl$logFC) >= thresholds$abs_logfc
  list(job = annotated$job, table = tbl[keep, , drop = FALSE], universe = unique(stats::na.omit(as.character(tbl$ENTREZID))))
}

rank_genes_lst <- function(annotated) {
  tbl <- annotated$table
  tbl <- tbl[!is.na(tbl$ENTREZID) & is.finite(tbl$logFC), , drop = FALSE]
  tbl <- tbl[order(-abs(tbl$logFC), tbl$Geneid), , drop = FALSE]
  tbl <- tbl[!duplicated(tbl$ENTREZID), , drop = FALSE]
  list(job = annotated$job, ranking = sort(stats::setNames(tbl$logFC, as.character(tbl$ENTREZID)), decreasing = TRUE))
}

calc_qc_lst <- function(count_data, metadata) {
  check_sample_alignment(metadata, count_data)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = count_data$counts))
  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  cpm <- edgeR::cpm(dge, log = FALSE)
  rpk <- sweep(count_data$counts, 1, count_data$lengths / 1000, "/")
  tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
  variable <- apply(log_cpm, 1, stats::sd) > 0
  if (sum(variable) < 2L) stop("可变表达基因不足，无法进行PCA")
  pca <- stats::prcomp(t(log_cpm[variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
  scores <- tibble::rownames_to_column(as.data.frame(pca$x[, 1:2, drop = FALSE]), "Sample")
  scores <- dplyr::left_join(scores, metadata, by = "Sample")
  list(log_cpm = log_cpm, cpm = cpm, tpm = tpm, scores = scores,
       variance = 100 * pca$sdev^2 / sum(pca$sdev^2), correlation = stats::cor(log_cpm),
       library = tibble::tibble(Sample = colnames(cpm), LibrarySize = colSums(count_data$counts), Detected = colSums(count_data$counts > 0)))
}
