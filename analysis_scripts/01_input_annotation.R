if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
input_dir <- module_dir("01_inputs")
anno_dir <- module_dir("02_annotation")
input_paths <- c(cfg$count_file, cfg$sample_file, cfg$contrast_file, cfg$gtf_file, cfg$gene2go_file)
signature <- input_signature(input_paths)
for (path in c(cfg$gtf_file, cfg$gene2go_file)) {
  if (grepl("\\.gz$", path) && nzchar(Sys.which("gzip"))) {
    status <- system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE)
    if (status != 0L) stop("压缩文件不完整，请确认下载完成：", path)
  }
}

# 适配既有分组文件
sample_info <- data.table::fread(cfg$sample_file, colClasses = "character", data.table = FALSE)
if (identical(names(sample_info)[1], "Sample")) names(sample_info)[1] <- "Samples"
if (!all(c("Samples", "Group") %in% names(sample_info))) stop("分组表必须包含Samples、Group")
if (anyDuplicated(sample_info$Samples) || anyNA(sample_info[, c("Samples", "Group")])) stop("分组表有重复样本或缺失值")
sample_info <- sample_info |>
  dplyr::mutate(Treatment = sub("_[^_]+$", "", Group), Time = sub("^.*_", "", Group),
                Replicate = as.integer(sub("^.*_", "", Samples))) |>
  restore_factors_tbl() |> dplyr::arrange(Time, Treatment, Replicate)
if (anyNA(sample_info)) stop("样本名或分组不符合Treatment_Time_Replicate规则")
expected_samples <- unlist(lapply(cfg$group_levels, function(group) paste0(group, "_", 1:3)), use.names = FALSE)
if (!setequal(sample_info$Samples, expected_samples)) stop("当前计划要求36个指定样本，实际样本集合不一致")
if (any(as.character(sample_info$Group) != sub("_[0-9]+$", "", sample_info$Samples))) stop("样本名与组合组名不对应")
save_table(sample_info, file.path(input_dir, "01_sample_metadata.tsv"))
sample_compatible <- sample_info[, c("Samples", "Group")]
save_table(sample_compatible, file.path(input_dir, "01_GroupSamples.tsv"))
sample_compatible_path <- file.path(input_dir, "01_GroupSamples.tsv")

# 基座比较文件不含表头
contrasts <- data.table::fread(cfg$contrast_file, header = FALSE, colClasses = "character", data.table = FALSE)
if (ncol(contrasts) != 2L) stop("比较文件必须为两列")
if (identical(as.character(unlist(contrasts[1, ])), c("Treatment", "Control"))) contrasts <- contrasts[-1, , drop = FALSE]
names(contrasts) <- c("GroupExperiment", "GroupControl")
contrasts$Comparison <- paste0(contrasts$GroupExperiment, "_vs_", contrasts$GroupControl)
if (anyDuplicated(contrasts$Comparison) || !setequal(contrasts$Comparison, cfg$contrast_levels)) stop("比较方案与计划的9组比较不一致")
contrasts <- contrasts[match(cfg$contrast_levels, contrasts$Comparison), , drop = FALSE]
if (any(!unlist(contrasts[, 1:2]) %in% cfg$group_levels)) stop("比较中有未知组名")
contrast_compatible_path <- file.path(input_dir, "03_GroupContrast.txt")
data.table::fwrite(contrasts[, 1:2], contrast_compatible_path, sep = "\t", col.names = FALSE)
save_table(contrasts, file.path(input_dir, "03_contrast_metadata.tsv"))

# 保留原始计数和长度
count_tbl <- data.table::fread(cfg$count_file, data.table = FALSE, check.names = FALSE)
if (!all(c("Geneid", "Length") %in% names(count_tbl))) stop("矩阵必须包含Geneid和Length")
if (anyDuplicated(count_tbl$Geneid) || anyNA(count_tbl$Geneid)) stop("Geneid存在重复或缺失，需先核对定量主键")
count_columns <- setdiff(names(count_tbl), c("Geneid", "Length"))
if (!setequal(count_columns, sample_info$Samples)) stop("count矩阵与分组文件的样本集合不同")
counts <- as.matrix(count_tbl[, sample_info$Samples, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- count_tbl$Geneid
if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) stop("count必须是无缺失的非负整数")
if (any(colSums(counts) == 0)) stop("存在总计数为0的样本")
lengths <- stats::setNames(as.numeric(count_tbl$Length), count_tbl$Geneid)
if (any(!is.finite(lengths)) || any(lengths <= 0)) stop("Length必须是正数")
design <- object_experiment_design$new(count_matrix = counts, gene_length = lengths,
                                       file_GroupSamples = sample_compatible_path, file_GroupContrast = contrast_compatible_path)
public_data <- PublicDataManager$new(kegg_organism_code = cfg$kegg_code, cache_dir = file.path(anno_dir, "cache"))
public_data$load_project_annotation(cfg$gtf_file, cfg$gene2go_file)
annotation <- public_data$get_gene_annotation()
annotation$ENTREZID <- as.character(annotation$ENTREZID)
save_table(annotation, file.path(anno_dir, "01_gene_annotation.tsv"))
save_table(public_data$get_go_annotation(), file.path(anno_dir, "02_gene2go_direct.tsv"))
mapping_n <- annotation |> dplyr::group_by(Geneid) |> dplyr::summarise(n_ids = dplyr::n_distinct(ENTREZID, na.rm = TRUE), .groups = "drop")
save_table(mapping_n[mapping_n$n_ids > 1, ], file.path(anno_dir, "03_ambiguous_gene_mapping.tsv"))
reverse_mapping <- annotation |> dplyr::filter(!is.na(ENTREZID)) |> dplyr::group_by(ENTREZID) |>
  dplyr::summarise(n_input_ids = dplyr::n_distinct(Geneid), .groups = "drop")
save_table(reverse_mapping[reverse_mapping$n_input_ids > 1, ], file.path(anno_dir, "03b_ambiguous_entrez_mapping.tsv"))
usable_annotation <- attach_annotation_tbl(data.frame(Geneid = rownames(counts)), annotation)
save_table(usable_annotation, file.path(anno_dir, "03c_unambiguous_count_gene_annotation.tsv"))
unmapped <- data.frame(Geneid = setdiff(rownames(counts), annotation$Geneid[!is.na(annotation$ENTREZID)]))
if (nrow(unmapped) == nrow(counts)) stop("矩阵Geneid与GTF GeneID映射完全不匹配，请核对定量使用的GTF版本与ID字段")
save_table(unmapped, file.path(anno_dir, "04_unmapped_genes.tsv"))
backgrounds <- list(GO = public_data$get_go_background(ontology = "BP", include_ancestors = cfg$go_ancestors), KEGG = NULL)
if (!nrow(backgrounds$GO$TERM2GENE)) stop("GO-BP背景为空，不能继续富集")
if (cfg$run_kegg) {
  backgrounds$KEGG <- tryCatch(public_data$get_kegg_background(), error = function(e) {
    record_status("01_annotation", "KEGG", "FAILED", conditionMessage(e))
    if (cfg$require_kegg) stop(e)
    warning("KEGG背景暂不可用：", conditionMessage(e), "；GO和差异分析仍可继续", call. = FALSE)
    NULL
  })
}
for (db in names(backgrounds)) {
  if (is.null(backgrounds[[db]])) next
  save_table(backgrounds[[db]]$TERM2GENE, file.path(anno_dir, paste0(db, "_TERM2GENE.tsv")))
  save_table(backgrounds[[db]]$TERM2NAME, file.path(anno_dir, paste0(db, "_TERM2NAME.tsv")))
  record_status("01_annotation", db, "SUCCESS", paste(nrow(backgrounds[[db]]$TERM2GENE), "gene-term associations"))
}
coverage <- data.frame(Source = c("NCBI_GeneID_any", "NCBI_GeneID_unambiguous", "GO_BP_unambiguous", "KEGG_unambiguous"),
                       Annotated = c(sum(rownames(counts) %in% annotation$Geneid[!is.na(annotation$ENTREZID)]),
                                     sum(!is.na(usable_annotation$ENTREZID)),
                                     sum(usable_annotation$ENTREZID %in% backgrounds$GO$TERM2GENE$GENE),
                                     if (is.null(backgrounds$KEGG)) NA_integer_ else
                                       sum(usable_annotation$ENTREZID %in% backgrounds$KEGG$TERM2GENE$GENE)))
coverage$Total <- nrow(counts)
coverage$Proportion <- coverage$Annotated / coverage$Total
save_table(coverage, file.path(anno_dir, "05_annotation_coverage.tsv"))
saveRDS(cfg, file.path(input_dir, "00_configuration.rds"))
inputs <- list(design = design, public_data = public_data, sample_info = sample_info, contrasts = contrasts, annotation = annotation,
               backgrounds = backgrounds, coverage = coverage, signature = signature,
               settings = cfg[c("tax_id", "kegg_code", "go_ontology", "go_ancestors", "group_levels", "contrast_levels")])
save_module(inputs, "01_inputs")
record_status("01_inputs", "36_samples_9_contrasts", "SUCCESS", paste(nrow(counts), "genes"))
