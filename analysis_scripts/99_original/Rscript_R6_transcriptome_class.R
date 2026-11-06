if (!require("R6", quietly = TRUE)) {
  # 如果 R6 未安装，则安装
  install.packages("R6")
  # 加载 R6 包
  library(R6)
  message("R6 已安装并加载")
} else {
  message("R6 已安装")
}

if (!require("tidyverse", quietly = TRUE)) {
  # 如果 tidyverse 未安装，则安装
  install.packages("tidyverse")
  # 加载 tidyverse 包
  library(tidyverse)
  message("tidyverse 已安装并加载")
} else {
  message("tidyverse 已安装")
}

# Bioconductor 包单独处理
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_pkgs <- c("SummarizedExperiment", "S4Vectors", "rtracklayer", "AnnotationDbi", "GO.db")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# 加载基础数据操作对象
.base_file <- "Rscript_R6_transcriptome_base.R"
if (!file.exists(.base_file)) {
  .ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(.ofile)) {
    .candidate <- file.path(dirname(normalizePath(.ofile, mustWork = FALSE)),
                            "Rscript_R6_transcriptome_base.R")
    if (file.exists(.candidate)) .base_file <- .candidate
  }
}
if (!file.exists(.base_file)) {
  .script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(.script_arg) > 0) {
    .script_file <- sub("^--file=", "", .script_arg[1])
    .candidate <- file.path(dirname(normalizePath(.script_file, mustWork = FALSE)),
                            "Rscript_R6_transcriptome_base.R")
    if (file.exists(.candidate)) .base_file <- .candidate
  }
}
if (!file.exists(.base_file)) {
  stop("无法找到基础脚本: Rscript_R6_transcriptome_base.R")
}
source(.base_file)

# ==================== run_DEG ====================
run_calculate_differential_gene <- function(design_obj, 
                    contrast_index, 
                    method = c("edgeR", "DESeq2", "limma", "t_test", "wilcox_test"),
                    out_dir = NULL) {
  
  # 参数检查
  method <- match.arg(method)
  if (!inherits(design_obj, "Object_experiment_design")) {
    stop("design_obj 必须是 'Object_experiment_design' 对象")
  }

  # ---------- 内部函数：格式化输出表格结果 ----------
  format_result <- function(res_table, contrast, method) {
    if (!"Geneid" %in% colnames(res_table)) {
      res_table <- tibble::rownames_to_column(res_table, "Geneid")
    }
    if (method == "edgeR") {
      res_table <- res_table %>%
        dplyr::rename(logFC = logFC, p.value = PValue, padj = FDR) %>%
        dplyr::mutate(padj = ifelse(is.na(padj), NA_real_, padj))
    } else if (method == "DESeq2") {
      res_table <- res_table %>%
        dplyr::rename(logFC = log2FoldChange, p.value = pvalue, padj = padj)
    } else if (method == "limma") {
      res_table <- res_table %>%
        dplyr::rename(logFC = logFC, p.value = P.Value, padj = adj.P.Val) %>%
        dplyr::mutate(padj = ifelse(is.na(padj), NA_real_, padj))
    } else if (method == "t_test") {
      res_table <- res_table %>%
        # 注意这里使用 padj = PValue ，使用 p value 硬替代
        dplyr::rename(logFC = logFC, p.value = PValue) %>%
        # dplyr::mutate(padj = NA_real_)
        # padj 直接使用 p.value 的值（因为 t_test 不提供调整后的 p 值）
        dplyr::mutate(padj = p.value)
    } else if (method == "wilcox_test") {
      res_table <- res_table |>
        dplyr::rename(p.value = PValue, padj = FDR)
    }
    
    res_table = dplyr::mutate(res_table, SampleA = contrast[1], SampleB = contrast[2], method = method) |> 
      dplyr::mutate(Result_expression = dplyr::case_when(
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down"
      ))
      
    
    base_cols <- c("Geneid", "logFC", "p.value", "padj", "method", "SampleA", "SampleB", "Result_expression")
    other_cols <- setdiff(colnames(res_table), base_cols)
    res_table <- res_table[, c(base_cols, other_cols)]
    return(tibble::as_tibble(res_table))
  }
  
  # ---------- 内部函数：edgeR 分析 ----------
  run_edgeR <- function(contrast_index, covariant = NULL, out_dir = NULL) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("需要安装 edgeR 包: BiocManager::install('edgeR')")
    }
    contrast <- design_obj$get_group_contrast_c(contrast_index)
    exp_group <- contrast[1]
    ctrl_group <- contrast[2]
    
    samples_exp <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == exp_group, 1]
    samples_ctrl <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == ctrl_group, 1]
    all_samples <- c(samples_ctrl, samples_exp)
    
    cnt <- design_obj$get_assay("counts")[, all_samples, drop = FALSE]
    cnt <- round(cnt)
    
    if (!is.null(covariant)) {
      covar_df <- design_obj$GroupSamples[match(all_samples, design_obj$GroupSamples[, 1]), 
                                          covariant, drop = FALSE]
      if (any(is.na(covar_df))) stop("协变量数据缺失")
      covar_df <- as.data.frame(covar_df)
      colnames(covar_df) <- covariant
    } else {
      covar_df <- NULL
    }
    
    keep <- rowSums(edgeR::cpm(cnt) > 1) >= min(length(samples_exp), length(samples_ctrl)) / 2
    cnt <- cnt[keep, , drop = FALSE]
    if (nrow(cnt) == 0) stop("过滤后无剩余基因")
    
    group <- factor(c(rep(ctrl_group, length(samples_ctrl)),
                      rep(exp_group, length(samples_exp))),
                    levels = c(ctrl_group, exp_group))
    
    if (!is.null(covariant) && length(covariant) > 0) {
      design_data <- cbind(data.frame(group = group), covar_df)
      form_str <- paste("~ group +", paste(covariant, collapse = "+"))
      design <- model.matrix(as.formula(form_str), data = design_data)
    } else {
      design <- model.matrix(~ group)
    }
    
    dge <- edgeR::DGEList(counts = cnt, group = group)
    dge <- edgeR::calcNormFactors(dge)
    
    if (length(samples_ctrl) > 1 && length(samples_exp) > 1) {
      dge <- edgeR::estimateDisp(dge, design)
      fit <- edgeR::glmFit(dge, design)
      coef_col <- grep("^group", colnames(design), value = TRUE)
      if (length(coef_col) == 0) stop("未找到组别系数")
      lrt <- edgeR::glmLRT(fit, coef = coef_col)
      res <- edgeR::topTags(lrt, n = Inf)$table
      cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
      res$logCPM <- rowMeans(cpm_mat[rownames(res), , drop = FALSE])
    } else {
      warning("缺少生物学重复，将忽略协变量并使用 exactTest")
      dge$common.dispersion <- 0.1
      et <- edgeR::exactTest(dge, pair = c(ctrl_group, exp_group))
      res <- edgeR::topTags(et, n = Inf)$table
      cpm_mat <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
      res$logCPM <- rowMeans(cpm_mat[rownames(res), , drop = FALSE])
    }
    
    # 添加read counts
    res_cnt = cnt[rownames(res), c(samples_exp, samples_ctrl)]
    # 添加 cpm
    res_cpm = cpm_mat[rownames(res), c(samples_exp, samples_ctrl), drop = FALSE]
    
    res <- dplyr::mutate(res,
                         mean_Sample_exp = rowMeans(res_cnt[rownames(res), samples_exp, drop = FALSE]),
                         mean_Sample_ctrl = rowMeans(res_cnt[rownames(res), samples_ctrl, drop = FALSE]),
                         mean_norm_Sample_exp = rowMeans(res_cpm[rownames(res), samples_exp, drop = FALSE]),
                         mean_norm_Sample_ctrl = rowMeans(res_cpm[rownames(res), samples_ctrl, drop = FALSE]))
    
    # colnames(res_cnt) = base::paste("read_counts", colnames(res_cnt), sep = "_")
    # res = base::cbind(res, res_cnt)
    # colnames(res_cpm) = base::paste("cpm", colnames(res_cpm), sep = "_")
    # res = base::cbind(res, res_cpm)

    res <- as.data.frame(res) %>%
      tibble::rownames_to_column("Geneid")
    
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      openxlsx::write.xlsx(as.data.frame(cnt),
                           file.path(out_dir, paste0("edgeR_filtered_counts_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
      openxlsx::write.xlsx(as.data.frame(edgeR::cpm(dge, normalized.lib.sizes = TRUE)),
                           file.path(out_dir, paste0("edgeR_normalized_CPM_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
    }
    
    return(format_result(res, contrast, method = "edgeR"))
  }
  
  # ---------- 内部函数：DESeq2 分析 ----------
  run_DESeq2 <- function(contrast_index, out_dir = NULL) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("需要安装 DESeq2 包: BiocManager::install('DESeq2')")
    }
    contrast <- design_obj$get_group_contrast_c(contrast_index)
    exp_group <- contrast[1]
    ctrl_group <- contrast[2]
    
    samples_exp <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == exp_group, 1]
    samples_ctrl <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == ctrl_group, 1]
    all_samples <- c(samples_ctrl, samples_exp)
    
    cnt <- design_obj$get_assay("counts")[, all_samples, drop = FALSE]
    cnt <- round(cnt)
    cnt <- cnt[rowMeans(cnt) > 1, , drop = FALSE]
    if (nrow(cnt) == 0) stop("过滤后无剩余基因")
    
    col_data <- data.frame(row.names = all_samples,
                           group = factor(c(rep(ctrl_group, length(samples_ctrl)),
                                            rep(exp_group, length(samples_exp))),
                                          levels = c(ctrl_group, exp_group)))
    design_formula <- ~ group
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt,
                                          colData = col_data,
                                          design = design_formula)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c("group", exp_group, ctrl_group))
    
    # 添加read counts
    res_cnt = cnt[rownames(res), c(samples_exp, samples_ctrl)]
    # 添加 cnt_norm
    cnt_norm <- DESeq2::counts(dds, normalized = TRUE)
    res_cnt_norm = cnt_norm[rownames(res), c(samples_exp, samples_ctrl), drop = FALSE]
    res <- dplyr::mutate(as.data.frame(res),
                         mean_cnt_Sample_exp = rowMeans(res_cnt[rownames(res), samples_exp, drop = FALSE]),
                         mean_cnt_Sample_ctrl = rowMeans(res_cnt[rownames(res), samples_ctrl, drop = FALSE]),
                         mean_norm_Sample_exp = rowMeans(cnt_norm[rownames(res), samples_exp, drop = FALSE]),
                         mean_norm_Sample_ctrl = rowMeans(cnt_norm[rownames(res), samples_ctrl, drop = FALSE])
    )
    
    # colnames(res_cnt) = base::paste("read_counts", colnames(res_cnt), sep = "_")
    # res = base::cbind(res, res_cnt)
    # colnames(res_cnt_norm) = base::paste("cnt_norm", colnames(res_cnt_norm), sep = "_")
    # res = base::cbind(res, res_cnt_norm)
    
    res <- as.data.frame(res) %>%
      tibble::rownames_to_column("Geneid")
    
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      openxlsx::write.xlsx(as.data.frame(cnt),
                           file.path(out_dir, paste0("DESeq2_filtered_counts_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
      openxlsx::write.xlsx(as.data.frame(cnt_norm),
                           file.path(out_dir, paste0("DESeq2_normalized_counts_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
    }
    
    return(format_result(res, contrast, method = "DESeq2"))
  }
  
  # ---------- 内部函数：limma-voom 分析 ----------
  run_limma <- function(contrast_index, out_dir = NULL) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("需要安装 limma 包: BiocManager::install('limma')")
    }
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("limma-voom 需要 edgeR 包用于 DGEList 和 TMM 标准化: BiocManager::install('edgeR')")
    }
    
    contrast <- design_obj$get_group_contrast_c(contrast_index)
    exp_group <- contrast[1]
    ctrl_group <- contrast[2]
    
    samples_exp <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == exp_group, 1]
    samples_ctrl <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == ctrl_group, 1]
    all_samples <- c(samples_ctrl, samples_exp)
    
    if (length(samples_exp) < 2 || length(samples_ctrl) < 2) {
      stop("limma-voom 需要每组至少 2 个生物学重复")
    }
    
    cnt <- design_obj$get_assay("counts")[, all_samples, drop = FALSE]
    cnt <- round(cnt)
    
    group <- factor(c(rep(ctrl_group, length(samples_ctrl)),
                      rep(exp_group, length(samples_exp))),
                    levels = c(ctrl_group, exp_group))
    
    dge <- edgeR::DGEList(counts = cnt, group = group)
    keep <- edgeR::filterByExpr(dge, group = group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge$counts) == 0) stop("过滤后无剩余基因")
    
    dge <- edgeR::calcNormFactors(dge)
    design <- model.matrix(~ 0 + group)
    design_colnames <- make.names(levels(group), unique = TRUE)
    names(design_colnames) <- levels(group)
    colnames(design) <- design_colnames
    exp_col <- design_colnames[exp_group]
    ctrl_col <- design_colnames[ctrl_group]
    
    v <- limma::voom(dge, design, plot = FALSE)
    contrast_matrix <- limma::makeContrasts(
      contrasts = paste0(exp_col, "-", ctrl_col),
      levels = design
    )
    fit <- limma::lmFit(v, design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "P")
    
    cnt_filtered <- dge$counts[rownames(res), c(samples_exp, samples_ctrl), drop = FALSE]
    norm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
    voom_log_cpm <- v$E
    
    res <- dplyr::mutate(as.data.frame(res),
                         mean_cnt_Sample_exp = rowMeans(cnt_filtered[rownames(res), samples_exp, drop = FALSE]),
                         mean_cnt_Sample_ctrl = rowMeans(cnt_filtered[rownames(res), samples_ctrl, drop = FALSE]),
                         mean_norm_Sample_exp = rowMeans(norm_cpm[rownames(res), samples_exp, drop = FALSE]),
                         mean_norm_Sample_ctrl = rowMeans(norm_cpm[rownames(res), samples_ctrl, drop = FALSE]),
                         mean_voom_logCPM_Sample_exp = rowMeans(voom_log_cpm[rownames(res), samples_exp, drop = FALSE]),
                         mean_voom_logCPM_Sample_ctrl = rowMeans(voom_log_cpm[rownames(res), samples_ctrl, drop = FALSE]))
    
    res <- as.data.frame(res) %>%
      tibble::rownames_to_column("Geneid")
    
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      openxlsx::write.xlsx(as.data.frame(dge$counts),
                           file.path(out_dir, paste0("limma_filtered_counts_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
      openxlsx::write.xlsx(as.data.frame(norm_cpm),
                           file.path(out_dir, paste0("limma_normalized_CPM_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
      openxlsx::write.xlsx(as.data.frame(voom_log_cpm),
                           file.path(out_dir, paste0("limma_voom_logCPM_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
    }
    
    return(format_result(res, contrast, method = "limma"))
  }
  
  # ---------- 内部函数：t 检验分析 ----------
  run_t_test <- function(contrast_index, filter_min_count = TRUE, out_dir = NULL) {
    contrast <- design_obj$get_group_contrast_c(contrast_index)
    exp_group <- contrast[1]
    ctrl_group <- contrast[2]
    
    samples_exp <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == exp_group, 1]
    samples_ctrl <- design_obj$GroupSamples[design_obj$GroupSamples[, 2] == ctrl_group, 1]
    all_samples <- c(samples_ctrl, samples_exp)
    
    cnt <- design_obj$get_assay("counts")[, all_samples, drop = FALSE]
    cnt <- round(cnt)
    if (filter_min_count) cnt <- cnt[rowSums(cnt) >= 1, , drop = FALSE]
    cnt <- cnt[apply(cnt, MARGIN = 1, FUN = \(x){sd(x) != 0}), , drop = FALSE]
    if (nrow(cnt) == 0) stop("过滤后无剩余基因")
    
    res <- data.frame(Geneid = rownames(cnt),
                      mean_SampleA = rowMeans(cnt[, samples_exp, drop = FALSE]),
                      mean_SampleB = rowMeans(cnt[, samples_ctrl, drop = FALSE]),
                      logFC = log2(rowMeans(cnt[, samples_exp, drop = FALSE]) /
                                     rowMeans(cnt[, samples_ctrl, drop = FALSE])),
                      PValue = NA,
                      stringsAsFactors = FALSE)
    
    for (i in seq_len(nrow(cnt))) {
      x <- as.numeric(cnt[i, samples_exp])
      y <- as.numeric(cnt[i, samples_ctrl])
      if (length(x) < 2 || length(y) < 2 || (sd(x) == 0 && sd(y) == 0)) {
        res$PValue[i] <- NA
      } else {
        t_res <- tryCatch(t.test(x, y, var.equal = FALSE), error = function(e) NULL)
        res$PValue[i] <- if (!is.null(t_res)) t_res$p.value else NA
      }
    }
    res <- res[!is.na(res$PValue), ]
    
    # 添加read counts
    res_cnt = cnt[res$Geneid, c(samples_exp, samples_ctrl)]
    colnames(res_cnt) = base::paste("read_counts", colnames(res_cnt), sep = "_")
    res = base::cbind(res, res_cnt)
    
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      openxlsx::write.xlsx(as.data.frame(cnt[res$Geneid, ]),
                           file.path(out_dir, paste0("t_test_counts_",
                                                     exp_group, "_vs_", ctrl_group, ".xlsx")),
                           rowNames = TRUE)
    }
    
    return(format_result(res, contrast, method = "t_test"))
  }
  
  # ---------- 内部函数： Wilcoxon 秩和检验 ----------
  run_wilcox_test <- function(contrast_index, out_dir = NULL) {
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("TMM 标准化需要 edgeR 包: BiocManager::install('edgeR')")
    contrast <- design_obj$get_group_contrast_c(contrast_index)
    exp_group <- contrast[1]
    ctrl_group <- contrast[2]
    if (exp_group == ctrl_group) stop("实验组和对照组不能相同")
    samples_exp <- design_obj$GroupSamples[design_obj$GroupSamples[, "Group"] == exp_group, "Samples"]
    samples_ctrl <- design_obj$GroupSamples[design_obj$GroupSamples[, "Group"] == ctrl_group, "Samples"]
    if (length(samples_exp) < 2 || length(samples_ctrl) < 2)
      stop("Wilcoxon 检验要求每组至少 2 个生物学重复")
    all_samples <- c(samples_ctrl, samples_exp)
    
    # 检查计数矩阵
    cnt <- design_obj$get_assay("counts")[, all_samples, drop = FALSE]
    if (!is.numeric(cnt)) stop("counts 必须是数值矩阵")
    if (any(!is.finite(cnt)) || any(cnt < 0))
      stop("counts 不能包含 NA、非有限值或负数")
    cnt <- round(cnt)
    if (nrow(cnt) == 0) stop("计数矩阵没有基因")
    if (any(colSums(cnt) == 0))
      stop("以下样本总计数为0: ", paste(colnames(cnt)[colSums(cnt) == 0], collapse = ", "))
    
    # TMM 标准化及过滤
    dge <- edgeR::DGEList(counts = cnt)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    keep <- rowSums(cnt) > 0
    cnt <- cnt[keep, , drop = FALSE]
    norm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)[keep, , drop = FALSE]
    
    # 逐基因检验
    p_values <- apply(norm_cpm, 1, function(values) {
      if (all(values == values[1])) return(1)
      stats::wilcox.test(
        x = values[samples_exp], y = values[samples_ctrl],
        alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE
      )$p.value
    })
    mean_exp <- rowMeans(norm_cpm[, samples_exp, drop = FALSE])
    mean_ctrl <- rowMeans(norm_cpm[, samples_ctrl, drop = FALSE])
    res <- tibble::tibble(
      Geneid = rownames(cnt),
      logFC = log2((mean_exp + 1) / (mean_ctrl + 1)),
      PValue = unname(p_values),
      FDR = stats::p.adjust(p_values, method = "BH"),
      mean_SampleA = rowMeans(cnt[, samples_exp, drop = FALSE]),
      mean_SampleB = rowMeans(cnt[, samples_ctrl, drop = FALSE]),
      mean_norm_Sample_exp = mean_exp,
      mean_norm_Sample_ctrl = mean_ctrl,
      median_norm_Sample_exp = apply(norm_cpm[, samples_exp, drop = FALSE], 1, median),
      median_norm_Sample_ctrl = apply(norm_cpm[, samples_ctrl, drop = FALSE], 1, median)
    )
    
    # 添加各样本计数
    res_cnt <- cnt[, c(samples_exp, samples_ctrl), drop = FALSE]
    colnames(res_cnt) <- paste0("read_counts_", colnames(res_cnt))
    res <- cbind(res, res_cnt)
    
    # 保存分析矩阵
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      contrast_name <- paste0(exp_group, "_vs_", ctrl_group)
      openxlsx::write.xlsx(as.data.frame(cnt),
                           file.path(out_dir, paste0("wilcox_test_filtered_counts_", contrast_name, ".xlsx")), rowNames = TRUE)
      openxlsx::write.xlsx(as.data.frame(norm_cpm),
                           file.path(out_dir, paste0("wilcox_test_normalized_CPM_", contrast_name, ".xlsx")), rowNames = TRUE)
    }
    return(format_result(res, contrast, method = "wilcox_test"))
  }
  
  # ---------- 主流程：批量运行 ----------
  ctr <- design_obj$get_group_contrast_c(contrast_index)
  contrast_names <- paste0(ctr[1], "_vs_", ctr[2])
  
  deg_func <- switch(method,
                     edgeR = run_edgeR,
                     DESeq2 = run_DESeq2,
                     limma = run_limma,
                     t_test = run_t_test,
                     wilcox_test = run_wilcox_test)
  
  deg_tbl = deg_func(contrast_index = contrast_index, out_dir = out_dir) |> 
    dplyr::mutate(group = contrast_names)
  
  
  return(list(
    # 组别
    contrast_names = contrast_names,
    # 差异基因表格
    deg_tbl        = deg_tbl,
    # 方法
    method         = method
  ))
}

# ==================== run_filter_zero_gene ====================
# 蛋白质组/蛋白丰度矩阵特殊过滤策略：
#   该函数用于过滤 0 值过多的词条/基因/蛋白。
#   典型场景是蛋白质组数据中存在大量未检出值或 0 值，需要先按检出比例筛选。
#   默认规则为：保留超过 50% 样本中 assay 值 > 0 的条目。
# 注意：
#   该策略不是所有转录组 counts 数据的默认预处理步骤。
#   转录组差异分析通常还需要结合 edgeR::filterByExpr() 等方法做表达过滤。
run_filter_zero_gene <- function(design_obj,
                                 assay_name = "counts",
                                 nonzero_ratio_threshold = 0.5,
                                 out_dir = NULL) {
  
  # 参数检查
  if (!inherits(design_obj, "read_counts")) {
    stop("design_obj 必须是 'read_counts' 或其子类对象")
  }
  if (!is.numeric(nonzero_ratio_threshold) || length(nonzero_ratio_threshold) != 1 ||
      !is.finite(nonzero_ratio_threshold) ||
      nonzero_ratio_threshold < 0 || nonzero_ratio_threshold > 1) {
    stop("nonzero_ratio_threshold 必须是 0 到 1 之间的单个数值")
  }
  
  # 计算每个基因/词条在样本中的非零比例
  assay_mat = design_obj$get_assay(assay_name)
  nonzero_ratio = apply(
    assay_mat,
    MARGIN = 1,
    FUN = function(x) {
      sum(x > 0, na.rm = TRUE) / length(x)
    }
  )
  keep_gene = nonzero_ratio > nonzero_ratio_threshold
  
  if (sum(keep_gene) == 0) {
    stop("0值过滤后没有剩余基因，请降低 nonzero_ratio_threshold")
  }
  
  filter_tbl = data.frame(
    Geneid = names(nonzero_ratio),
    nonzero_ratio = as.numeric(nonzero_ratio),
    keep_gene = as.logical(keep_gene),
    stringsAsFactors = FALSE
  )
  
  filter_stats = data.frame(
    filter_strategy = "proteomics_nonzero_ratio_filter",
    filter_note = "蛋白质组/蛋白丰度矩阵特殊过滤策略：保留超过阈值比例样本中 assay 值 > 0 的条目",
    assay_name = assay_name,
    nonzero_ratio_threshold = nonzero_ratio_threshold,
    gene_number_before = nrow(assay_mat),
    gene_number_keep = sum(keep_gene),
    gene_number_filter = sum(!keep_gene),
    stringsAsFactors = FALSE
  )
  
  message(
    "[", Sys.time(), "] 0值过滤完成: ",
    "过滤前 ", filter_stats$gene_number_before, " 个基因; ",
    "保留 ", filter_stats$gene_number_keep, " 个基因; ",
    "过滤 ", filter_stats$gene_number_filter, " 个基因; ",
    "阈值 nonzero_ratio > ", nonzero_ratio_threshold
  )
  
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    write.table(
      filter_tbl,
      file = file.path(out_dir, "gene_nonzero_ratio.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    write.table(
      filter_stats,
      file = file.path(out_dir, "zero_filter_stats.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  filter_design_obj = design_obj$filter_genes(names(nonzero_ratio)[keep_gene])
  
  return(list(
    # 过滤后的设计对象
    filter_design_obj = filter_design_obj,
    # 过滤策略名称
    filter_strategy = "proteomics_nonzero_ratio_filter",
    # 过滤策略说明
    filter_note = "蛋白质组/蛋白丰度矩阵特殊过滤策略：保留超过阈值比例样本中 assay 值 > 0 的条目",
    # 每个基因/词条的非零比例表
    filter_tbl = filter_tbl,
    # 过滤统计信息
    filter_stats = filter_stats,
    # 使用的表达矩阵名称
    assay_name = assay_name,
    # 非零比例阈值
    nonzero_ratio_threshold = nonzero_ratio_threshold
  ))
}

run_calculate_significant_DEG <- function(deg_list, padj_threshold = NULL, abs_logFC_threshold = NULL, out_dir = NULL){
  different_expression_gene = deg_list$deg_tbl
  program = deg_list$method
  
  
  # 自适应计算 logFC 和 padj 阈值
  calculate_threshold = function(tbl){
    p_col <- if ("padj" %in% colnames(tbl) && any(!is.na(tbl$padj))) "padj" else "p.value"
    adaptive_padj_threshold <- quantile(tbl[[p_col]], probs = 0.05, na.rm = TRUE)
    if(!is.finite(adaptive_padj_threshold) || adaptive_padj_threshold > 0.05) {
      adaptive_padj_threshold = 0.05
    }
    adaptive_logFC_threshold = quantile(abs(tbl$logFC), probs = 0.20, na.rm = TRUE)
    if(!is.finite(adaptive_logFC_threshold)) {
      adaptive_logFC_threshold = 1
    }
    
    return(list(
      logFC_threshold = adaptive_logFC_threshold,
      padj_threshold = adaptive_padj_threshold
    ))
  }
  
  
  if(is.null(padj_threshold) || is.null(abs_logFC_threshold)){
    # 调用自适应计算 logFC 和 padj 阈值
    res = calculate_threshold(tbl = different_expression_gene)
    if(is.null(padj_threshold)){
      padj_threshold = res[["padj_threshold"]]
    }
    if(is.null(abs_logFC_threshold)){
      abs_logFC_threshold = res[["logFC_threshold"]]
    }
  }else{
    message(paste0("[", Sys.time(), "] ", "padj_threshold 设置为", padj_threshold))
    message(paste0("[", Sys.time(), "] ", "abs_logFC_threshold 设置为", abs_logFC_threshold))
  }
  
  # 筛选符合阈值的显著差异基因
  if(program %in% c("t", "t_test")){
    Significant_DEG <- dplyr::filter(
      different_expression_gene, 
      p.value < padj_threshold,
      abs(logFC) > abs_logFC_threshold
    )
  }else{
    Significant_DEG <- dplyr::filter(
      different_expression_gene, 
      padj < padj_threshold, 
      abs(logFC) > abs_logFC_threshold
    )
  }
  
  if(dim(Significant_DEG)[1]==0) {
    warning(paste0("[", Sys.time(), "] ", "无显著差异基因表格"))
    return(NULL)
  }
  
  if(!is.null(out_dir)){
    out_prefix = "significant_differential_gene"
    file_prefix <- stringr::str_glue("{out_dir}/{program}_{out_prefix}_{base::unique(Significant_DEG$SampleA)}_versus_{base::unique(Significant_DEG$SampleB)}")
    DEG_list <- split(x = Significant_DEG, f = Significant_DEG$Result_expression)
    # 上调基因和下调基因分为两个子表
    message(paste0("[", Sys.time(), "] ", "输出Excel文件 ", stringr::str_glue("{file_prefix}.xlsx")))
    openxlsx::write.xlsx(x = DEG_list, file = stringr::str_glue("{file_prefix}.xlsx"))
  }
  
  return(list(
    contrast_names      = deg_list$contrast_names,
    DEG_tbl             = Significant_DEG,
    method              = deg_list$method,
    padj_threshold      = padj_threshold,
    abs_logFC_threshold = abs_logFC_threshold
  ))
}

# ==================== run_deg_visualization ====================
#' 差异表达结果可视化（面向过程）
#'
#' 该函数替代 DEGVisualizer R6 类，直接处理 run_DEG 返回的合并数据框。
#' 自动根据 group 列拆分为多个对比组，并生成相应的可视化图形。
#'
#' @param deg_df run_DEG 返回的数据框，必须包含列：Geneid, logFC, p.value, padj, group。
#'       若缺少 SYMBOL 列，则会通过 public_data 自动补充。
#' @param public_data PublicDataManager 对象，用于获取基因注释（SYMBOL, ENTREZID）。
#' @param type 可视化类型，可选 "volcano", "rank", "ma", "bar", "all"（默认）。
#' @param pvalue_cutoff 显著性阈值（默认 0.05）。
#' @param logFC_cutoff log2 倍数变化阈值（默认 1）。
#' @param out_dir 输出目录，若为 NULL 则不保存图形，否则保存 PDF/PNG。
#' @param ... 其他参数传递给具体绘图函数（如 volcano 的 label_top_n 等）。
#'
#' @return 若 type != "all"，返回对应的 ggplot 或 patchwork 对象；若 type = "all"，返回包含所有图形的列表。
#'
#' @examples
#' \dontrun{
#' deg_all <- run_DEG(design, contrast_range = 1:2, method = "edgeR")
#' run_deg_visualization(deg_all, public, type = "volcano", out_dir = "DEG_plots")
#' }
run_deg_visualization <- function(deg_list,
                                  public_data = NULL,
                                  type = c("volcano", "rank", "ma", "bar"),
                                  label_top_n = 10,
                                  out_dir = NULL) {
  
  type <- match.arg(type)
  
  # ---------- 内部辅助函数 ----------
  .save_plot <- function(plot, file_prefix, width, height) {
    if (is.null(plot)) return(invisible())
    message(sprintf("[%s] 输出图像 %s.pdf", Sys.time(), file_prefix))
    message(sprintf("[%s] 输出图像 %s.png", Sys.time(), file_prefix))
    if (!dir.exists(dirname(file_prefix))) dir.create(dirname(file_prefix), recursive = TRUE)
    ggplot2::ggsave(paste0(file_prefix, ".pdf"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
    ggplot2::ggsave(paste0(file_prefix, ".png"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
  }
  
  # 确保差异表包含 SYMBOL 列
  .add_symbol_column <- function(df) {
    if (!"SYMBOL" %in% colnames(df)) {
      gene_anno <- if (!is.null(public_data)) public_data$get_gene_annotation() else NULL
      if (!is.null(gene_anno) && "Geneid" %in% colnames(df)) {
        df$Geneid <- as.character(df$Geneid)
        gene_anno$Geneid <- as.character(gene_anno$Geneid)
        df <- merge(df, gene_anno[, c("Geneid", "SYMBOL")], by = "Geneid", all.x = TRUE)
        df$SYMBOL <- ifelse(is.na(df$SYMBOL), df$Geneid, df$SYMBOL)
      } else {
        df$SYMBOL <- if ("Geneid" %in% colnames(df)) df$Geneid else rownames(df)
      }
    }
    if (!"padj" %in% colnames(df)) df$padj <- df$p.value
    if (!"Result_expression" %in% colnames(df)) {
      df$Result_expression <- ifelse(df$logFC > 0, "Up", "Down")
    }
    return(df)
  }
  
  # ---------- 1. 火山图（单组）----------
  plot_volcano <- function(df, contrast_name, pvalue_cutoff, logFC_cutoff, label_top_n = label_top_n) {
    df <- df %>%
      dplyr::mutate(
        threshold = dplyr::case_when(
          padj < pvalue_cutoff & logFC > logFC_cutoff ~ "Up",
          padj < pvalue_cutoff & logFC < -logFC_cutoff ~ "Down",
          TRUE ~ "Stale"
        )
      )
    
    top_genes <- df %>%
      dplyr::filter(threshold != "Stale") %>%
      dplyr::arrange(dplyr::desc(-log10(padj))) %>%
      dplyr::pull(SYMBOL)
    top_genes = top_genes[1:label_top_n]
    
    title <- paste("Volcano Plot:", contrast_name)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(padj))) +
      ggplot2::geom_point(alpha = 0.4, size = 3.5, ggplot2::aes(color = threshold)) +
      ggrepel::geom_text_repel(data = dplyr::filter(df, SYMBOL %in% top_genes),
                               ggplot2::aes(label = SYMBOL), size = 3) +
      ggplot2::ylab("-log10(padj)") +
      ggplot2::labs(title = title) +
      ggplot2::scale_color_manual(values = c("Down" = "blue4", "Stale" = "grey", "Up" = "red3")) +
      ggplot2::geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
                          lty = 4, col = "black", lwd = 0.8) +
      ggplot2::geom_hline(yintercept = -log10(pvalue_cutoff),
                          lty = 4, col = "black", lwd = 0.8) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.key.height = ggplot2::unit(0.6, "cm"))
    return(p)
  }
  
  # ---------- 2. 排序图（单组）----------
  plot_rank <- function(df, contrast_name, label_top_n = label_top_n) {
    rank_df <- df %>%
      dplyr::arrange(dplyr::desc(logFC)) %>%
      dplyr::mutate(rank = dplyr::row_number(),
                    color = ifelse(logFC > 0, "Up", "Down"))
    top_genes <- df %>%
      dplyr::mutate(abs_logFC = abs(logFC)) %>%
      dplyr::arrange(dplyr::desc(abs_logFC)) %>%
      dplyr::slice_head(n = label_top_n) %>%
      dplyr::pull(SYMBOL)
    rank_df <- rank_df %>%
      dplyr::mutate(show_label = SYMBOL %in% top_genes)
    
    title <- paste("Rank Plot:", contrast_name)
    p <- ggplot2::ggplot(rank_df, ggplot2::aes(x = rank, y = logFC, color = color)) +
      ggplot2::geom_point(size = 1, alpha = 0.6) +
      ggplot2::scale_color_manual(values = c("Up" = "red3", "Down" = "blue4")) +
      ggrepel::geom_text_repel(data = dplyr::filter(rank_df, show_label),
                               ggplot2::aes(label = SYMBOL), size = 3) +
      ggplot2::labs(title = title, x = "Rank", y = "log2 Fold Change") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")
    return(p)
  }
  
  # ---------- 3. MA 图（单组）----------
  plot_ma <- function(df, contrast_name, label_top_n = label_top_n) {
    if ("baseMean" %in% colnames(df)) {
      baseMean <- df$baseMean
    } else if ("logCPM" %in% colnames(df)) {
      if (is.matrix(df$logCPM) || is.data.frame(df$logCPM)) {
        baseMean <- 2 ^ rowMeans(df$logCPM, na.rm = TRUE)
      } else {
        baseMean <- 2 ^ df$logCPM
      }
    } else if ("mean_SampleA" %in% colnames(df) && "mean_SampleB" %in% colnames(df)) {
      baseMean <- (df$mean_SampleA + df$mean_SampleB) / 2
    } else {
      stop("无法确定平均表达量（baseMean），请确保差异表包含 logCPM、baseMean 或 mean_SampleA/mean_SampleB 列")
    }
    
    ma_df <- data.frame(
      name = df$SYMBOL,
      log2FoldChange = df$logFC,
      padj = df$padj,
      baseMean = baseMean,
      stringsAsFactors = FALSE
    )
    ma_df <- ma_df[is.finite(ma_df$log2FoldChange) & is.finite(ma_df$padj) & is.finite(ma_df$baseMean), ]
    if (nrow(ma_df) == 0) stop("没有有效的基因用于 MA 图")
    
    p <- ggpubr::ggmaplot(
      ma_df,
      fdr = pvalue_cutoff,
      fc = logFC_cutoff,
      genenames = ma_df$name,
      size = 2.5,
      alpha = 0.5,
      palette = c("#E64B35", "#4DBBD5", "grey70"),
      top = label_top_n,
      select.top.method = "padj",
      font.label = c(8, "bold", "black"),
      label.rectangle = TRUE,
      main = contrast_name,
      xlab = "Log2 mean expression",
      ylab = "Log2 fold change"
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "right"
      )
    return(p)
  }
  
  # ---------- 主逻辑 ----------
  diff_tbl <- .add_symbol_column(deg_list$deg_tbl)
  
  # 根据 type 分发
  if (type == "volcano") {
    p <- plot_volcano(diff_tbl, deg_list$contrast_names, pvalue_cutoff = 0.05, logFC_cutoff = 1, label_top_n = label_top_n)
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, paste0("volcano_", deg_list$contrast_names)), width = 14, height = 9)
    }
    return(p)
    
  } else if (type == "rank") {
    p <- plot_rank(diff_tbl, deg_list$contrast_names, label_top_n = label_top_n)
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, paste0("volcano_", deg_list$contrast_names)), width = 14, height = 9)
    }
    return(p)
    
  } else if (type == "ma") {
    p <- plot_ma(diff_tbl, deg_list$contrast_names, label_top_n = label_top_n)
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, paste0("volcano_", deg_list$contrast_names)), width = 14, height = 9)
    }
    return(p)
    
  } else {
    stop("选项错误")
  }
}

run_DEG_visualization <- function(DEG_list,
                                  public_data = NULL,
                                  type = c("bar", "set_venn", "set_flower", "set_upset"),
                                  direction = c("Up", "Down", "both"),
                                  out_dir = NULL,
                                  # 集合图专用参数
                                  venn_label_show = "counts",
                                  flower_angle = 90,
                                  flower_ellipse_col_pal = "Spectral") {
  
  type <- match.arg(type)
  direction <- match.arg(direction)
  
  # ---------- 内部辅助函数 ----------
  # 保存图像
  .save_plot <- function(plot, file_prefix, width, height) {
    if (is.null(plot)) return(invisible())
    message(sprintf("[%s] 输出图像 %s.pdf", Sys.time(), file_prefix))
    message(sprintf("[%s] 输出图像 %s.png", Sys.time(), file_prefix))
    if (!dir.exists(dirname(file_prefix))) dir.create(dirname(file_prefix), recursive = TRUE)
    ggplot2::ggsave(paste0(file_prefix, ".pdf"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
    ggplot2::ggsave(paste0(file_prefix, ".png"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
  }
  
  # 确保差异表包含 SYMBOL 列
  .add_symbol_column <- function(df) {
    if (!"SYMBOL" %in% colnames(df)) {
      gene_anno <- if (!is.null(public_data)) public_data$get_gene_annotation() else NULL
      if (!is.null(gene_anno) && "Geneid" %in% colnames(df)) {
        df$Geneid <- as.character(df$Geneid)
        gene_anno$Geneid <- as.character(gene_anno$Geneid)
        df <- merge(df, gene_anno[, c("Geneid", "SYMBOL")], by = "Geneid", all.x = TRUE)
        df$SYMBOL <- ifelse(is.na(df$SYMBOL), df$Geneid, df$SYMBOL)
      } else {
        df$SYMBOL <- if ("Geneid" %in% colnames(df)) df$Geneid else rownames(df)
      }
    }
    if (!"padj" %in% colnames(df)) df$padj <- df$p.value
    if (!"Result_expression" %in% colnames(df)) {
      df$Result_expression <- ifelse(df$logFC > 0, "Up", "Down")
    }
    return(df)
  }

  # 差异基因数量柱状图（多组）
  plot_DEG_bar <- function(diff_group_tbl) {
    sig_df = diff_group_tbl |> 
      dplyr::group_by(group, Result_expression) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop")
    
    if (nrow(sig_df) == 0) {
      warning("没有显著差异基因，无法绘制柱状图")
      return(invisible(NULL))
    }
    
    sig_df <- sig_df %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(Total = sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(Total), group, Result_expression) %>%
      dplyr::mutate(Result_expression = factor(Result_expression, levels = c("Up", "Down")))
    
    p_vert <- ggplot2::ggplot(sig_df, ggplot2::aes(x = group, y = Count, fill = Result_expression)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.2,
                        position = ggplot2::position_dodge2(preserve = "single")) +
      ggplot2::geom_text(ggplot2::aes(label = Count),
                         position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
                         vjust = -0.5, size = 3, fontface = "bold") +
      ggplot2::scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5"),
                                 name = "表达方向") +
      ggplot2::labs(title = "差异基因数量统计", x = "对比组", y = "基因数量") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title = ggplot2::element_text(size = 12, face = "bold"),
                     legend.title = ggplot2::element_text(face = "bold"),
                     panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8))
    
    ratio_data <- sig_df %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(Total_DEG = sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Ratio = Count / Total_DEG * 100)
    
    p_ratio <- ggplot2::ggplot(ratio_data, ggplot2::aes(x = group, y = Ratio, fill = Result_expression)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.2) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Ratio)),
                         position = ggplot2::position_stack(vjust = 0.5),
                         size = 3, fontface = "bold", color = "white") +
      ggplot2::scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5"),
                                 name = "表达方向") +
      ggplot2::labs(title = "差异基因上下调比例", x = "对比组", y = "百分比 (%)") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title = ggplot2::element_text(size = 12, face = "bold"),
                     legend.title = ggplot2::element_text(face = "bold"),
                     panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8))
    
    combined <- patchwork::wrap_plots(p_vert, p_ratio, ncol = 2, nrow = 2)
    return(combined)
  }
  
  # Venn 图
  plot_venn <- function(gene_sets, title = NULL) {
    gene_sets <- gene_sets[lengths(gene_sets) > 0]
    if (length(gene_sets) == 0) {
      warning("没有可用于绘制韦恩图的基因集合")
      return(NULL)
    }
    p <- venn_plot(
      data = gene_sets,
      title_size = 1,
      label_show = venn_label_show,
      label_size = 0.8,
      border_show = TRUE,
      line_type = "longdash",
      ellipse_shape = "circle",
      sci_fill_color = "Sci_AAAS",
      sci_fill_alpha = 0.65
    )
    if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
    return(p)
  }

  # 根据 direction = "both" 绘制上下调两个图形并拼合
  plot_venn_patchwork <- function(diff_group_tbl, plot_func, title_up = "Up-regulated", title_down = "Down-regulated") {
    gene_sets_up <- dplyr::filter(diff_group_tbl, Result_expression == "Up")
    gene_sets_up <- split(gene_sets_up$Geneid, gene_sets_up$group)
    gene_sets_up <- gene_sets_up[lengths(gene_sets_up) > 0]
    
    gene_sets_down <- dplyr::filter(diff_group_tbl, Result_expression == "Down")
    gene_sets_down <- split(gene_sets_down$Geneid, gene_sets_down$group)
    gene_sets_down <- gene_sets_down[lengths(gene_sets_down) > 0]
    
    
    p_up <- plot_func(gene_sets_up, title = title_up)
    p_down <- plot_func(gene_sets_down, title = title_down)
    
    # 如果某个方向无有效集，则只返回另一个
    if (is.null(p_up) && is.null(p_down)) return(NULL)
    if (is.null(p_up)) return(p_down)
    if (is.null(p_down)) return(p_up)
    
    if (inherits(p_up, "ggplot") && inherits(p_down, "ggplot")) {
      combined <- patchwork::wrap_plots(p_up, p_down, ncol = 2)
      return(combined)
    } else {
      warning("韦恩图对象不是 ggplot，无法拼合")
      return(NULL)
    }
    return(combined)
  }
  
  # Venn 图 - UpSet
  plot_venn_upset <- function(gene_sets, title = NULL) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
      stop("需要安装 UpSetR 包")
    }
    gene_sets <- gene_sets[lengths(gene_sets) > 0]
    if (length(gene_sets) == 0) {
      warning("没有可用于绘制 UpSet 图的基因集合")
      return(NULL)
    }
    sets_num <- length(gene_sets)
    p <- upsetr_plot(
      data = gene_sets,
      sets_num = sets_num,
      keep_order = FALSE,
      order_by = "freq",
      decrease = TRUE,
      mainbar_color = "#006600",
      number_angle = 45,
      matrix_color = "#cc0000",
      point_size = 4.5,
      point_alpha = 0.50,
      line_size = 0.8,
      shade_color = "#cdcdcd",
      shade_alpha = 0.50,
      setsbar_color = "#000066",
      setsnum_size = 6,
      text_scale = 1.2
    )
    if (!is.null(title)) message("UpSet图不支持添加标题，标题已忽略: ", title)
    return(p)
  }

  # Venn 图 - 花朵图- 由于原始 flower_plot 不返回 ggplot，给出警告并返回 NULL
  plot_venn_flower <- function(gene_sets, title = NULL) {
    warning("flower_plot 不返回 ggplot 对象，无法在拼图中使用。建议使用 set_venn 或 set_upset。")
    return(NULL)
  }
  
  # ---------- 主逻辑 ----------
  diff_group_tbl = map(DEG_list, possibly(~ .x$DEG_tbl, otherwise = NULL)) |> 
    dplyr::bind_rows()
  if (nrow(diff_group_tbl) == 0) {
    stop("DEG_list 中没有可用于可视化的显著差异表格")
  }
  diff_group_tbl <- .add_symbol_column(diff_group_tbl)
  
  # 根据 type 分发
  if (type == "bar") {
    p <- plot_DEG_bar(diff_group_tbl)
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, "DEG_bar"), width = 14, height = 9)
    }
    return(p)
    
  } else if (type == "set_venn") {
    if (direction == "both") {
      p <- plot_venn_patchwork(diff_group_tbl, plot_venn, "Up-regulated genes", "Down-regulated genes")
    } else {
      gene_sets <- dplyr::filter(diff_group_tbl, Result_expression == direction)
      gene_sets <- split(gene_sets$Geneid, gene_sets$group)
      gene_sets <- gene_sets[lengths(gene_sets) > 0]
      p <- plot_venn(gene_sets, title = paste(direction, "-regulated genes"))
    }
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, paste0("venn_", direction)), width = 10, height = 8)
    }
    return(p)
    
  } else if (type == "set_upset") {
    if (direction == "both") {
      p <- plot_venn_patchwork(diff_group_tbl, plot_venn_upset, "Up-regulated genes", "Down-regulated genes")
    } else {
      gene_sets <- dplyr::filter(diff_group_tbl, Result_expression == direction)
      gene_sets <- split(gene_sets$Geneid, gene_sets$group)
      gene_sets <- gene_sets[lengths(gene_sets) > 0]
      p <- plot_venn_upset(gene_sets, title = paste(direction, "-regulated genes"))
    }
    if (!is.null(out_dir)) {
      .save_plot(p, file.path(out_dir, paste0("venn_upset_", direction)), width = 12, height = 8)
    }
    return(p)
    
  } else if (type == "set_flower") {
    warning("flower_plot 不返回 ggplot 对象，无法保存或拼图。")
    # 如果需要，可以调用 flower_plot 直接保存文件（不返回图形）
    if (!is.null(out_dir)) {
      if (direction == "both") {
        gene_sets_up <- dplyr::filter(diff_group_tbl, Result_expression == "Up")
        gene_sets_down <- dplyr::filter(diff_group_tbl, Result_expression == "Down")
        df_up <- .gene_sets_to_df(gene_sets_up)
        df_down <- .gene_sets_to_df(gene_sets_down)
        flower_plot(df_up, angle = flower_angle, ellipse_col_pal = flower_ellipse_col_pal,
                    dir_path = out_dir, file_name = "flower_Up.png")
        flower_plot(df_down, angle = flower_angle, ellipse_col_pal = flower_ellipse_col_pal,
                    dir_path = out_dir, file_name = "flower_Down.png")
        message("花朵图已保存至 ", file.path(out_dir, "flower_Up.png"))
        message("花朵图已保存至 ", file.path(out_dir, "flower_Down.png"))
      } else {
        gene_sets <- .extract_gene_sets(diff_group_tbl, direction)
        df <- .gene_sets_to_df(gene_sets)
        flower_plot(df, angle = flower_angle, ellipse_col_pal = flower_ellipse_col_pal,
                    dir_path = out_dir, file_name = paste0("venn_flower_", direction, ".png"))
        message("花朵图已保存至 ", file.path(out_dir, paste0("venn_flower_", direction, ".png")))
      }
    }
    return(invisible(NULL))
  } 
}

# ==================== run_enrichment_analysis ====================

#' 面向过程的富集分析函数（替代 EnrichmentAnalysis R6 类）
#'
#' 该函数支持 ORA（Over-Representation Analysis）和 GSEA（Gene Set Enrichment Analysis），
#' 基于差异表达结果和 PublicDataManager 对象（KEGG/GO 背景数据）进行通路富集分析。
#'
#' @param diff_list 差异表达结果，可以是单个数据框或数据框列表。每个数据框至少应包含
#'        Geneid / SYMBOL / ENTREZID / logFC 等列（若缺少 SYMBOL/ENTREZID，将通过 public_data 补充）。
#' @param public_data PublicDataManager 对象，提供基因注释及 KEGG/GO 背景数据。
#' @param analysis_type 字符向量，指定要执行的分析类型，可选 "ORA"、"GSEA" 或 c("ORA","GSEA")（默认两者都做）。
#' @param direction ORA 分析的方向，可选 "Up"、"Down"、"All"（默认 "All"）。
#' @param db 使用的数据库，可选 "KEGG" 或 "GO"（默认 "KEGG"）。
#' @param pvalue_cutoff 显著性阈值（默认 0.05）。
#' @param qvalue_cutoff ORA 的 q 值阈值（默认 0.2）。
#' @param minGSSize 通路最小基因数（默认 10）。
#' @param maxGSSize 通路最大基因数（默认 500）。
#' @param out_dir 输出目录，若不为 NULL 则保存结果（RDS 和 Excel 文件）；若为 NULL 则仅返回结果对象。
#'
#' @return 返回一个列表，包含 ORA 和/或 GSEA 的结果（每个结果均为 enrichResult/gseaResult 对象或 NULL）。
#'         列表结构：list(ORA = list(...), GSEA = list(...))
#'
#' @examples
#' \dontrun{
#' res <- run_enrichment_analysis(
#'   diff_list = deg_result,          # 差异表或差异表列表
#'   public_data = pub_data,
#'   analysis_type = "ORA",
#'   direction = "Up",
#'   db = "KEGG",
#'   out_dir = "Enrichment/KEGG_ORA_Up"
#' )
#' }
run_enrichment_analysis <- function(
    diff_list,
    public_data,
    analysis_type = c("ORA", "GSEA"),
    direction = c("Up", "Down"),
    db = c("KEGG", "GO"),
    out_dir = NULL
) {
  
  # ------------------------- 1. 参数检查与准备 -------------------------
  analysis_type <- match.arg(analysis_type, several.ok = TRUE)
  direction <- match.arg(direction)
  db <- match.arg(db)
  
  if (missing(diff_list)) stop("必须提供 diff_list")
  if (length(unique(diff_list$group)) != 1) stop("diff_list group 组别存在多个分组")
  if (missing(public_data)) stop("必须提供 public_data (PublicDataManager 对象)")
  if (!inherits(public_data, "PublicDataManager")) {
    stop("public_data 必须是 PublicDataManager 对象")
  }
  
  # ------------------------- 2. 内部辅助函数 -------------------------
  # 2.1 获取背景数据（TERM2GENE, TERM2NAME）
  get_background <- function(db) {
    if (db == "KEGG") {
      bg <- public_data$get_kegg_background()
    } else {
      bg <- public_data$get_go_background(ontology = "BP")
    }
    if (is.null(bg$TERM2GENE) || nrow(bg$TERM2GENE) == 0) {
      warning(db, " 背景数据为空，无法进行富集分析")
      return(NULL)
    }
    return(bg)
  }
  
  # 2.2 确保差异表包含 SYMBOL 和 ENTREZID 列
  prepare_diff <- function(df) {
    df <- as.data.frame(df)
    # 添加 SYMBOL
    if (!"SYMBOL" %in% colnames(df)) {
      gene_anno <- public_data$get_gene_annotation()
      if (!is.null(gene_anno) && "Geneid" %in% colnames(df)) {
        df$Geneid <- as.character(df$Geneid)
        gene_anno$Geneid <- as.character(gene_anno$Geneid)
        df <- merge(df, gene_anno[, c("Geneid", "SYMBOL")], by = "Geneid", all.x = TRUE)
        df$SYMBOL <- ifelse(is.na(df$SYMBOL), df$Geneid, df$SYMBOL)
      } else {
        df$SYMBOL <- if ("Geneid" %in% colnames(df)) df$Geneid else rownames(df)
      }
    }
    
    # 添加 ENTREZID
    if (!"ENTREZID" %in% colnames(df)) {
      gene_anno <- public_data$get_gene_annotation()
      if (!is.null(gene_anno) && "Geneid" %in% colnames(df)) {
        df <- merge(df, gene_anno[, c("Geneid", "ENTREZID")], by = "Geneid", all.x = TRUE)
      } else {
        df$ENTREZID <- NA_character_
      }
    }
    # 添加上下调方向列（基于 logFC）
    if (!"Result_expression" %in% colnames(df) && "logFC" %in% colnames(df)) {
      df$Result_expression <- ifelse(df$logFC > 0, "Up", "Down")
    }
    return(df)
  }
  
  # 2.3 规范化差异列表（确保为带名称的列表，且每个表都含有必要列）
  normalize_diff <- function(diff_list) {
    if (!is.data.frame(diff_list)) stop("diff_list 必须是数据框")
    if (length(diff_list) == 0) stop("差异列表为空")
    diff_list = prepare_diff(diff_list)
    return(diff_list)
  }
  
  # 2.3 保存结果到文件（若 out_dir 非空）
  save_result <- function(res, file_prefix) {
    if (!is.null(res) && inherits(res, "enrichResult") && nrow(res@result) > 0) {
      dir.create(dirname(file_prefix), recursive = TRUE, showWarnings = FALSE)
      saveRDS(res, paste0(file_prefix, ".RDS"))
      openxlsx::write.xlsx(as.data.frame(res), paste0(file_prefix, ".xlsx"), rowNames = FALSE)
      message("已保存 ", file_prefix, ".{RDS,xlsx}")
    }
  }
  
  # ------------------------- 3. ORA 分析 -------------------------
  run_ora <- function(df, bg, contrast_name_label) {
    # 提取基因
    if (direction == "All") {
      genes <- df$ENTREZID
    } else {
      genes <- df$ENTREZID[df$Result_expression == direction]
    }
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0) {
      message(contrast_name, ": 方向 ", direction, " 无有效 ENTREZID 基因，ORA 跳过")
      return(NULL)
    }
    res <- clusterProfiler::enricher(
      gene = genes,
      TERM2GENE = bg$TERM2GENE,
      TERM2NAME = bg$TERM2NAME
    )
    if (!is.null(res) && nrow(res@result) > 0) {
      if (!is.null(out_dir)) {
        file_prefix <- file.path(out_dir, paste0(db, "_ORA_", contrast_name, "_", direction))
        save_result(res, file_prefix)
      }
      return(res)
    } else {
      message(contrast_name, ": ", db, " ORA 无显著结果 (", direction, ")")
      return(NULL)
    }
  }
  
  # ------------------------- 4. GSEA 分析 -------------------------
  run_gsea <- function(df, bg, contrast_name_label) {
    # 构建排序基因列表 (logFC)
    ranking_tbl <- df[, c("ENTREZID", "logFC")]
    ranking_tbl <- ranking_tbl[!is.na(ranking_tbl$ENTREZID), ]
    ranking_tbl <- ranking_tbl[!duplicated(ranking_tbl$ENTREZID), ]
    if (nrow(ranking_tbl) == 0) {
      message(contrast_name, ": 无有效 ENTREZID 基因，GSEA 跳过")
      return(NULL)
    }
    gene_ranking <- setNames(ranking_tbl$logFC, ranking_tbl$ENTREZID)
    gene_ranking <- sort(gene_ranking, decreasing = TRUE)
    
    res <- clusterProfiler::GSEA(
      geneList = gene_ranking,
      TERM2GENE = bg$TERM2GENE,
      TERM2NAME = bg$TERM2NAME,
      eps = 0,
    )
    if (!is.null(res) && nrow(res@result) > 0) {
      if (!is.null(out_dir)) {
        file_prefix <- file.path(out_dir, paste0(db, "_GSEA_", contrast_name))
        save_result(res, file_prefix)
      }
      return(res)
    } else {
      message(contrast_name, ": ", db, " GSEA 无显著结果")
      return(NULL)
    }
  }
  
  # ------------------------- 5. 主流程 -------------------------
  # 规范化差异列表
  diff_list <- normalize_diff(diff_list)
  
  # 获取背景数据
  bg_data <- get_background(db)
  if (is.null(bg_data)) return(NULL)
  
  # contrast_name 需要将组别信息标注到图像里面，便于区分
  unique(diff_list$group)
  
  
  results <- list(contrast_name = unique(diff_list$group), ORA = NULL, GSEA = NULL)
  
  # 执行 ORA
  if ("ORA" %in% analysis_type) {
    ora_res <- run_ora(diff_list, bg_data, unique(diff_list$group))
    results$ORA <- ora_res
  }
  
  # 执行 GSEA
  if ("GSEA" %in% analysis_type) {
    gsea_res <- run_gsea(diff_list, bg_data, unique(diff_list$group))
    results$GSEA <- gsea_res
  }
  
  return(results)
}


# ==================== 富集结果可视化函数集 ====================
# 将 enrichResult/gseaResult 中的 ENTREZID 转换为 SYMBOL
convert_result_to_symbol <- function(res, public_data) {
  if (is.null(res)) return(NULL)
  if (!inherits(res, "enrichResult") && !inherits(res, "gseaResult")) return(res)
  
  gene_anno <- public_data$get_gene_annotation()
  id2sym <- setNames(gene_anno$SYMBOL, gene_anno$ENTREZID)
  id2sym <- id2sym[!is.na(names(id2sym)) & !duplicated(names(id2sym))]
  
  convert_geneid <- function(geneid_str) {
    if (is.na(geneid_str) || geneid_str == "" || geneid_str == "NA") return(NA_character_)
    ids <- strsplit(as.character(geneid_str), "/")[[1]]
    syms <- id2sym[ids]
    syms <- syms[!is.na(syms)]
    if (length(syms) == 0) return(NA_character_)
    paste(syms, collapse = "/")
  }

  # enrich 分析结果
  if (nrow(res@result) > 0 & class(res) == "enrichResult") {
    new_geneID <- sapply(res@result$geneID, convert_geneid)
    res@result$geneID <- new_geneID
    
    new_geneID_2 <- sapply(names(res@gene), convert_geneid)
    names(res@gene) = new_geneID_2
  }

  # gsea 分析结果
  if (length(res@geneSets) > 0 & class(res) == "gseaResult") {
    new_geneSets <- lapply(res@result$core_enrichment, convert_geneid)
    res@result$core_enrichment <- new_geneSets
    
    new_geneID_2 <- sapply(names(res@geneList), convert_geneid)
    names(res@geneList) = new_geneID_2
  }
  return(res)
}

# 保存 ggplot 对象为 PDF 和 PNG
save_plot <- function(plot, file_prefix, width, height) {
  if (is.null(plot)) return(invisible())
  message(sprintf("[%s] 输出图像 %s.pdf", Sys.time(), file_prefix))
  message(sprintf("[%s] 输出图像 %s.png", Sys.time(), file_prefix))
  if (!dir.exists(dirname(file_prefix))) dir.create(dirname(file_prefix), recursive = TRUE)
  ggplot2::ggsave(paste0(file_prefix, ".pdf"), plot, width = width, height = height,
                  limitsize = FALSE, units = "in")
  ggplot2::ggsave(paste0(file_prefix, ".png"), plot, width = width, height = height,
                  limitsize = FALSE, units = "in")
}

# 计算图形尺寸（基于通路数量）
calc_size <- function(n_terms, max_label_len = 20) {
  width <- min(12, max(6, n_terms * 0.3 + 3))
  height <- min(12, max(5, n_terms * 0.25 + 2))
  if (!is.null(max_label_len) && max_label_len > 30) width <- width + 2
  return(list(width = width, height = height))
}

# ========================= ORA 可视化-- 网络图 ===========================
# plot_ora_cnet <- function(enrich_result, public_data, showCategory = 5, circular = FALSE,
#                           node_label = "all", out_dir = NULL,
#                           file_name = NULL) {
#   res <- convert_result_to_symbol(enrich_result, public_data)
#   p <- enrichplot::cnetplot(res, showCategory = showCategory, circular = circular,
#                             node_label = node_label, colorEdge = TRUE,
#                             foldChange = foldChange)
#   
#   if (!is.null(out_dir)) {
#     if (is.null(file_name)) file_name <- "ORA_cnet"
#     file_prefix <- file.path(out_dir, file_name)
#     private$save_plot(p, file_prefix, width = 10, height = 8)
#   }
#   return(p)
# }

# 需要重新修改，示例
# enrichplot::cnetplot(
#   convert_result_to_symbol(res$ORA, transcriptome_obj$public_data), 
#   showCategory = 10, circular = F,
#   node_label = "all", colorEdge = TRUE,
#   foldChange = setNames(sig_res$DEG_tbl$logFC, sig_res$DEG_tbl$Geneid))


plot_ora_dot <- function(enrich_result, public_data, showCategory = 10, title = NULL,
                         out_dir = NULL, file_name = NULL) {
  res <- convert_result_to_symbol(enrich_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  if (is.null(title)) title <- "ORA Dotplot"
  p <- enrichplot::dotplot(res, showCategory = showCategory, title = title)
  
  if (!is.null(out_dir)) {
    if (is.null(file_name)) file_name <- "ORA_dot"
    file_prefix <- file.path(out_dir, file_name)
    n_terms <- min(nrow(res@result), showCategory)
    sz <- calc_size(n_terms)
    save_plot(p, file_prefix, width = sz$width, height = sz$height)
  }
  return(p)
}

plot_ora_geom_col <- function(enrich_result, public_data, showCategory = 5, 
                              direction = c("Up", "Down"), title = NULL,
                              out_dir = NULL, file_name = NULL) {
  direction <- match.arg(direction)
  res <- convert_result_to_symbol(enrich_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  
  df <- as.data.frame(res)
  df <- df[order(df$p.adjust), ][1:min(showCategory, nrow(df)), ]
  if (nrow(df) == 0) {
    warning("No significant terms to plot")
    return(invisible(NULL))
  }
  df$Description <- factor(df$Description, levels = rev(df$Description))
  fill_color <- if (direction == "Up") "#E64B35" else "#4DBBD5"
    if (is.null(title)) title <- paste("Top", showCategory, "Terms")
  p <- ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
    geom_col(fill = fill_color) +
    labs(title = title, x = "-log10(p.adjust)", y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  if (!is.null(out_dir)) {
    if (is.null(file_name)) file_name <- "ORA_geom_col"
    file_prefix <- file.path(out_dir, file_name)
    sz <- calc_size(nrow(df))
    save_plot(p, file_prefix, width = sz$width, height = sz$height)
  }
  return(p)
}

plot_ora_upset <- function(enrich_result, public_data, n = 10, title = NULL,
                           out_dir = NULL, file_name = NULL) {
  res <- convert_result_to_symbol(enrich_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  if (is.null(title)) title <- "ORA UpSet Plot"
  p <- enrichplot::upsetplot(res, n = n) + labs(title = title)
  
  if (!is.null(out_dir)) {
    if (is.null(file_name)) file_name <- "ORA_upset"
    file_prefix <- file.path(out_dir, file_name)
    save_plot(p, file_prefix, width = 12, height = 7)
  }
  return(p)
}

plot_ora_pathview <- function(enrich_result, public_data, logFC, species,
                              out_dir = NULL, kegg.dir = "./") {
  if (!requireNamespace("pathview", quietly = TRUE)) {
    stop("需要安装 pathview 包: BiocManager::install('pathview')")
  }
  res <- convert_result_to_symbol(enrich_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  
  # 提取显著通路 ID
  pathway_ids <- rownames(res@result)[res@result$p.adjust < 0.05]
  if (length(pathway_ids) == 0) {
    warning("No significant pathways found")
    return(invisible(NULL))
  }
  
  # 准备基因数据（需要 ENTREZID）
  gene_anno <- public_data$get_gene_annotation()
  id_map <- setNames(gene_anno$ENTREZID, gene_anno$SYMBOL)
  if (is.null(names(logFC))) stop("logFC 向量必须含有名称（SYMBOL）")
  entrez_vec <- logFC
  if (all(names(logFC) %in% gene_anno$SYMBOL)) {
    names(entrez_vec) <- id_map[names(logFC)]
    entrez_vec <- entrez_vec[!is.na(names(entrez_vec))]
  } else {
    warning("logFC 名称不是 SYMBOL，尝试直接用作 ENTREZID")
  }
  
  for (pid in pathway_ids) {
    pv <- tryCatch({
      pathview::pathview(gene.data = entrez_vec,
                         pathway.id = pid,
                         species = species,
                         kegg.dir = kegg.dir,
                         limit = list(gene = max(abs(entrez_vec), na.rm = TRUE)),
                         mid = list(gene = "yellow"),
                         same.layer = FALSE)
    }, error = function(e) NULL)
    if (!is.null(pv) && !is.null(out_dir)) {
      src_png <- file.path(kegg.dir, paste0(pid, ".pathview.png"))
      if (file.exists(src_png)) {
        dest_dir <- file.path(out_dir, "pathview")
        if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
        file.rename(src_png, file.path(dest_dir, paste0(pid, ".pathview.png")))
        message("已保存 ", file.path(dest_dir, paste0(pid, ".pathview.png")))
      }
    }
  }
}

# ========================= GSEA 可视化 ===========================

plot_gsea_double_bar <- function(gsea_result,
                                 public_data = NULL,
                                 contrast_name = NULL,
                                 pvalue_col = "pvalue",
                                 pvalue_signif = 0.001,
                                 color_up = "#36638a",
                                 color_down = "#7bcd7b",
                                 color_not = "#cccccc",
                                 out_dir = NULL,
                                 width = NULL,
                                 height = NULL) {
  if (!class(gsea_result) == "gseaResult" || length(gsea_result) == 0) {
    stop("gsea_result 必须是非空列表或 gseaResult 对象")
  }
  
  res <- gsea_result
  if (is.null(res) || !inherits(res, "gseaResult")) {
    warning("对比组 ", contrast_name, " 不是有效的 gseaResult 对象，跳过")
  }
  # 提取数据框
  df <- as.data.frame(res)
  required_cols <- c("Description", "NES", pvalue_col)
  if (!all(required_cols %in% colnames(df))) {
    warning("对比组 ", contrast_name, " 缺少所需列，跳过")
    next
  }
  
  # 构建标题
  title <- paste("GSEA enrichment:", contrast_name)
  
  # 准备数据
  df$lg_p <- -log10(df[[pvalue_col]]) * sign(df$NES)
  df$threshold <- ifelse(
    df[[pvalue_col]] < pvalue_signif,
    ifelse(df$NES > 0, "Up", "Down"),
    "Not"
  )
  df <- df[order(df$lg_p), ]
  df$Description <- factor(df$Description, levels = df$Description)
  
  # 绘图
  p <- ggplot(df, aes(x = Description, y = lg_p, fill = threshold)) +
    geom_col() +
    geom_hline(yintercept = c(-log10(pvalue_signif), log10(pvalue_signif)),
               color = "white", linewidth = 0.5, linetype = "dashed") +
    geom_text(data = subset(df, threshold == "Up"),
              aes(x = Description, y = -0.1, label = Description),
              hjust = 1, size = 3, color = "#000000") +
    geom_text(data = subset(df, threshold == "Not" & NES > 0),
              aes(x = Description, y = -0.1, label = Description),
              hjust = 1, size = 3, color = color_not) +
    geom_text(data = subset(df, threshold == "Not" & NES < 0),
              aes(x = Description, y = 0.1, label = Description),
              hjust = 0, size = 3, color = color_not) +
    geom_text(data = subset(df, threshold == "Down"),
              aes(x = Description, y = 0.1, label = Description),
              hjust = 0, size = 3, color = "#000000") +
    scale_y_continuous(breaks = scales::extended_breaks(n = 10)) +
    coord_flip() +
    xlab("") + ylab("-log10(pvalue) of GSEA") +
    labs(title = title) +
    guides(fill = "none") +
    scale_fill_manual(values = c("Up" = color_up, "Not" = color_not, "Down" = color_down)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 1)
    )
  
  # 保存图像
  if(!is.null(out_dir)){
    contrast_dir <- file.path(out_dir, contrast_name)
    if (!dir.exists(contrast_dir)) dir.create(contrast_dir, recursive = TRUE)
    file_prefix <- file.path(contrast_dir, "GSEA_double_bar")
    # 动态计算宽高
    if (is.null(width)) {
      max_desc_len <- max(nchar(as.character(df$Description)), na.rm = TRUE)
      w <- max(6, max_desc_len * 0.1 + 4)
    } else {
      w <- width
    }
    if (is.null(height)) {
      h <- max(4, nrow(df) * 0.2 + 1)
    } else {
      h <- height
    }
    save_plot(p, file_prefix, width = w, height = h)
    message("已保存对比组 ", nm, " 的 GSEA 双向柱状图")
  }
  
  
  return(p)
}

plot_gsea_gseaplot2 <- function(gsea_result,
                                public_data = NULL,
                                pvalue_cutoff = 0.05,
                                color = "firebrick",
                                rel_heights = c(1, 0.2, 0.6),
                                out_dir = NULL,
                                width = 6.45,
                                height = 4.52) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("需要安装 enrichplot 包: BiocManager::install('enrichplot')")
  }
  
  # 辅助函数：处理单个 gseaResult 对象
  process_single <- function(res, contrast_name) {
    if (is.null(res) || !inherits(res, "gseaResult")) {
      warning(contrast_name, " 不是有效的 gseaResult 对象，跳过")
      return(NULL)
    }
    # 显著通路索引
    sig_idx <- which(res@result$p.adjust < pvalue_cutoff & !is.na(res@result$NES))
    if (length(sig_idx) == 0) {
      message(contrast_name, " 没有显著富集的通路")
      return(NULL)
    }
    plots <- list()
    for (idx in sig_idx) {
      pathway_id <- rownames(res@result)[idx]
      pathway_desc <- res@result$Description[idx]
      title <- paste0(pathway_id, " ", pathway_desc)
      p <- tryCatch({
        enrichplot::gseaplot2(res, geneSetID = idx, color = color,
                              title = title, rel_heights = rel_heights)
      }, error = function(e) {
        warning("绘制 ", contrast_name, " 通路 ", pathway_id, " 失败: ", e$message)
        return(NULL)
      })
      if (!is.null(p)) {
        plots[[pathway_id]] <- p
        if (!is.null(out_dir)) {
          # 创建对比组子目录
          sub_dir <- file.path(out_dir, contrast_name)
          if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
          # 文件名去除特殊字符
          safe_id <- gsub("/", "_", pathway_id)
          file_prefix <- file.path(sub_dir, paste0("gseaplot2_", safe_id))
          save_plot(p, file_prefix, width = width, height = height)
        }
      }
    }
    return(plots)
  }
  
  # 判断输入类型
  if (inherits(gsea_result, "gseaResult")) {
    # 单个对象
    plots <- process_single(gsea_result, "Default")
    if (is.null(out_dir)) return(plots)
    invisible(plots)
  } else if (is.list(gsea_result)) {
    # 列表（可能来自 run_gsea）
    if (length(gsea_result) == 0) {
      warning("输入列表为空")
      return(invisible(NULL))
    }
    if (is.null(names(gsea_result))) {
      names(gsea_result) <- paste0("contrast_", seq_along(gsea_result))
    }
    all_plots <- list()
    for (nm in names(gsea_result)) {
      res <- gsea_result[[nm]]
      all_plots[[nm]] <- process_single(res, nm)
    }
    if (is.null(out_dir)) return(all_plots)
    invisible(all_plots)
  } else {
    stop("gsea_result 必须是 gseaResult 对象或由 enrich$run_gsea 返回的列表")
  }
}



plot_gsea_cnet <- function(gsea_result, public_data, contrast_name = NULL, showCategory = 5, circular = FALSE,
                           node_label = "all", out_dir = NULL,
                           file_name = NULL) {
  res <- convert_result_to_symbol(gsea_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  
  # 构建标题
  title <- paste("GSEA enrichment :", contrast_name)
  
  p <- enrichplot::cnetplot(res, showCategory = showCategory, circular = circular,
                            node_label = node_label, colorEdge = TRUE,
                            foldChange = res@geneList) + 
    ggplot2::ggtitle(label = title)
  
  if (!is.null(out_dir)) {
    if (is.null(file_name)) file_name <- "GSEA_cnet"
    file_prefix <- file.path(out_dir, file_name)
    save_plot(p, file_prefix, width = 10, height = 8)
  }
  return(p)
}

plot_gsea_pathview <- function(gsea_result, public_data, logFC, species,
                               out_dir = NULL, kegg.dir = "./") {
  if (!requireNamespace("pathview", quietly = TRUE)) {
    stop("需要安装 pathview 包: BiocManager::install('pathview')")
  }
  res <- convert_result_to_symbol(gsea_result, public_data)
  if (is.null(res)) return(invisible(NULL))
  
  pathway_ids <- rownames(res@result)[res@result$p.adjust < 0.05]
  if (length(pathway_ids) == 0) {
    warning("No significant pathways found")
    return(invisible(NULL))
  }
  
  gene_anno <- public_data$get_gene_annotation()
  id_map <- setNames(gene_anno$ENTREZID, gene_anno$SYMBOL)
  if (is.null(names(logFC))) stop("logFC 向量必须含有名称（SYMBOL）")
  entrez_vec <- logFC
  if (all(names(logFC) %in% gene_anno$SYMBOL)) {
    names(entrez_vec) <- id_map[names(logFC)]
    entrez_vec <- entrez_vec[!is.na(names(entrez_vec))]
  } else {
    warning("logFC 名称不是 SYMBOL，尝试直接用作 ENTREZID")
  }
  
  for (pid in pathway_ids) {
    pv <- tryCatch({
      pathview::pathview(gene.data = entrez_vec,
                         pathway.id = pid,
                         species = species,
                         kegg.dir = kegg.dir,
                         limit = list(gene = max(abs(entrez_vec), na.rm = TRUE)),
                         mid = list(gene = "yellow"),
                         same.layer = FALSE)
    }, error = function(e) NULL)
    if (!is.null(pv) && !is.null(out_dir)) {
      src_png <- file.path(kegg.dir, paste0(pid, ".pathview.png"))
      if (file.exists(src_png)) {
        dest_dir <- file.path(out_dir, "pathview")
        if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
        file.rename(src_png, file.path(dest_dir, paste0(pid, ".pathview.png")))
        message("已保存 ", file.path(dest_dir, paste0(pid, ".pathview.png")))
      }
    }
  }
}

# ====================  辅助构造函数 ====================
create_from_files <- function(
    counts_file,            # 表达矩阵文件（含/不含 Length 列）
    sample_file,            # 样本分组文件（Samples, Group）
    contrast_file = NULL,   # 对比文件，NULL 则使用构造函数默认路径
    kegg_organism_code = NULL,
    gtf_file_ncbi = NULL,
    gtf_file_ensembl = NULL,
    use_orgdb = TRUE,
    gene2go.gz = NULL,
    gene2ensembl.gz = NULL,
    eggnog = NULL,
    cache_dir = NULL,
    force_update = FALSE,
    annotate_rowData = FALSE,
    ...
) {

  # ---------- 1. 读取表达矩阵 ----------
  message("[", Sys.time(), "] 读取表达矩阵: ", counts_file)
  mat <- data.table::fread(counts_file, data.table = FALSE)
  rownames(mat) <- mat[[1]]
  mat <- mat[, -1, drop = FALSE]
  
  # 提取 Length 列（如果有）
  if ("Length" %in% colnames(mat)) {
    gene_length <- setNames(as.numeric(mat$Length), rownames(mat))
    mat <- mat[, colnames(mat) != "Length", drop = FALSE]
  } else {
    gene_length <- rep(1, nrow(mat))
    names(gene_length) <- rownames(mat)
    message("未检测到 Length 列，基因长度默认设为 1")
  }
  
  count_matrix <- as.matrix(mat)
  mode(count_matrix) <- "numeric"
  
  # ---------- 2. 读取设计文件并验证样本一致性 ----------
  message("[", Sys.time(), "] 读取样本分组: ", sample_file)
  design_raw <- read.table(sample_file, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, row.names = NULL)
  if (!all(c("Samples", "Group") %in% colnames(design_raw))) {
    stop("样本分组文件必须包含 Samples 和 Group 列")
  }
  design_samples <- design_raw$Samples
  missing_samples <- setdiff(design_samples, colnames(count_matrix))
  if (length(missing_samples) > 0) {
    stop("设计文件中的以下样本在表达矩阵中不存在: ",
         paste(missing_samples, collapse = ", "))
  }
  # 按设计文件中的样本顺序重排矩阵
  count_matrix <- count_matrix[, design_samples, drop = FALSE]
  
  # ---------- 3. 创建设计对象（继承自 read_counts，表达矩阵已包含在内） ----------
  # 如果用户未提供 contrast_file，使用构造函数的默认路径（可能不存在，故仍建议检查）
  if (is.null(contrast_file)) {
    contrast_file <- "02.Quant_定量结果/GroupContrast.txt"
    message("未提供 contrast_file，使用默认路径: ", contrast_file)
  }
  
  design_obj <- object_experiment_design$new(
    count_matrix = count_matrix,
    gene_length = gene_length,
    col_data = NULL,                     # 可后续通过 colData 添加
    file_GroupSamples = sample_file,     # 构造函数会重新读取，但保持一致
    file_GroupContrast = contrast_file
  )
  
  # ---------- 4. 创建 PublicDataManager 并加载注释 ----------
  public_data_manager <- NULL
  if (!is.null(kegg_organism_code)) {
    message("[", Sys.time(), "] 初始化公共数据管理器，物种: ", kegg_organism_code)
    public_data_manager <- PublicDataManager$new(
      kegg_organism_code = kegg_organism_code,
      cache_dir = cache_dir,
      force_update = force_update
    )

        # 组装 load_combined_gtf 的参数
    load_args <- list(
      gtf_file_ncbi = gtf_file_ncbi,
      gtf_file_ensembl = gtf_file_ensembl,
      use_orgdb = use_orgdb,
      gene2go.gz = gene2go.gz,
      gene2ensembl.gz = gene2ensembl.gz,
      eggnog = eggnog,
      use_cache = TRUE,
      ...
    )
    load_args <- load_args[!sapply(load_args, is.null)]
    
    has_gtf <- !is.null(gtf_file_ncbi) || !is.null(gtf_file_ensembl)
    has_supp <- !is.null(gene2go.gz) || !is.null(gene2ensembl.gz) || !is.null(eggnog)
    
    if (has_gtf || has_supp) {
      do.call(public_data_manager$load_combined_gtf, load_args)
    } else {
      message("未提供 GTF 或补充注释文件，尝试从 OrgDb 加载基因注释...")
      tryCatch({
        public_data_manager$load_annotation_from_orgdb(use_cache = TRUE, load_go = TRUE, load_kegg = FALSE)
      }, error = function(e) {
        warning("从 OrgDb 加载基因注释失败: ", e$message,
                "\n后续 ORA/GSEA 可能无法使用基因符号（将用 Geneid 代替）。")
      })
    }
    
    if (annotate_rowData) {
      message("尝试将公共注释补充到 design_obj$rowData ...")
      tryCatch({
        public_data_manager$annotate_rowData(read_obj = design_obj, overwrite = FALSE)
      }, error = function(e) {
        warning("rowData 注释补充失败: ", e$message)
      })
    }
    
    # 加载 KEGG / GO 背景
    message("加载 KEGG 背景数据...")
    tryCatch({
      public_data_manager$get_kegg_background()
    }, error = function(e) {
      warning("KEGG 数据加载失败: ", e$message)
    })
    
    message("加载 GO 背景数据...")
    tryCatch({
      public_data_manager$get_go_background(ontology = "BP")
    }, error = function(e) {
      warning("GO 数据加载失败: ", e$message)
    })
  } else {
    message("未提供 kegg_organism_code，将不创建公共数据管理器。")
  }
  
  # ---------- 5. 返回结果（reads 与 design 指向同一对象，保持兼容） ----------
  result_list <- list(
    design = design_obj,      # 设计对象自身
    public_data = public_data_manager
  )
  
  message("[", Sys.time(), "] 数据对象创建完成，返回 list（reads/design 为同一对象，含表达矩阵+设计）")
  return(result_list)
}


# featureCounts 文件读取
read_featureCounts = function(file_path){
  tbl = data.table::fread(file = file_path, skip = 1) |> 
    dplyr::select(1, 6, 7:last_col()) |> 
    dplyr::as_tibble()
  colnames(tbl) = colnames(tbl) |> basename() |> gsub(pattern = ".sorted.bam", replacement = "")
  return(tbl)
}


# tranexp_object --------------------------------------------------------------------------------------------------


# 试验设计对象
tranexp_object = R6Class("tranexp_object",
                         public = list(
                           design = NULL,
                           public_data = NULL,
                           compareCluster_res = list(
                             res_kegg      = NULL,
                             res_go        = NULL,
                             res_kegg_GSEA   = NULL,
                             res_go_GSEA   = NULL
                           ),
                           
                           initialize = function(transcript_object) {
                             suppressMessages(suppressWarnings(require(tidyverse)))
                             suppressMessages(suppressWarnings(require(clusterProfiler)))
                             
                             self$design = transcript_object$design
                             self$public_data = transcript_object$public_data
                           },
                           
                           go_item = function(){
                             cat("# 免疫反应热点基因集合\n")
                             cat("# GO:0051607 defense response to virus\n")
                             cat("# GO:0006955 immune response\n")
                             cat("# GO:0006954 inflammatory response\n")
                             cat("# GO:0005125 cytokine activity\n")
                             cat("# GO:0050896 response to stimulus\n")
                           },
                           
                           # 比较转录组学 经典富集分析
                           compareCluster_GO_enrich      = function(contrast_range = NULL){
                             
                             dataset_DEG_tbl = purrr::map_dfr(.x = contrast_range, .f = \(x){
                               run_calculate_differential_gene(design_obj = self$design,
                                                               contrast_index = x,
                                                               method = "edgeR")$deg_tbl
                             }) |> 
                               dplyr::left_join(self$public_data$get_gene_annotation(), by = join_by(Geneid))
                             
                             res = clusterProfiler::compareCluster(geneClusters = as.formula(paste("SYMBOL ~", "group + Result_expression")),
                                                                   fun = "enricher",
                                                                   TERM2GENE = self$public_data$get_go_background()[["TERM2GENE"]],
                                                                   TERM2NAME = self$public_data$get_go_background()[["TERM2NAME"]],
                                                                   data =  dataset_DEG_tbl)
                             res = dplyr::as_tibble(res) |> 
                               dplyr::mutate(group = gsub(pattern = "(.*)\\.(Up|Down)", replacement = "\\1", x = Cluster), 
                                             Result_expression = gsub(pattern = "(.*)\\.(Up|Down)", replacement = "\\2", x = Cluster))
                             
                             res = dplyr::group_by(res, Cluster) |> 
                               dplyr::arrange(p.adjust, .by_group = TRUE)
                             
                             self$compareCluster_res$res_go = res
                             return(res)
                           },
                           compareCluster_GO_enrich_plot = function(contrast_range = NULL, 
                                                                    GO_term = NULL, 
                                                                    string_pattern_include = NULL, 
                                                                    string_pattern_exclude = NULL, 
                                                                    out_dir = NULL){
                             self$go_item()
                             
                             # 更新数据
                             res = self$compareCluster_GO_enrich(contrast_range)
                             
                             if(!is.null(GO_term) || (!is.null(string_pattern_include) || !is.null(string_pattern_exclude))){
                               limited_go_term = self$go_offspring(GO_term = GO_term, 
                                                                   string_pattern_include = string_pattern_include, 
                                                                   string_pattern_exclude = string_pattern_exclude)
                               
                               res = dplyr::filter(res, ID %in% limited_go_term)
                             }
                             
                             res = dplyr::mutate(res, GeneRatio = map_dbl(.x = GeneRatio, .f = \(x){
                               eval(parse(text = x))
                             }))
                             
                             res_plot = ggplot(res, aes(x = group, y = Description)) +
                               geom_point(aes(fill = -log10(p.adjust), size = GeneRatio, shape = Result_expression), color = "white") +
                               scale_fill_gradientn(colours = c("#000000", "#009900", "#FFFF00")) +
                               scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
                               ggtitle(label = "GO Comparative Enrichment Analysis") +
                               theme_bw() +
                               theme(
                                 panel.grid = element_blank(),
                                 axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
                                 strip.background = element_rect(fill = NA)
                               )
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/go_enrich_plot")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               return(res_plot)
                             }
                           },
                           
                           compareCluster_KEGG_enrich      = function(contrast_range = NULL){
                             dataset_DEG_tbl = purrr::map_dfr(.x = contrast_range, .f = \(x){
                               run_calculate_differential_gene(design_obj = self$design,
                                                               contrast_index = x,
                                                               method = "edgeR")$deg_tbl
                             }) |> 
                               dplyr::left_join(self$public_data$get_gene_annotation(), by = join_by(Geneid))
                             
                             res = clusterProfiler::compareCluster(geneClusters = as.formula(paste("SYMBOL ~", "group + Result_expression")),
                                                                   fun = "enricher",
                                                                   TERM2GENE = self$public_data$get_kegg_background()[["TERM2GENE"]],
                                                                   TERM2NAME = self$public_data$get_kegg_background()[["TERM2NAME"]],
                                                                   data =  dataset_DEG_tbl)
                             # 添加通路分级
                             kegg_background = readRDS("/data/shenlab_data/00.Database/kegg_db.RDS")
                             kegg_background = dplyr::select(kegg_background, MapID, Class_level_1, Class_level_2) |> dplyr::distinct()
                             
                             res =  dplyr::as_tibble(res) |> 
                               dplyr::mutate(group = gsub(pattern = "(.*)\\.(Up|Down)", replacement = "\\1", x = Cluster), 
                                             Result_expression = gsub(pattern = "(.*)\\.(Up|Down)", replacement = "\\2", x = Cluster)) |> 
                               dplyr::mutate(ID = sub(pattern = "[a-z]{2,3}", replacement = "\\1", x = ID)) |> 
                               dplyr::left_join(dplyr::mutate(kegg_background, MapID = sub(pattern = "map", replacement = "", x = MapID)), by = c("ID" = "MapID"))
                             
                             res = dplyr::group_by(res, Cluster) |> 
                               dplyr::arrange(p.adjust, .by_group = TRUE)
                             
                             self$compareCluster_res$res_kegg = res
                             return(res)
                           },
                           compareCluster_KEGG_enrich_plot = function(contrast_range = NULL, out_dir = NULL){
                             # 更新数据
                             res = self$compareCluster_KEGG_enrich(contrast_range)
                             self$compareCluster_res$res_kegg = res
                             
                             res = dplyr::mutate(res, GeneRatio = map_dbl(.x = GeneRatio, .f = \(x){
                               eval(parse(text = x))
                             }))
                             
                             res_plot = ggplot(res, aes(x = group, y = Description)) +
                               geom_point(aes(fill = -log10(p.adjust), size = GeneRatio, shape = Result_expression), color = "white") +
                               scale_fill_gradientn(colours = c("#000000", "#009900", "#FFFF00")) +
                               scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
                               ggtitle(label = "KEGG Comparative Enrichment Analysis") +
                               theme_bw() +
                               theme(
                                 panel.grid = element_blank(),
                                 axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
                                 strip.background = element_rect(fill = NA)
                               )
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/GSEA_kegg_enrich_plot")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               return(res_plot)
                             }
                             
                             # 手动转换
                             # compareCluster_KEGG_enrich_formated = dplyr::mutate(compareCluster_KEGG_enrich_formated, geneID = map_chr(.x = geneID, .f = \(x){
                             #   tibble(ENTREZID = str_split(x, pattern = "/")[[1]]) |>
                             #     dplyr::left_join(dplyr::select(trans$gtf_annotation$tbl_gtf, SYMBOL, ENTREZID) |> dplyr::distinct(), by = join_by(ENTREZID)) |>
                             #     dplyr::pull(SYMBOL) |>
                             #     base::paste(collapse = "/")
                             # }))
                           },
                           
                           # 比较转录组学 GSEA富集分析
                           compareCluster_GSEA_GO_enrich      = function(contrast_range = NULL){
                             dataset_DEG_tbl = purrr::map_dfr(.x = contrast_range, .f = \(x){
                               run_calculate_differential_gene(design_obj = self$design,
                                                               contrast_index = x,
                                                               method = "edgeR")$deg_tbl
                             }) |> 
                               dplyr::left_join(self$public_data$get_gene_annotation(), by = join_by(Geneid))
                             
                             if(is.null(self$compareCluster_res$res_go_GSEA)){
                               res = clusterProfiler::compareCluster(geneClusters = as.formula(paste("ENTREZID | logFC ~", "group")), 
                                                                     fun = "GSEA",
                                                                     TERM2GENE = self$public_data$get_go_background()[["TERM2GENE"]],
                                                                     TERM2NAME = self$public_data$get_go_background()[["TERM2NAME"]],
                                                                     data =  dataset_DEG_tbl)
                               
                               self$compareCluster_res$res_go_GSEA = res
                               return(res)
                             }else{
                               return(self$compareCluster_res$res_go_GSEA)
                             }
                           },
                           compareCluster_GSEA_GO_enrich_plot = function(contrast_range = NULL,
                                                                         GO_term = NULL, 
                                                                         string_pattern_include = NULL, 
                                                                         string_pattern_exclude = NULL, 
                                                                         out_dir = NULL){
                             
                             res = self$compareCluster_GSEA_GO_enrich(contrast_range = contrast_range)
                             
                             if(!is.null(GO_term) || (!is.null(string_pattern_include) || !is.null(string_pattern_exclude))){
                               limited_go_term = self$go_offspring(GO_term = GO_term, 
                                                                   string_pattern_include = string_pattern_include, 
                                                                   string_pattern_exclude = string_pattern_exclude)
                               
                               res = dplyr::filter(res, ID %in% limited_go_term)
                             }
                             
                             # 数据格式化
                             res = dplyr::as_tibble(res) %>% 
                               dplyr::mutate(Expression = ifelse(NES > 0, "Up", "Down")) |> 
                               dplyr::mutate(leading_edge = as.numeric(gsub("tags=([0-9]{1,2})%.*", "\\1", x = leading_edge))) %>%
                               dplyr::mutate(Cluster = gsub(pattern = "\\.", replacement = " vs. ", x = Cluster))
                             
                             
                             res_plot = ggplot(res, aes(x = Cluster, y = Description, size = abs(NES), color = -log10(p.adjust))) +
                               geom_point(aes(fill = -log10(p.adjust), shape = Expression), color = "white") +
                               scale_fill_gradientn(colours = c("#000000", "#009900", "#FFFF00")) +
                               scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
                               theme_bw()+
                               theme(axis.text.y = element_text(size = 10),
                                     axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
                                     axis.title = element_blank(),
                                     panel.grid.major = element_blank())
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/GSEA_go_enrich_plot")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               return(res_plot)
                             }
                           },
                           
                           compareCluster_GSEA_KEGG_enrich      = function(contrast_range = NULL){
                             dataset_DEG_tbl = purrr::map_dfr(.x = contrast_range, .f = \(x){
                               run_calculate_differential_gene(design_obj = self$design,
                                                               contrast_index = x,
                                                               method = "edgeR")$deg_tbl
                             }) |> 
                               dplyr::left_join(self$public_data$get_gene_annotation(), by = join_by(Geneid))
                             
                             res = clusterProfiler::compareCluster(geneClusters = as.formula(paste("ENTREZID | logFC ~", "group")), 
                                                                   fun = "GSEA",
                                                                   TERM2GENE = self$public_data$get_kegg_background()[["TERM2GENE"]],
                                                                   TERM2NAME = self$public_data$get_kegg_background()[["TERM2NAME"]],
                                                                   data =  dataset_DEG_tbl)
                             
                             self$compareCluster_res$res_kegg_GSEA = res
                             return(res)
                           },
                           compareCluster_GSEA_KEGG_enrich_plot = function(contrast_range = NULL,
                                                                           filter_term_pattern = NULL,
                                                                           out_dir = NULL){
                             
                             res = self$compareCluster_GSEA_KEGG_enrich(contrast_range = contrast_range)
                             
                             # 添加通路分级
                             kegg_background = readRDS("/data/shenlab_data/00.Database/kegg_db.RDS")
                             kegg_background = dplyr::select(kegg_background, MapID, Class_level_1, Class_level_2) |> dplyr::distinct()
                             
                             # 数据格式化
                             res = dplyr::as_tibble(res) |> 
                               dplyr::mutate(Expression = ifelse(NES > 0, "Up", "Down")) |> 
                               dplyr::mutate(ID = sub(pattern = "[a-z]{2,3}", replacement = "\\1", x = ID)) |> 
                               dplyr::left_join(dplyr::mutate(kegg_background, ID = sub(pattern = "map", replacement = "", x = MapID)), by = join_by(ID)) |> 
                               dplyr::mutate(leading_edge = as.numeric(gsub("tags=([0-9]{1,2})%.*", "\\1", x = leading_edge))) |> 
                               dplyr::mutate(Cluster = gsub(pattern = "\\.", replacement = " vs. ", x = Cluster))
                             
                             res = dplyr::filter(res, 
                                                 if_any(.cols = everything(), 
                                                        .fns = \(x){grepl(pattern = paste0("(", paste(filter_term_pattern, collapse = "|"), ")"), x = x)}))
                             
                             
                             # 检查res是否剩余
                             if(dim(res)[1] == 0){
                               warning(paste0("[", Sys.time(), "] ", "无结果"))
                               return(NULL)
                             }
                             
                             res_plot = ggplot(res, aes(x = Cluster, y = Description, size = abs(NES), color = -log10(p.adjust))) +
                               geom_point(aes(fill = -log10(p.adjust), shape = Expression), color = "white") +
                               scale_fill_gradientn(colours = c("#000000", "#009900", "#FFFF00")) +
                               scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
                               theme_bw()+
                               theme(axis.text.y = element_text(size = 10),
                                     axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
                                     # panel.grid = element_blank(),
                                     panel.grid = element_blank(),
                                     axis.title = element_blank())
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/GSEA_kegg_enrich_plot")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               return(res_plot)
                             }
                           },
                           
                           # 网络图分析
                           GSEA_cnetplot = function(contrast_numbers = NULL){
                             res = self$compareCluster_GSEA_GO_enrich(contrast_numbers = 1:7)
                             res_plot = enrichplot::cnetplot(res)
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/GSEA_cnetplot")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               return(res_plot)
                             }
                           },
                           
                           # 可视化
                           # 根据GOID和关键字搜索和返回，分管下属GOID和输入GOIDID
                           # GO term收缩限制
                           go_offspring = function(GO_term = NULL, 
                                                   string_pattern_include = NULL,
                                                   string_pattern_exclude = NULL) {
                             
                             # ========== 参数检查 ==========
                             if (!is.null(GO_term) && !is.null(string_pattern_include)) {
                               stop("GO_term 和 string_pattern_include 互斥，只能指定其中一个")
                             }
                             
                             if (is.null(GO_term) && is.null(string_pattern_include)) {
                               stop("必须指定 GO_term 或 string_pattern_include 其中一个")
                             }
                             
                             if (!require(GO.db, quietly = TRUE)) {
                               stop("需要安装 GO.db 包")
                             }
                             
                             # ========== 模式 1：按 GO_term 取 offspring ==========
                             if (!is.null(GO_term)) {
                               if (any(!GO_term %in% names(as.list(GO.db::GOBPOFFSPRING)))) {
                                 no_index = which(!GO_term %in% names(as.list(GO.db::GOBPOFFSPRING)))
                                 stop("GO_term '", GO_term[no_index], "' 不存在于 GOBPOFFSPRING 中")
                               }
                               
                               offspring <- unlist(as.list(GO.db::GOBPOFFSPRING)[GO_term])
                               result <- unique(c(GO_term, offspring))
                               
                               message("GO_term 模式：", GO_term)
                               message("  - Offspring 数量：", length(offspring))
                               message("  - 返回总数：", length(result))
                               
                               return(result)
                             }
                             
                             # 所有 GO term
                             all_terms <- as.list(GO.db::GOTERM)
                             
                             # 只保留 BP
                             bp_terms <- all_terms[sapply(all_terms, function(term) term@Ontology == "BP")]
                             
                             bp_ids   <- names(bp_terms)
                             bp_names <- sapply(bp_terms, function(term) term@Term)
                             
                             message("BP term 总数：", length(bp_ids))
                             
                             # ---------- include 过滤 ----------
                             include_hit <- rep(FALSE, length(bp_names))
                             
                             # ========== 模式 2：按字符串在 BP 中搜索 ==========
                             if (!is.null(string_pattern_include)) {
                               for (pat in string_pattern_include) {
                                 include_hit <- include_hit | grepl(pat, bp_names, ignore.case = TRUE)
                               }
                               
                               matched_ids = bp_ids[include_hit]
                               matched_names = bp_names[include_hit]
                               
                               message("Include 匹配关键词：", paste(string_pattern_include, collapse = ", "))
                               message("  - 匹配到：", length(matched_ids), " 个 BP term")
                             }
                             
                             # ---------- exclude 过滤 ----------
                             exclude_hit <- rep(FALSE, length(bp_names))
                             
                             if (!is.null(string_pattern_exclude) && length(matched_ids) > 0) {
                               for (pat in string_pattern_exclude) {
                                 exclude_hit <- exclude_hit | grepl(pat, bp_names, ignore.case = TRUE)
                               }
                               
                               matched_ids = bp_ids[!exclude_hit]
                               matched_names = bp_names[!exclude_hit]
                               
                               message("Exclude 排除关键词：", paste(string_pattern_exclude, collapse = ", "))
                               message("  - 排除掉：", length(matched_ids), " 个 BP term")
                             }
                             
                             if (!is.null(string_pattern_include) && !is.null(string_pattern_exclude)) {
                               message(" string_pattern_include - string_pattern_exclude 排除掉：", sum(include_hit & exclude_hit), " 个 BP term")
                             }
                             
                             res_matched_ids <- bp_ids[include_hit & (!exclude_hit)]
                             res_matched_names <- bp_names[include_hit & (!exclude_hit)]
                             
                             if (length(res_matched_ids) == 0) {
                               warning("未找到符合 include/exclude 条件的 BP GO term")
                               return(character(0))
                             }
                             
                             message("最终返回 BP term 数量：", length(res_matched_ids))
                             
                             return(res_matched_ids)
                           },
                           
                           # 促炎因子 抗炎因子 调节因子 趋化因子
                           cytokines_chemokines_gene = function() {
                             # 编码细胞因子的基因
                             gene_cytokines_Pro_inflammatory <- c(
                               "IL1A", "IL1B",
                               "IL6",
                               "IL12A", "IL12B",  # IL-12 (p35 + p40)
                               "IL15",
                               "IL17A",
                               "IL18",
                               "IL23A",
                               "IL31",
                               "TNF",      # TNF-alpha
                               "LTA",      # TNF-beta
                               "IFNA1",    # IFN-alpha
                               "IFNG"      # IFN-gamma
                             )
                             gene_cytokines_Anti_inflammatory = c(
                               "IL1RN",    # IL-1 receptor antagonist
                               "IL10",
                               "IL13",
                               "IL22"
                             )
                             gene_cytokines_versatility_inflammatory <- c(
                               "CSF2",     # GM-CSF
                               "IL2",
                               "IL4",
                               "IL5",
                               "IL7",
                               "IL9",
                               "IL21",
                               "IL27"
                             )
                             gene_chemokines <- c(
                               "CXCL1",    # GRO-alpha
                               "CXCL8",    # IL-8
                               "CXCL10",   # IP-10
                               "CXCL12",   # SDF-1alpha
                               "CCL2",     # MCP-1
                               "CCL3",     # MIP-1alpha
                               "CCL4",     # MIP-1beta
                               "CCL5",     # RANTES
                               "CCL11"     # Eotaxin
                             )
                             
                             tbl_cytokines_chemokines_gene <- dplyr::bind_rows(
                               tibble(fun = "cytokines_Pro_inflammatory", gene = gene_cytokines_Pro_inflammatory),
                               tibble(fun = "cytokines_Anti_inflammatory", gene = gene_cytokines_Anti_inflammatory),
                               tibble(fun = "cytokines_versatility_inflammatory", gene = gene_cytokines_versatility_inflammatory),
                               tibble(fun = "chemokines", gene = gene_chemokines)
                               # tibble(fun = "immune_response_pathway", gene = gene_immune_response_pathway),
                             ) |>
                               dplyr::mutate(gene = fct_inorder(gene), 
                                             fun  = fct_inorder(fun))
                             
                             return(tbl_cytokines_chemokines_gene)
                           },
                           
                           # 热图, 组与组相比较
                           pheatmap_go_group2group = function(contrast_numbers = NULL, 
                                                              GO_term = NULL, 
                                                              string_pattern_include = NULL,
                                                              string_pattern_exclude = NULL,
                                                              cytokines_chemokines_gene = FALSE, 
                                                              out_dir = NULL){
                             self$go_item()
                             
                             # 限制在免疫反应热图
                             res              = self$transcript_object$merge_all_GEG_list(Contrast_number_range = contrast_numbers)
                             res_significant  = self$transcript_object$merge_significant_GEG_list(Contrast_number_range = contrast_numbers)
                             
                             # 限制输出
                             res = dplyr::filter(res, Geneid %in% res_significant$Geneid)
                             
                             if(!is.null(GO_term) || (!is.null(string_pattern_include) || !is.null(string_pattern_exclude)) && !cytokines_chemokines_gene){
                               limited_go_term = self$public_data$get_go_background()[["TERM2GENE"]] |> 
                                 dplyr::filter(TERM %in% tranexp$go_offspring(GO_term = GO_term, 
                                                                              string_pattern_include = string_pattern_include, 
                                                                              string_pattern_exclude = string_pattern_exclude)) |> 
                                 dplyr::pull(GENE)
                               
                               res = dplyr::filter(res, Geneid %in% limited_go_term)
                             } else if(cytokines_chemokines_gene){
                               limited_go_term = self$cytokines_chemokines_gene()
                               
                               res = dplyr::filter(res, Geneid %in% limited_go_term$gene)
                             }
                             
                             # 按照contrast出现顺序进行排序
                             res = dplyr::mutate(res, group = forcats::fct_inorder(group))
                             
                             expression_matrix = tidyr::pivot_wider(res, id_cols = Geneid, names_from = group, names_sort = TRUE, values_from = logFC, values_fn = mean, values_fill = 0)
                             
                             expression_matrix = as.data.frame(expression_matrix)
                             rownames(expression_matrix) = expression_matrix$Geneid
                             expression_matrix$Geneid = NULL
                             
                             # colorRampPalette 颜色插值
                             res_plot = pheatmap::pheatmap(
                               expression_matrix, 
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(8,"PiYG")))(100), # 自定义颜色渐变
                               cluster_cols = FALSE, 
                               border_color = NA, 
                               breaks = seq(-max(abs(expression_matrix), na.rm = T), max(abs(expression_matrix), na.rm = T), length.out = 101),
                               silent = TRUE)
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/pheatmap_go")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               print(res_plot)
                               return(expression_matrix)
                             }
                           },
                           pheatmap_go_sample2sample = function(contrast_numbers = NULL, 
                                                                GO_term = NULL, 
                                                                string_pattern_include = NULL,
                                                                string_pattern_exclude = NULL,
                                                                cytokines_chemokines_gene = FALSE, 
                                                                out_dir = NULL){
                             self$go_item()
                             
                             # 限制在免疫反应热图
                             res              = self$transcript_object$TPM
                             res_significant  = self$transcript_object$merge_significant_GEG_list(Contrast_number_range = contrast_numbers)
                             
                             # 限制输出
                             res = res[res_significant$Geneid, ]
                             
                             if(!is.null(GO_term) || (!is.null(string_pattern_include) || !is.null(string_pattern_exclude)) && !cytokines_chemokines_gene){
                               limited_go_term = self$public_data$get_go_background()[["TERM2GENE"]] |> 
                                 dplyr::filter(TERM %in% tranexp$go_offspring(GO_term = GO_term, 
                                                                              string_pattern_include = string_pattern_include, 
                                                                              string_pattern_exclude = string_pattern_exclude)) |> 
                                 dplyr::pull(GENE)
                               res = res[rownames(res) %in% limited_go_term, ]
                             } else if(cytokines_chemokines_gene){
                               limited_go_term = self$cytokines_chemokines_gene()
                               res = res[rownames(res) %in% limited_go_term, ]
                             }
                             
                             expression_matrix = scale(t(res)) |> t()
                             
                             expression_matrix = as.data.frame(expression_matrix)
                             rownames(expression_matrix) = rownames(res)
                             expression_matrix$Geneid = NULL
                             
                             # colorRampPalette 颜色插值
                             res_plot = pheatmap::pheatmap(
                               expression_matrix, 
                               color = colorRampPalette(c("blue", "white", "red"))(100), # colorRampPalette(rev(RColorBrewer::brewer.pal(8,"Spectral")))(100), # 自定义颜色渐变
                               cluster_cols = FALSE, 
                               border_color = NA, 
                               breaks = seq(-max(abs(expression_matrix), na.rm = T), max(abs(expression_matrix), na.rm = T), length.out = 101),
                               silent = TRUE)
                             
                             if(!is.null(out_dir)){
                               file_prefix = stringr::str_glue("{out_dir}/pheatmap_go")
                               private$figure_output_ggplot(res_plot = res_plot, file_prefix = file_prefix)
                             }else{
                               print(res_plot)
                               return(expression_matrix)
                             }
                           }
                           
                         ), private = list(
                           figure_output_ggplot = function(res_plot, file_prefix = NULL, header_title = FALSE){
                             
                             width = (max(nchar(res$Description)) + length(unique(res$Cluster)) + 5)/5
                             if(header_title){
                               height = (max(unique(length(res$Description))))/10 + 1
                             }else{
                               height = max(unique(length(res$Description)))/10
                             }
                             private$figure_output(
                               plot = res_plot,
                               file_prefix = file_prefix,
                               width = width, height = height
                             )
                           }, 
                           figure_output = function(plot, file_prefix, width, height){
                             message(paste0("[", Sys.time(), "] ", "输出图像 ", stringr::str_glue("{file_prefix}.pdf")))
                             message(paste0("[", Sys.time(), "] ", "输出图像 ", stringr::str_glue("{file_prefix}.png")))
                             ggsave(stringr::str_glue("{file_prefix}.pdf"), plot, width = width, height = height, limitsize = FALSE, units = "in")
                             ggsave(stringr::str_glue("{file_prefix}.png"), plot, width = width, height = height, limitsize = FALSE, units = "in")
                           }
                         ))
