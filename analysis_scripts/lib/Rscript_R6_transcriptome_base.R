# 基础类对象
.require_cran_package <- function(pkg, install_missing = FALSE, load = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing) {
      stop("需要安装 CRAN 包: install.packages('", pkg, "')")
    }
    install.packages(pkg)
  }
  if (load) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}

.require_bioc_package <- function(pkg, install_missing = FALSE, load = TRUE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (!install_missing) {
      stop("需要安装 BiocManager: install.packages('BiocManager')")
    }
    install.packages("BiocManager")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing) {
      stop("需要安装 Bioconductor 包: BiocManager::install('", pkg, "')")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  if (load) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}

check_transcriptome_package_environment <- function(install_missing = FALSE,
                                                   include_optional = TRUE) {
  cran_pkgs <- data.frame(
    package = c("R6", "tidyverse", "data.table", "openxlsx", "ggplot2"),
    source = "CRAN",
    role = c("R6 class system", "data wrangling", "fast table IO", "Excel output", "plotting"),
    stringsAsFactors = FALSE
  )
  bioc_pkgs <- data.frame(
    package = c("SummarizedExperiment", "S4Vectors", "AnnotationDbi", "GO.db"),
    source = "Bioconductor",
    role = c("assay container", "rowData/colData", "annotation select API", "GO term metadata"),
    stringsAsFactors = FALSE
  )
  optional_pkgs <- data.frame(
    package = c("rtracklayer", "clusterProfiler", "edgeR", "DESeq2", "limma",
                "enrichplot", "pathview"),
    source = "Bioconductor",
    role = c("GTF import", "KEGG/GO enrichment", "edgeR DEG", "DESeq2 DEG",
             "limma-voom DEG", "enrichment plots", "KEGG pathway plots"),
    stringsAsFactors = FALSE
  )
  pkg_tbl <- rbind(cran_pkgs, bioc_pkgs, if (include_optional) optional_pkgs else NULL)
  pkg_tbl$installed <- vapply(pkg_tbl$package, requireNamespace, logical(1), quietly = TRUE)
  pkg_tbl$version <- vapply(pkg_tbl$package, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1))
  
  if (install_missing && any(!pkg_tbl$installed)) {
    for (i in which(!pkg_tbl$installed)) {
      if (pkg_tbl$source[i] == "CRAN") {
        .require_cran_package(pkg_tbl$package[i], install_missing = TRUE, load = FALSE)
      } else {
        .require_bioc_package(pkg_tbl$package[i], install_missing = TRUE, load = FALSE)
      }
    }
    pkg_tbl$installed <- vapply(pkg_tbl$package, requireNamespace, logical(1), quietly = TRUE)
    pkg_tbl$version <- vapply(pkg_tbl$package, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
    }, character(1))
  }
  
  pkg_tbl
}

.require_cran_package("R6")
.require_cran_package("tidyverse")
.require_bioc_package("SummarizedExperiment")
.require_bioc_package("S4Vectors")


# ==================== read_counts  原始矩阵 基础数据预处理和储存类 ======================
# 基于 SummarizedExperiment 存储计数矩阵和基因长度
read_counts <- R6::R6Class(
  classname = "read_counts",
  private = list(
    .se = NULL,
    calc_tpm = function() {
      len <- self$rowData$Length
      if (any(len == 0)) warning("存在长度为0的基因，TPM 计算可能产生 Inf")
      cnt <- self$get_assay("counts")
      rpk <- sweep(cnt, 1, len, FUN = "/")
      tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
      message("TPM 计算完成，已保存为 assay 'tpm'")
      return(tpm)
    },
    
    calc_fpkm = function() {
      len <- self$rowData$Length
      cnt <- self$get_assay("counts")
      fpkm <- sweep(cnt, 1, len, FUN = "/")
      fpkm <- sweep(fpkm, 2, colSums(cnt), FUN = "/") * 1e9
      message("FPKM 计算完成，已保存为 assay 'fpkm'")
      return(fpkm)
    },
    
    add_assay = function(name, matrix) {
      if (!identical(dim(matrix), dim(private$.se))) stop("矩阵维度不一致")
      matrix <- matrix[rownames(private$.se), colnames(private$.se), drop = FALSE]
      assays_list <- assays(private$.se)
      assays_list[[name]] <- matrix
      assays(private$.se) <- assays_list
      invisible(self)
    }
  ),
  active = list(
    assays  = function() SummarizedExperiment::assays(private$.se),
    rowData = function() SummarizedExperiment::rowData(private$.se),
    colData = function() SummarizedExperiment::colData(private$.se),
    dims    = function() dim(private$.se)
  ),
  public = list(
    initialize = function(count_matrix, gene_length = NULL, col_data = NULL) {
      if (missing(count_matrix)) stop("必须提供 count_matrix")
      if (!is.matrix(count_matrix)) count_matrix <- as.matrix(count_matrix)
      
      # 处理 gene_length
      if (is.null(gene_length)) {
        gene_length <- rep(1, nrow(count_matrix))
        names(gene_length) <- rownames(count_matrix)
      } else {
        if (length(gene_length) != nrow(count_matrix)) {
          stop("gene_length 长度必须与 count_matrix 行数一致")
        }
        if (is.null(names(gene_length))) {
          names(gene_length) <- rownames(count_matrix)
        }
      }
      
      # 构建基础 rowData：包含 Length
      row_df <- S4Vectors::DataFrame(
        Length = gene_length,
        row.names = rownames(count_matrix)
      )
      
      # 处理 colData（不变）
      if (is.null(col_data)) {
        col_df <- S4Vectors::DataFrame(
          sample = colnames(count_matrix),
          row.names = colnames(count_matrix)
        )
      } else {
        if (!is(col_data, "DataFrame")) col_data <- S4Vectors::DataFrame(col_data)
        col_df <- col_data
        if (!identical(rownames(col_df), colnames(count_matrix))) {
          warning("col_data 行名与矩阵列名不一致，将重新对齐")
          col_df <- col_df[colnames(count_matrix), , drop = FALSE]
        }
      }
      
      # 创建 SummarizedExperiment（注意行名已经同步）
      private$.se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = count_matrix[rownames(row_df), , drop = FALSE]),
        rowData = row_df,
        colData = col_df
      )
      
      message(sprintf("[%s] read_counts 对象创建: %d 基因 × %d 样本",
                      Sys.time(), nrow(private$.se), ncol(private$.se)))
    },
    
    keep_samples = function(samples) {
      if (!all(samples %in% colnames(private$.se))) {
        missing <- setdiff(samples, colnames(private$.se))
        stop("以下样本不存在: ", paste(missing, collapse = ", "))
      }
      private$.se <- private$.se[, samples, drop = FALSE]
      invisible(self)
    },
    
    get_assay = function(name = "counts") {
      if (!name %in% names(assays(private$.se))) {
        message("assay ", name, " 不存在，重新计算")
        if(name == "tpm") {
          private$add_assay("tpm", private$calc_tpm())
          assay(private$.se, name)
        }else if(name == "fpkm"){
          private$add_assay("fpkm", private$calc_fpkm())
          assay(private$.se, name)
        }else{
          stop("assay ", name, " 计算方式不存在")
        }
      }
      assay(private$.se, name)
    },
    
    filter_genes = function(genes) {
      if (!all(genes %in% rownames(private$.se))) {
        missing <- setdiff(genes, rownames(private$.se))
        warning("以下基因不存在: ", paste(missing, collapse = ", "))
        genes <- intersect(genes, rownames(private$.se))
      }
      new_obj <- self$clone()
      new_obj$.__enclos_env__$private$.se <- private$.se[genes, ]
      new_obj
    },
    
    filter_samples = function(samples) {
      if (!all(samples %in% colnames(private$.se))) {
        missing <- setdiff(samples, colnames(private$.se))
        warning("以下样本不存在: ", paste(missing, collapse = ", "))
        samples <- intersect(samples, colnames(private$.se))
      }
      new_obj <- self$clone()
      new_obj$.__enclos_env__$private$.se <- private$.se[, samples]
      new_obj
    },
    
    get_rowData_table = function(add_Geneid = TRUE) {
      row_tbl <- as.data.frame(SummarizedExperiment::rowData(private$.se))
      if (add_Geneid && !"Geneid" %in% colnames(row_tbl)) {
        row_tbl <- data.frame(Geneid = rownames(private$.se), row_tbl,
                              check.names = FALSE, stringsAsFactors = FALSE)
      }
      row_tbl
    },
    
    add_rowData = function(row_data,
                           by = "Geneid",
                           overwrite = FALSE,
                           keep_unmatched = TRUE) {
      if (missing(row_data) || is.null(row_data)) {
        stop("必须提供 row_data")
      }
      row_data <- as.data.frame(row_data, stringsAsFactors = FALSE)
      if (nrow(row_data) == 0) {
        stop("row_data 为空")
      }
      
      current <- as.data.frame(SummarizedExperiment::rowData(private$.se))
      current_rownames <- rownames(private$.se)
      
      if (by %in% colnames(row_data)) {
        match_key <- row_data[[by]]
      } else if (!is.null(rownames(row_data))) {
        match_key <- rownames(row_data)
      } else {
        stop("row_data 必须包含匹配列 ", by, " 或有效行名")
      }
      
      if (anyDuplicated(match_key)) {
        warning("row_data 中存在重复匹配键，将保留第一条记录")
        row_data <- row_data[!duplicated(match_key), , drop = FALSE]
        match_key <- match_key[!duplicated(match_key)]
      }
      
      idx <- match(current_rownames, match_key)
      if (!keep_unmatched && any(is.na(idx))) {
        stop("存在无法匹配注释的行: ", paste(head(current_rownames[is.na(idx)], 10), collapse = ", "))
      }
      
      add_cols <- setdiff(colnames(row_data), by)
      if (length(add_cols) == 0) {
        warning("row_data 没有可添加到 rowData 的列")
        return(invisible(self))
      }
      
      for (col in add_cols) {
        if (col %in% colnames(current) && !overwrite) {
          warning("rowData 已存在列 ", col, "，overwrite = FALSE，已跳过")
          next
        }
        new_col <- rep(NA, length(current_rownames))
        new_col[!is.na(idx)] <- row_data[[col]][idx[!is.na(idx)]]
        current[[col]] <- new_col
      }
      
      SummarizedExperiment::rowData(private$.se) <- S4Vectors::DataFrame(
        current,
        row.names = current_rownames
      )
      message("rowData 更新完成，当前列: ", paste(colnames(current), collapse = ", "))
      invisible(self)
    },
    
    print = function() {
      cat("<read_counts>\n")
      cat(sprintf("  Dimensions: %d genes × %d samples\n", 
                  nrow(private$.se), ncol(private$.se)))
      cat(sprintf("  Assays: %s\n", paste(names(assays(private$.se)), collapse = ", ")))
      cat(sprintf("  rowData columns: %s\n", paste(names(rowData(private$.se)), collapse = ", ")))
      cat(sprintf("  colData columns: %s\n", paste(names(colData(private$.se)), collapse = ", ")))
    }
  )
)

`[.read_counts` <- function(x, i, j, drop = TRUE) {
  if (missing(i) && missing(j)) return(x)
  # 先按行筛选
  if (!missing(i)) {
    x <- x$filter_genes(rownames(x$get_assay("counts"))[i])
  }
  # 再按列筛选
  if (!missing(j)) {
    x <- x$filter_samples(colnames(x$get_assay("counts"))[j])
  }
  x
}

# 2. 中缀运算符：基因/样本保留/排除（返回新对象）

# 保留基因
`%keep_genes%` <- function(x, genes) {
  x$filter_genes(genes)
}

# 排除基因（取反）
`%exclude_genes%` <- function(x, genes) {
  all_genes <- rownames(x$get_assay("counts"))
  keep <- setdiff(all_genes, genes)
  x$filter_genes(keep)
}

# 保留样本
`%keep_samples%` <- function(x, samples) {
  x$filter_samples(samples)
}

# 排除样本
`%exclude_samples%` <- function(x, samples) {
  all_samples <- colnames(x$get_assay("counts"))
  keep <- setdiff(all_samples, samples)
  x$filter_samples(keep)
}

# 3. subset 泛型方法
subset.read_counts <- function(x, genes = NULL, samples = NULL,
                               exclude_genes = NULL, exclude_samples = NULL, ...) {
  if (!is.null(genes)) x <- x$filter_genes(genes)
  if (!is.null(exclude_genes)) {
    all_genes <- rownames(x$get_assay("counts"))
    keep <- setdiff(all_genes, exclude_genes)
    x <- x$filter_genes(keep)
  }
  if (!is.null(samples)) x <- x$filter_samples(samples)
  if (!is.null(exclude_samples)) {
    all_samples <- colnames(x$get_assay("counts"))
    keep <- setdiff(all_samples, exclude_samples)
    x <- x$filter_samples(keep)
  }
  x
}

# ==================== object_experiment_design 实验设计功能 基础数据预处理和储存类 ======================
object_experiment_design <- R6::R6Class(
  classname = "Object_experiment_design",
  # 继承关系 ！！！
  inherit = read_counts,
  public = list(
    # 数据存储位置
    GroupSamples = NULL,
    GroupContrast = NULL,
    
    # 校验和验证 GroupSamples 和 GroupContrast 数据对应一致
    validate_data = function() {
      # 检查分组数据
      if (!is.matrix(self$GroupSamples) || !is.matrix(self$GroupContrast))
        stop("GroupSamples 和 GroupContrast 必须是矩阵")
      if (ncol(self$GroupSamples) != 2 || ncol(self$GroupContrast) != 2)
        stop("GroupSamples 和 GroupContrast 必须为2列矩阵")
      if (!is.character(self$GroupSamples) || !is.character(self$GroupContrast))
        stop("GroupSamples 和 GroupContrast 必须是字符型")
      if (anyNA(self$GroupSamples) || anyNA(self$GroupContrast))
        stop("GroupSamples 和 GroupContrast 不能包含 NA")
      if (!all(c("Samples", "Group") %in% colnames(self$GroupSamples)))
        stop("GroupSamples 必须包含 Samples 和 Group 列")
      if (any(!nzchar(trimws(self$GroupSamples))) || any(!nzchar(trimws(self$GroupContrast))))
        stop("GroupSamples 和 GroupContrast 不能包含空值")
      
      # 检查样本名称
      count_samples <- colnames(self$get_assay("counts"))
      group_samples <- self$GroupSamples[, "Samples"]
      if (is.null(count_samples) || length(count_samples) == 0)
        stop("read_counts 必须包含样本列名")
      if (anyNA(count_samples) || any(!nzchar(trimws(count_samples))))
        stop("read_counts 列样本名不能包含 NA 或空值")
      if (anyDuplicated(count_samples))
        stop("read_counts 列样本名重复: ", paste(unique(count_samples[duplicated(count_samples)]), collapse = ", "))
      if (anyDuplicated(group_samples))
        stop("GroupSamples 样本名重复: ", paste(unique(group_samples[duplicated(group_samples)]), collapse = ", "))
      
      # 检查双向对应
      missing_counts <- setdiff(group_samples, count_samples)
      missing_groups <- setdiff(count_samples, group_samples)
      errors <- character()
      if (length(missing_counts))
        errors <- c(errors, paste0("GroupSamples 中的样本未在 read_counts 中找到: ", paste(missing_counts, collapse = ", ")))
      if (length(missing_groups))
        errors <- c(errors, paste0("read_counts 中的样本未在 GroupSamples 中找到: ", paste(missing_groups, collapse = ", ")))
      if (length(errors)) stop(paste(errors, collapse = "\n"), call. = FALSE)
      
      # 检查比较分组
      groups_in_samples <- unique(self$GroupSamples[, "Group"])
      groups_in_contrast <- unique(c(self$GroupContrast[, 1], self$GroupContrast[, 2]))
      if (length(miss <- setdiff(groups_in_samples, groups_in_contrast)))
        stop("GroupSamples 中的组未在 GroupContrast 中出现: ", paste(miss, collapse = ", "))
      if (length(miss <- setdiff(groups_in_contrast, groups_in_samples)))
        stop("GroupContrast 中的组未在 GroupSamples 中找到: ", paste(miss, collapse = ", "))
      
      # 按分组表顺序对齐
      if (!identical(unname(count_samples), unname(group_samples))) {
        self$keep_samples(group_samples)
        message("read_counts 已按 GroupSamples 样本顺序对齐")
      }
      invisible(TRUE)
    },
    
    # 初始化
    initialize = function(
    count_matrix,
    gene_length = NULL,
    col_data = NULL,
    file_GroupSamples = "02.Quant_定量结果/GroupSamples.txt",
    file_GroupContrast = "02.Quant_定量结果/GroupContrast.txt"
    ) {
      
      # 1. 初始化父类
      super$initialize(count_matrix, gene_length, col_data)
      
      # 2. 读取 GroupSamples 文件
      self$read_GroupSamples(file_GroupSamples)
      
      # 3. 读取 GroupContrast 文件
      self$read_GroupContrast(file_GroupContrast)
      
      # 4. 检查设计样本与表达矩阵的一致性
      self$validate_data()
    },
    
    # 手动读取和替换分组和试验分组分析方案
    read_GroupSamples  = function(file_path) {
      if (!file.exists(file_path)) stop("文件不存在: ", file_path)
      
      res <- read.table(file = file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if (!all(colnames(res)[1:2] == c("Samples", "Group"))) {
        stop(paste0(file_path, "表格文件头读取错误，需要 Samples 和 Group 列"))
      }
      self$GroupSamples <- as.matrix(res[, 1:2])
      invisible(self)
    },
    read_GroupContrast = function(file_path) {
      if (!file.exists(file_path)) stop("文件不存在: ", file_path)
      res <- read.table(file = file_path, header = FALSE, sep = "\t",
                        col.names = c("GroupExperiment", "GroupControl"),
                        stringsAsFactors = FALSE)
      self$GroupContrast <- as.matrix(res[, 1:2])
      invisible(self)
    },
  
    
    # variable 索引序列构建
    validate_variable = function(variables) {
      if (!is.numeric(variables) || length(variables) == 0)
        stop("variables 必须是非空数值向量")
      if (any(!is.finite(variables))) stop("variables 包含非有限值")
      if (any(variables < 1 | variables > nrow(self$GroupContrast)))
        stop("variables 超出 1 到 ", nrow(self$GroupContrast), " 的范围")
      TRUE
    },
    
    # 获取样本分组信息
    get_all_samples = function() {
      self$GroupSamples[, "Samples"]
    },
    get_all_groups  = function() {
      unique(self$GroupSamples[, "Group"])
    },
    
    # 获得分组比较信息
    get_group_contrast_c = function(contrast_index) {
      self$validate_variable(contrast_index)
      if (length(contrast_index) != 1L)
        stop("contrast_index 必须是单个正整数")
      return(c(self$GroupContrast[contrast_index, 1], self$GroupContrast[contrast_index, 2]))
    },
    get_group_contrast_list   = function(contrast_range) {
      self$validate_variable(contrast_range)
      res_list <- lapply(contrast_range, function(i) {
        c(self$GroupContrast[i, 1], self$GroupContrast[i, 2])
      })
      names(res_list) <- contrast_range
      return(res_list)
    },
    
    # 四种筛选操作 深度拷贝关系
    # 选取样本子集
    operation_Samplessubset = function(GroupSamples_Samplessubset) {
      new_obj <- self$clone(deep = TRUE)
      GroupSamplesSubset <- new_obj$GroupSamples[
        new_obj$GroupSamples[, "Samples"] %in% GroupSamples_Samplessubset, , drop = FALSE
      ]
      if (nrow(GroupSamplesSubset) == 0) stop("筛选后没有剩余样本")
      new_obj$keep_samples(GroupSamples_Samplessubset)
      remaining <- unique(GroupSamplesSubset[, "Group"])
      new_obj$GroupContrast <- new_obj$GroupContrast[apply(
        new_obj$GroupContrast, 1, function(x) all(x %in% remaining)
      ), , drop = FALSE]
      new_obj$GroupSamples <- GroupSamplesSubset
      new_obj$validate_data()
      new_obj
    },
    # 剔除样本子集
    operation_Samplesexclusion = function(GroupSamples_Samplesexclusion) {
      new_obj <- self$clone(deep = TRUE)
      GroupSamplesSubset <- new_obj$GroupSamples[
        !new_obj$GroupSamples[, "Samples"] %in% GroupSamples_Samplesexclusion, , drop = FALSE
      ]
      if (nrow(GroupSamplesSubset) == 0) stop("筛选后没有剩余样本")
      new_obj$keep_samples(GroupSamplesSubset[, "Samples"])
      remaining <- unique(GroupSamplesSubset[, "Group"])
      new_obj$GroupContrast <- new_obj$GroupContrast[apply(
        new_obj$GroupContrast, 1, function(x) all(x %in% remaining)
      ), , drop = FALSE]
      new_obj$GroupSamples <- GroupSamplesSubset
      new_obj$validate_data()
      new_obj
    },
    # 选取分组子集
    operation_Groupsubset = function(GroupSamples_Groupsubset) {
      new_obj <- self$clone(deep = TRUE)
      GroupSamplesSubset <- new_obj$GroupSamples[
        new_obj$GroupSamples[, "Group"] %in% GroupSamples_Groupsubset, , drop = FALSE
      ]
      if (nrow(GroupSamplesSubset) == 0) stop("筛选后没有剩余样本")
      new_obj$keep_samples(GroupSamplesSubset[, "Samples"])
      new_obj$GroupContrast <- new_obj$GroupContrast[apply(
        new_obj$GroupContrast, 1, function(x) all(x %in% GroupSamples_Groupsubset)
      ), , drop = FALSE]
      new_obj$GroupSamples <- GroupSamplesSubset
      new_obj$validate_data()
      return(new_obj)
    },
    # 剔除分组子集
    operation_Groupexclusion = function(GroupSamples_Groupexclusion) {
      new_obj <- self$clone(deep = TRUE)
      GroupSamplesSubset <- new_obj$GroupSamples[
        !new_obj$GroupSamples[, "Group"] %in% GroupSamples_Groupexclusion, , drop = FALSE
      ]
      if (nrow(GroupSamplesSubset) == 0) stop("筛选后没有剩余样本")
      new_obj$keep_samples(GroupSamplesSubset[, "Samples"])
      remaining <- unique(GroupSamplesSubset[, "Group"])
      new_obj$GroupContrast <- new_obj$GroupContrast[apply(
        new_obj$GroupContrast, 1, function(x) all(x %in% remaining)
      ), , drop = FALSE]
      new_obj$GroupSamples <- GroupSamplesSubset
      new_obj$validate_data()
      new_obj
    },

    # 打印方法
    print = function() {
      cat("<object_experiment_design> (inherits from read_counts)\n")
      cat("  Samples count:", length(self$get_all_samples()), "\n")
      cat("  Samples:", paste(self$get_all_samples(), collapse = ", "), "\n")
      cat("  Groups count:", length(self$get_all_groups()), "\n")
      cat("  Groups:", paste(self$get_all_groups(), collapse = ", "), "\n")
      cat("  Contrasts number:", nrow(self$GroupContrast), "\n")
      cat("  Contrasts (first 6):\n")
      for (i in 1:min(6, nrow(self$GroupContrast))) {
        cat("    ", i, ":", self$GroupContrast[i, 1], " vs ", self$GroupContrast[i, 2], "\n")
      }
      if (nrow(self$GroupContrast) > 6) cat("    ...\n")
      cat("  Underlying expression data:\n")
      super$print()
    }
  )
)

# 符号重载 [
`[.object_experiment_design` <- function(x, i, type = c("group", "sample"), exclude = FALSE) {
  type <- match.arg(type)
  if (missing(i)) {
    return(x)
  }
  
  if (type == "group") {
    if (exclude) {
      x$operation_Groupexclusion(i)
    } else {
      x$operation_Groupsubset(i)
    }
  } else if (type == "sample") {
    if (exclude) {
      x$operation_Samplesexclusion(i)
    } else {
      x$operation_Samplessubset(i)
    }
  }
  
  invisible(x)
}

# design <- object_experiment_design$new("GroupSamples.txt", "GroupContrast.txt")
# design[c("A", "B"), type = "group"]                # 保留组 A 和 B
# design[c("C", "D"), type = "group", exclude = TRUE] # 排除组 C 和 D
# design[c("S1", "S2"), type = "sample"]             # 保留样本 S1, S2
# design[c("S3"), type = "sample", exclude = TRUE]   # 排除样本 S3

# 自定义中缀运算符

`%keep_groups%` <- function(design, groups) {
  design$operation_Groupsubset(groups)
  invisible(design)
}

`%exclude_groups%` <- function(design, groups) {
  design$operation_Groupexclusion(groups)
  invisible(design)
}

`%keep_samples%` <- function(design, samples) {
  design$operation_Samplessubset(samples)
  invisible(design)
}

`%exclude_samples%` <- function(design, samples) {
  design$operation_Samplesexclusion(samples)
  invisible(design)
}

# design <- object_experiment_design$new(...)
# design %keep_groups% c("Group1", "Group2")
# design %exclude_groups% c("Group3")
# design %keep_samples% c("Sample_ctrl", "Sample_exp")
# design %exclude_samples% c("SampleC")
# 
# # 支持链式操作
# design %keep_groups% c("G1","G2") %exclude_samples% c("bad1","bad2")

# subset 泛型方法
subset.object_experiment_design <- function(x, groups = NULL, samples = NULL, 
                                            exclude_groups = NULL, exclude_samples = NULL, ...) {
  if (!is.null(groups)) x$operation_Groupsubset(groups)
  if (!is.null(samples)) x$operation_Samplessubset(samples)
  if (!is.null(exclude_groups)) x$operation_Groupexclusion(exclude_groups)
  if (!is.null(exclude_samples)) x$operation_Samplesexclusion(exclude_samples)
  invisible(x)
}

# design <- object_experiment_design$new(...)
# subset(design, groups = c("A","B"))                    # 保留组
# subset(design, exclude_groups = "C")                  # 排除组
# subset(design, samples = c("S1","S2"))                # 保留样本
# subset(design, exclude_samples = "bad")               # 排除样本
# # 组合使用
# subset(design, groups = c("A","B"), exclude_samples = "badSample")

# ==================== PublicDataManager 公共数据解析 基础数据预处理和储存类 ====================
#' 公共数据管理器：负责 GTF、KEGG、GO 数据库的本地缓存、网络测试与解析
#' 
#' 提供以下核心功能：
#' 1. 根据物种代码（KEGG 三字母码）自动下载/缓存 KEGG 通路信息
#' 2. 解析 GTF 文件（支持 NCBI GTF、OrgDb、AnnotationHub）提取基因注释（SYMBOL, ENTREZID, GENENAME）
#' 3. 解析 GO 注释（支持 OrgDb、NCBI 的 gene2go 文件、AnnotationHub）
#' 4. 生成 clusterProfiler 所需的 TERM2GENE 和 TERM2NAME 背景数据框
#' 5. 测试网络连通性（KEGG REST API）
#' 
#' @field kegg_organism_code KEGG 物种代码，如 "hsa" 人, "mmu" 小鼠
#' @field cache_dir 缓存目录，默认为 "~/.cache/PublicDataManager"
#' @field force_update 是否强制重新下载（忽略已有缓存）
PublicDataManager <- R6::R6Class(
  classname = "PublicDataManager",
  private = list(
    .kegg_code = NULL,
    .cache_dir = NULL,
    .force_update = FALSE,
    
    # 缓存数据存储
    .gtf_tbl = NULL,          # 基因  注释表  (Geneid, SYMBOL, ENTREZID, GENENAME)
    # 表格信息解释
    # Geneid    featureCounts 主键
    # SYMBOL    基因的SYMBOL英文缩写
    # ENTREZID  NCBI的基因ENTREZID
    # GENENAME  NCBI的基因description全称
    
    .go_tbl = NULL,           # GO    注释表  (SYMBOL, ENTREZID, GO)
    # 表格信息解释
    # SYMBOL    基因的SYMBOL英文缩写
    # ENTREZID  NCBI的基因ENTREZID
    # GO        GO数据库的GO：XXXX 信息
    
    
    .kegg_tbl = NULL,         # KEGG  注释表  (Keggpathid, ENTREZID, Keggpathname)
    # 表格信息解释
    # Keggpathid    KEGG数据库 koXXXXX 数据ID
    # ENTREZID      NCBI的基因ENTREZID
    # Keggpathname  KEGG数据库 koXXXXX 全称描述信息
    
    
    .orgdb = NULL,            # OrgDb 对象    (如 org.Hs.eg.db)
    # Bioconduct 模式生物注释信息
    
    # ========== 1. 私有函数 KEGG 数据处理 ==========
    # 下载/缓存 KEGG 通路信息，返回数据框 (Keggpathid, ENTREZID, Keggpathname)
    load_kegg = function(use_cache = TRUE) {
      if (use_cache) {
        kegg <- private$read_cache("kegg")
        if (!is.null(kegg)) {
          private$.kegg_tbl <- kegg
          return(invisible(kegg))
        }
      }
      
      message("正在从 KEGG 下载 ", private$.kegg_code, " 通路数据...")
      .require_bioc_package("clusterProfiler", load = FALSE)
      kegg_db <- tryCatch({
        clusterProfiler::download_KEGG(species = private$.kegg_code)
      }, error = function(e) {
        stop("KEGG 下载失败: ", e$message)
      })
      
      # 整理格式
      path2ext <- kegg_db$KEGGPATHID2EXTID
      path2name <- kegg_db$KEGGPATHID2NAME
      colnames(path2ext) <- c("Keggpathid", "ENTREZID")
      colnames(path2name) <- c("Keggpathid", "Keggpathname")
      kegg_tbl <- dplyr::left_join(path2ext, path2name, by = "Keggpathid")
      kegg_tbl$Keggpathname <- gsub(" - .*", "", kegg_tbl$Keggpathname)
      
      private$.kegg_tbl <- kegg_tbl
      if (use_cache) private$write_cache(kegg_tbl, "kegg")
      message("KEGG 数据加载完成，包含 ", nrow(kegg_tbl), " 条记录")
      invisible(kegg_tbl)
    },
    
    # ========== 2. 私有函数 GTF 注释处理==========
    # 解析 NCBI 格式的 GTF 文件，提取 Geneid, SYMBOL, ENTREZID, GENENAME
    # @param gtf_file GTF 文件路径
    # @param use_cache 是否使用缓存
    parse_gtf_ncbi = function(gtf_file, use_cache = TRUE) {
      if (use_cache) {
        cached <- private$read_cache("gtf_ncbi")
        if (!is.null(cached)) {
          private$.gtf_tbl <- cached
          return(invisible(cached))
        }
      }
      if (!file.exists(gtf_file)) stop("GTF 文件不存在: ", gtf_file)
      
      # 预处理：删除非 GeneID 的 db_xref（可选，提升速度）
      modified_file <- paste0(gtf_file, ".modify")
      if (!file.exists(modified_file)) {
        cmd <- sprintf('perl -pe \'s/db_xref "(?!GeneID)[^"]+"; //g\' %s > %s', 
                       gtf_file, modified_file)
        system(cmd)
      }
      
      if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        .require_bioc_package("rtracklayer", load = FALSE)
      }
      gtf <- rtracklayer::import(modified_file, format = "gtf")
      
      # 提取注释信息
      gtf_df <- data.frame(
        gene_id = gtf$gene_id,
        gene_name = gtf$gene,
        db_xref = gtf$db_xref,
        description = if ("description" %in% names(gtf)) gtf$description else gtf$product,
        stringsAsFactors = FALSE
      )
      
      # 聚合到基因级别
      library(dplyr)
      gtf_tbl <- gtf_df %>%
        dplyr::distinct() %>%
        dplyr::group_by(gene_id, gene_name) %>%
        dplyr::summarise(
          ENTREZID = {
            x <- db_xref[grepl("GeneID", db_xref)]
            if (length(x) == 0) NA else gsub("GeneID:(.*)", "\\1", x[1])
          },
          Description = {
            d <- description[!is.na(description)]
            if (length(d) == 0) NA else strsplit(d[1], split = ", ")[[1]][1]
          },
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Geneid = gene_id,
          SYMBOL = ifelse(is.na(gene_name), gene_id, gene_name),
          GENENAME = Description
        ) %>%
        dplyr::select(Geneid, SYMBOL, ENTREZID, GENENAME)
      
      private$.gtf_tbl <- gtf_tbl
      if (use_cache) private$write_cache(gtf_tbl, "gtf_ncbi")
      message("GTF 解析完成，共 ", nrow(gtf_tbl), " 个基因")
      invisible(gtf_tbl)
    },
    
    # 在 PublicDataManager 类的 public 列表中添加以下方法
    
    # 解析 Ensembl 格式的 GTF 文件，提取基因注释
    # @param gtf_file GTF 文件路径
    # @param use_cache 是否使用缓存
    parse_gtf_ensembl = function(gtf_file, use_cache = TRUE) {
      if (use_cache) {
        cached <- private$read_cache("gtf_ensembl")
        if (!is.null(cached)) {
          private$.gtf_tbl <- cached
          return(invisible(cached))
        }
      }
      if (!file.exists(gtf_file)) stop("GTF 文件不存在: ", gtf_file)
      
      if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        .require_bioc_package("rtracklayer", load = FALSE)
      }
      
      # 导入 GTF 文件
      gtf <- rtracklayer::import(gtf_file, format = "gtf")
      
      # 只保留 gene 类型的记录（每个基因一条）
      gtf_gene <- gtf[gtf$type == "gene"]
      
      # 提取基本信息
      gene_ids <- mcols(gtf_gene)$gene_id
      gene_names <- mcols(gtf_gene)$gene_name
      # 若有 description 字段，提取；否则用 gene_name 代替
      descriptions <- if ("description" %in% names(mcols(gtf_gene))) {
        mcols(gtf_gene)$description
      } else {
        gene_names
      }
      
      # 尝试提取 ENTREZID（可能在 db_xref 或 entrez_id 中）
      entrez_ids <- rep(NA_character_, length(gtf_gene))
      if ("db_xref" %in% names(mcols(gtf_gene))) {
        db_xref <- mcols(gtf_gene)$db_xref
        # 查找 "GeneID:" 模式
        entrez_ids <- sapply(db_xref, function(x) {
          if (is.na(x) || x == "") return(NA_character_)
          match <- regmatches(x, regexpr("GeneID:[0-9]+", x))
          if (length(match) > 0) gsub("GeneID:", "", match[1]) else NA_character_
        })
      } else if ("entrez_id" %in% names(mcols(gtf_gene))) {
        entrez_ids <- mcols(gtf_gene)$entrez_id
      }
      
      # 构建数据框
      gtf_tbl <- data.frame(
        Geneid = gene_ids,
        SYMBOL = ifelse(is.na(gene_names), gene_ids, gene_names),
        ENTREZID = entrez_ids,
        GENENAME = descriptions,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::distinct(Geneid, .keep_all = TRUE)  # 确保基因 ID 唯一
      
      private$.gtf_tbl <- gtf_tbl
      if (use_cache) private$write_cache(gtf_tbl, "gtf_ensembl")
      message(sprintf("[%s] 解析 Ensembl GTF 完成，共 %d 个基因", Sys.time(), nrow(gtf_tbl)))
      return(gtf_tbl)
    },
    
    # 基于 OrgDb Bioconductor注释函数包获取基因注释（无需 GTF 文件）
    # @param use_cache 是否使用缓存
    load_gtf_from_orgdb = function(use_cache = TRUE) {
      if (use_cache) {
        cached <- private$read_cache("gtf_orgdb")
        if (!is.null(cached)) {
          private$.gtf_tbl <- cached
          return(invisible(cached))
        }
      }
      .require_bioc_package("AnnotationDbi", load = FALSE)
      gtf_tbl <- AnnotationDbi::select(
        private$.orgdb,
        keys = AnnotationDbi::keys(private$.orgdb, keytype = "SYMBOL"),
        columns = c("ENTREZID", "GENENAME"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::rename(SYMBOL = !!sym("SYMBOL")) %>%
        dplyr::mutate(Geneid = SYMBOL) %>%
        dplyr::select(Geneid, SYMBOL, ENTREZID, GENENAME) %>%
        dplyr::distinct()
      
      private$.gtf_tbl <- gtf_tbl
      if (use_cache) private$write_cache(gtf_tbl, "gtf_orgdb")
      message("从 OrgDb 加载基因注释，共 ", nrow(gtf_tbl), " 个基因")
      invisible(gtf_tbl)
    },
    
    # 通过 NCBI 表合并 ENTREZID（通用，支持 Geneid / SYMBOL 匹配）
    supplement_entrez_by_orgdb = function(gtf_tbl) {
      if (is.null(private$.orgdb)) stop("OrgDb 未加载")
      .require_bioc_package("AnnotationDbi", load = FALSE)
      message("通过 OrgDb 补充 ENTREZID...")
      # 获取 SYMBOL -> ENTREZID 映射
      org_map <- AnnotationDbi::select(
        private$.orgdb,
        keys = AnnotationDbi::keys(private$.orgdb, keytype = "SYMBOL"),
        columns = "ENTREZID",
        keytype = "SYMBOL"
      ) %>% dplyr::rename(ENTREZID = ENTREZID) %>% dplyr::distinct() |> 
        dplyr::mutate(SYMBOL = toupper(SYMBOL))
      # 确保映射表列名匹配 merge_entrez_from_ncbi 要求
      colnames(org_map)[colnames(org_map) == "SYMBOL"] <- "SYMBOL"
      # 调用通用合并函数（匹配 SYMBOL 字段）
      private$merge_entrez_from_ncbi(gtf_tbl, org_map, match_by = "SYMBOL")
    },
    
    # 用 gene2ensembl 文件补充 ENTREZID（通过 Ensembl Gene ID 匹配）
    supplement_entrez_by_gene2ensembl = function(gtf_tbl, gene2ensembl) {
      message("通过 gene2ensembl 文件补充 ENTREZID...")
      # gene2ensembl 应包含 GeneID (ENTREZID) 和 Ensembl_gene_identifier
      # 重命名使匹配字段为 Geneid
      map_tbl <- gene2ensembl %>%
        dplyr::select(ENTREZID, Geneid = Ensembl_gene_identifier) %>%
        dplyr::distinct()
      private$merge_entrez_from_ncbi(gtf_tbl, map_tbl, match_by = "Geneid")
    },
    
    # ========== 3. 私有函数 GO 注释处理==========
    # 加载 GO 注释（优先使用 OrgDb，备选 gene2go 文件）
    # @param gene2go_file NCBI gene2go 文件路径（可选）
    load_go_org = function(use_cache = TRUE){
      .require_bioc_package("AnnotationDbi", load = FALSE)
      go_tbl <- AnnotationDbi::select(
        private$.orgdb,
        keys = AnnotationDbi::keys(private$.orgdb, keytype = "SYMBOL"),
        columns = c("ENTREZID", "GO"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::select(SYMBOL, ENTREZID, GO) %>%
        dplyr::filter(!is.na(GO)) %>%
        dplyr::distinct()
      private$.go_tbl <- go_tbl
      if (use_cache) private$write_cache(go_tbl, "go")
      message("GO 注释加载完成，共 ", nrow(go_tbl), " 条记录")
      return(go_tbl)
    },
    load_go_gene2go.gz = function(gene2go.gz = NULL, use_cache = TRUE) {
      if (is.null(gene2go.gz) || !file.exists(gene2go.gz)) {
        stop("gene2go.gz 文件不存在: ", gene2go.gz)
      }
      rds_cache <- gsub(pattern = "\\.gz$", replacement = ".RDS", x = gene2go.gz)
      if (use_cache && file.exists(rds_cache)) {
        go_tbl = readRDS(gsub(pattern = "\\.gz", replacement = ".RDS", x = gene2go.gz)) |> 
          dplyr::as_tibble()
      }else{
        go_tbl <- data.table::fread(gene2go.gz, data.table = FALSE) %>%
          dplyr::select(SYMBOL = Symbol, ENTREZID = GeneID, GO = GO_ID) %>%
          dplyr::filter(!is.na(GO)) %>%
          dplyr::distinct()
      }
      private$.go_tbl <- go_tbl
      message("GO 注释加载完成，共 ", nrow(go_tbl), " 条记录")
      return(go_tbl)
    },
    
    load_gene2ensembl.gz = function(gene2ensembl.gz = NULL, use_cache = TRUE) {
      
      
      # 检查文件存在性
      if (is.null(gene2ensembl.gz) || !file.exists(gene2ensembl.gz)) 
        stop("gene2ensembl.gz 文件不存在: ", gene2ensembl.gz)
      
      # 读取并整理为两列：ENTREZID (GeneID) 和 Ensembl_gene_identifier
      # 尝试从缓存读取
      rds_cache <- gsub(pattern = "\\.gz$", replacement = ".RDS", x = gene2ensembl.gz)
      if (use_cache && file.exists(rds_cache)) {
        gene2ensembl_tbl = readRDS(gsub(pattern = "\\.gz", replacement = ".RDS", x = gene2ensembl.gz)) |> 
          dplyr::as_tibble()
      }else{
        gene2ensembl_tbl <- data.table::fread(gene2ensembl.gz, data.table = FALSE) %>%
          dplyr::select(ENTREZID = GeneID, Ensembl_gene_identifier) %>%
          dplyr::distinct()
      }
      message("gene2ensembl 映射加载完成，共 ", nrow(gene2ensembl_tbl), " 条记录")
      return(gene2ensembl_tbl)
    },
    merge_entrez_from_ncbi = function(ensembl_tbl, ncbi_tbl, match_by = c("Geneid", "SYMBOL")) {
      merged <- ensembl_tbl
      has_entrez <- !is.na(merged$ENTREZID)
      message("Ensembl 表中已有 ENTREZID 的基因数: ", sum(has_entrez))
      for (field in match_by) {
        if (field %in% colnames(merged) && field %in% colnames(ncbi_tbl)) {
          idx_need <- which(is.na(merged$ENTREZID) & !is.na(merged[[field]]))
          if (length(idx_need) == 0) break
          ncbi_sub <- ncbi_tbl[!is.na(ncbi_tbl[[field]]), c(field, "ENTREZID")]
          ncbi_sub <- ncbi_sub[!duplicated(ncbi_sub[[field]]), ]
          matched <- match(merged[[field]][idx_need], ncbi_sub[[field]])
          new_entrez <- ncbi_sub$ENTREZID[matched]
          merged$ENTREZID[idx_need] <- ifelse(is.na(merged$ENTREZID[idx_need]), new_entrez, merged$ENTREZID[idx_need])
          message("通过 ", field, " 匹配补充了 ", sum(!is.na(new_entrez)), " 个 ENTREZID")
        }
      }
      final_has <- sum(!is.na(merged$ENTREZID))
      message("合并后共有 ", final_has, " 个基因具有 ENTREZID（占比 ", round(100 * final_has / nrow(merged), 1), "%）")
      return(merged)
    },
    
    # 辅助函数：创建缓存目录
    ensure_cache_dir = function() {
      if (!dir.exists(private$.cache_dir)) {
        dir.create(private$.cache_dir, recursive = TRUE)
        message("创建缓存目录: ", private$.cache_dir)
      }
    },
    
    # 内部缓存文件路径
    cache_path = function(name) {
      file.path(private$.cache_dir, paste0(private$.kegg_code, "_", name, ".rds"))
    },
    
    # 读取或写入 RDS 缓存
    read_cache = function(name) {
      path <- private$cache_path(name)
      if (file.exists(path) && !private$.force_update) {
        message("从缓存加载 ", name, ": ", path)
        return(readRDS(path))
      }
      return(NULL)
    },
    
    write_cache = function(obj, name) {
      path <- private$cache_path(name)
      saveRDS(obj, path)
      message("写入缓存 ", name, ": ", path)
    },
    
    # 测试网络连通性（KEGG）
    test_kegg_connection = function() {
      # 使用更合适的 /info/kegg 端点进行连接测试
      test_url <- "https://rest.kegg.jp/info/kegg"
      # 设置一个较短超时，避免长时间等待
      res <- tryCatch({
        # 使用 readLines 并设置超时来测试连接
        con <- url(test_url, open = "rb")
        on.exit(close(con), add = TRUE)
        # 尝试读取少量数据来确认连接
        status <- suppressWarnings(readLines(con, n = 1, warn = FALSE))
        # 如果能读到以 "kegg" 开头的内容，表明成功
        if (length(status) > 0 && grepl("^kegg", status, ignore.case = TRUE)) {
          TRUE
        } else {
          FALSE
        }
      }, error = function(e) {
        # 连接出错时返回 FALSE
        FALSE
      })
      
      if (res) {
        message("KEGG API 连接成功。")
      } else {
        warning("无法连接到 KEGG API (", test_url, ")，请检查网络或使用本地缓存。")
      }
      return(res)
    },
    
    # 获取 OrgDb（基于 KEGG 代码）
    get_orgdb = function() {
      mapping <- list(
        hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db", rno = "org.Rn.eg.db",
        dre = "org.Dr.eg.db", dme = "org.Dm.eg.db", cel = "org.Ce.eg.db",
        sce = "org.Sc.sgd.db", ath = "org.At.tair.db", eco = "org.EcK12.eg.db",
        bsu = "org.Bs.eg.db", pfa = "org.Pf.plasmo.db", gga = "org.Gg.eg.db"
      )
      if (is.null(private$.kegg_code) || !private$.kegg_code %in% names(mapping)) {
        stop("不支持的 KEGG 物种代码: ", private$.kegg_code)
      }
      pkg <- mapping[[private$.kegg_code]]
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("正在安装 ", pkg, " ...")
        .require_bioc_package(pkg, load = FALSE)
      }
      library(pkg, character.only = TRUE)
      return(get(pkg))
    }
  ),
  public = list(
    # ========== 1. 初始化函数 ==========
    # 构造函数
    # @param kegg_organism_code KEGG 物种代码（如 "hsa"）
    # @param cache_dir 缓存目录
    # @param force_update 是否强制重新下载
    initialize = function(kegg_organism_code, cache_dir = NULL, force_update = FALSE) {
      if (missing(kegg_organism_code)) stop("必须提供 kegg_organism_code")
      private$.kegg_code <- kegg_organism_code
      private$.force_update <- force_update
      if (is.null(cache_dir)) {
        private$.cache_dir <- file.path("~/.cache/PublicDataManager", kegg_organism_code)
      } else {
        private$.cache_dir <- cache_dir
      }
      private$ensure_cache_dir()
      
      # 网络测试
      if (!private$test_kegg_connection()) {
        warning("KEGG 网络不可达，将仅使用本地缓存数据")
      }
    },
    
    # ========== 2. 公共接口 构造函数 ==========
    # 支持模式不同输入选择模式
    # 配置RDS缓存
    # 模式1：使用 NCBI         参考基因组 模式生物   输入形参：NCBI GTF文件 Orgdb 人、小鼠模式生物
    # 模式2：使用 NCBI         参考基因组 非模式生物 输入形参：GTF文件 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    # ENTREZID 配合 gene2go.gz 映射 GO:XXXX 等信息
    
    # 模式3：使用 Ensembl      参考基因组 模式生物   输入形参：NCBI GTF文件 
    
    # 模式4：使用 Ensembl      参考基因组 模式生物和非模式生物   输入形参：ENSEMBL GTF文件、NCBI GTF文件、https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz、https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    # 最麻烦，需要进行 ID 转换，因为 KEGG 需要 ENTREZID 进行输入
    # ENTREZID 配合 gene2go.gz 映射 GO:XXXX 等信息
    
    # 模式5：使用 三方基因组   模式生物和非模式生物 输入形参：三方参考基因组GTF文件、eggnog 综合注释文件 GTF KEGG GO 注释
    
    # @param gtf_file_ncbi    第一个 NCBI GTF 文件路径
    # @param gtf_file_ensembl 第二个 Ensembl GTF 文件路径（可选，若提供则视为 NCBI）
    # @param use_orgdb        是否使用模式生物 Bioconductor 内置注释函数包补充 ENTREZID
    # @param gene2go.gz       是否使用 NCBI Gene go 补充 goid
    # @param gene2ensembl.gz  是否使用 NCBI Gene ensemble 补充 ENTREZID
    # @param use_cache        是否使用缓存
    # @param match_by         字段匹配顺序
    load_combined_gtf = function(gtf_file_ncbi = NULL, 
                                 gtf_file_ensembl = NULL, 
                                 # Orgdb 已经提前内置
                                 use_orgdb  = FALSE,
                                 gene2go.gz = NULL,
                                 gene2ensembl.gz = NULL, 
                                 eggnog = NULL,
                                 use_cache = TRUE, 
                                 match_by = c("Geneid", "SYMBOL")) {
      
      # KEGG 解析，使用 缓存数据存储 特性
      # .kegg_tbl = NULL,         # KEGG 注释表 (Keggpathid, ENTREZID, Keggpathname)
      private$load_kegg(use_cache = use_cache)
      
      #  解析其它，使用 缓存数据存储 特性
      # .gtf_tbl = NULL,          # 基因注释表  (Geneid, SYMBOL, ENTREZID, GENENAME)
      # .go_tbl = NULL,           # GO 注释表   (SYMBOL, ENTREZID, GO)
      # .orgdb = NULL,            # OrgDb 对象  (如 org.Hs.eg.db)
      
      # 模式1：使用 NCBI 参考基因组 模式生物 输入形参：GTF文件 Orgdb 人、小鼠
      if (
        !is.null(gtf_file_ncbi) && 
        is.null(gtf_file_ensembl) && 
        use_orgdb && 
        is.null(gene2go.gz) && 
        is.null(gene2ensembl.gz) && 
        is.null(eggnog)
      ) {
        message("模式1：使用 NCBI 参考基因组 模式生物 输入形参：GTF文件 Orgdb 人、小鼠")
        private$.gtf_tbl = private$parse_gtf_ncbi(gtf_file_ncbi, use_cache = use_cache)
        private$.go_tbl  = private$load_go_org()
      }
      # 模式2：使用 NCBI 参考基因组 非模式生物 输入形参：GTF文件 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
      else if (
        !is.null(gtf_file_ncbi) && 
        is.null(gtf_file_ensembl) && 
        !use_orgdb && 
        !is.null(gene2go.gz) && 
        is.null(gene2ensembl.gz) && 
        is.null(eggnog)
      ) {
        message("模式2：使用 NCBI 参考基因组 非模式生物 输入形参：GTF文件 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
        private$.gtf_tbl = private$parse_gtf_ncbi(gtf_file_ncbi, use_cache = use_cache)
        private$.go_tbl  = private$load_go_gene2go.gz(gene2go.gz = gene2go.gz, use_cache = use_cache)
      }
      # 模式4：使用 Ensembl 参考基因组 模式生物和非模式生物   输入形参：GTF文件 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
      else if (
        !is.null(gtf_file_ncbi) && 
        !is.null(gtf_file_ensembl) && 
        use_orgdb && 
        !is.null(gene2go.gz) && 
        !is.null(gene2ensembl.gz) && 
        is.null(eggnog)
      ) {
        message("模式4：使用 Ensembl 参考基因组 模式生物和非模式生物   输入形参：GTF文件 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
        # 读取 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
        gene2ensembl = private$load_gene2ensembl.gz(gene2ensembl.gz = gene2ensembl.gz)
        # 读取 https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
        gene2go =  private$load_go_gene2go.gz(gene2go.gz = gene2go.gz)
        
        private$.gtf_tbl <- private$parse_gtf_ensembl(gtf_file_ensembl, use_cache = use_cache) |> 
          # 通过 gene2ensembl.gz 补充添加 ENTREZID 信息
          private$supplement_entrez_by_gene2ensembl(gene2ensembl = gene2ensembl)
        
        # 修改 GENENAME 列信息，使用 gtf_file_ncbi 补充基因的全称
        private$.gtf_tbl$GENENAME = NULL
        ncbi_gtf_tbl = private$parse_gtf_ncbi(gtf_file_ncbi, use_cache = use_cache)
        ncbi_gtf_tbl = ncbi_gtf_tbl[,c("ENTREZID", "GENENAME")]
        
        private$.gtf_tbl = dplyr::left_join(private$.gtf_tbl, ncbi_gtf_tbl, by = c("ENTREZID" = "ENTREZID"))
        
        private$.go_tbl <- gene2go[gene2go$ENTREZID %in% private$.gtf_tbl$ENTREZID, ]
      }
      # 模式5：使用 三方基因组 模式生物和非模式生物 输入形参：三方参考基因组GTF文件、eggnog 综合注释文件 GTF KEGG GO 注释
      else if(
        is.null(gtf_file_ncbi) && 
        !is.null(gtf_file_ensembl) && 
        !use_orgdb && 
        !is.null(gene2go.gz) && 
        !is.null(gene2ensembl.gz) && 
        is.null(eggnog)
      ){
        message("模式5：使用 三方基因组 模式生物和非模式生物 输入形参：三方参考基因组GTF文件、eggnog 综合注释文件 GTF KEGG GO 注释")
        stop("")
      }
      else{
        stop("形参配置错误")
      }
    },
    
    # 使用 OrgDb 作为注释来源，适合模式生物且没有 GTF 文件时快速构建注释和 GO 背景。
    load_annotation_from_orgdb = function(use_cache = TRUE,
                                          load_go = TRUE,
                                          load_kegg = TRUE) {
      if (is.null(private$.orgdb)) {
        private$.orgdb <- private$get_orgdb()
      }
      private$.gtf_tbl <- private$load_gtf_from_orgdb(use_cache = use_cache)
      if (load_go) {
        private$.go_tbl <- private$load_go_org(use_cache = use_cache)
      }
      if (load_kegg) {
        private$load_kegg(use_cache = use_cache)
      }
      invisible(self)
    },
    
    # 生成与 read_counts/object_experiment_design 行顺序一致的注释表。
    prepare_rowData_annotation = function(read_obj,
                                          match_by = c("Geneid", "SYMBOL"),
                                          annotation_tbl = NULL) {
      if (!inherits(read_obj, "read_counts")) {
        stop("read_obj 必须是 read_counts 或其子类对象")
      }
      if (is.null(annotation_tbl)) {
        annotation_tbl <- self$get_gene_annotation()
      }
      annotation_tbl <- as.data.frame(annotation_tbl, stringsAsFactors = FALSE)
      gene_ids <- rownames(read_obj$get_assay("counts"))
      match_by <- match_by[match_by %in% colnames(annotation_tbl)]
      if (length(match_by) == 0) {
        stop("annotation_tbl 不包含可匹配字段: Geneid / SYMBOL")
      }
      
      overlap_n <- vapply(match_by, function(field) {
        sum(gene_ids %in% annotation_tbl[[field]], na.rm = TRUE)
      }, numeric(1))
      best_field <- match_by[which.max(overlap_n)]
      if (max(overlap_n) == 0) {
        stop("注释表与表达矩阵行名无法匹配，请检查 Geneid/SYMBOL 体系")
      }
      
      annotation_tbl <- annotation_tbl[!duplicated(annotation_tbl[[best_field]]), , drop = FALSE]
      idx <- match(gene_ids, annotation_tbl[[best_field]])
      aligned <- annotation_tbl[idx, , drop = FALSE]
      aligned$Geneid_input <- gene_ids
      rownames(aligned) <- gene_ids
      aligned
    },
    
    # 将 PublicDataManager 中的基因注释补充到 read_counts 的 rowData。
    annotate_rowData = function(read_obj,
                                match_by = c("Geneid", "SYMBOL"),
                                overwrite = FALSE) {
      row_anno <- self$prepare_rowData_annotation(read_obj = read_obj, match_by = match_by)
      read_obj$add_rowData(row_anno, by = "Geneid_input", overwrite = overwrite)
      invisible(read_obj)
    },
    
    # 生成背景基因 universe，用于 ORA/GSEA 的 universe 或结果追踪。
    get_background_genes = function(read_obj = NULL,
                                    id_col = c("ENTREZID", "SYMBOL", "Geneid"),
                                    annotated_only = TRUE) {
      id_col <- match.arg(id_col)
      if (!is.null(read_obj)) {
        if (!inherits(read_obj, "read_counts")) {
          stop("read_obj 必须是 read_counts 或其子类对象")
        }
        row_tbl <- read_obj$get_rowData_table(add_Geneid = TRUE)
        if (!id_col %in% colnames(row_tbl)) {
          row_tbl <- self$prepare_rowData_annotation(read_obj = read_obj)
        }
      } else {
        row_tbl <- self$get_gene_annotation()
      }
      
      if (!id_col %in% colnames(row_tbl)) {
        stop("背景基因表中不存在字段: ", id_col)
      }
      genes <- as.character(row_tbl[[id_col]])
      if (annotated_only) {
        genes <- genes[!is.na(genes) & genes != ""]
      }
      unique(genes)
    },
    
    # 统一富集背景接口，返回 clusterProfiler::enricher/GSEA 可用的 TERM2GENE/TERM2NAME。
    get_enrichment_background = function(db = c("KEGG", "GO"), ontology = "BP") {
      db <- match.arg(db)
      if (db == "KEGG") {
        self$get_kegg_background()
      } else {
        self$get_go_background(ontology = ontology)
      }
    },
    
    # 检查当前运行环境中的 CRAN/Bioconductor 包状态。
    check_package_environment = function(install_missing = FALSE,
                                         include_optional = TRUE) {
      check_transcriptome_package_environment(
        install_missing = install_missing,
        include_optional = include_optional
      )
    },
    
    # ========== 3. 获取富集分析背景数据获取 ==========
    # 生成 KEGG 背景：TERM2GENE 和 TERM2NAME
    get_kegg_background = function() {
      if (is.null(private$.kegg_tbl)) {
        private$load_kegg()
      }
      kegg <- private$.kegg_tbl
      TERM2GENE <- data.frame(TERM = kegg$Keggpathid, GENE = kegg$ENTREZID) %>% na.omit()
      TERM2NAME <- data.frame(TERM = kegg$Keggpathid, NAME = kegg$Keggpathname) %>% dplyr::distinct()
      return(list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME))
    },
    
    # 生成 GO 背景（仅 BP 分支）
    get_go_background = function(ontology = "BP") {
      if (!requireNamespace("GO.db", quietly = TRUE)) {
        .require_bioc_package("GO.db", load = FALSE)
      }
      .require_bioc_package("AnnotationDbi", load = FALSE)
      # GO.db 数据库获取 GO:XXXX 和 Ontology
      go_all <- AnnotationDbi::select(GO.db::GO.db, 
                                      keys = AnnotationDbi::keys(GO.db::GO.db),
                                      columns = c("TERM", "ONTOLOGY"))
      go_bp <- go_all[go_all$ONTOLOGY == ontology, ]
      TERM2NAME <- go_bp[, c("GOID", "TERM")]
      colnames(TERM2NAME) <- c("TERM", "NAME")
      
     
      if (is.null(private$.go_tbl)) {
        stop("GO注释尚未加载")
      }
      # 只保留 ontology = "BP" 对应的基因关联
      go_tbl <- private$.go_tbl
      TERM2GENE <- data.frame(TERM = go_tbl$GO, GENE = go_tbl$ENTREZID) %>%
        dplyr::filter(TERM %in% TERM2NAME$TERM) %>%
        na.omit() %>%
        distinct()
      
      return(list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME))
    },
    
    # 获取基因注释表（供差异表达结果注释使用）
    get_gene_annotation = function() {
      if (is.null(private$.gtf_tbl)) {
        stop("基因注释尚未加载")
      }
      return(private$.gtf_tbl)
    },
    
    # ========== 4. 工具方法 ==========
    # 测试网络连通性（对外暴露）
    check_connection = function() {
      private$test_kegg_connection()
    },
    
    # 清空所有缓存
    clear_cache = function() {
      files <- list.files(private$.cache_dir, pattern = "\\.rds$", full.names = TRUE)
      if (length(files) > 0) {
        unlink(files)
        message("已清除缓存目录 ", private$.cache_dir, " 中的所有 .rds 文件")
      } else {
        message("没有缓存文件")
      }
    },
    
    # 打印对象信息
    print = function() {
      cat("<PublicDataManager>\n")
      cat("  KEGG species: ", private$.kegg_code, "\n")
      cat("  Cache dir: ", private$.cache_dir, "\n")
      cat("  Force update: ", private$.force_update, "\n")
      cat("  Loaded data:\n")
      cat("    - Orgdb: ", ifelse(is.null(private$.orgdb), "no", "yes"), "\n")
      cat("    - GTF: ", ifelse(is.null(private$.gtf_tbl), "no", "yes"), "\n")
      cat("    - GO:  ", ifelse(is.null(private$.go_tbl), "no", "yes"), "\n")
      cat("    - KEGG:", ifelse(is.null(private$.kegg_tbl), "no", "yes"), "\n")
    }
  )
)

# ==================== run_project_visualization ======================
#' 面向过程的样本可视化函数
#'
#' 该函数替代 ProjectVisualizer R6 类，根据提供的实验设计对象生成多种样本级别可视化图形。
#'
#' @param design_obj object_experiment_design 对象（需包含表达矩阵和分组信息）
#' @param type 可视化类型，可选 "boxplot", "pca", "heatmap", "dendro" 或 "all"（默认）
#' @param assay_name 使用的表达矩阵名称（如 "tpm", "fpkm", "counts"），默认为 "tpm"
#' @param out_dir 输出目录，若为 NULL 则不保存图形，默认为 "03.Samples_统计结果"
#' @param ... 其他参数传递给具体的绘图函数（如 log2_transform, scale, method 等）
#'
#' @return 若 type != "all"，返回对应的 ggplot 或 pheatmap 对象；若 type = "all"，返回包含所有图形的列表
#'
#' @examples
#' \dontrun{
#' res <- create_from_files(...)
#' run_visualization(res$design, type = "pca", assay_name = "tpm")
#' run_visualization(res$design, type = "all", out_dir = "QC")
#' }
#' 
run_project_visualization <- function(design_obj,
                                      type = c("boxplot", "pca", "heatmap", "dendro", "normalized_matrix"),
                                      assay_name = "tpm",
                                      out_dir = "03.Samples_统计结果",
                                      ...) {
  
  # 参数匹配
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
  
  # 获取表达矩阵（自动计算 TPM/FPKM 若所需 assay 不存在）
  .get_expr_matrix <- function(design_obj, assay_name) {
    if (assay_name %in% names(design_obj$assays)) {
      mat <- design_obj$get_assay(assay_name)
    } else if (assay_name == "tpm") {
      message("TPM 矩阵不存在，正在计算...")
      mat <- design_obj$get_assay("tpm")
    } else if (assay_name == "fpkm") {
      message("FPKM 矩阵不存在，正在计算...")
      mat <- design_obj$get_assay("fpkm")
    } else if (assay_name == "counts") {
      mat <- design_obj$get_assay("counts")
    } else {
      stop("表达矩阵 ", assay_name, " 不存在，且无法自动计算")
    }
    return(mat)
  }
  
  # 获取分组向量（样本名 -> 组名）
  .get_group_vec <- function(design_obj, samples) {
    group_vec <- design_obj$GroupSamples[, 2]
    names(group_vec) <- design_obj$GroupSamples[, 1]
    return(group_vec[samples])
  }
  
  # 1. 箱线图
  .plot_boxplot <- function(expr_mat, group_vec, log2_transform = TRUE, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装 ggplot2")
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("需要安装 reshape2")
    
    mat <- expr_mat
    if (log2_transform) mat <- log2(mat + 1)
    samples <- colnames(mat)
    group <- group_vec[samples]
    
    df <- reshape2::melt(mat)
    colnames(df) <- c("Gene", "Sample", "Expression")
    df$Group <- group[as.character(df$Sample)]
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Expression, fill = Group)) +
      ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
      ggplot2::labs(x = "Samples", 
                    y = ifelse(log2_transform, "log2(Expression + 1)", "Expression"),
                    title = paste("Sample Expression Boxplot (", assay_name, ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    if (!is.null(out_dir)) {
      file_prefix <- file.path(out_dir, paste0("Samples_boxplot_", assay_name))
      width <- max(6, ncol(mat) * 0.3)
      height <- ifelse("height" %in% names(list(...)), list(...)$height, 6)
      .save_plot(p, file_prefix, width = width, height = height)
    }
    return(p)
  }
  
  # 2. PCA 图
  .plot_pca <- function(expr_mat, group_vec, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装 ggplot2")
    if (!suppressMessages(require(FactoMineR))) stop("需要安装FactoMineR包")
    if (!suppressMessages(require(factoextra))) stop("需要安装factoextra包")
    
    # 转置，剔除方差变化为零的基因
    mat_t <- t(expr_mat)
    keep <- apply(mat_t, 2, function(x) length(unique(x)) > 1)
    if (sum(keep) == 0) stop("所有基因方差为零，无法进行 PCA")
    mat_t <- mat_t[, keep, drop = FALSE]
    
    data.pca <- FactoMineR::PCA(mat_t, graph = FALSE)
    
    # 3D PCA
    if (!suppressMessages(require(scatterplot3d))) stop("需要安装 scatterplot3d 包")
    
    pca.result = data.pca$ind$contrib
    
    Color_self = setNames(randomcoloR::distinctColorPalette(length(unique(group_vec))), unique(group_vec))
    
    # 计算每个主成分的方差解释比例
    explained_variance <- data.pca$eig[, 2]  # 获取解释比例（百分比）
    
    # 提取前3个主成分的解释比例
    explained_variance_pct <- explained_variance[1:3]
    explained_variance_pct <- round(explained_variance_pct, 2)  # 四舍五入为小数点后两位
    # PDF打印输出
    {
      width = 10
      height = 9
      pdf(stringr::str_glue("{out_dir}/Sample_PCA_3Dplot.pdf"), width = width, height = height)
      
      scatterplot3d(pca.result[, 1:3],
                    color = Color_self[group_vec],
                    pch = 16, angle = 30,
                    cex.symbols = 2,
                    box = T, type = "p",
                    main = "3D PCA Plot",
                    lty.hide = 2, lty.grid = 2,
                    xlab = paste0("Dim (", explained_variance_pct[1], "%)"),
                    ylab = paste0("Dim (", explained_variance_pct[2], "%)"),
                    zlab = paste0("Dim (", explained_variance_pct[3], "%)"),
      )
      legend("topright", unique(group_vec),
             fill = Color_self,
             horiz = F
      )
      
      dev.off()
    }
    message(stringr::str_glue("{Sys.time()} 输出图像 {out_dir}/Sample_PCA_3Dplot.pdf"))
    
    # 2D PCA
    p <- factoextra::fviz_pca_ind(data.pca,
                 col.ind = group_vec,
                 mean.point = F,
                 addEllipses = TRUE,
                 legend.title = "Groups",
                 ellipse.type = "confidence",
                 ellipse.level = 0.95
    ) +
      ggplot2::labs(x = paste0("PC1 (", round(data.pca$eig[1,2], 1), "%)"),
                    y = paste0("PC2 (", round(data.pca$eig[2,2], 1), "%)"),
                    title = paste("PCA of Samples (", assay_name, ")")) +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1, 
        panel.grid = element_blank()
      )
    
    if (!is.null(out_dir)) {
      file_prefix <- file.path(out_dir, paste0("Samples_pca_", assay_name))
      width <- ifelse("width" %in% names(list(...)), list(...)$width, 7)
      height <- ifelse("height" %in% names(list(...)), list(...)$height, 6)
      .save_plot(p, file_prefix, width = width, height = height)
    }
    return(p)
  }
  
  # 3. 相关性热图
  .plot_heatmap <- function(expr_mat, group_vec, method = "pearson",
                            cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE, ...) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) stop("需要安装 pheatmap")
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("需要安装 RColorBrewer")
    
    cor_mat <- cor(expr_mat, method = method, use = "pairwise.complete.obs")
    annotation <- data.frame(Group = group_vec[colnames(expr_mat)])
    rownames(annotation) <- colnames(expr_mat)
    
    p <- pheatmap::pheatmap(cor_mat,
                            annotation_col = annotation,
                            annotation_row = annotation,
                            color = colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(100),
                            cluster_rows = cluster_rows,
                            cluster_cols = cluster_cols,
                            display_numbers = display_numbers,
                            number_format = "%.2f",
                            main = paste("Sample Correlation (", method, ")"),
                            silent = TRUE,
                            ...)
    
    if (!is.null(out_dir)) {
      file_prefix <- file.path(out_dir, paste0("Samples_heatmap_", assay_name))
      width <- max(6, ncol(cor_mat) * 0.4)
      height <- width
      .save_plot(p$gtable, file_prefix, width = width, height = height)
    }
    return(p)
  }
  
  # 4. 聚类树状图
  .plot_dendro <- function(expr_mat, group_vec, method = "average", ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装 ggplot2")
    if (!requireNamespace("ggdendro", quietly = TRUE)) stop("需要安装 ggdendro")
    
    mat_t <- t(expr_mat)
    dist_mat <- dist(mat_t, method = "euclidean")
    hc <- hclust(dist_mat, method = method)
    
    library(ggdendro)
    dend <- ggdendro::dendro_data(hc)
    labels <- hc$labels[hc$order]
    leaf_groups <- group_vec[labels]
    
    p <- ggplot2::ggplot(ggdendro::segment(dend)) +
      ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::geom_text(data = ggdendro::label(dend),
                         ggplot2::aes(x = x, y = y, label = label, color = leaf_groups),
                         hjust = -0.1, size = 3) +
      ggplot2::coord_flip() + 
      ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
      ggplot2::labs(title = paste("Sample Dendrogram (", assay_name, ")"),
                    x = "Samples", y = "Height") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(), 
                     aspect.ratio = 1)
    
    if (!is.null(out_dir)) {
      file_prefix <- file.path(out_dir, paste0("Samples_dendro_", assay_name))
      width <- max(6, ncol(expr_mat) * 0.3)
      height <- ifelse("height" %in% names(list(...)), list(...)$height, 6)
      .save_plot(p, file_prefix, width = width, height = height)
    }
    return(p)
  }
  
  # 5. 标准化的表达矩阵
  calculate_normalized_matrix = function(design_obj, out_dir){
    if (!suppressMessages(require(openxlsx))) stop("需要安装 openxlsx 包")
    # 输出表格
    message("输出表格", stringr::str_glue("{out_dir}/reads_count_FPKM.xlsx"))
    openxlsx::write.xlsx(design_obj$get_assay(name = "fpkm"), stringr::str_glue("{out_dir}/reads_count_FPKM.xlsx"), rowNames =TRUE)
    message("输出表格", stringr::str_glue("{out_dir}/reads_count_TPM.xlsx"))
    openxlsx::write.xlsx(design_obj$get_assay(name = "tpm"), stringr::str_glue("{out_dir}/reads_count_TPM.xlsx"), rowNames =TRUE)
  }
  
  # z.主逻辑
  # 提取数据（如果 design_obj 是 Object_experiment_design，则自动处理）
  if (!inherits(design_obj, "Object_experiment_design")) {
    stop("design_obj 必须是 Object_experiment_design 对象")
  }
  
  expr_mat <- .get_expr_matrix(design_obj, assay_name)
  samples <- colnames(expr_mat)
  group_vec <- .get_group_vec(design_obj, samples)
  
  # 根据 type 执行对应的绘图函数
  if (type == "boxplot") {
    return(.plot_boxplot(expr_mat, group_vec, ...))
  } else if (type == "pca") {
    return(.plot_pca(expr_mat, group_vec, ...))
  } else if (type == "heatmap") {
    return(.plot_heatmap(expr_mat, group_vec, ...))
  } else if (type == "dendro") {
    return(.plot_dendro(expr_mat, group_vec, ...))
  } else if (type == "normalized_matrix") {
    return(calculate_normalized_matrix(design_obj = design_obj, out_dir = out_dir))
  } 
  else{
    stop("选项错误")
  }
}

# ==================== run_pca_style_visualization ====================
# 多种 PCA 图像风格测试/输出函数。
# 该函数属于样本级基础可视化工具，复用 design_obj 中的 assay 和分组信息。
run_pca_style_visualization <- function(design_obj,
                                        assay_name = "counts",
                                        styles = c("all", "points_ellipse", "sample_labels",
                                                   "group_centroid", "minimal_points"),
                                        log2_transform = TRUE,
                                        center = TRUE,
                                        scale. = TRUE,
                                        out_dir = NULL,
                                        width = 7,
                                        height = 6) {
  
  # 参数检查
  if (!inherits(design_obj, "read_counts")) {
    stop("design_obj 必须是 'read_counts' 或其子类对象")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要安装 ggplot2 包")
  }
  
  available_styles = c("points_ellipse", "sample_labels", "group_centroid", "minimal_points")
  styles = match.arg(styles, choices = c("all", available_styles), several.ok = TRUE)
  if ("all" %in% styles) styles = available_styles
  
  if (!is.null(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # ---------- 1. 准备 PCA 输入矩阵 ----------
  expr_mat = design_obj$get_assay(assay_name)
  if (log2_transform) {
    expr_mat = log2(expr_mat + 1)
  }
  
  # 转置后行为样本、列为基因/词条；去除零方差特征，避免 PCA 报错。
  pca_input = t(expr_mat)
  keep_feature = apply(pca_input, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  pca_input = pca_input[, keep_feature, drop = FALSE]
  if (ncol(pca_input) == 0) {
    stop("PCA 绘图失败：过滤零方差特征后没有剩余特征")
  }
  
  pca_res = stats::prcomp(pca_input, center = center, scale. = scale.)
  if (ncol(pca_res$x) < 2) {
    stop("PCA 绘图失败：主成分数量少于 2")
  }
  
  variance_pct = round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)
  
  group_vec = design_obj$GroupSamples[, "Group"]
  names(group_vec) = design_obj$GroupSamples[, "Samples"]
  
  pca_df = data.frame(
    Sample = rownames(pca_res$x),
    Group = group_vec[rownames(pca_res$x)],
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(out_dir)) {
    write.table(
      pca_df,
      file = file.path(out_dir, "pca_coordinates.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    write.table(
      data.frame(PC = paste0("PC", seq_along(variance_pct)), variance_pct = variance_pct),
      file = file.path(out_dir, "pca_variance_pct.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  # ---------- 2. PCA 绘图风格轮子 ----------
  pca_axis_labs = list(
    x = paste0("PC1 (", variance_pct[1], "%)"),
    y = paste0("PC2 (", variance_pct[2], "%)")
  )
  
  pca_style_plotters = list()
  
  pca_style_plotters[["points_ellipse"]] = function(pca_df, variance_pct) {
    p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::labs(x = pca_axis_labs$x, y = pca_axis_labs$y, title = "PCA: group points / eligible ellipses") +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
    eligible <- names(which(table(pca_df$Group) >= 4))
    if (length(eligible)) p <- p + ggplot2::stat_ellipse(data = pca_df[pca_df$Group %in% eligible, ], linewidth = 0.6, linetype = 2)
    if (any(table(pca_df$Group) < 4)) p <- p + ggplot2::labs(caption = "Ellipses omitted for groups with fewer than 4 samples")
    p
  }
  
  pca_style_plotters[["sample_labels"]] = function(pca_df, variance_pct) {
    p = ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_point(size = 2.8, alpha = 0.9) +
      ggplot2::labs(x = pca_axis_labs$x, y = pca_axis_labs$y,
                    title = "PCA - Sample Labels") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.grid.minor = ggplot2::element_blank(),
                     aspect.ratio = 1)
    
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p = p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 3,
                                       max.overlaps = Inf)
    } else {
      p = p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = 3,
                                 vjust = -0.8, show.legend = FALSE)
    }
    p
  }
  
  pca_style_plotters[["group_centroid"]] = function(pca_df, variance_pct) {
    centroid_df = stats::aggregate(cbind(PC1, PC2) ~ Group, data = pca_df, FUN = mean)
    plot_df = merge(pca_df, centroid_df, by = "Group", suffixes = c("", "_center"))
    
    ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_segment(ggplot2::aes(xend = PC1_center, yend = PC2_center),
                            alpha = 0.35, linewidth = 0.4) +
      ggplot2::geom_point(size = 2.8, alpha = 0.9) +
      ggplot2::geom_point(data = centroid_df,
                          ggplot2::aes(x = PC1, y = PC2, fill = Group),
                          shape = 23, color = "black", size = 4, inherit.aes = FALSE) +
      ggplot2::labs(x = pca_axis_labs$x, y = pca_axis_labs$y,
                    title = "PCA - Group Centroid") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     aspect.ratio = 1)
  }
  
  pca_style_plotters[["minimal_points"]] = function(pca_df, variance_pct) {
    ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
      ggplot2::geom_point(size = 3.2, alpha = 0.9) +
      ggplot2::labs(x = pca_axis_labs$x, y = pca_axis_labs$y,
                    title = "PCA - Minimal Points") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.grid.minor = ggplot2::element_blank(),
                     aspect.ratio = 1)
  }
  
  save_pca_style_plot = function(plot, style_name) {
    if (is.null(out_dir)) return(invisible(NULL))
    file_prefix = file.path(out_dir, paste0("pca_", style_name))
    ggplot2::ggsave(paste0(file_prefix, ".pdf"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
    ggplot2::ggsave(paste0(file_prefix, ".png"), plot, width = width, height = height,
                    limitsize = FALSE, units = "in")
    saveRDS(plot, file = paste0(file_prefix, ".rds"))
  }
  
  # ---------- 3. 批量测试并输出 ----------
  pca_style_results = list()
  for (style_name in styles) {
    message("[", Sys.time(), "] 开始测试 PCA 风格：", style_name)
    pca_style_results[[style_name]] = tryCatch({
      p = pca_style_plotters[[style_name]](pca_df = pca_df, variance_pct = variance_pct)
      save_pca_style_plot(p, style_name = style_name)
      list(status = "success", plot = p, error = NA_character_)
    }, error = function(e) {
      warning("PCA 风格测试失败：", style_name, "；原因：", e$message)
      list(status = "failed", plot = NULL, error = e$message)
    })
  }
  
  pca_style_summary = data.frame(
    style_name = names(pca_style_results),
    status = vapply(pca_style_results, function(x) x$status, character(1)),
    error = vapply(pca_style_results, function(x) x$error, character(1)),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(out_dir)) {
    write.table(
      pca_style_summary,
      file = file.path(out_dir, "pca_style_test_summary.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    saveRDS(pca_style_results, file = file.path(out_dir, "pca_style_results.rds"))
  }
  
  return(list(
    pca_res = pca_res,
    pca_df = pca_df,
    variance_pct = variance_pct,
    style_results = pca_style_results,
    style_summary = pca_style_summary,
    assay_name = assay_name,
    styles = styles
  ))
}
