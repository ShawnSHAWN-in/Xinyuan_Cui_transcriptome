cast_variable <- function(x, spec, name) {
  raw <- as.character(x)
  if (anyNA(raw) || any(!nzchar(trimws(raw)))) stop("变量缺失：", name)
  if (spec$type == "factor") {
    levels <- as.character(unlist(spec$levels))
    if (!length(levels) || anyDuplicated(levels)) stop("factor必须提供唯一levels：", name)
    if (length(setdiff(raw, levels))) stop("变量有未知水平：", name, "：", paste(setdiff(raw, levels), collapse = ", "))
    return(factor(raw, levels = levels))
  }
  if (spec$type %in% c("numeric", "integer")) {
    value <- suppressWarnings(as.numeric(raw))
    if (any(!is.finite(value))) stop("变量无法转为数值：", name)
    if (spec$type == "integer") {
      if (any(value != round(value))) stop("变量不是整数：", name)
      value <- as.integer(value)
    }
    return(value)
  }
  if (spec$type == "character") return(raw)
  stop("不支持变量类型：", spec$type)
}

read_sample_metadata_tbl <- function(file, config) {
  variable_names <- c(names(config$variables), names(config$derived), names(config$combinations))
  if (anyDuplicated(variable_names) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", variable_names))) stop("变量名须唯一且使用字母、数字、下划线")
  tbl <- data.table::fread(file, colClasses = "character", data.table = FALSE)
  sample_col <- config$sample_column
  if (!sample_col %in% names(tbl) && "Samples" %in% names(tbl)) sample_col <- "Samples"
  if (!sample_col %in% names(tbl)) stop("样本表缺少样本ID列：", config$sample_column)
  names(tbl)[names(tbl) == sample_col] <- "Sample"
  if (anyNA(tbl$Sample) || any(!nzchar(tbl$Sample)) || anyDuplicated(tbl$Sample)) stop("样本ID缺失或重复")
  if (isTRUE(config$legacy_parse) && config$legacy_group_column %in% names(tbl)) {
    group <- tbl[[config$legacy_group_column]]
    if (!"Condition" %in% names(tbl)) tbl$Condition <- sub("_[^_]+$", "", group)
    if (!"Time_f" %in% names(tbl)) tbl$Time_f <- sub("^.*_", "", group)
    if (!"Replicate" %in% names(tbl)) tbl$Replicate <- sub("^.*_", "", tbl$Sample)
  }
  for (name in names(config$variables)) {
    if (!name %in% names(tbl)) stop("样本表缺少变量：", name)
    tbl[[name]] <- cast_variable(tbl[[name]], config$variables[[name]], name)
  }
  for (name in names(config$derived)) {
    spec <- config$derived[[name]]
    if (!is.null(spec$constant)) value <- rep(spec$constant, nrow(tbl)) else {
      if (!spec$from %in% names(tbl)) stop("派生变量来源不存在：", spec$from)
      value <- tbl[[spec$from]]
    }
    tbl[[name]] <- cast_variable(value, spec, name)
  }
  for (name in names(config$combinations)) {
    spec <- config$combinations[[name]]
    variables <- as.character(unlist(spec$variables))
    order_vars <- as.character(unlist(spec$order))
    if (!setequal(variables, order_vars) || !all(variables %in% names(tbl))) stop("组合变量及排序声明不一致：", name)
    if (any(vapply(tbl[variables], is.numeric, logical(1)))) stop("组合分组请使用factor或character变量：", name)
    labels <- do.call(paste, c(lapply(tbl[variables], as.character), sep = spec$separator))
    tuples <- unique(as.data.frame(lapply(tbl[variables], as.character)))
    tuple_labels <- do.call(paste, c(tuples, sep = spec$separator))
    if (anyDuplicated(tuple_labels)) stop("组合标签冲突，请更换separator：", name)
    sort_keys <- lapply(tbl[order_vars], function(x) if (is.factor(x)) as.integer(x) else x)
    label_levels <- unique(labels[do.call(order, sort_keys)])
    tbl[[name]] <- factor(labels, levels = label_levels)
  }
  tibble::as_tibble(tbl)
}

check_sample_alignment <- function(metadata, count_data) {
  if (!setequal(metadata$Sample, colnames(count_data$counts))) stop("样本表与count矩阵样本集合不一致")
  invisible(TRUE)
}

sort_metadata_tbl <- function(metadata, variables) {
  variables <- as.character(unlist(variables))
  if (!all(variables %in% names(metadata))) stop("绘图排序含未知变量")
  metadata[do.call(order, lapply(metadata[variables], function(x) if (is.factor(x)) as.integer(x) else x)), , drop = FALSE]
}

order_comparison_jobs_lst <- function(jobs, variables) {
  tbl <- comparison_manifest_tbl(jobs)
  variables <- intersect(as.character(unlist(variables)), names(tbl))
  tbl$Scheme <- factor(tbl$Scheme, levels = unique(tbl$Scheme))
  tbl <- sort_metadata_tbl(tbl, c("Scheme", variables))
  jobs[match(tbl$Job, vapply(jobs, `[[`, character(1), "id"))]
}

order_comparison_results_lst <- function(results, variables) {
  jobs <- lapply(results, `[[`, "job")
  ordered <- order_comparison_jobs_lst(jobs, variables)
  results[match(vapply(ordered, `[[`, character(1), "id"), vapply(jobs, `[[`, character(1), "id"))]
}
