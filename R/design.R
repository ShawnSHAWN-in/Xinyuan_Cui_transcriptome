make_comparison_jobs_lst <- function(metadata, config) {
  jobs <- list()
  for (scheme in config$schemes) {
    if (!isTRUE(scheme$enabled)) next
    if (!scheme$type %in% c("pairwise", "explicit")) stop("首版支持pairwise和explicit方案；交互效应需独立模型")
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", scheme$id)) stop("方案ID仅允许字母、数字和下划线")
    if (!scheme$group %in% names(metadata)) stop("不存在分组变量：", scheme$group)
    strata <- as.character(unlist(scheme$strata))
    if (!all(strata %in% names(metadata))) stop("分层变量不存在")
    if (scheme$type == "pairwise" && (!scheme$variable %in% names(metadata) || scheme$variable %in% strata)) stop("比较变量不存在或与分层变量相同")
    strata_tbl <- if (length(strata)) unique(metadata[strata]) else data.frame(.all = 1)
    if (length(strata)) strata_tbl <- sort_metadata_tbl(strata_tbl, strata)
    for (s in seq_len(nrow(strata_tbl))) {
      use <- rep(TRUE, nrow(metadata))
      for (name in strata) use <- use & as.character(metadata[[name]]) == as.character(strata_tbl[[name]][s])
      block <- metadata[use, , drop = FALSE]
      pairs <- if (scheme$type == "explicit") scheme$pairs else lapply(scheme$treatments, function(x) {
        list(numerator = as.character(x), denominator = as.character(scheme$reference))
      })
      for (pair in pairs) {
        variable <- if (scheme$type == "explicit") scheme$group else scheme$variable
        exp <- block[as.character(block[[variable]]) == pair$numerator, , drop = FALSE]
        ctrl <- block[as.character(block[[variable]]) == pair$denominator, , drop = FALSE]
        if (!nrow(exp) || !nrow(ctrl)) stop("比较缺少实验或对照样本：", scheme$id, " ", pair$numerator, " vs ", pair$denominator)
        a <- unique(as.character(exp[[scheme$group]]))
        b <- unique(as.character(ctrl[[scheme$group]]))
        if (length(a) != 1L || length(b) != 1L || a == b) stop("每侧必须对应唯一且不同的组合分组，请增加strata")
        comparison <- paste0(a, "_vs_", b)
        id <- paste(scheme$id, comparison, sep = "__")
        job <- list(id = id, scheme = scheme$id, comparison = comparison, numerator = a, denominator = b,
                    strata = strata, samples = tibble::tibble(Sample = c(exp$Sample, ctrl$Sample), Group = c(rep(a, nrow(exp)), rep(b, nrow(ctrl)))))
        job$labels <- as.list(strata_tbl[s, strata, drop = FALSE])
        job$labels[[variable]] <- exp[[variable]][1]
        if (id %in% names(jobs)) stop("比较任务ID重复：", id)
        jobs[[id]] <- job
      }
    }
  }
  if (!length(jobs)) stop("没有启用的比较方案")
  ids <- vapply(jobs, `[[`, character(1), "id")
  if (anyDuplicated(ids) || anyDuplicated(safe_id(ids))) stop("比较任务ID重复或文件名冲突")
  jobs
}

comparison_manifest_tbl <- function(jobs) {
  dplyr::bind_rows(lapply(jobs, function(job) {
    data.frame(Job = job$id, Scheme = job$scheme, Comparison = job$comparison, Numerator = job$numerator,
               Denominator = job$denominator, N_experiment = sum(job$samples$Group == job$numerator),
               N_control = sum(job$samples$Group == job$denominator), as.data.frame(job$labels), check.names = FALSE)
  }))
}
