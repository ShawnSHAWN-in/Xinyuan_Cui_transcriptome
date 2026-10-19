read_config_lst <- function(file) {
  value <- yaml::read_yaml(file, eval.expr = FALSE)
  if (!is.list(value)) stop("配置必须为映射：", file)
  value
}

project_path <- function(path) {
  if (!nzchar(path)) return(character())
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(Sys.getenv("TRANSCRIPTOME_PROJECT_DIR", unset = getwd()), path)
}

write_table_path <- function(tbl, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tbl <- as.data.frame(tbl)
  if (!ncol(tbl)) tbl <- data.frame(Status = character())
  for (field in names(tbl)) if (is.character(tbl[[field]])) tbl[[field]] <- gsub("[\r\n\t]", " ", tbl[[field]])
  data.table::fwrite(tbl, path, sep = "\t", quote = FALSE, na = "NA")
  path
}

write_object_path <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  path
}

safe_id <- function(x) {
  raw <- enc2utf8(as.character(x))
  value <- gsub("[^A-Za-z0-9_.-]", "_", raw)
  changed <- value != raw | nchar(value) > 180
  if (any(changed)) {
    hash <- vapply(raw[changed], function(item) substr(digest::digest(item, algo = "xxhash64", serialize = FALSE), 1, 8), character(1))
    value[changed] <- paste0(substr(value[changed], 1, 170), "_", hash)
  }
  if (any(!nzchar(value))) stop("输出标识不能为空")
  value
}

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("缺少包：", paste(missing, collapse = ", "), "；运行scripts/check_environment.R查看安装说明")
  invisible(TRUE)
}

read_count_data_lst <- function(file) {
  tbl <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  if (!all(c("Geneid", "Length") %in% names(tbl))) stop("count表必须包含Geneid和Length")
  if (anyDuplicated(tbl$Geneid) || anyNA(tbl$Geneid) || any(!nzchar(tbl$Geneid))) stop("Geneid重复或缺失")
  samples <- setdiff(names(tbl), c("Geneid", "Length"))
  if (!length(samples) || anyDuplicated(samples)) stop("样本列为空或重复")
  counts <- as.matrix(tbl[, samples, drop = FALSE])
  storage.mode(counts) <- "double"
  rownames(counts) <- tbl$Geneid
  lengths <- stats::setNames(as.numeric(tbl$Length), tbl$Geneid)
  if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) stop("count必须为非负整数")
  if (any(colSums(counts) == 0)) stop("存在总计数为0的样本")
  if (any(!is.finite(lengths)) || any(lengths <= 0)) stop("Length必须为正数")
  list(counts = counts, lengths = lengths)
}
