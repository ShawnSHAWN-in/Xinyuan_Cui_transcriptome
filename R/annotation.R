annotation_source_chr <- function(reference) {
  source <- reference$annotation_source
  if (is.null(source)) source <- "ncbi"
  if (!source %in% c("ncbi", "eggnog")) stop("annotation_source仅支持ncbi或eggnog")
  source
}

annotation_input_paths <- function(paths, reference) {
  fields <- if (annotation_source_chr(reference) == "ncbi") c("gtf_file", "gene2go_file") else "eggnog_file"
  if (annotation_source_chr(reference) == "eggnog") {
    if (!isTRUE(reference$query_is_gene_id)) fields <- c(fields, "query_gene_file")
    if (!is.null(paths$kegg_names_file) && nzchar(paths$kegg_names_file)) fields <- c(fields, "kegg_names_file")
  }
  if (any(vapply(fields, function(name) is.null(paths[[name]]) || !nzchar(paths[[name]]), logical(1)))) {
    stop("注释输入路径缺失：", paste(fields, collapse = ", "))
  }
  stats::setNames(vapply(fields, function(name) project_path(paths[[name]]), character(1)), fields)
}

read_project_annotation_lst <- function(source_files, reference, engine_files) {
  if (!identical(reference$go_ontology, "BP")) stop("本项目GO图像仅使用BP")
  if (annotation_source_chr(reference) == "ncbi") {
    return(read_go_annotation_lst(source_files[["gtf_file"]], source_files[["gene2go_file"]], reference, engine_files))
  }
  context <- new.env(parent = globalenv())
  helper <- engine_files[basename(engine_files) == "eggnog_annotation.R"]
  if (length(helper) != 1L) stop("engine_files需要包含eggnog_annotation.R")
  sys.source(helper, context)
  optional <- function(name) if (name %in% names(source_files)) source_files[[name]] else NULL
  context$build_eggnog_annotation_lst(source_files[["eggnog_file"]], optional("query_gene_file"), isTRUE(reference$query_is_gene_id),
    reference$reference_species, reference$go_ontology, reference$go_ancestors, optional("kegg_names_file"), isTRUE(reference$run_kegg))
}

read_go_annotation_lst <- function(gtf_file, gene2go_file, reference, engine_files) {
  if (!identical(reference$go_ontology, "BP")) stop("本项目GO图像仅使用BP")
  context <- new_legacy_context(engine_files, reference)
  cache <- context$cfg$output_dir
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  manager <- context$PublicDataManager$new(reference$kegg_code, cache_dir = file.path(cache, "cache"))
  manager$load_project_annotation(gtf_file, gene2go_file)
  annotation <- manager$get_gene_annotation()
  annotation$ENTREZID <- as.character(annotation$ENTREZID)
  background <- manager$get_go_background("BP", include_ancestors = reference$go_ancestors)
  if (!nrow(background$TERM2GENE)) stop("GO-BP背景为空")
  list(annotation = tibble::as_tibble(annotation), direct_go = tibble::as_tibble(manager$get_go_annotation()),
       background = background, database = "GO", status = "SUCCESS")
}

kegg_input_paths <- function(config, config_file) {
  paths <- c(config$kegg_term2gene, config$kegg_term2name)
  paths <- paths[nzchar(paths)]
  if (length(paths) == 1L) stop("KEGG本地背景必须同时提供TERM2GENE及TERM2NAME")
  if (length(paths)) vapply(paths, project_path, character(1)) else config_file
}

build_kegg_snapshot_paths <- function(config, source_files, output_dir, annotation = NULL) {
  dir <- file.path(output_dir, "02_annotation", "KEGG")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  status <- "DISABLED"
  detail <- "run_kegg=FALSE"
  bg <- list(TERM2GENE = data.frame(TERM = character(), GENE = character()), TERM2NAME = data.frame(TERM = character(), NAME = character()))
  if (isTRUE(config$run_kegg)) {
    value <- tryCatch({
      if (annotation_source_chr(config) == "eggnog") {
        if (is.null(annotation$backgrounds$KEGG)) stop("eggNOG KEGG背景尚未构建")
        annotation$backgrounds$KEGG
      } else if (nzchar(config$kegg_term2gene)) {
        assoc <- data.table::fread(source_files[1], colClasses = "character", data.table = FALSE)
        names <- data.table::fread(source_files[2], colClasses = "character", data.table = FALSE)
        if (!all(c("TERM", "GENE") %in% colnames(assoc)) || !all(c("TERM", "NAME") %in% colnames(names))) stop("KEGG本地背景列名错误")
        list(TERM2GENE = unique(assoc[, c("TERM", "GENE")]), TERM2NAME = unique(names[, c("TERM", "NAME")]))
      } else {
        check_packages("KEGGREST")
        links <- KEGGREST::keggLink("pathway", config$kegg_code)
        terms <- KEGGREST::keggList("pathway", config$kegg_code)
        list(TERM2GENE = data.frame(TERM = sub("^path:", "", unname(links)), GENE = sub("^[^:]+:", "", names(links))),
             TERM2NAME = data.frame(TERM = sub("^path:", "", names(terms)), NAME = unname(terms)))
      }
    }, error = function(e) e)
    if (inherits(value, "error") || !nrow(value$TERM2GENE) || !nrow(value$TERM2NAME)) {
      detail <- if (inherits(value, "error")) conditionMessage(value) else "KEGG背景为空"
      if (isTRUE(config$require_kegg)) stop(detail)
      warning("KEGG不可用：", detail, call. = FALSE)
      status <- "UNAVAILABLE"
    } else {
      bg <- value
      status <- "SUCCESS"
      detail <- paste("snapshot", config$kegg_refresh)
    }
  }
  c(write_table_path(bg$TERM2GENE, file.path(dir, "TERM2GENE.tsv")), write_table_path(bg$TERM2NAME, file.path(dir, "TERM2NAME.tsv")),
    write_table_path(data.frame(Database = "KEGG", Status = status, Detail = detail), file.path(dir, "status.tsv")))
}

make_background_jobs_lst <- function(go, kegg_files) {
  status <- data.table::fread(kegg_files[3], data.table = FALSE)
  list(GO = list(database = "GO", status = "SUCCESS", background = go$background),
       KEGG = list(database = "KEGG", status = status$Status[1], background = list(
         TERM2GENE = data.table::fread(kegg_files[1], colClasses = "character", data.table = FALSE),
         TERM2NAME = data.table::fread(kegg_files[2], colClasses = "character", data.table = FALSE))))
}

