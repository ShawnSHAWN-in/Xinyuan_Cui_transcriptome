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

build_kegg_snapshot_paths <- function(config, source_files, output_dir) {
  dir <- file.path(output_dir, "02_annotation", "KEGG")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  status <- "DISABLED"
  detail <- "run_kegg=FALSE"
  bg <- list(TERM2GENE = data.frame(TERM = character(), GENE = character()), TERM2NAME = data.frame(TERM = character(), NAME = character()))
  if (isTRUE(config$run_kegg)) {
    value <- tryCatch({
      if (nzchar(config$kegg_term2gene)) {
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
