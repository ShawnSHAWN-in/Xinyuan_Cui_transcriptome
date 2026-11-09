# 读取注释主表
read_eggnog_annotation_tbl <- function(file) {
  if (!file.exists(file)) stop("eggNOG注释文件不存在：", file)
  input <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(input), add = TRUE)
  clean_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(clean_file), add = TRUE)
  output <- base::file(clean_file, "wt")
  header <- NULL
  comments <- character()
  tryCatch({
    repeat {
      lines <- readLines(input, n = 10000L, warn = FALSE)
      if (!length(lines)) break
      if (is.null(header)) {
        index <- grep("^#?(query|query_name)\\t", lines)
        if (!length(index)) {
          comments <- utils::head(c(comments, lines[startsWith(lines, "##")]), 50L)
          if (any(nzchar(lines) & !startsWith(lines, "#"))) stop("未找到eggNOG列名；请提供.emapper.annotations文件")
          next
        }
        index <- index[1]
        comments <- utils::head(c(comments, lines[seq_len(index - 1L)]), 50L)
        header <- strsplit(sub("^#", "", lines[index]), "\t", fixed = TRUE)[[1]]
        header[header == "query_name"] <- "query"
        if (anyDuplicated(header) || !all(c("query", "GOs", "KEGG_Pathway") %in% header)) {
          stop("eggNOG列名必须唯一且包含query、GOs、KEGG_Pathway；不能使用hits或seed_orthologs文件")
        }
        writeLines(paste(header, collapse = "\t"), output)
        lines <- if (index < length(lines)) lines[seq.int(index + 1L, length(lines))] else character()
      }
      lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
      if (length(lines)) {
        fields <- nchar(lines, type = "bytes") - nchar(gsub("\t", "", lines, fixed = TRUE), type = "bytes") + 1L
        if (any(fields != length(header))) stop("eggNOG数据行列数不一致，请检查TSV格式")
        writeLines(lines, output)
      }
    }
  }, finally = close(output))
  if (is.null(header)) stop("未找到eggNOG表头")
  fields <- intersect(c("query", "GOs", "KEGG_Pathway", "KEGG_ko", "Preferred_name", "Description"), header)
  tbl <- data.table::fread(clean_file, select = fields, colClasses = "character", na.strings = c("", "-", "NA"), data.table = FALSE)
  if (!nrow(tbl) || anyNA(tbl$query) || any(!nzchar(trimws(tbl$query))) || anyDuplicated(tbl$query)) stop("eggNOG query为空或重复")
  attr(tbl, "eggnog_header") <- comments
  message("读取eggNOG：", nrow(tbl), "个query")
  tibble::as_tibble(tbl)
}

# 显式连接蛋白与基因
read_eggnog_gene_map_tbl <- function(query_gene_file, queries, query_is_gene_id = FALSE) {
  if (is.null(query_gene_file) || !length(query_gene_file) || !nzchar(query_gene_file)) {
    if (!isTRUE(query_is_gene_id)) stop("请提供query、Geneid两列映射表；仅query确为分析基因ID时设置query_is_gene_id=TRUE")
    return(tibble::tibble(query = queries, Geneid = queries))
  }
  if (isTRUE(query_is_gene_id)) stop("映射表与query_is_gene_id=TRUE不能同时使用")
  tbl <- data.table::fread(query_gene_file, colClasses = "character", na.strings = c("", "-", "NA"), data.table = FALSE)
  if (!all(c("query", "Geneid") %in% names(tbl))) stop("蛋白-基因映射表需要query、Geneid两列")
  tbl <- unique(tbl[, c("query", "Geneid"), drop = FALSE])
  if (!nrow(tbl) || anyNA(tbl) || any(!nzchar(trimws(unlist(tbl)))) || anyDuplicated(tbl$query)) {
    stop("映射表不能为空；每个query只能对应一个Geneid，多个query可以对应同一Geneid")
  }
  if (!any(queries %in% tbl$query)) stop("映射表query与eggNOG结果完全不匹配；不要自动删前缀或版本号")
  tibble::as_tibble(tbl)
}

# 展开逗号分隔关联
expand_eggnog_terms_tbl <- function(tbl, field, mapping) {
  values <- strsplit(ifelse(is.na(tbl[[field]]), "", tbl[[field]]), ",", fixed = TRUE)
  assoc <- data.frame(query = rep(tbl$query, lengths(values)), TERM = trimws(unlist(values)), stringsAsFactors = FALSE)
  assoc <- assoc[nzchar(assoc$TERM) & assoc$TERM != "-", , drop = FALSE]
  assoc$GENE <- mapping$Geneid[match(assoc$query, mapping$query)]
  unique(assoc[!is.na(assoc$GENE), c("TERM", "GENE"), drop = FALSE])
}

# 复用GO祖先展开
build_go_background_lst <- function(direct_go, ontology = "BP", include_ancestors = TRUE) {
  if (!ontology %in% c("BP", "CC", "MF")) stop("GO ontology必须是BP/CC/MF")
  metadata <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), columns = c("TERM", "ONTOLOGY")))
  names <- metadata[!is.na(metadata$ONTOLOGY) & metadata$ONTOLOGY == ontology, c("GOID", "TERM")]
  colnames(names) <- c("TERM", "NAME")
  names <- unique(names)
  assoc <- unique(data.frame(TERM = direct_go$GO, GENE = as.character(direct_go$ENTREZID)))
  assoc <- assoc[assoc$TERM %in% names$TERM & !is.na(assoc$GENE), , drop = FALSE]
  if (isTRUE(include_ancestors) && nrow(assoc)) {
    ancestors <- as.list(get(paste0("GO", ontology, "ANCESTOR"), envir = asNamespace("GO.db")))
    expansion <- dplyr::bind_rows(lapply(unique(assoc$TERM), function(id) {
      values <- unique(c(id, ancestors[[id]]))
      data.frame(Source = id, TERM = values[!is.na(values) & values %in% names$TERM])
    }))
    assoc <- dplyr::inner_join(assoc, expansion, by = c("TERM" = "Source"), relationship = "many-to-many") |>
      dplyr::select(TERM = TERM.y, GENE) |> dplyr::distinct()
  }
  list(TERM2GENE = assoc, TERM2NAME = names)
}

# 统一参考通路ID
normalize_eggnog_pathway_vec <- function(ids) {
  ids <- sub("^path:", "", ids)
  if (any(!grepl("^(map|ko)[0-9]{5}$", ids))) stop("KEGG_Pathway或名称表包含非map/ko通路ID")
  sub("^ko", "map", ids)
}

# 获取通路名称
read_eggnog_kegg_names_tbl <- function(file = NULL) {
  if (is.null(file) || !length(file) || !nzchar(file)) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("需要KEGGREST获取通路名称，或提供本地TERM、NAME两列名称表")
    terms <- tryCatch(KEGGREST::keggList("pathway"), error = function(e) stop("KEGG名称下载失败；请提供本地名称表：", conditionMessage(e)))
    tbl <- data.frame(TERM = names(terms), NAME = unname(terms))
  } else {
    tbl <- data.table::fread(file, colClasses = "character", data.table = FALSE)
    if (!all(c("TERM", "NAME") %in% names(tbl))) stop("KEGG名称表需要TERM、NAME两列")
    tbl <- tbl[, c("TERM", "NAME"), drop = FALSE]
  }
  if (!nrow(tbl) || anyNA(tbl) || any(!nzchar(trimws(unlist(tbl))))) stop("KEGG名称表为空或存在缺失值")
  tbl$TERM <- normalize_eggnog_pathway_vec(tbl$TERM)
  tbl <- unique(tbl)
  if (anyDuplicated(tbl$TERM)) stop("同一KEGG通路对应多个名称，请核对名称表版本")
  tibble::as_tibble(tbl)
}

# 构建标准富集背景
build_eggnog_annotation_lst <- function(eggnog_file, query_gene_file = NULL, query_is_gene_id = FALSE, species,
                                        ontology = "BP", include_ancestors = TRUE, kegg_term2name_file = NULL, run_kegg = TRUE) {
  if (missing(species) || length(species) != 1L || is.na(species) || !nzchar(trimws(species))) stop("必须提供实际参考物种species")
  tbl <- read_eggnog_annotation_tbl(eggnog_file)
  header <- attr(tbl, "eggnog_header")
  mapping <- read_eggnog_gene_map_tbl(query_gene_file, tbl$query, query_is_gene_id)
  go <- expand_eggnog_terms_tbl(tbl, "GOs", mapping)
  if (any(!grepl("^GO:[0-9]{7}$", go$TERM))) stop("GOs存在非法GO编号")
  direct_go <- tibble::tibble(ENTREZID = go$GENE, GO = go$TERM, SYMBOL = go$GENE)
  go_background <- build_go_background_lst(direct_go, ontology, include_ancestors)
  if (!nrow(go_background$TERM2GENE)) stop("映射后的GO-", ontology, "背景为空，请核对注释及ID映射")
  kegg <- expand_eggnog_terms_tbl(tbl, "KEGG_Pathway", mapping)
  kegg$TERM <- normalize_eggnog_pathway_vec(kegg$TERM)
  kegg <- unique(kegg)
  rownames(kegg) <- NULL
  kegg_names <- tibble::tibble(TERM = character(), NAME = character())
  if (isTRUE(run_kegg) && nrow(kegg)) {
    kegg_names <- read_eggnog_kegg_names_tbl(kegg_term2name_file)
    missing_names <- setdiff(kegg$TERM, kegg_names$TERM)
    if (length(missing_names)) stop("KEGG名称表缺少通路：", paste(utils::head(missing_names, 10), collapse = ", "))
    kegg_names <- kegg_names[kegg_names$TERM %in% kegg$TERM, , drop = FALSE]
  }
  kegg_background <- list(TERM2GENE = if (isTRUE(run_kegg)) kegg else kegg[FALSE, ], TERM2NAME = kegg_names)
  genes <- unique(mapping$Geneid)
  annotation <- tibble::tibble(Geneid = genes, SYMBOL = genes, ENTREZID = genes, GENENAME = genes)
  audit <- dplyr::full_join(mapping, tibble::tibble(query = tbl$query, HasAnnotationRow = TRUE), by = "query")
  audit$HasAnnotationRow[is.na(audit$HasAnnotationRow)] <- FALSE
  audit$HasGO <- audit$query %in% tbl$query[!is.na(tbl$GOs)]
  audit$HasKEGG <- audit$query %in% tbl$query[!is.na(tbl$KEGG_Pathway)]
  coverage <- tibble::tibble(Database = c(paste0("GO_", ontology), "KEGG"),
    Annotated = c(length(intersect(genes, go_background$TERM2GENE$GENE)), length(intersect(genes, kegg_background$TERM2GENE$GENE))),
    Total = length(genes), Proportion = Annotated / Total)
  known_go <- AnnotationDbi::keys(GO.db::GO.db)
  unknown_go <- tibble::tibble(GO = setdiff(unique(direct_go$GO), known_go))
  if (nrow(unknown_go)) warning("部分GO编号不在当前GO.db中，已记录unknown_go；请核对数据库版本", call. = FALSE)
  if (anyNA(audit$Geneid)) warning("部分eggNOG query缺少基因映射，已记录query_mapping且不进入背景", call. = FALSE)
  paths <- c(eggnog_file, query_gene_file, kegg_term2name_file)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  metadata <- list(source = "eggNOG-mapper", species = species, id_type = "Geneid", compatibility_id_column = "ENTREZID",
    ontology = ontology, include_ancestors = include_ancestors, run_kegg = run_kegg, header = header,
    files = stats::setNames(unname(tools::md5sum(paths)), normalizePath(paths)), go_db_version = as.character(utils::packageVersion("GO.db")),
    kegg_names = if (!isTRUE(run_kegg) || !nrow(kegg)) "not_used" else if (is.null(kegg_term2name_file)) "KEGGREST::keggList(pathway)" else kegg_term2name_file,
    created = as.character(Sys.time()))
  message("eggNOG背景完成：", length(genes), "个基因；GO关联", nrow(go_background$TERM2GENE), "条；KEGG关联", nrow(kegg_background$TERM2GENE), "条")
  list(annotation = annotation, direct_go = direct_go, background = go_background, backgrounds = list(GO = go_background, KEGG = kegg_background),
    database = "GO", status = "SUCCESS", metadata = metadata, query_mapping = audit, coverage = coverage, unknown_go = unknown_go)
}

# 输出可复用参考文件
save_eggnog_reference_paths <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(gene_annotation = result$annotation, gene2go_direct = result$direct_go, query_mapping = result$query_mapping,
    coverage = result$coverage, unknown_go = result$unknown_go,
    gene_map = data.frame(gene = result$annotation$Geneid, enrichment_id = result$annotation$ENTREZID))
  for (db in names(result$backgrounds)) {
    tables[[paste0(db, "_TERM2GENE")]] <- result$backgrounds[[db]]$TERM2GENE
    tables[[paste0(db, "_TERM2NAME")]] <- result$backgrounds[[db]]$TERM2NAME
  }
  paths <- vapply(names(tables), function(name) {
    path <- file.path(output_dir, paste0(name, ".tsv"))
    data.table::fwrite(tables[[name]], path, sep = "\t", quote = FALSE, na = "NA")
    path
  }, character(1))
  gene_sets <- lapply(names(result$backgrounds), function(db) {
    bg <- result$backgrounds[[db]]
    if (!nrow(bg$TERM2GENE)) return(NULL)
    list(term2gene = stats::setNames(as.data.frame(bg$TERM2GENE), c("term_id", "gene_id")),
      term2name = stats::setNames(as.data.frame(bg$TERM2NAME), c("term_id", "term_name")),
      species = result$metadata$species, id_type = result$metadata$id_type,
      source = paste(result$metadata$source, paste(result$metadata$header[grepl("^##.*(emapper-|version)", result$metadata$header)], collapse = "; "),
        paste0("GO.db=", result$metadata$go_db_version), result$metadata$created, sep = "; "))
  })
  names(gene_sets) <- c(paste0("GO_", result$metadata$ontology), "KEGG")
  gene_sets <- Filter(Negate(is.null), gene_sets)
  saveRDS(gene_sets, file.path(output_dir, "gene_sets.rds"))
  saveRDS(result$metadata, file.path(output_dir, "annotation_metadata.rds"))
  c(unname(paths), file.path(output_dir, c("gene_sets.rds", "annotation_metadata.rds")))
}
