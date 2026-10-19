# 增加标准化矩阵入口
read_counts$set("public", "set_assay", function(name, matrix) {
  private$add_assay(name, matrix)
  invisible(self)
})

# 压缩GTF与缓存
PublicDataManager$set("private", "parse_gtf_ncbi", function(gtf_file, use_cache = TRUE) {
  signature <- input_signature(gtf_file)
  cache_file <- file.path(cfg$output_dir, "02_annotation", "gtf_annotation_cache.rds")
  if (use_cache && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$signature, signature)) {
      private$.gtf_tbl <- cached$annotation
      return(invisible(private$.gtf_tbl))
    }
  }
  message("解析NCBI GTF：", gtf_file)
  input <- gtf_file
  if (grepl("\\.gz$", gtf_file)) {
    input <- tempfile(fileext = ".gtf")
    on.exit(unlink(input), add = TRUE)
    con_in <- gzfile(gtf_file, "rb")
    con_out <- file(input, "wb")
    tryCatch({
      repeat {
        bytes <- readBin(con_in, "raw", n = 1024L * 1024L)
        if (!length(bytes)) break
        writeBin(bytes, con_out)
      }
    }, finally = { close(con_in); close(con_out) })
  }
  gtf <- rtracklayer::import(input, format = "gtf")
  if (any(gtf$type == "gene")) gtf <- gtf[gtf$type == "gene"]
  tbl <- as.data.frame(gtf)
  if (!"gene_id" %in% names(tbl)) stop("GTF缺少gene_id")
  flatten <- function(x) vapply(as.list(x), function(value) paste(as.character(value), collapse = ";"), character(1))
  if (!"db_xref" %in% names(tbl)) stop("GTF缺少db_xref，无法获得NCBI GeneID")
  xrefs <- flatten(tbl$db_xref)
  gene_ids <- regmatches(xrefs, gregexpr("GeneID:[0-9]+", xrefs))
  symbols <- if ("gene" %in% names(tbl)) flatten(tbl$gene) else as.character(tbl$gene_id)
  descriptions <- if ("description" %in% names(tbl)) flatten(tbl$description) else if ("product" %in% names(tbl)) flatten(tbl$product) else symbols
  annotation <- dplyr::bind_rows(lapply(seq_len(nrow(tbl)), function(i) {
    ids <- unique(sub("^GeneID:", "", gene_ids[[i]]))
    if (!length(ids)) ids <- NA_character_
    data.frame(Geneid = as.character(tbl$gene_id[i]), SYMBOL = symbols[i], ENTREZID = ids, GENENAME = descriptions[i])
  })) |>
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", Geneid, SYMBOL)) |>
    dplyr::distinct()
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(signature = signature, annotation = annotation), cache_file)
  private$.gtf_tbl <- annotation
  invisible(annotation)
}, overwrite = TRUE)

# 按物种流式读取GO
PublicDataManager$set("private", "load_go_gene2go.gz", function(gene2go.gz = NULL, use_cache = TRUE) {
  if (is.null(private$.gtf_tbl)) stop("必须先加载GTF")
  signature <- input_signature(gene2go.gz)
  gene_ids <- sort(unique(stats::na.omit(as.character(private$.gtf_tbl$ENTREZID))))
  cache_file <- file.path(cfg$output_dir, "02_annotation", "gene2go_annotation_cache.rds")
  if (use_cache && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$signature, signature) && identical(cached$gene_ids, gene_ids) && identical(cached$tax_id, cfg$tax_id)) {
      private$.go_tbl <- cached$annotation
      return(private$.go_tbl)
    }
  }
  con <- if (grepl("\\.gz$", gene2go.gz)) gzfile(gene2go.gz, "rt") else file(gene2go.gz, "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  header <- sub("^#", "", header)
  header <- gsub(" ", "_", header, fixed = TRUE)
  required <- c("tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier", "GO_term", "Category")
  if (!all(required %in% header)) stop("gene2go列不符合标准格式：", paste(setdiff(required, header), collapse = ", "))
  chunks <- list()
  index <- 0L
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[startsWith(lines, paste0(cfg$tax_id, "\t"))]
    if (!length(lines)) next
    tbl <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE, col.names = header, colClasses = "character", data.table = FALSE)
    tbl <- tbl[tbl$tax_id == as.character(cfg$tax_id) & tbl$GeneID %in% gene_ids, , drop = FALSE]
    tbl <- tbl[!grepl("(^|[|,])NOT([|,]|$)", tbl$Qualifier) & grepl("^GO:[0-9]{7}$", tbl$GO_ID), , drop = FALSE]
    if (nrow(tbl)) {
      index <- index + 1L
      chunks[[index]] <- tbl
    }
  }
  if (!length(chunks)) stop("gene2go中没有可匹配的参考物种注释；请核对文件完整性、tax_id和GTF GeneID")
  annotation <- dplyr::bind_rows(chunks) |>
    dplyr::transmute(ENTREZID = GeneID, GO = GO_ID, Evidence, Qualifier, GO_term, Category, tax_id) |>
    dplyr::distinct()
  symbols <- private$.gtf_tbl |> dplyr::distinct(ENTREZID, SYMBOL)
  annotation <- dplyr::left_join(annotation, symbols, by = "ENTREZID")
  saveRDS(list(signature = signature, gene_ids = gene_ids, tax_id = cfg$tax_id, annotation = annotation), cache_file)
  private$.go_tbl <- annotation
  annotation
}, overwrite = TRUE)

# GO祖先关联统一处理
PublicDataManager$set("public", "get_go_background", function(ontology = "BP", include_ancestors = cfg$go_ancestors) {
  if (!ontology %in% c("BP", "CC", "MF")) stop("GO ontology必须是BP/CC/MF")
  if (is.null(private$.go_tbl)) stop("GO注释尚未加载")
  metadata <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), columns = c("TERM", "ONTOLOGY")))
  names <- metadata[metadata$ONTOLOGY == ontology, c("GOID", "TERM")]
  colnames(names) <- c("TERM", "NAME")
  names <- unique(names)
  assoc <- unique(data.frame(TERM = private$.go_tbl$GO, GENE = as.character(private$.go_tbl$ENTREZID)))
  assoc <- assoc[assoc$TERM %in% names$TERM & !is.na(assoc$GENE), , drop = FALSE]
  if (include_ancestors && nrow(assoc)) {
    ancestors <- as.list(get(paste0("GO", ontology, "ANCESTOR"), envir = asNamespace("GO.db")))
    expansion <- dplyr::bind_rows(lapply(unique(assoc$TERM), function(id) {
      values <- unique(c(id, ancestors[[id]]))
      data.frame(Source = id, TERM = values[!is.na(values) & values %in% names$TERM])
    }))
    assoc <- dplyr::inner_join(assoc, expansion, by = c("TERM" = "Source"), relationship = "many-to-many") |>
      dplyr::select(TERM = TERM.y, GENE) |> dplyr::distinct()
  }
  list(TERM2GENE = assoc, TERM2NAME = names)
}, overwrite = TRUE)

# KEGG下载仅依赖KEGGREST
PublicDataManager$set("private", "load_kegg", function(use_cache = TRUE) {
  if (use_cache && !private$.force_update) {
    cached <- private$read_cache("kegg")
    if (!is.null(cached)) {
      private$.kegg_tbl <- cached
      return(invisible(cached))
    }
  }
  if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("KEGG模块需要KEGGREST")
  links <- KEGGREST::keggLink("pathway", private$.kegg_code)
  terms <- KEGGREST::keggList("pathway", private$.kegg_code)
  if (!length(links) || !length(terms)) stop("KEGG返回空背景")
  assoc <- data.frame(Keggpathid = sub("^path:", "", unname(links)), ENTREZID = sub("^[^:]+:", "", names(links)))
  names <- data.frame(Keggpathid = sub("^path:", "", names(terms)), Keggpathname = unname(terms))
  tbl <- dplyr::left_join(assoc, names, by = "Keggpathid") |> dplyr::distinct()
  private$.kegg_tbl <- tbl
  private$write_cache(tbl, "kegg")
  invisible(tbl)
}, overwrite = TRUE)

# 本地注释不依赖网络预检
PublicDataManager$set("public", "initialize", function(kegg_organism_code, cache_dir = NULL, force_update = FALSE) {
  private$.kegg_code <- kegg_organism_code
  private$.force_update <- force_update
  private$.cache_dir <- if (is.null(cache_dir)) file.path(cfg$output_dir, "02_annotation", "cache") else cache_dir
  private$ensure_cache_dir()
}, overwrite = TRUE)
PublicDataManager$set("public", "load_project_annotation", function(gtf_file, gene2go_file) {
  private$parse_gtf_ncbi(gtf_file)
  private$load_go_gene2go.gz(gene2go_file)
  invisible(self)
})
PublicDataManager$set("public", "get_go_annotation", function() private$.go_tbl)

# ID转换只用于图上标签
convert_result_to_symbol <- function(res, public_data) res

# 保留既有富集接口
run_enrichment_analysis <- function(diff_list, public_data, analysis_type = c("ORA", "GSEA"), direction = c("Up", "Down"),
                                    db = c("KEGG", "GO"), out_dir = NULL, background = NULL, universe = NULL,
                                    minGSSize = cfg$min_gs_size, maxGSSize = cfg$max_gs_size) {
  require_enrichment_packages()
  analysis_type <- match.arg(analysis_type, several.ok = TRUE)
  direction <- match.arg(direction)
  db <- match.arg(db)
  df <- attach_annotation_tbl(as.data.frame(diff_list), public_data$get_gene_annotation())
  if (length(unique(df$group)) != 1L) stop("单比较富集要求唯一group")
  if (is.null(background)) background <- if (db == "GO") public_data$get_go_background("BP") else public_data$get_kegg_background()
  bg <- background
  df <- df[!is.na(df$ENTREZID) & is.finite(df$logFC), , drop = FALSE]
  if (is.null(universe) && "ORA" %in% analysis_type) stop("ORA必须显式提供该比较受检GeneID背景")
  results <- list(contrast_name = unique(as.character(diff_list$group)), ORA = NULL, GSEA = NULL, ranking = NULL)
  if ("ORA" %in% analysis_type) {
    genes <- unique(df$ENTREZID[df$Result_expression == direction])
    effective <- intersect(as.character(universe), as.character(bg$TERM2GENE$GENE))
    genes <- intersect(genes, effective)
    if (length(genes)) {
      results$ORA <- clusterProfiler::enricher(gene = genes, universe = effective, TERM2GENE = bg$TERM2GENE, TERM2NAME = bg$TERM2NAME,
                                             pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = minGSSize, maxGSSize = maxGSSize)
    }
    results$foreground <- genes
    results$universe <- effective
  }
  if ("GSEA" %in% analysis_type) {
    df <- df[order(-abs(df$logFC), df$Geneid), , drop = FALSE]
    df <- df[!duplicated(df$ENTREZID), , drop = FALSE]
    ranks <- sort(stats::setNames(df$logFC, df$ENTREZID), decreasing = TRUE)
    results$ranking <- ranks
    if (length(ranks) >= minGSSize && any(bg$TERM2GENE$GENE %in% names(ranks))) {
      set.seed(cfg$seed)
      results$GSEA <- clusterProfiler::GSEA(geneList = ranks, TERM2GENE = bg$TERM2GENE, TERM2NAME = bg$TERM2NAME, eps = 0,
                                          pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = minGSSize, maxGSSize = maxGSSize,
                                          seed = TRUE, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    }
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(results, file.path(out_dir, "enrichment_objects.rds"))
    for (method in analysis_type) save_table(result_table(results[[method]]), file.path(out_dir, paste0(method, "_all.tsv")))
  }
  results
}

# 双向柱图统一使用FDR
plot_gsea_double_bar <- function(gsea_result, public_data = NULL, contrast_name = NULL, pvalue_col = "p.adjust",
                                 pvalue_signif = cfg$fdr_cutoff, color_up = cfg$direction_colors[["Up"]],
                                 color_down = cfg$direction_colors[["Down"]], color_not = cfg$direction_colors[["Not"]],
                                 out_dir = NULL, width = NULL, height = NULL) {
  df <- result_table(gsea_result, significant = TRUE)
  if (!nrow(df)) return(NULL)
  pos <- utils::head(df[df$NES > 0, ][order(df$p.adjust[df$NES > 0]), , drop = FALSE], cfg$top_terms)
  neg <- utils::head(df[df$NES < 0, ][order(df$p.adjust[df$NES < 0]), , drop = FALSE], cfg$top_terms)
  df <- rbind(pos, neg)
  df$score <- sign(df$NES) * -log10(pmax(df[[pvalue_col]], .Machine$double.xmin))
  df$Direction <- ifelse(df$NES > 0, "Up", "Down")
  df <- df[order(df$score), , drop = FALSE]
  df$TermLabel <- make.unique(paste0(df$Description, " [", df$ID, "]"))
  df$TermLabel <- factor(df$TermLabel, levels = df$TermLabel)
  p <- ggplot2::ggplot(df, ggplot2::aes(score, TermLabel, fill = Direction)) + ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) + ggplot2::scale_fill_manual(values = c(Up = color_up, Down = color_down)) +
    ggplot2::labs(x = "sign(NES) × −log10(FDR)", y = NULL, title = contrast_name) + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (!is.null(out_dir)) save_plot_pair(p, file.path(out_dir, "GSEA_double_bar"), width = width, height = height, n_terms = nrow(df))
  p
}

# 比较图复用已算结果
tranexp_object$set("public", "enrichment_results", NULL)
tranexp_object$set("public", "deg_results", NULL)
tranexp_object$set("public", "sample_info", NULL)
tranexp_object$set("public", "normalized_expression", NULL)
tranexp_object$set("public", "initialize", function(transcript_object) {
  self$design <- transcript_object$design
  self$public_data <- transcript_object$public_data
  self$enrichment_results <- transcript_object$enrichment
  self$deg_results <- transcript_object$deg
  self$sample_info <- transcript_object$sample_info
  self$normalized_expression <- transcript_object$normalized_expression
}, overwrite = TRUE)
tranexp_object$set("public", "plot_comparative", function(db, method, contrast_range = NULL, GO_term = NULL,
                                                         string_pattern_include = NULL, string_pattern_exclude = NULL,
                                                         out_dir = NULL, filter_term_pattern = NULL, selected_ids = NULL) {
  available <- any(vapply(self$enrichment_results, function(item) !is.null(item[[db]]), logical(1)))
  if (!available) return(empty_panel_plot(paste(db, method, "comparative enrichment"), "Background unavailable or enrichment not run"))
  tbl <- collect_comparative_tbl(self$enrichment_results, db, method, contrast_range)
  if (db == "GO" && (!is.null(GO_term) || !is.null(string_pattern_include))) {
    ids <- self$go_offspring(GO_term = GO_term, string_pattern_include = string_pattern_include,
                             string_pattern_exclude = string_pattern_exclude)
    tbl <- tbl[tbl$ID %in% ids, , drop = FALSE]
  } else if (!is.null(string_pattern_exclude) && nrow(tbl)) {
    tbl <- tbl[!grepl(paste(string_pattern_exclude, collapse = "|"), tbl$Description, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(filter_term_pattern) && nrow(tbl)) {
    tbl <- tbl[grepl(paste(filter_term_pattern, collapse = "|"), tbl$Description, ignore.case = TRUE), , drop = FALSE]
  }
  p <- plot_comparative_enrichment_plot(tbl, db, method, selected_ids)
  if (!is.null(out_dir) && !is.null(p)) save_plot_pair(p, file.path(out_dir, paste0(db, "_", method, "_comparative")),
                                                     n_terms = length(unique(p$data$ID)))
  p
})
tranexp_object$set("public", "compareCluster_GO_enrich", function(contrast_range = NULL) {
  collect_comparative_tbl(self$enrichment_results, "GO", "ORA", contrast_range)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_KEGG_enrich", function(contrast_range = NULL) {
  collect_comparative_tbl(self$enrichment_results, "KEGG", "ORA", contrast_range)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_GSEA_GO_enrich", function(contrast_range = NULL) {
  collect_comparative_tbl(self$enrichment_results, "GO", "GSEA", contrast_range)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_GSEA_KEGG_enrich", function(contrast_range = NULL) {
  collect_comparative_tbl(self$enrichment_results, "KEGG", "GSEA", contrast_range)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_GO_enrich_plot", function(contrast_range = NULL, ...) {
  self$plot_comparative("GO", "ORA", contrast_range, ...)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_KEGG_enrich_plot", function(contrast_range = NULL, ...) {
  self$plot_comparative("KEGG", "ORA", contrast_range, ...)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_GSEA_GO_enrich_plot", function(contrast_range = NULL, ...) {
  self$plot_comparative("GO", "GSEA", contrast_range, ...)
}, overwrite = TRUE)
tranexp_object$set("public", "compareCluster_GSEA_KEGG_enrich_plot", function(contrast_range = NULL, ...) {
  self$plot_comparative("KEGG", "GSEA", contrast_range, ...)
}, overwrite = TRUE)

# 旧热图改接当前结果
tranexp_object$set("public", "pheatmap_go_group2group", function(contrast_numbers = NULL, GO_term = NULL, string_pattern_include = NULL,
                                                               string_pattern_exclude = NULL, cytokines_chemokines_gene = FALSE, out_dir = NULL) {
  ids <- if (!is.null(GO_term) || !is.null(string_pattern_include)) self$go_offspring(GO_term, string_pattern_include, string_pattern_exclude) else NULL
  matrix <- build_effect_heatmap_mat(self$deg_results, self$public_data, contrast_numbers, ids)
  p <- plot_matrix_heatmap(matrix, "Gene × comparison: log2FC", scale_rows = FALSE)
  if (!is.null(out_dir)) save_plot_pair(p, file.path(out_dir, "GO_gene_logFC"), n_terms = nrow(matrix))
  list(matrix = matrix, plot = p)
}, overwrite = TRUE)
tranexp_object$set("public", "pheatmap_go_sample2sample", function(contrast_numbers = NULL, GO_term = NULL, string_pattern_include = NULL,
                                                                 string_pattern_exclude = NULL, cytokines_chemokines_gene = FALSE, out_dir = NULL) {
  effects <- self$pheatmap_go_group2group(contrast_numbers, GO_term, string_pattern_include, string_pattern_exclude)
  matrix <- self$normalized_expression[intersect(rownames(effects$matrix), rownames(self$normalized_expression)), , drop = FALSE]
  p <- plot_matrix_heatmap(matrix, "Sample expression: row Z score", scale_rows = TRUE, sample_info = self$sample_info)
  if (!is.null(out_dir)) save_plot_pair(p, file.path(out_dir, "GO_gene_expression"), n_terms = nrow(matrix))
  list(matrix = matrix, plot = p)
}, overwrite = TRUE)
