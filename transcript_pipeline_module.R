# transcript_pipeline_module.R
# transcript_pipeline_module.R
# 完整实现：把你原始 transcript_object 中的功能按模块化服务实现，
# 使用 SummarizedExperiment 作为数据模型 (DataModelService)，并提供
# Annotation / Normalization / DE / Enrichment / Visualization 服务，以及
# 一个轻量 PipelineHelper 供交互式使用。
#
# 保存后在 R 会话中使用： source("transcript_pipeline_module.R")
#
# 说明：
# - 我尽量照搬并调整你原脚本中的逻辑：读取、TPM/FPKM、edgeR/DESeq2 差异、注释解析、
#   KEGG/GO 背景构建、GSEA/enricher、常用绘图（PCA/heatmap/volcano/MA/venn）。
# - 方法签名偏向“纯函数式”：大多数计算方法返回对象或矩阵，调用者决定是否 attach 到 SE。
# - 你可以逐步替换/修改某个服务内的实现（符合你要求的“层依赖、易替换”）。
#
# 依赖（按需安装）: R6, SummarizedExperiment, data.table, stringr, dplyr, tidyr, ggplot2,
# edgeR, DESeq2, clusterProfiler, enrichplot, pathview, AnnotationDbi, AnnotationHub, rtracklayer,
# GO.db, pheatmap, ComplexUpset, openxlsx, ggrepel, ggsci, scales, forcats, patchwork
#
# 作者: based on user's original script, refactored into services
# Date: 2025-10-30

# --- Helpers -----------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# minimal pkg-check helper
.ensure_pkg <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) stop("缺少 R 包: ", paste(missing, collapse = ", "), ". 请先安装。")
}

# ---------------- IO Service -------------------------------------------------
IOService <- R6::R6Class("IOService",
  public = list(
    read_counts = function(path) stop("not implemented"),
    read_groups = function(path) stop("not implemented"),
    read_contrasts = function(path) stop("not implemented"),
    write_table = function(df, path) stop("not implemented"),
    save_rds = function(obj, path) saveRDS(obj, path)
  )
)

DefaultIOService <- R6::R6Class("DefaultIOService",
  inherit = IOService,
  public = list(
    read_counts = function(path) {
      .ensure_pkg(c("data.table"))
      df <- data.table::fread(path, data.table = FALSE)
      # expect first col geneid, a Length column (named Length) or second column numeric
      df
    },
    read_groups = function(path) {
      df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("Samples", "Group"))
      df
    },
    read_contrasts = function(path) {
      df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("GroupExperiment", "GroupControl"))
      df
    },
    write_table = function(df, path) {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(df, path, rowNames = TRUE)
      } else {
        write.table(df, file = path, sep = "\t", quote = FALSE, row.names = TRUE)
      }
    }
  )
)

# ---------------- Data Model Service ----------------------------------------
DataModelService <- R6::R6Class("DataModelService",
  public = list(
    se = NULL,
    build_se = function(counts_df, groups_df, length_col = "Length") {
      if (is.matrix(counts_df)) counts_df <- as.data.frame(counts_df, stringsAsFactors = FALSE)
      if (!is.data.frame(counts_df)) stop("counts_df must be data.frame or matrix-like")
      geneid <- as.character(counts_df[[1]])
      if (length_col %in% colnames(counts_df)) {
        length_vec <- as.numeric(counts_df[[length_col]])
        expr <- counts_df[, setdiff(colnames(counts_df), c(colnames(counts_df)[1], length_col)), drop = FALSE]
      } else {
        # assume second column is length if numeric
        if (ncol(counts_df) >= 3 && is.numeric(counts_df[[2]])) {
          length_vec <- as.numeric(counts_df[[2]])
          expr <- counts_df[, -(1:2), drop = FALSE]
        } else stop("找不到 Length 列，请保证 counts 文件第一列基因 ID，第二列为 Length 或指定 length_col")
      }
      rownames(expr) <- geneid;

      # sanitize sample names both in counts and groups to match featureCounts/perl cleaning
      original_colnames <- colnames(expr)
      sanitized_colnames <- make.names(original_colnames, unique = TRUE)
      if (!identical(original_colnames, sanitized_colnames)) {
        colnames(expr) <- sanitized_colnames;
      }

      # ensure sample order matches groups_df
      # sanitize groups_df$Samples similarly but keep original mapping
      groups_df <- as.data.frame(groups_df, stringsAsFactors = FALSE)
      groups_df$Samples_original <- groups_df$Samples;
      groups_df$Samples <- make.names(groups_df$Samples, unique = TRUE);

      missing_in_counts <- setdiff(groups_df$Samples, colnames(expr))
      missing_in_groups <- setdiff(colnames(expr), groups_df$Samples);
      if (length(missing_in_counts) > 0 || length(missing_in_groups) > 0) {
        stop("样本名称不匹配:\n",
             if (length(missing_in_counts) > 0) paste0("  在 GroupSamples 中但不在 counts 中: ", paste(missing_in_counts, collapse = ", "), "\n") else "",
             if (length(missing_in_groups) > 0) paste0("  在 counts 中但不在 GroupSamples 中: ", paste(missing_in_groups, collapse = ", "), "\n") else "",
             "建议检查特征计数文件和 GroupSamples 文件中的样本命名一致性（可以使用 make.names 规范化）")
      }

      expr <- expr[, groups_df$Samples, drop = FALSE];

      # build colData and set rownames to sample names
      coldf <- groups_df[match(colnames(expr), groups_df$Samples), , drop = FALSE];
      rownames(coldf) <- coldf$Samples;
      # restore original sample names in a column to keep mapping
      if ("Samples_original" %in% colnames(coldf)) {
        coldf$Samples_pretty <- coldf$Samples_original;
      }

      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(expr)),
        rowData = S4Vectors::DataFrame(Length = length_vec),
        colData = S4Vectors::DataFrame(coldf)
      )
      self$se <- se;
      invisible(se);
    },
    attach_assay = function(name, mat) {
      SummarizedExperiment::assay(self$se, name) <- mat;
      invisible(self$se);
    },
    get_assay = function(name) {
      SummarizedExperiment::assay(self$se, name);
    },
    save_se = function(path) {
      saveRDS(self$se, path);
    }
  )
);

# ---------------- Annotation Service ---------------------------------------
AnnotationService <- R6::R6Class("AnnotationService",
  public = list(
    gtf_annotation = list(),
    kegg_organism_code = NULL,
    parse_GTF_NCBI_file = function(file_path, do_perl_clean = FALSE) {
      if (!file.exists(file_path)) stop(file_path, " 不存在");
      .ensure_pkg(c("rtracklayer", "tibble", "dplyr", "stringr"));
      # Native R processing instead of perl system call
      gtf_source <- file_path;
      if (do_perl_clean) {
        mod <- paste0(file_path, ".modify");
        if (!file.exists(mod)) {
          # Read and clean in R instead of using perl system call
          lines <- readLines(file_path);
          # Remove db_xref entries that don't contain GeneID
          lines <- gsub('db_xref "(?!GeneID)[^"]+"; ', '', lines, perl = TRUE);
          writeLines(lines, mod);
        }
        gtf_source <- mod;
      }
      # Parse GTF with rtracklayer
      gtf <- rtracklayer::import(gtf_source, format = "gtf");
      df <- tibble::tibble(
        gene_id = as.character(gtf$gene_id),
        gene_name = as.character(gtf$gene),
        db_xref = as.list(gtf$db_xref),
        description = as.character(gtf$description),
        product = as.character(gtf$product)
      );
      ENTREZID <- vapply(df$db_xref, FUN.VALUE = character(1), FUN = function(x) {
        if (is.null(x)) return(NA_character_);
        v <- grep("GeneID:", x, value = TRUE);
        if (length(v) == 0) return(NA_character_);
        sub(".*GeneID:(\\d+).*", "\\1", v[1]);
      });
      out <- tibble::tibble(Geneid = df$gene_id, SYMBOL = df$gene_name, ENTREZID = ENTREZID, GENENAME = ifelse(is.na(df$description) | df$description == "", df$product, df$description));
      out <- dplyr::distinct(out);
      self$gtf_annotation$source <- "NCBI_GTF";
      self$gtf_annotation$tbl_gtf <- out;
      invisible(out);
    },
    parse_GTF_AnnotationHub_query = function(query = NULL, update = FALSE) {
      .ensure_pkg(c("AnnotationHub", "AnnotationDbi", "dplyr", "tibble"));
      hub <- if (update) AnnotationHub::AnnotationHub(localHub = FALSE) else AnnotationHub::AnnotationHub(localHub = TRUE);
      orgDB <- hub[[query]];
      df <- AnnotationDbi::select(orgDB, keys = AnnotationDbi::keys(orgDB, keytype = "SYMBOL"), columns = c("ENTREZID", "GENENAME", "GO"), keytype = "SYMBOL");
      df <- dplyr::select(df, SYMBOL, ENTREZID, GENENAME, GO);
      df$Geneid <- df$SYMBOL;
      df <- dplyr::distinct(df);
      self$gtf_annotation$source <- "AnnotationHub";
      self$gtf_annotation$tbl_gtf <- df;
      self$gtf_annotation$tbl_GO_tbl <- df %>% dplyr::select(SYMBOL, ENTREZID, GO);
      invisible(df);
    },
    parse_GTF_Orgdb = function(kegg_organism_code) {
      .ensure_pkg(c("AnnotationDbi"));
      Orgdb <- private$get_orgdb(kegg_organism_code);
      if (is.null(Orgdb)) stop("无法获取 OrgDb");
      df <- AnnotationDbi::select(Orgdb, keys = AnnotationDbi::keys(Orgdb, keytype = "SYMBOL"),
                                  columns = c("ENTREZID", "GENENAME"), keytype = "SYMBOL");
      df <- dplyr::select(df, SYMBOL, ENTREZID, GENENAME);
      df$Geneid <- df$SYMBOL;
      df <- dplyr::distinct(df);
      self$gtf_annotation$source <- "OrgDb";
      self$gtf_annotation$tbl_gtf <- df;
      invisible(df);
    },
    parse_kegg = function(kegg_organism_code) {
      .ensure_pkg(c("clusterProfiler", "dplyr"));
      kegg_database <- clusterProfiler::download_KEGG(species = kegg_organism_code);
      KEGGPATHID2EXTID <- kegg_database$KEGGPATHID2EXTID;
      KEGGPATHID2NAME <- kegg_database$KEGGPATHID2NAME;
      colnames(KEGGPATHID2EXTID) <- c("Keggpathid", "ENTREZID");
      colnames(KEGGPATHID2NAME) <- c("Keggpathid", "Keggpathname");
      merged <- dplyr::left_join(KEGGPATHID2EXTID, KEGGPATHID2NAME, by = "Keggpathid");
      merged$Keggpathname <- gsub(" - .*", "", merged$Keggpathname);
      self$gtf_annotation$tbl_kegg_tbl <- merged;
      invisible(merged);
    },
    parse_go = function(kegg_organism_code = NULL, gene_ontology_file_path = NULL) {
      .ensure_pkg(c("dplyr"));
      Orgdb <- NULL;
      if (!is.null(kegg_organism_code)) Orgdb <- private$get_orgdb(kegg_organism_code);
      if (!is.null(Orgdb)) {
        df <- AnnotationDbi::select(Orgdb, keys = AnnotationDbi::keys(Orgdb, keytype = "SYMBOL"),
                                    columns = c("ENTREZID", "GO"), keytype = "SYMBOL");
        df <- dplyr::select(df, SYMBOL, ENTREZID, GO);
      } else if (!is.null(gene_ontology_file_path) && file.exists(gene_ontology_file_path)) {
        .ensure_pkg("data.table");
        df <- data.table::fread(gene_ontology_file_path, data.table = FALSE);
        df <- dplyr::rename(df, SYMBOL = Symbol, ENTREZID = GeneID, GO = GO_ID);
      } else if (file.exists("/data/shenlab_data/cuixinyuan/00.Database/14.Human_genome/gene2go.RDS")) {
        tmp <- readRDS("/data/shenlab_data/cuixinyuan/00.Database/14.Human_genome/gene2go.RDS");
        df <- dplyr::select(tmp, ENTREZID = GeneID, GO = GO_ID) %>% dplyr::mutate(ENTREZID = as.character(ENTREZID));
        df <- dplyr::left_join(self$gtf_annotation$tbl_gtf %>% dplyr::select(ENTREZID, SYMBOL), df, by = "ENTREZID");
        df <- dplyr::select(df, SYMBOL, ENTREZID, GO);
      } else stop("无法获得 GO 注释，请指定 kegg_organism_code 或 gene_ontology_file_path 或提供本地缓存");
      df <- dplyr::mutate(df, ENTREZID = as.character(ENTREZID));
      self$gtf_annotation$tbl_GO_tbl <- df;
      invisible(df);
    },
    inject_annotation_to_se = function(data_model_service, gtf_tbl) {
      if (is.null(data_model_service$se)) stop("SE 为空");
      rd <- as.data.frame(SummarizedExperiment::rowData(data_model_service$se));
      ann <- gtf_tbl[match(rownames(data_model_service$se), gtf_tbl$Geneid), , drop = FALSE];
      ann2 <- ann[, setdiff(colnames(ann), "Geneid"), drop = FALSE];
      rd <- cbind(rd, ann2);
      SummarizedExperiment::rowData(data_model_service$se) <- S4Vectors::DataFrame(rd, check.names = FALSE);
      invisible(data_model_service$se);
    },
    update_SYMBOL2ORTHOLOGY = function() {
      .ensure_pkg(c("data.table", "tidyr"));
      tbl <- data.table::fread("https://rest.kegg.jp/list/ko", col.names = c("ORTHOLOGY", "SYMBOL"), sep = "\t");
      # split and tidy
      if (requireNamespace("tidyr", quietly = TRUE)) {
        tmp <- tidyr::separate_wider_delim(tbl, cols = SYMBOL, delim = "; ", names = c("SYMBOL", "SYMBOLNAME"), too_many = "drop", too_few = "align_end");
        tmp <- tidyr::separate_rows(tmp, SYMBOL, sep = ", ");
        res <- setNames(tmp$ORTHOLOGY, tmp$SYMBOL);
      } else {
        res <- setNames(tbl$ORTHOLOGY, tbl$SYMBOL);
      }
      private$SYMBOL2ORTHOLOGY <- res;
      invisible(res);
    }
  ),
  private = list(
    get_orgdb = function(kegg_organism_code) {
      mapping <- list(hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db", rno = "org.Rn.eg.db", dre = "org.Dr.eg.db",
                      dme = "org.Dm.eg.db", cel = "org.Ce.eg.db", sce = "org.Sc.sgd.db", ath = "org.At.tair.db",
                      eco = "org.EcK12.eg.db", bsu = "org.Bs.eg.db", pfa = "org.Pf.plasmo.db");
      if (!kegg_organism_code %in% names(mapping)) {
        message("OrgDb 不支持该 KEGG 代码: ", kegg_organism_code);
        return(NULL);
      }
      pkg <- mapping[[kegg_organism_code]];
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) stop("需要 BiocManager 安装 OrgDb");
        BiocManager::install(pkg);
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE));
      get(pkg);
    },
    SYMBOL2ORTHOLOGY = NULL
  )
);

# ---------------- Normalization Service ------------------------------------
NormalizationService <- R6::R6Class("NormalizationService",
  public = list(
    calc_FPKM = function(se) {
      counts <- SummarizedExperiment::assay(se, "counts");
      len <- SummarizedExperiment::rowData(se)$Length;
      lib <- colSums(counts);
      fpkm <- sweep(counts * 1e9, 2, lib, FUN = "/") / len;
      fpkm;
    },
    calc_TPM = function(se) {
      counts <- SummarizedExperiment::assay(se, "counts");
      len <- SummarizedExperiment::rowData(se)$Length;
      rpk <- sweep(counts, 1, len, FUN = "/");
      tpm <- sweep(rpk, 2, colSums(rpk, na.rm = TRUE), FUN = "/") * 1e6;
      tpm;
    },
    normalize_edgeR = function(mat_counts) {
      .ensure_pkg("edgeR");
      mat_counts <- round(mat_counts);
      dge <- edgeR::DGEList(counts = mat_counts);
      dge <- edgeR::calcNormFactors(dge);
      edgeR::cpm(dge, normalized.lib.sizes = TRUE);
    },
    normalize_DESeq2 = function(mat_counts, coldata) {
      .ensure_pkg("DESeq2");
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(mat_counts), colData = coldata, design = ~ conditions);
      dds <- DESeq2::DESeq(dds, quiet = TRUE);
      DESeq2::counts(dds, normalized = TRUE);
    }
  )
);

# ---------------- DE Analysis Service --------------------------------------
DEAnalysisService <- R6::R6Class("DEAnalysisService",
  public = list(
    run_edgeR_contrast = function(se, group_exp, group_ctrl, out_dir = NULL, filter_cpm = 1) {
      .ensure_pkg("edgeR");
      gs <- as.data.frame(SummarizedExperiment::colData(se));
      samples_exp <- rownames(gs)[gs$Group == group_exp];
      samples_ctrl <- rownames(gs)[gs$Group == group_ctrl];
      if (length(samples_exp) == 0 || length(samples_ctrl) == 0) stop("某组样本为空");
      counts <- round(SummarizedExperiment::assay(se, "counts")[, c(samples_ctrl, samples_exp), drop = FALSE]);
      keep <- rowSums(edgeR::cpm(counts) > filter_cpm) >= min(length(samples_exp), length(samples_ctrl)) / 2;
      counts_f <- counts[keep, , drop = FALSE];
      if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(counts_f, file.path(out_dir, paste0("read_counts_edgeR_", group_exp, "_versus_", group_ctrl, ".xlsx")), rowNames = TRUE);
      conditions <- factor(c(rep(group_ctrl, length(samples_ctrl)), rep(group_exp, length(samples_exp))));
      exp_study <- edgeR::DGEList(counts = counts_f, group = conditions);
      exp_study <- edgeR::calcNormFactors(exp_study);
      if (length(samples_ctrl) > 1 && length(samples_exp) > 1) {
        exp_study <- edgeR::estimateDisp(exp_study);
      } else {
        exp_study$common.dispersion <- 0.1;
      }
      et <- edgeR::exactTest(exp_study, pair = c(group_ctrl, group_exp));
      tTags <- edgeR::topTags(et, n = NULL);
      result_table <- tTags$table;
      norm_cpm <- edgeR::cpm(exp_study);
      # Compute mean samples efficiently
      mean_sample <- data.frame(
        mean_SampleA = rowMeans(norm_cpm[, samples_exp, drop = FALSE]),
        mean_SampleB = rowMeans(norm_cpm[, samples_ctrl, drop = FALSE])
      );
      result_table <- base::merge(mean_sample, result_table, by = "row.names");
      rownames(result_table) <- result_table[, 1];
      result_table <- result_table[, -1];
      result_table <- data.frame(sampleA = group_exp, sampleB = group_ctrl, result_table);
      # save normalized cpm if requested
      if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(norm_cpm, file.path(out_dir, paste0("read_counts_edgeR_normalized_", group_exp, "_versus_", group_ctrl, ".xlsx")), rowNames = TRUE);
      invisible(list(tbl = result_table, norm = norm_cpm));
    },
    run_DESeq2_contrast = function(se, group_exp, group_ctrl, out_dir = NULL, mean_filter = 1) {
      .ensure_pkg("DESeq2");
      gs <- as.data.frame(SummarizedExperiment::colData(se));
      samples_exp <- rownames(gs)[gs$Group == group_exp];
      samples_ctrl <- rownames(gs)[gs$Group == group_ctrl];
      counts <- round(SummarizedExperiment::assay(se, "counts")[, c(samples_ctrl, samples_exp), drop = FALSE]);
      keep <- rowMeans(counts) > mean_filter;
      counts_f <- counts[keep, , drop = FALSE];
      if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(counts_f, file.path(out_dir, paste0("read_counts_DESeq2_", group_exp, "_versus_", group_ctrl, ".xlsx")), rowNames = TRUE);
      colData <- data.frame(row.names = c(colnames(counts_f)), conditions = factor(c(rep(group_ctrl, length(samples_ctrl)), rep(group_exp, length(samples_exp))), levels = c(group_ctrl, group_exp)));
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_f, colData = colData, design = ~ conditions);
      dds <- DESeq2::DESeq(dds, quiet = TRUE);
      res <- DESeq2::results(dds);
      norm_counts <- DESeq2::counts(dds, normalized = TRUE);
      # Compute means efficiently
      meanA <- rowMeans(norm_counts[, samples_exp, drop = FALSE]);
      meanB <- rowMeans(norm_counts[, samples_ctrl, drop = FALSE]);
      res_df <- as.data.frame(res);
      res_df$Geneid <- rownames(res_df);
      meanA_sub <- as.numeric(meanA[rownames(res_df)]);
      meanB_sub <- as.numeric(meanB[rownames(res_df)]);
      result_table <- cbind(sampleA = group_exp, sampleB = group_ctrl,
                            mean_SampleA = meanA_sub, mean_SampleB = meanB_sub,
                            res_df);
      rownames(result_table) <- result_table$Geneid;
      result_table$Geneid <- NULL;
      if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(norm_counts, file.path(out_dir, paste0("read_counts_DESeq2_normalized_", group_exp, "_versus_", group_ctrl, ".xlsx")), rowNames = TRUE);
      invisible(list(tbl = result_table, norm = norm_counts));
    }
  )
);

# ---------------- Enrichment Service ---------------------------------------
EnrichmentService <- R6::R6Class("EnrichmentService",
  public = list(
    run_GSEA = function(ranks, TERM2GENE, TERM2NAME = NULL, out_dir = NULL, minGSSize = 10, maxGSSize = 500) {
      .ensure_pkg("clusterProfiler");
      res <- clusterProfiler::GSEA(geneList = ranks, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize, eps = 0);
      if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE);
      if (!is.null(out_dir)) saveRDS(res, file = file.path(out_dir, "GSEA_result.rds"));
      invisible(res);
    },
    run_enricher = function(gene_vector, TERM2GENE, TERM2NAME = NULL, out_dir = NULL) {
      .ensure_pkg("clusterProfiler");
      res <- clusterProfiler::enricher(gene = gene_vector, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME);
      if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE);
      if (!is.null(out_dir)) saveRDS(res, file = file.path(out_dir, "enricher_result.rds"));
      invisible(res);
    },
    run_pathview_for_kegg = function(gene_list_named_logFC, kegg_id, species = NULL, out_dir = ".") {
      .ensure_pkg("pathview");
      tmpd <- tempfile("pathview_");
      dir.create(tmpd);
      oldwd <- getwd();
      on.exit({
        try(setwd(oldwd), silent = TRUE);
      }, add = TRUE);
      setwd(tmpd);
      pv <- try(pathview::pathview(gene.data = gene_list_named_logFC,
                                   pathway.id = kegg_id,
                                   species = species,
                                   kegg.dir = tmpd,
                                   same.layer = FALSE), silent = TRUE);
      if (inherits(pv, "try-error")) {
        message("pathview failed for ", kegg_id, ": ", conditionMessage(attr(pv, "condition")));
        return(NULL);
      }
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE);
      generated <- list.files(tmpd, pattern = "pathview|\.pathview|\.kegg", full.names = TRUE, ignore.case = TRUE);
      if (length(generated) == 0) {
        generated <- list.files(tmpd, pattern = paste0(kegg_id, ".*(png|svg)$"), full.names = TRUE);
      }
      moved <- character(0);
      for (f in generated) {
        dest <- file.path(out_dir, basename(f));
        tryCatch({
          file.rename(f, dest);
          moved <- c(moved, dest);
        }, error = function(e) {
          message("failed to move ", f, " -> ", dest, ": ", e$message);
        });
      }
      invisible(list(pathview_result = pv, moved_files = moved));
    }
  )
);

# ---------------- Visualization Service ------------------------------------
VisualizationService <- R6::R6Class("VisualizationService",
  public = list(
    plot_volcano = function(de_tbl, threshold_padj = 0.05, threshold_Log2FC = 1, label_top = 20) {
      .ensure_pkg(c("ggplot2", "ggrepel", "dplyr"));
      df <- as.data.frame(de_tbl);
      if (!"padj" %in% colnames(df)) df$padj <- stats::p.adjust(df$p.value);
      df$threshold <- "Stale";
      df$threshold[df$padj < threshold_padj & df$logFC > threshold_Log2FC] <- "Up";
      df$threshold[df$padj < threshold_padj & df$logFC < -threshold_Log2FC] <- "Down";
      p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(padj), color = threshold)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c("Down" = "blue4", "Stale" = "grey", "Up" = "red3")) +
        ggplot2::geom_vline(xintercept = c(-threshold_Log2FC, threshold_Log2FC), lty = 2) +
        ggplot2::geom_hline(yintercept = -log10(threshold_padj), lty = 2) +
        ggplot2::theme_minimal();
      // label top by -log10(padj)
      if (label_top > 0) {
        lab <- df[df$threshold != "Stale", ];
        if (nrow(lab) > 0) {
          lab <- lab[order(-abs(lab$logFC) * -log10(lab$padj)), ];
          lab <- head(lab, label_top);
          if ("SYMBOL" %in% colnames(df)) {
            p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = SYMBOL), size = 3);
          } else {
            p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = rownames(lab)), size = 3);
          }
        }
      }
      p;
    },
    plot_MA = function(de_tbl) {
      .ensure_pkg("ggplot2");
      df <- as.data.frame(de_tbl);
      # Efficiently compute A (average expression)
      if (!("baseMeanLog2" %in% colnames(df))) {
        if (all(c("mean_SampleA", "mean_SampleB") %in% colnames(df))) {
          df$A <- (df$mean_SampleA + df$mean_SampleB) * 0.5;
        } else df$A <- 0;
      } else {
        df$A <- df$baseMeanLog2;
      }
      p <- ggplot2::ggplot(df, ggplot2::aes(x = A, y = logFC)) + ggplot2::geom_point(alpha = 0.4) + ggplot2::theme_minimal() + ggplot2::xlab("Average expression") + ggplot2::ylab("log2 Fold Change");
      p;
    },
    plot_pca = function(se, assay = "TPM", color_by = "Group") {
      .ensure_pkg("ggplot2");
      if (!(assay %in% SummarizedExperiment::assayNames(se))) stop("Assay not found");
      mat <- SummarizedExperiment::assay(se, assay);
      mat_log <- t(log2(mat + 1));
      pca <- stats::prcomp(mat_log, center = TRUE, scale. = TRUE);
      df <- as.data.frame(pca$x[, 1:2]);
      df$Sample <- rownames(df);
      coldata <- as.data.frame(SummarizedExperiment::colData(se));
      df <- merge(df, coldata, by.x = "Sample", by.y = "row.names", all.x = TRUE);
      ggplot2::ggplot(df, ggplot2::aes_string(x = "PC1", y = "PC2", color = color_by)) + ggplot2::geom_point(size = 3) + ggplot2::theme_classic();
    },
    plot_sample_corr_heatmap = function(se, assay = "TPM", out_file = NULL) {
      .ensure_pkg("pheatmap");
      mat <- SummarizedExperiment::assay(se, assay);
      cor_m <- cor(mat, method = "pearson", use = "pairwise.complete.obs");
      ph <- pheatmap::pheatmap(cor_m, main = paste0(assay, " sample correlation"));
      if (!is.null(out_file)) {
        png(out_file, width = 800, height = 800); pheatmap::pheatmap(cor_m, main = paste0(assay, " sample correlation")); dev.off();
      }
      invisible(ph);
    },
    plot_samples_box_violin = function(tpm_mat, groups_df, out_dir = NULL) {
      .ensure_pkg(c("reshape2", "ggplot2", "ggsci", "dplyr"));
      dfm <- reshape2::melt(log2(tpm_mat + 1));
      colnames(dfm) <- c("Geneid", "Sample", "value");
      dfm <- dplyr::left_join(dfm, groups_df, by = c("Sample" = "Samples"));
      boxplot <- ggplot2::ggplot(dfm, ggplot2::aes(x = Sample, y = value, color = Group)) +
        ggplot2::geom_boxplot() + ggsci::scale_color_aaas() + ggplot2::theme_classic();
      violin <- ggplot2::ggplot(dfm, ggplot2::aes(x = Sample, y = value, fill = Group)) + ggplot2::geom_violin() + ggsci::scale_fill_frontiers() + ggplot2::theme_classic();
      if (!is.null(out_dir)) {
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE);
        try({ ggplot2::ggsave(file.path(out_dir, "Samples_boxplot.png"), boxplot, width = 8, height = 4) }, silent = TRUE);
        try({ ggplot2::ggsave(file.path(out_dir, "Samples_violin.png"), violin, width = 8, height = 4) }, silent = TRUE);
      }
      list(boxplot = boxplot, violin = violin);
    },
    plot_deg_venn_upset = function(significant_list, out_dir = NULL) {
      .ensure_pkg("ComplexUpset");
      // significant_list: named list of data.frames each containing Geneid, Result_expression
      df_all <- lapply(names(significant_list), function(nm) {
        d <- significant_list[[nm]];
        if (!"Geneid" %in% colnames(d)) stop("significant table must contain Geneid");
        d <- d[d$Result_expression == "Up", , drop = FALSE];
        data.frame(Geneid = d$Geneid, group = nm, present = TRUE, stringsAsFactors = FALSE);
      }) |> do.call(rbind, .);
      if (nrow(df_all) == 0) return(NULL);
      library(tidyr); library(dplyr);
      wide <- tidyr::pivot_wider(df_all, id_cols = Geneid, names_from = group, values_from = present, values_fill = FALSE);
      p <- ComplexUpset::upset(wide, colnames(wide)[-1]);
      if (!is.null(out_dir)) {
        try({ ggplot2::ggsave(file.path(out_dir, "DEG_upset.png"), p, width = 8, height = 6) }, silent = TRUE);
      }
      p;
    }
  )
);

# ---------------- Pipeline Helper ------------------------------------------
PipelineHelper <- R6::R6Class("PipelineHelper",
  public = list(
    io = NULL, model = NULL, annot = NULL, norm = NULL, de = NULL, enrich = NULL, viz = NULL,
    initialize = function(io = NULL, model = NULL, annot = NULL, norm = NULL, de = NULL, enrich = NULL, viz = NULL) {
      self$io <- io %||% DefaultIOService$new();
      self$model <- model %||% DataModelService$new();
      self$annot <- annot %||% AnnotationService$new();
      self$norm <- norm %||% NormalizationService$new();
      self$de <- de %||% DEAnalysisService$new();
      self$enrich <- enrich %||% EnrichmentService$new();
      self$viz <- viz %||% VisualizationService$new();
    },
    build_se = function(counts_path = NULL, groups_path = NULL, length_col = "Length") {
      // If paths are not provided, try to use My_QUANT_DIR environment variable (set by upstream bash scripts)
      quant_dir <- Sys.getenv("My_QUANT_DIR", unset = "02.Quant_定量结果");
      if (is.null(counts_path) || counts_path == "") counts_path <- file.path(quant_dir, "result_count.with_Length.tsv");
      if (is.null(groups_path) || groups_path == "") groups_path <- file.path(quant_dir, "GroupSamples.txt");

      counts_df <- self$io$read_counts(counts_path);
      groups_df <- self$io$read_groups(groups_path);
      se <- self$model$build_se(counts_df, groups_df, length_col = length_col);
      invisible(se);
    },
    compute_and_attach_TPM = function() {
      se <- self$model$se;
      tpm <- self$norm$calc_TPM(se);
      self$model$attach_assay("TPM", tpm);
      invisible(tpm);
    },
    compute_and_attach_FPKM = function() {
      se <- self$model$se;
      fpkm <- self$norm$calc_FPKM(se);
      self$model$attach_assay("FPKM", fpkm);
      invisible(fpkm);
    },
    run_edgeR_by_contrast_index = function(contrast_index) {
      if (is.null(self$model$se)) stop("请先 build_se()");
      gc <- self$io$read_contrasts; // not used; expecting user to have set GroupContrast externally if needed
      // We assume user has attached GroupContrast as colData groups in model$se
      // Here we accept group names explicitly
      stop("请使用 run_edgeR_for_groups(group_exp, group_ctrl) 方法以传入两个组名");
    },
    run_edgeR_for_groups = function(group_exp, group_ctrl, out_dir = NULL) {
      res <- self$de$run_edgeR_contrast(self$model$se, group_exp, group_ctrl, out_dir = out_dir);
      // attach normalized matrix
      self$model$attach_assay(paste0("edgeR_CPM_", group_exp, "_vs_", group_ctrl), res$norm);
      invisible(res);
    },
    run_DESeq2_for_groups = function(group_exp, group_ctrl, out_dir = NULL) {
      res <- self$de$run_DESeq2_contrast(self$model$se, group_exp, group_ctrl, out_dir = out_dir);
      self$model$attach_assay(paste0("DESeq2_norm_", group_exp, "_vs_", group_ctrl), res$norm);
      invisible(res);
    },
    annotate_with_gtf = function(gtf_file, method = c("ncbi", "orgdb", "annotationhub"), ...) {
      method <- match.arg(method);
      if (method == "ncbi") {
        tbl <- self$annot$parse_GTF_NCBI_file(gtf_file);
      } else if (method == "orgdb") {
        tbl <- self$annot$parse_GTF_Orgdb(...);
      } else {
        tbl <- self$annot$parse_GTF_AnnotationHub_query(...);
      }
      self$annot$inject_annotation_to_se(self$model, tbl);
      invisible(tbl);
    },
    run_GSEA_from_deg = function(de_tbl, kind = c("KEGG", "GO"), out_dir = NULL) {
      kind <- match.arg(kind);
      // de_tbl should contain SYMBOL and logFC
      ranks <- de_tbl[order(de_tbl$logFC, decreasing = TRUE), ];
      rnk <- structure(ranks$logFC, names = ranks$SYMBOL);
      if (kind == "KEGG") {
        bg <- self$model; // user should have self$annot$gtf_annotation populated
        TERM2 <- self$annot$gtf_annotation$tbl_kegg_tbl;
        if (is.null(TERM2)) stop("缺少 KEGG 背景，请先 parse_kegg()");
        t2g <- data.frame(TERM = TERM2$Keggpathid, GENE = TERM2$ENTREZID);
        res <- self$enrich$run_GSEA(rnk, TERM2GENE = t2g, out_dir = out_dir);
      } else {
        gobg <- self$annot$gtf_annotation$tbl_GO_tbl;
        if (is.null(gobg)) stop("缺少 GO 注释，请先 parse_go()");
        t2g <- data.frame(TERM = gobg$GO, GENE = gobg$SYMBOL);
        res <- self$enrich$run_GSEA(rnk, TERM2GENE = t2g, out_dir = out_dir);
      }
      invisible(res);
    },
    plot_volcano = function(de_tbl, out_file = NULL) {
      p <- self$viz$plot_volcano(de_tbl);
      if (!is.null(out_file)) ggplot2::ggsave(out_file, p, width = 8, height = 6);
      invisible(p);
    },
    plot_pca = function(assay = "TPM", out_file = NULL, color_by = "Group") {
      p <- self$viz$plot_pca(self$model$se, assay = assay, color_by = color_by);
      if (!is.null(out_file)) ggplot2::ggsave(out_file, p, width = 8, height = 6);
      invisible(p);
    }
  )
);

# ---------------- Example Usage --------------------------------------------
# Example (interactive):
# source("Xinyuan_Cui_transcript_pipeline_module.R")
# ph <- PipelineHelper$new();
# se <- ph$build_se();
# ph$compute_and_attach_TPM();
# res <- ph$run_edgeR_for_groups("G1","G2", out_dir = "04.Differential_差异分析结果");
