# 仅在本地生成合成数据回归基线
project <- normalizePath(getwd())
fixtures <- file.path(project, "tests", "fixtures")
dir.create(fixtures, recursive = TRUE, showWarnings = FALSE)
set.seed(21)
groups <- as.vector(outer(c("Cat", "Leopard", "Tiger", "NC"), c("06", "12", "24"), paste, sep = "_"))
samples <- unlist(lapply(groups, function(group) paste0(group, "_", 1:3)))
gene_ids <- sprintf("Gene%04d", seq_len(400))
mu <- matrix(runif(400, 60, 400), nrow = 400, ncol = 36)
for (treatment in c("Cat", "Leopard", "Tiger")) {
  columns <- startsWith(samples, paste0(treatment, "_"))
  offset <- match(treatment, c("Cat", "Leopard", "Tiger")) * 30
  mu[offset + seq_len(15), columns] <- mu[offset + seq_len(15), columns] * 5
  mu[offset + 15 + seq_len(15), columns] <- mu[offset + 15 + seq_len(15), columns] / 5
}
counts <- matrix(rnbinom(length(mu), mu = as.vector(mu), size = 30), nrow = 400, dimnames = list(gene_ids, samples))
write.table(data.frame(Geneid = gene_ids, Length = 1000, counts, check.names = FALSE), file.path(fixtures, "counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Samples = samples, Group = rep(groups, each = 3)), file.path(fixtures, "samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gtf <- vapply(seq_along(gene_ids), function(i) paste("chr1", "synthetic", "gene", sprintf("%d", i * 2000), sprintf("%d", i * 2000 + 999), ".", "+", ".",
  sprintf('gene_id "%s"; gene "%s"; db_xref "GeneID:%d"; description "Synthetic test gene";',
          gene_ids[i], gene_ids[i], 9000000 + i), sep = "\t"), character(1))
writeLines(gtf, file.path(fixtures, "reference.gtf"))
go_ids <- c("GO:0006260", "GO:0006281", "GO:0006412", "GO:0006096", "GO:0005975", "GO:0006810")
go <- do.call(rbind, lapply(1:6, function(i) data.frame(tax_id = 9685, GeneID = 9000000 + 30 + (i - 1) * 15 + 1:15,
  GO_ID = go_ids[i], Evidence = "IEA", Qualifier = "involved_in", GO_term = paste("Synthetic BP", i), PubMed = "-", Category = "Process")))
con <- file(file.path(fixtures, "gene2go.tsv"), "wt")
writeLines(paste0("#", paste(names(go), collapse = "\t")), con)
write.table(go, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)
kegg_ids <- paste0("fca", c("00010", "00020", "00030", "00040", "00051", "00052"))
kegg <- data.frame(TERM = kegg_ids[match(go$GO_ID, go_ids)], GENE = as.character(go$GeneID))
write.table(kegg, file.path(fixtures, "KEGG_TERM2GENE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(TERM = kegg_ids, NAME = paste("Synthetic KEGG", 1:6)), file.path(fixtures, "KEGG_TERM2NAME.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
required <- file.path(fixtures, c("counts.tsv", "samples.tsv", "reference.gtf", "gene2go.tsv", "KEGG_TERM2GENE.tsv", "KEGG_TERM2NAME.tsv"))
if (any(!file.exists(required))) stop("请在包含tests/fixtures的代码根目录执行")
root <- tempfile("legacy_baseline_")
dir.create(root)
old_root <- Sys.getenv("TRANSCRIPTOME_PROJECT_DIR", unset = NA_character_)
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = root)
source("analysis_scripts/00_config.R")
cfg$count_file <- required[1]
cfg$sample_file <- required[2]
cfg$gtf_file <- required[3]
cfg$gene2go_file <- required[4]
cfg$contrast_file <- file.path(project, "analysis_scripts", "01_input", "03_GroupContrast.txt")
cfg$run_kegg <- FALSE
cfg$run_pathview <- FALSE
cfg$go_ancestors <- FALSE
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(cfg$seed)
# 直接运行旧模块，跳过与数值基线无关的绘图
attempt_plot <- function(...) NULL
save_plot_pair <- function(...) invisible(NULL)
combine_plots <- function(...) invisible(NULL)
source("analysis_scripts/01_input_annotation.R")
inputs <- read_module("01_inputs")
inputs$backgrounds$KEGG <- list(
  TERM2GENE = data.table::fread(required[5], colClasses = "character", data.table = FALSE),
  TERM2NAME = data.table::fread(required[6], colClasses = "character", data.table = FALSE))
save_module(inputs, "01_inputs")
source("analysis_scripts/03_differential_expression.R")
source("analysis_scripts/06_enrichment.R")
deg <- read_module("03_deg")
rows <- dplyr::bind_rows(lapply(names(deg$comparisons), function(name) {
  tbl <- deg$comparisons[[name]]$all
  data.frame(Comparison = name, tbl[, c("Geneid", "logFC", "p.value", "padj", "logCPM", "baseMean")])
}))
write.table(rows, file.path(fixtures, "legacy_edgeR_expected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
enrichment <- read_module("06_enrichment")
rows <- list()
for (comparison in names(enrichment$results)) {
  for (db in c("GO", "KEGG")) {
    item <- enrichment$results[[comparison]][[db]]
    objects <- list(ORA_Up = item$ORA$Up, ORA_Down = item$ORA$Down, GSEA = item$GSEA)
    for (kind in names(objects)) {
      object <- objects[[kind]]
      if (is.null(object) || !nrow(object@result)) next
      tbl <- as.data.frame(object@result)
      rows[[length(rows) + 1L]] <- data.frame(Comparison = comparison, Database = db, Kind = kind,
        ID = tbl$ID, pvalue = tbl$pvalue, p.adjust = tbl$p.adjust, NES = if (kind == "GSEA") tbl$NES else NA_real_)
    }
  }
}
write.table(dplyr::bind_rows(rows), file.path(fixtures, "legacy_enrichment_expected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (is.na(old_root)) Sys.unsetenv("TRANSCRIPTOME_PROJECT_DIR") else Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = old_root)
unlink(root, recursive = TRUE)
cat("SYNTHETIC_BASELINES_GENERATED\n")
