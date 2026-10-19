# 使用合成数据验证接口
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- normalizePath(file.path(dirname(sub("^--file=", "", arg[1])), ".."))
test_root <- tempfile("transcriptome_smoke_", tmpdir = Sys.getenv("TEST_WORK_DIR", tempdir()))
dir.create(test_root)
dir.create(file.path(test_root, "00.database"))
dir.create(file.path(test_root, "02.Quant_定量结果_Calicivirus"))
Sys.setenv(TRANSCRIPTOME_PROJECT_DIR = test_root)
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
input <- data.frame(Geneid = gene_ids, Length = 1000, counts, check.names = FALSE)
quant_dir <- file.path(test_root, "02.Quant_定量结果_Calicivirus")
write.table(input, file.path(quant_dir, "result_count.with_Length.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
metadata <- data.frame(Sample = samples, Group = rep(groups, each = 3))
write.table(metadata, file.path(quant_dir, "01_GroupSamples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
contrast <- do.call(rbind, lapply(c("06", "12", "24"), function(time) {
  data.frame(Treatment = paste0(c("Cat", "Leopard", "Tiger"), "_", time), Control = paste0("NC_", time))
}))
write.table(contrast, file.path(quant_dir, "03_GroupContrast.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
gtf <- vapply(seq_along(gene_ids), function(i) paste("chr1", "synthetic", "gene", sprintf("%d", i * 2000), sprintf("%d", i * 2000 + 999), ".", "+", ".",
  sprintf('gene_id "%s"; gene "%s"; db_xref "GeneID:%d"; description "Synthetic test gene";', gene_ids[i], gene_ids[i], 9000000 + i), sep = "\t"), character(1))
gtf_file <- file.path(test_root, "00.database", "GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.gtf.gz")
con <- gzfile(gtf_file, "wt"); writeLines(gtf, con); close(con)
go <- data.frame(tax_id = 9685, GeneID = 9000000 + seq_along(gene_ids), GO_ID = rep(c("GO:0008152", "GO:0009987"), 200),
                  Evidence = "IEA", Qualifier = "involved_in", GO_term = rep(c("metabolic process", "cellular process"), 200),
                  PubMed = "-", Category = "Process")
go <- rbind(go, transform(go[1, ], Qualifier = "NOT"), transform(go[1, ], tax_id = 9606))
con <- gzfile(file.path(test_root, "gene2go.gz"), "wt")
writeLines(paste0("#", paste(names(go), collapse = "\t")), con)
write.table(go, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE); close(con)
source(file.path(script_dir, "00_config.R"))
cfg$run_kegg <- FALSE
cfg$dpi <- 60L
cfg$top_genes_each <- 3L
for (file in c("01_input_annotation.R", "02_sample_qc.R", "03_differential_expression.R", "04_gene_sets.R", "05_expression_heatmaps.R")) {
  message("TEST MODULE: ", file)
  source(file.path(script_dir, file))
}
inputs <- read_module("01_inputs")
deg <- read_module("03_deg")
stopifnot(ncol(inputs$design$get_assay("counts")) == 36L, nrow(inputs$contrasts) == 9L)
stopifnot(identical(levels(inputs$sample_info$Group), cfg$group_levels))
stopifnot(identical(colnames(inputs$design$get_assay("counts")), as.character(inputs$sample_info$Samples)))
stopifnot(nrow(inputs$public_data$get_go_annotation()) == 400L)
stopifnot(all(abs(colSums(inputs$design$get_assay("tpm")) - 1e6) < 1e-6))
stopifnot(length(deg$comparisons) == 9L, all(deg$summary$Up > 0), all(deg$summary$Down > 0))
for (comparison in names(deg$comparisons)) {
  offset <- match(sub("_.*$", "", comparison), c("Cat", "Leopard", "Tiger")) * 30
  tbl <- deg$comparisons[[comparison]]$all
  positive <- sprintf("Gene%04d", offset + seq_len(15))
  negative <- sprintf("Gene%04d", offset + 15 + seq_len(15))
  stopifnot(all(tbl$logFC[match(positive, tbl$Geneid)] > 0), all(tbl$logFC[match(negative, tbl$Geneid)] < 0))
}
stopifnot(all(vapply(deg$comparisons, function(x) "baseMean" %in% names(x$all), logical(1))))
statuses <- data.table::fread(file.path(cfg$output_dir, "10_objects_and_summary", "module_status.tsv"))
if (any(statuses$Status == "FAILED")) stop("存在失败的核心图：", paste(statuses$Item[statuses$Status == "FAILED"], collapse = ", "))
capture.output(sessionInfo(), file = file.path(test_root, "core_test_sessionInfo.txt"))
capture.output(warnings(), file = file.path(test_root, "core_test_warnings.txt"))
cat("CORE_SMOKE_PASS\n", "TEST_ROOT=", test_root, "\n", sep = "")
