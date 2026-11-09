# 合成数据验证注释端口
source("analysis_scripts/lib/eggnog_annotation.R")
source("R/io.R")
source("R/annotation.R")
source("R/statistics.R")
source("R/export.R")
root <- tempfile("eggnog_test_")
dir.create(root)
map_file <- file.path(root, "protein_gene.tsv")
names_file <- file.path(root, "pathway_names.tsv")
data.table::fwrite(data.frame(query = c("p1", "p2", "p3", "p4"), Geneid = c("gene-A", "gene-A", "gene-B", "gene-C")), map_file, sep = "\t")
data.table::fwrite(data.frame(TERM = c("map00010", "map00020"), NAME = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")),
  names_file, sep = "\t")
make_input <- function(file, version = 2L) {
  tbl <- data.frame(query = c("p1", "p2", "p3"), GOs = c("GO:0008152,GO:0003674", "GO:0008152", "GO:0008152"),
    KEGG_Pathway = c("map00010,ko00010", "ko00010,map00020", "-"), KEGG_ko = c("ko:K00844", "ko:K00844", "-"),
    Preferred_name = c("same_symbol", "same_symbol", "same_symbol"), stringsAsFactors = FALSE)
  if (version == 3L) tbl <- tbl[, c("query", "Preferred_name", "KEGG_ko", "GOs", "KEGG_Pathway")]
  con <- if (grepl("gz$", file)) gzfile(file, "wt") else base::file(file, "wt")
  writeLines(c(paste0("## emapper-", version, ".test"), paste0("#", paste(names(tbl), collapse = "\t"))), con)
  write.table(tbl, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines("## 3 queries scanned", con)
  close(con)
  file
}
v2 <- make_input(file.path(root, "v2.emapper.annotations"))
v3 <- make_input(file.path(root, "v3.emapper.annotations.gz"), 3L)
build <- function(file = v2, ...) build_eggnog_annotation_lst(file, map_file, species = "Ailurus fulgens",
  include_ancestors = FALSE, kegg_term2name_file = names_file, ...)

testthat::test_that("v2/v3按列名读取并保留来源", {
  a <- build()
  b <- build(v3)
  testthat::expect_equal(a$backgrounds, b$backgrounds)
  testthat::expect_true(any(grepl("emapper-2.test", a$metadata$header)))
  testthat::expect_true(any(grepl("emapper-3.test", b$metadata$header)))
})

testthat::test_that("蛋白合并为基因且map/ko关联不重复", {
  value <- build()
  testthat::expect_equal(value$backgrounds$KEGG$TERM2GENE, data.frame(TERM = c("map00010", "map00020"), GENE = "gene-A"))
  testthat::expect_setequal(value$background$TERM2GENE$GENE, c("gene-A", "gene-B"))
  testthat::expect_false("GO:0003674" %in% value$background$TERM2GENE$TERM)
  testthat::expect_setequal(value$annotation$ENTREZID, c("gene-A", "gene-B", "gene-C"))
  testthat::expect_false(value$query_mapping$HasAnnotationRow[value$query_mapping$query == "p4"])
  testthat::expect_equal(value$coverage$Annotated, c(2L, 1L))
  testthat::expect_equal(value$coverage$Total, c(3L, 3L))
  testthat::expect_identical(value$metadata$id_type, "Geneid")
  testthat::expect_false(any(grepl("^K[0-9]", value$backgrounds$KEGG$TERM2GENE$GENE)))
})

testthat::test_that("GO祖先展开不丢失直接注释", {
  value <- build_eggnog_annotation_lst(v2, map_file, species = "Ailurus fulgens", include_ancestors = TRUE,
    kegg_term2name_file = names_file)
  testthat::expect_true(all(c("GO:0008152", "GO:0008150") %in% value$background$TERM2GENE$TERM))
})

testthat::test_that("缺失映射和歧义映射明确停止", {
  testthat::expect_error(build_eggnog_annotation_lst(v2, species = "Ailurus fulgens"), "映射表")
  testthat::expect_error(read_eggnog_gene_map_tbl(map_file, c("x", "y")), "完全不匹配")
  bad <- file.path(root, "bad_map.tsv")
  data.table::fwrite(data.frame(query = c("p1", "p1"), Geneid = c("a", "b")), bad, sep = "\t")
  testthat::expect_error(read_eggnog_gene_map_tbl(bad, "p1"), "只能对应一个")
  identity <- build_eggnog_annotation_lst(v2, query_is_gene_id = TRUE, species = "Ailurus fulgens",
    kegg_term2name_file = names_file, include_ancestors = FALSE)
  testthat::expect_setequal(identity$background$TERM2GENE$GENE, c("p1", "p2", "p3"))
  partial <- file.path(root, "partial_map.tsv")
  data.table::fwrite(data.frame(query = "p1", Geneid = "a"), partial, sep = "\t")
  testthat::expect_warning(build_eggnog_annotation_lst(v2, partial, species = "Ailurus fulgens", kegg_term2name_file = names_file), "缺少基因映射")
})

testthat::test_that("错误文件或通路名称不伪造背景", {
  bad <- file.path(root, "bad.tsv")
  writeLines(c("#query\tseed_ortholog", "p1\tref"), bad)
  testthat::expect_error(read_eggnog_annotation_tbl(bad), "列名")
  writeLines(c("#query\tGOs\tKEGG_Pathway", "p1\tBAD\tmap00010"), bad)
  testthat::expect_error(build(bad), "非法GO")
  writeLines(c("#query\tGOs\tKEGG_Pathway", "p1\tGO:0008152\tK00844"), bad)
  testthat::expect_error(build(bad), "非map/ko")
  names_bad <- file.path(root, "bad_names.tsv")
  writeLines(c("TERM\tNAME", "map00010\tGlycolysis"), names_bad)
  testthat::expect_error(build_eggnog_annotation_lst(v2, map_file, species = "Ailurus fulgens", kegg_term2name_file = names_bad), "缺少通路")
  testthat::expect_equal(nrow(build(run_kegg = FALSE)$backgrounds$KEGG$TERM2GENE), 0L)
})

testthat::test_that("输出适配单细胞端口与既有bulk富集ID", {
  value <- build()
  paths <- save_eggnog_reference_paths(value, file.path(root, "out"))
  testthat::expect_true(all(file.exists(paths)))
  sets <- readRDS(file.path(root, "out/gene_sets.rds"))
  testthat::expect_named(sets, c("GO_BP", "KEGG"))
  testthat::expect_named(sets$GO_BP, c("term2gene", "term2name", "species", "id_type", "source"))
  testthat::expect_named(sets$KEGG$term2gene, c("term_id", "gene_id"))
  fit <- list(job = list(id = "test"), method = "edgeR", table = data.frame(Geneid = c("gene-A", "gene-B"), logFC = c(2, -2), padj = 0.01))
  annotated <- attach_gene_annotation_tbl(fit, value$annotation)
  selected <- select_deg_lst(annotated, list(fdr = 0.05, abs_logfc = 1))
  testthat::expect_setequal(selected$universe, c("gene-A", "gene-B"))
  snapshot <- build_kegg_snapshot_paths(list(annotation_source = "eggnog", run_kegg = TRUE, require_kegg = TRUE, kegg_refresh = "test"),
    character(), file.path(root, "snapshot"), value)
  testthat::expect_true(all(file.exists(snapshot)))
  testthat::expect_setequal(data.table::fread(snapshot[1])$GENE, "gene-A")
  testthat::expect_true(all(file.exists(export_annotation_files(value, file.path(root, "export")))))
})

testthat::test_that("PublicDataManager旧入口可分派eggNOG而不访问物种KEGG", {
  env <- new.env(parent = globalenv())
  env$`%>%` <- magrittr::`%>%`
  env$R6Class <- R6::R6Class
  definitions <- parse("analysis_scripts/lib/Rscript_R6_transcriptome_base.R")
  for (expr in definitions) {
    if (is.call(expr) && as.character(expr[[1]])[1] %in% c("<-", "=")) eval(expr, env)
  }
  manager <- env$PublicDataManager$new(cache_dir = file.path(root, "cache"))
  manager$load_combined_gtf(eggnog = v2, query_gene_file = map_file, species = "Ailurus fulgens",
    kegg_term2name_file = names_file, include_ancestors = FALSE)
  testthat::expect_setequal(manager$get_gene_annotation()$ENTREZID, c("gene-A", "gene-B", "gene-C"))
  testthat::expect_equal(lapply(manager$get_kegg_background(), as.data.frame), lapply(build()$backgrounds$KEGG, as.data.frame))
  testthat::expect_identical(manager$get_annotation_metadata()$id_type, "Geneid")
  testthat::expect_equal(manager$get_go_background()$TERM2GENE, build()$background$TERM2GENE)
  definitions <- parse("analysis_scripts/lib/Rscript_R6_transcriptome_class.R")
  for (expr in definitions) {
    if (is.call(expr) && as.character(expr[[1]])[1] %in% c("<-", "=")) eval(expr, env)
  }
  env$cfg <- list(output_dir = root, go_ancestors = FALSE)
  sys.source("analysis_scripts/lib/compatibility.R", env)
  manager <- env$PublicDataManager$new(cache_dir = file.path(root, "cache_compat"))
  manager$load_project_annotation(eggnog = v3, query_gene_file = map_file, species = "Ailurus fulgens",
    kegg_term2name_file = names_file, include_ancestors = FALSE)
  testthat::expect_equal(manager$get_go_background()$TERM2GENE, build()$background$TERM2GENE)
  testthat::expect_equal(lapply(manager$get_kegg_background(), as.data.frame), lapply(build()$backgrounds$KEGG, as.data.frame))
  value <- read_project_annotation_lst(c(eggnog_file = v3, query_gene_file = map_file, kegg_names_file = names_file),
    list(annotation_source = "eggnog", reference_species = "Ailurus fulgens", query_is_gene_id = FALSE,
      go_ontology = "BP", go_ancestors = FALSE, run_kegg = TRUE), "analysis_scripts/lib/eggnog_annotation.R")
  testthat::expect_equal(value$backgrounds, build()$backgrounds)
})
unlink(root, recursive = TRUE)
message("eggNOG合成验证完成；未运行真实物种数据")
