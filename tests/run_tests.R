if (!all(file.exists(file.path("tests", "fixtures", c("legacy_edgeR_expected.tsv", "legacy_enrichment_expected.tsv"))))) {
  stop("请先在代码根目录执行Rscript tests/generate_baselines.R，生成本地合成回归基线")
}
Sys.setenv(ENGINEERING_PROJECT_ROOT = normalizePath(getwd()))
for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(file)
testthat::test_dir("tests/unit", reporter = "summary", stop_on_failure = TRUE)
