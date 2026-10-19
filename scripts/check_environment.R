packages <- c("targets", "yaml", "renv", "R6", "tidyverse", "data.table", "edgeR", "SummarizedExperiment", "S4Vectors", "BiocManager",
              "rtracklayer", "AnnotationDbi", "GO.db", "clusterProfiler", "enrichplot", "BiocParallel", "ggplot2", "patchwork",
              "pheatmap", "ggdendro", "ggupset", "openxlsx", "digest", "KEGGREST", "png")
installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(Package = packages, Installed = installed, Version = vapply(packages, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
}, character(1))), row.names = FALSE)
if (any(!installed)) stop("缺少依赖；按README安装或使用renv::restore()，流程不会自动安装包")
