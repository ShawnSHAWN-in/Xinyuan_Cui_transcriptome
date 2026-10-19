# 在项目顶层运行
cran_packages <- c("R6", "tidyverse", "data.table", "openxlsx", "ggplot2", "patchwork", "pheatmap", "ggrepel", "R.utils", "digest", "ggpubr", "ggupset")
bioc_packages <- c("SummarizedExperiment", "S4Vectors", "rtracklayer", "AnnotationDbi", "GO.db", "edgeR", "clusterProfiler", "enrichplot", "KEGGREST")
optional_packages <- c("UpSetR", "ggdendro", "scatterplot3d", "pathview", "png")
environment_table <- data.frame(Package = c(cran_packages, bioc_packages, optional_packages),
                                Role = c(rep("CRAN", length(cran_packages)), rep("Bioconductor", length(bioc_packages)),
                                         rep("Optional", length(optional_packages))))
environment_table$Installed <- vapply(environment_table$Package, requireNamespace, logical(1), quietly = TRUE)
environment_table$Version <- vapply(environment_table$Package, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
}, character(1))
print(environment_table, row.names = FALSE)
if (any(!environment_table$Installed & environment_table$Role != "Optional")) {
  message("请先按README安装缺少的基础或富集包；此文件不会自动安装软件")
}
if (requireNamespace("dplyr", quietly = TRUE) && utils::packageVersion("dplyr") < "1.1.1") stop("需要dplyr >= 1.1.1")
