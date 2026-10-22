# transcript_pipeline_module.R 包依赖清单

本文档总结了 `transcript_pipeline_module.R` 文件所依赖的所有 R 包及其用途。

## 依赖包概览

该模块共依赖 **26** 个 R 包，分为以下几个类别：

---

## 1. 核心框架包

### R6
- **用途**: 提供面向对象编程（OOP）的 R6 类系统
- **使用场景**: 所有服务类的基础框架（IOService, DataModelService, AnnotationService 等）
- **安装**: `install.packages("R6")`

---

## 2. 数据结构与生物信息学基础包

### SummarizedExperiment
- **用途**: 提供 SummarizedExperiment 对象，用于存储基因表达数据及其元数据
- **使用场景**: DataModelService 中构建和管理转录组数据模型
- **安装**: `BiocManager::install("SummarizedExperiment")`

### S4Vectors
- **用途**: 提供 S4 向量和 DataFrame 类
- **使用场景**: 构建 SummarizedExperiment 的 rowData 和 colData
- **安装**: `BiocManager::install("S4Vectors")`

---

## 3. 数据处理与操作包

### data.table
- **用途**: 高效的数据表读取和处理
- **使用场景**: 
  - IOService: 读取 counts 文件
  - AnnotationService: 读取 GO 注释文件和 KEGG 数据
- **安装**: `install.packages("data.table")`

### dplyr
- **用途**: 数据操作（选择、过滤、连接、变换等）
- **使用场景**: 
  - AnnotationService: 注释数据处理
  - VisualizationService: 数据准备
- **安装**: `install.packages("dplyr")`

### tidyr
- **用途**: 数据整理和重塑
- **使用场景**: 
  - AnnotationService: 处理 KEGG ORTHOLOGY 数据
  - VisualizationService: pivot_wider 用于 upset 图数据准备
- **安装**: `install.packages("tidyr")`

### stringr
- **用途**: 字符串处理
- **使用场景**: AnnotationService 中处理 GTF 注释
- **安装**: `install.packages("stringr")`

### tibble
- **用途**: 现代化的数据框
- **使用场景**: AnnotationService 中构建注释表
- **安装**: `install.packages("tibble")`

### reshape2
- **用途**: 数据重塑（宽格式和长格式转换）
- **使用场景**: VisualizationService 中准备箱线图和小提琴图数据
- **安装**: `install.packages("reshape2")`

---

## 4. 统计分析包

### edgeR
- **用途**: RNA-seq 差异表达分析（基于负二项分布）
- **使用场景**: DEAnalysisService 中进行差异表达分析和标准化
- **安装**: `BiocManager::install("edgeR")`

### DESeq2
- **用途**: RNA-seq 差异表达分析（基于负二项分布）
- **使用场景**: DEAnalysisService 中进行差异表达分析和标准化
- **安装**: `BiocManager::install("DESeq2")`

---

## 5. 注释与数据库包

### AnnotationDbi
- **用途**: 访问生物学注释数据库的通用接口
- **使用场景**: AnnotationService 中查询 OrgDb 数据库
- **安装**: `BiocManager::install("AnnotationDbi")`

### AnnotationHub
- **用途**: 访问 Bioconductor AnnotationHub 资源
- **使用场景**: AnnotationService 中通过 Hub 查询注释
- **安装**: `BiocManager::install("AnnotationHub")`

### rtracklayer
- **用途**: 导入和导出基因组注释文件（GTF/GFF）
- **使用场景**: AnnotationService 中解析 NCBI GTF 文件
- **安装**: `BiocManager::install("rtracklayer")`

### GO.db
- **用途**: Gene Ontology (GO) 注释数据库
- **使用场景**: GO 富集分析背景数据
- **安装**: `BiocManager::install("GO.db")`

---

## 6. 富集分析包

### clusterProfiler
- **用途**: 基因富集分析（GSEA, ORA）
- **使用场景**: 
  - EnrichmentService: GSEA 和 enricher 分析
  - AnnotationService: 下载 KEGG 数据库
- **安装**: `BiocManager::install("clusterProfiler")`

### enrichplot
- **用途**: 富集分析结果可视化
- **使用场景**: 富集分析结果的绘图支持
- **安装**: `BiocManager::install("enrichplot")`

### pathview
- **用途**: KEGG 通路可视化
- **使用场景**: EnrichmentService 中可视化 KEGG 通路
- **安装**: `BiocManager::install("pathview")`

---

## 7. 可视化包

### ggplot2
- **用途**: 强大的图形绘制系统
- **使用场景**: 
  - VisualizationService: 火山图、MA 图、PCA 图、箱线图、小提琴图
- **安装**: `install.packages("ggplot2")`

### ggrepel
- **用途**: 避免标签重叠的文本标注
- **使用场景**: VisualizationService 中火山图的基因标签
- **安装**: `install.packages("ggrepel")`

### ggsci
- **用途**: 科研期刊配色方案
- **使用场景**: VisualizationService 中箱线图和小提琴图的颜色
- **安装**: `install.packages("ggsci")`

### pheatmap
- **用途**: 绘制热图
- **使用场景**: VisualizationService 中样本相关性热图
- **安装**: `install.packages("pheatmap")`

### ComplexUpset
- **用途**: 绘制 UpSet 图（集合交集可视化）
- **使用场景**: VisualizationService 中差异基因集合的交集展示
- **安装**: `install.packages("ComplexUpset")`

### scales
- **用途**: 图形坐标轴刻度和颜色映射
- **使用场景**: 图形美化支持
- **安装**: `install.packages("scales")`

### forcats
- **用途**: 因子变量处理
- **使用场景**: 图形中因子级别管理
- **安装**: `install.packages("forcats")`

### patchwork
- **用途**: 组合多个 ggplot2 图形
- **使用场景**: 图形排版和组合
- **安装**: `install.packages("patchwork")`

---

## 8. 文件输入输出包

### openxlsx
- **用途**: 读写 Excel 文件（.xlsx）
- **使用场景**: 
  - IOService: 导出表格
  - DEAnalysisService: 保存差异分析结果
- **安装**: `install.packages("openxlsx")`

---

## 物种特异性注释包（按需安装）

根据分析的物种，可能需要安装以下 OrgDb 包之一：

- **org.Hs.eg.db** (人类, hsa): `BiocManager::install("org.Hs.eg.db")`
- **org.Mm.eg.db** (小鼠, mmu): `BiocManager::install("org.Mm.eg.db")`
- **org.Rn.eg.db** (大鼠, rno): `BiocManager::install("org.Rn.eg.db")`
- **org.Dr.eg.db** (斑马鱼, dre): `BiocManager::install("org.Dr.eg.db")`
- **org.Dm.eg.db** (果蝇, dme): `BiocManager::install("org.Dm.eg.db")`
- **org.Ce.eg.db** (线虫, cel): `BiocManager::install("org.Ce.eg.db")`
- **org.Sc.sgd.db** (酵母, sce): `BiocManager::install("org.Sc.sgd.db")`
- **org.At.tair.db** (拟南芥, ath): `BiocManager::install("org.At.tair.db")`
- **org.EcK12.eg.db** (大肠杆菌, eco): `BiocManager::install("org.EcK12.eg.db")`
- **org.Bs.eg.db** (枯草芽孢杆菌, bsu): `BiocManager::install("org.Bs.eg.db")`
- **org.Pf.plasmo.db** (疟原虫, pfa): `BiocManager::install("org.Pf.plasmo.db")`

---

## 批量安装脚本

### 安装 CRAN 包

```r
cran_packages <- c(
  "R6", "data.table", "stringr", "dplyr", "tidyr", "tibble", "reshape2",
  "ggplot2", "ggrepel", "ggsci", "pheatmap", "ComplexUpset", 
  "scales", "forcats", "patchwork", "openxlsx"
)

install.packages(cran_packages)
```

### 安装 Bioconductor 包

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c(
  "SummarizedExperiment", "S4Vectors", "edgeR", "DESeq2",
  "clusterProfiler", "enrichplot", "pathview",
  "AnnotationDbi", "AnnotationHub", "rtracklayer", "GO.db"
)

BiocManager::install(bioc_packages)
```

### 安装特定物种的注释包（示例：人类）

```r
BiocManager::install("org.Hs.eg.db")
```

---

## 依赖检查

在使用模块前，可以运行以下代码检查所有依赖是否已安装：

```r
required_packages <- c(
  # CRAN
  "R6", "data.table", "stringr", "dplyr", "tidyr", "tibble", "reshape2",
  "ggplot2", "ggrepel", "ggsci", "pheatmap", "ComplexUpset", 
  "scales", "forcats", "patchwork", "openxlsx",
  # Bioconductor
  "SummarizedExperiment", "S4Vectors", "edgeR", "DESeq2",
  "clusterProfiler", "enrichplot", "pathview",
  "AnnotationDbi", "AnnotationHub", "rtracklayer", "GO.db"
)

missing_packages <- required_packages[!sapply(required_packages, 
  function(pkg) requireNamespace(pkg, quietly = TRUE))]

if (length(missing_packages) > 0) {
  cat("缺少以下包，请安装：\n")
  cat(paste(missing_packages, collapse = "\n"), "\n")
} else {
  cat("所有必需的包都已安装！\n")
}
```

---

## 包分类汇总表

| 类别 | 包数量 | 包名称 |
|------|--------|--------|
| 核心框架 | 1 | R6 |
| 数据结构 | 2 | SummarizedExperiment, S4Vectors |
| 数据处理 | 6 | data.table, dplyr, tidyr, stringr, tibble, reshape2 |
| 统计分析 | 2 | edgeR, DESeq2 |
| 注释数据库 | 4 | AnnotationDbi, AnnotationHub, rtracklayer, GO.db |
| 富集分析 | 3 | clusterProfiler, enrichplot, pathview |
| 可视化 | 7 | ggplot2, ggrepel, ggsci, pheatmap, ComplexUpset, scales, forcats, patchwork |
| 文件IO | 1 | openxlsx |
| **总计** | **26** | |

---

## 版本要求

- **R 版本**: >= 4.0.0（推荐 4.2.0 或更高）
- **Bioconductor 版本**: >= 3.14（推荐最新版本）

---

## 注意事项

1. **Bioconductor 包**: 必须通过 `BiocManager::install()` 安装，不能使用 `install.packages()`
2. **物种特异性包**: 根据实际分析的物种选择安装对应的 OrgDb 包
3. **网络依赖**: AnnotationHub 和 clusterProfiler 的某些功能需要网络访问
4. **依赖关系**: 某些包会自动安装其依赖包，实际安装的包可能更多

---

**最后更新**: 2025-11-06  
**文件版本**: 1.0
