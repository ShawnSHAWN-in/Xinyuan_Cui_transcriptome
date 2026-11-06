# 杯状病毒项目转录组分析

本目录按模块复用原有 R6 数据容器、差异分析和绘图函数。原始脚本保存在 `99_original`；`lib` 是项目适配版本，输入文件和服务器上的原有轮子不会被覆盖。

## 放置位置

把 `analysis_scripts` 整个目录放到项目顶层，与 `00.database`、`02.Quant_定量结果_Calicivirus` 和 `gene2go.gz` 并列。路径默认相对于这个位置计算。

已有文件：

- `00.database/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.gtf.gz`
- `gene2go.gz`
- `02.Quant_定量结果_Calicivirus/result_count.with_Length.tsv`
- `02.Quant_定量结果_Calicivirus/01_GroupSamples.tsv`
- `02.Quant_定量结果_Calicivirus/03_GroupContrast.txt`

样本来自狮子的说明与家猫参考注释分别记录。参考注释使用 NCBI tax_id 9685 和 KEGG fca；KEGG 的 fca 页面对应 GCF_018350175.1：<https://www.kegg.jp/entry/fca>。

## 依赖检查

从项目顶层进入 R/RStudio，首先运行：

```r
source("analysis_scripts/00_check_environment.R")
```

代码不会自动安装包。缺少时可自行运行以下安装语句；Bioconductor 使用与当前 R 匹配的版本，避免混用不同版本的 ggtree/enrichplot/ggplot2。

```r
install.packages(c("R6", "tidyverse", "data.table", "openxlsx", "ggplot2", "patchwork", "pheatmap", "ggrepel", "R.utils", "digest", "ggpubr", "ggupset"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "S4Vectors", "rtracklayer", "AnnotationDbi", "GO.db", "edgeR", "clusterProfiler", "enrichplot", "KEGGREST"))
# 特殊图的可选依赖
install.packages(c("ggdendro", "ggpubr", "scatterplot3d", "png"))
BiocManager::install("pathview")
```

R >= 4.1；dplyr >= 1.1.1。编译依赖的安装权限和系统库由服务器环境决定。`clusterProfiler` 与 `enrichplot` 在进入06模块前必须可载入；缺少时会停止，不会用其他算法悄悄替代。

## 按文件逐步执行

```r
source("analysis_scripts/00_config.R")
source("analysis_scripts/01_input_annotation.R")
```

第一步核对 `03_analysis/01_inputs` 下的样本信息和比较信息，以及 `02_annotation` 下的映射覆盖率、未匹配基因和背景。确认36样本、12组合组、9比较及注释对应后继续：

```r
source("analysis_scripts/02_sample_qc.R")
source("analysis_scripts/03_differential_expression.R")
source("analysis_scripts/04_gene_sets.R")
source("analysis_scripts/05_expression_heatmaps.R")
source("analysis_scripts/06_enrichment.R")
source("analysis_scripts/07_enrichment_plots.R")
source("analysis_scripts/08_comparative_enrichment.R")
source("analysis_scripts/09_combine_export.R")
```

每个模块读取前置RDS，不要求把所有中间变量一直留在同一R会话。重启会话后先source配置，再source需要的模块即可。批量入口：

```bash
Rscript analysis_scripts/run_all.R > analysis_run.log 2>&1
```

如果脚本目录不放在项目顶层，可在启动R前设置 `TRANSCRIPTOME_PROJECT_DIR` 为项目绝对路径；配置文件本身不改变工作目录。

## 文件与职责

| 文件 | 内容 |
|---|---|
| `00_config.R` | 全部路径、factor、阈值、配色、GO/KEGG与绘图参数 |
| `01_input_annotation.R` | 输入适配、R6对象、GTF/gene2go/KEGG、覆盖率 |
| `02_sample_qc.R` | TMM/logCPM、TPM/FPKM导出、表达箱线图、四种PCA、相关性及聚类 |
| `03_differential_expression.R` | 9组差异分析、火山/排序/MA单图与3×3组合 |
| `04_gene_sets.R` | 数量/比例柱图、两类三集合Venn、UpSet、花朵图 |
| `05_expression_heatmaps.R` | 样本行Z热图、比较logFC热图及可选BP主题热图 |
| `06_enrichment.R` | GO-BP/KEGG各自ORA和GSEA；保存完整对象和表格 |
| `07_enrichment_plots.R` | 富集气泡、柱图、UpSet、网络、曲线及Pathview |
| `08_comparative_enrichment.R` | 四类上下三角比较图、ORA/GSEA配对、NES矩阵 |
| `09_combine_export.R` | 跨比较曲线、Pathview栅格拼图、图像索引与运行汇总 |
| `lib/compatibility.R` | 对原有R6接口的项目适配与错误修正 |
| `lib/plot_helpers.R` | 特殊图适配、图对象组合及统一PDF/PNG输出 |

## 输入兼容处理

01模块接受当前的 `Sample/Group` 样本表，并另存为基座要求的 `Samples/Group`。比较表的 `Treatment/Control` 表头会去掉，另存无表头的两列文件。原始文件保留不动。

压缩文件先检查完整性。GTF使用压缩流解压到临时文件后解析；标准gene2go按块读取，筛选参考tax_id与GTF GeneID，排除NOT注释，保留证据类型，默认保留IEA。

GO仅BP，默认向BP祖先条目传播关联；参数 `go_ancestors` 可关闭，但改变后需从01模块重跑。多义映射不自动取第一条：富集只使用可明确对应的Geneid/GeneID，保留未匹配和多义记录。原始count行不因此被删除。

## 统计和方向

- 当前模型按两组独立重复处理，未包含个体或批次协变量。
- 首轮默认edgeR，沿用原轮子的TMM + GLM/LRT路线与比较内CPM过滤；不是新建另一套差异分析框架。
- 默认FDR < 0.05且|log2FC| >= 1；显著DEG、火山图、MA和ORA统一参数。
- 同时间处理组相对NC；正logFC为处理组更高。
- 样本质控的PCA使用统一全样本TMM/logCPM及表达过滤，默认不进行基因单位方差缩放。该输入选择与参数会保留在结果中。
- ORA使用显著上调/下调基因分别计算，受检GeneID作为背景；有效注释背景另存。
- GSEA使用全部有效受检基因的logFC排序；不先筛显著DEG，不先拆方向。保存实际rank列表。
- 富集计算保留返回的全部有效条目，绘图按FDR筛选；BH校正的实际范围按每个分析结果解释，不称为跨所有比较的全局校正。
- GO/KEGG ORA三角图的方向来自DEG集合；GSEA方向来自NES。不能直接写成过程激活/抑制。
- 输入或阈值发生变化时，需重跑对应前置模块；GO/KEGG背景改变需重跑01和06之后的模块。富集检查点以输入、差异结果校验和及背景内容区分。

## 绘图与可选项

- 图像保持固定factor顺序；聚类版相关图按聚类排列，另有固定设计顺序版本。
- 每组3重复，四风格PCA的椭圆版会省略不满足绘图最低样本数的椭圆，并明确说明。
- 默认PDF和300dpi PNG；组合图按时间/处理拆面板。原始KEGG通路图是栅格，嵌入PDF后仍是栅格。
- ORA气泡/柱图默认Top12；GSEA曲线每比较默认Top3，可修改 `top_terms`、`curve_top_each`。结果表保留完整内容。
- 主题热图/比较图默认不预设主题；在 `theme_go_ids` 填BP父项后启用。
- Pathview默认尝试，每对象最多Top2通路，统一所有比较的logFC色阶；缺包或网络错误单独记录。
- 花朵图明确展示各集合独有计数和全部集合交集；UpSet保存完整交集计数和成员表。
- 山脊图、GSVA/ssGSEA、条目相似性网络、其他差异算法敏感性、时间交互模型属于规划的可选扩展，本批未启用。

## 结果与失败信息

结果位于 `03_analysis`。各模块统计对象保存在 `10_objects_and_summary`，组合图集中在 `09_combined_figures`。

重点查看：`module_status.tsv`、`figure_inventory.tsv`、`analysis_parameters.tsv`、`run_summary.md`。无显著结果、背景不可用、绘图失败、未配置专题分别记录；不能仅以目录里有图片判断全部分析成功。

如果必须GO和KEGG都完成，把 `require_kegg` 设为TRUE。默认FALSE时，KEGG网络失败会记录并允许GO及差异模块继续，后续可恢复KEGG背景后重跑。

## 验证范围

`tests/smoke_core.R` 使用合成36样本数据测试真实R6接口、压缩GTF、无Symbol的标准gene2go、NOT/物种过滤、9个edgeR比较及核心图像；不含本项目真实测序数据。测试结果和可用包版本以附带 `VALIDATION.md` 为准。

本代码包不意味着已在用户服务器运行；真实数据结果需要在该服务器按模块执行后检查。
