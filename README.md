# Xinyuan_Cui_transcriptome

> A blog-style introduction and usage guide for the transcriptome upstream pipeline  
> (download → QC → align → quantify → DE / optional de novo assembly)

---

![banner](https://img.shields.io/badge/Transcriptome-Pipeline-blue)

作者：Xinyuan Cui  
最后更新：2025-10-30

---

## 写在前面（博客开放式简介）

欢迎来到 Xinyuan_Cui_transcriptome —— 这是一个轻量级、模块化的转录组上游分析流水线，覆盖从数据下载、质控、比对到 featureCounts 计数的常见步骤。它的设计灵感来自于真实项目工作流：清晰、可复用，并允许你在不同步骤中灵活替换工具或参数。

本文档以“博客风”方式呈现：先讲为什么这样做、再讲如何使用、最后给出进阶与部署建议。

---

## 核心特点（为什么用它）

- 模块化：三个主要脚本拆分职责（环境、工具库、主流水线），便于维护和扩展。
- 并行友好：使用 GNU parallel 加速样本级任务。
- 可配置：支持按物种自动选择索引 / 可手动指定 GTF；featureCounts 的 -t/-g 可自由指定。
- 便于集成：可放入 HPC、CI 或容器环境中运行。

---

## 快速开始（3 分钟上手）

1. 在本地克隆仓库（假设仓库已创建为 `ShawnSHAWN-in/Xinyuan_Cui_transcriptome`）：
   - gh 方式（推荐）：
     - gh repo clone ShawnSHAWN-in/Xinyuan_Cui_transcriptome
   - 或 git 方式：
     - git clone https://github.com/ShawnSHAWN-in/Xinyuan_Cui_transcriptome.git

2. 进入目录并赋予执行权限：
   - chmod +x transcriptome_pipeline.sh
   - chmod +x lib_transcriptome.sh
   - chmod +x Script_Transcript_Upstream_0_source_environment.sh

3. 运行示例（对一个样本列表 list_file）：
   - ./transcriptome_pipeline.sh all -l list_file -o human -t 40

4. 如果只希望运行计数并自定义 featureCounts 参数：
   - ./transcriptome_pipeline.sh count -l list_file -a path/to/genome.gtf -t 12 -F exon -G gene_id

---

## 文件说明（仓库文件速览）

- Script_Transcript_Upstream_0_source_environment.sh — 环境与目录变量
- lib_transcriptome.sh — 常用函数（下载、hisat2 流程、featureCounts、Trinity）
- transcriptome_pipeline.sh — 主流水线，子命令化组织（download/preprocess/align/count/de/trinity/all）
- transcript_pipeline_module.R — R 语言模块化转录组分析框架（DE/富集/可视化）
- transcript_pipeline_module_dependencies.md — R 模块的包依赖清单与安装指南
- README.md — 你正在阅读的博客风说明

---

## R 语言分析模块（下游分析）

本仓库还提供了一个模块化的 R 语言转录组分析框架 `transcript_pipeline_module.R`，用于下游分析：

### 主要功能
- 数据模型：基于 SummarizedExperiment 管理表达数据
- 标准化：支持 TPM, FPKM, edgeR, DESeq2 多种方法
- 差异分析：edgeR 和 DESeq2 差异表达分析
- 注释服务：支持 NCBI GTF, OrgDb, AnnotationHub
- 富集分析：GSEA, ORA (KEGG/GO)
- 可视化：火山图、MA图、PCA、热图、箱线图、UpSet图等

### 快速使用示例

```r
# 加载模块
source("transcript_pipeline_module.R")

# 创建 pipeline helper
ph <- PipelineHelper$new()

# 构建数据模型
se <- ph$build_se(
  counts_path = "02.Quant_定量结果/result_count.with_Length.tsv",
  groups_path = "02.Quant_定量结果/GroupSamples.txt"
)

# 计算 TPM 和 FPKM
ph$compute_and_attach_TPM()
ph$compute_and_attach_FPKM()

# 差异分析
res <- ph$run_edgeR_for_groups("Treatment", "Control", 
                                out_dir = "04.Differential_差异分析结果")

# 可视化
ph$plot_pca(assay = "TPM", out_file = "PCA_plot.png")
ph$plot_volcano(res$tbl, out_file = "volcano_plot.png")
```

### 依赖包安装

该模块依赖 26 个 R 包。详细的包列表、用途说明和安装脚本请参见：
- **[transcript_pipeline_module_dependencies.md](transcript_pipeline_module_dependencies.md)** — 完整的包依赖清单与安装指南

---

## 常见操作举例

### Bash 流水线操作

- 仅下载 SRR 列表：
  - ./transcriptome_pipeline.sh download -l list_file -t 6
- 并行比对（已存在 fastq 目录）：
  - ./transcriptome_pipeline.sh align -l list_file -g /path/to/hisat_index -t 40
- 计数并使用 transcript 级别注释：
  - ./transcriptome_pipeline.sh count -l list_file -a genome.gtf -F transcript -G transcript_id

---

## 进阶与运维建议

- 在 HPC 上：用 Slurm/SGE 调度 parallel，限制并行任务数，避免 I/O 瓶颈。
- 容器化：将环境（fastp/hisat2/samtools/featureCounts）打包到一个 Docker 或 Singularity 镜像里，保证可重现性。
- 日志与监控：定期检查各样本 logs/ 下的工具日志，特别是 hisat2 和 featureCounts 的输出。
- 权限与安全：不要把敏感凭证（PAT）写入仓库。CI 中将凭证保存为 secrets。

---

## 许可证与致谢

- 许可证：MIT（可按需替换）
- 致谢：参考了常见转录组上游流程的实现与最佳实践。

---

## 联系方式 / 贡献

有问题、建议或想贡献特性？打开 Issue 或发 PR。也可以在仓库主页联系我（ShawnSHAWN-in）。