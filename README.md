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
- README.md — 你正在阅读的博客风说明

---

## 常见操作举例

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

## 性能优化说明

本流水线已经过性能优化，针对大规模样本处理进行了以下改进：

- **批量文件处理**：Trinity 组装时使用哈希映射批量查找文件，避免重复调用 find 命令
- **高效数据格式化**：featureCounts 输出使用单个 awk 命令替代多个 sed/cut/perl 管道
- **动态资源分配**：预处理步骤的并行任务数根据可用线程自动调整
- **原生 R 处理**：GTF 解析使用原生 R 函数替代 perl 系统调用，提高可移植性
- **优化计算逻辑**：简化了差异表达分析中的均值计算和条件判断

这些优化在大样本集（50+ 样本）上可显著提升运行速度，同时保持结果的准确性。

---

## 许可证与致谢

- 许可证：MIT（可按需替换）
- 致谢：参考了常见转录组上游流程的实现与最佳实践。

---

## 联系方式 / 贡献

有问题、建议或想贡献特性？打开 Issue 或发 PR。也可以在仓库主页联系我（ShawnSHAWN-in）。