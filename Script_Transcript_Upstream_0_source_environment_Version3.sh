#!/usr/bin/env bash
# 环境变量与默认目录（可被外部脚本 source）
set -euo pipefail

[ ! -d "00.database" ] && mkdir -p "00.database"
[ ! -d "01.Data_质量控制结果" ] && mkdir -p "01.Data_质量控制结果"
[ ! -d "02.Quant_定量结果" ] && mkdir -p "02.Quant_定量结果"
[ ! -d "03.Samples_统计结果" ] && mkdir -p "03.Samples_统计结果"
[ ! -d "04.Differential_差异分析结果" ] && mkdir -p "04.Differential_差异分析结果"
[ ! -d "05.Enrichment_富集分析结果" ] && mkdir -p "05.Enrichment_富集分析结果"

export My_DATABASE_DIR="$(pwd)/00.database"
export My_DATA_DIR="$(pwd)/01.Data_质量控制结果"
export My_QUANT_DIR="$(pwd)/02.Quant_定量结果"
export My_SAMPLES_STAT_DIR="$(pwd)/03.Samples_统计结果"
export My_DIFFERENTIAL_DIR="$(pwd)/04.Differential_差异分析结果"
export My_ENRICHMENT_DIR="$(pwd)/05.Enrichment_富集分析结果"

# 可选：默认 thread 设置（可在命令行覆盖）
export My_DEFAULT_THREADS=40

# 可选：如果需要自动激活 conda base，请在此添加（按环境调整）
# source ~/miniconda3/etc/profile.d/conda.sh

# 提示（可删除）
echo "Loaded transcriptome environment:"
echo "  My_DATA_DIR=${My_DATA_DIR}"
echo "  My_QUANT_DIR=${My_QUANT_DIR}"