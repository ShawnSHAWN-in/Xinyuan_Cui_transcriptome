#!/usr/bin/env bash
# 常用函数库：被主脚本 source
set -euo pipefail

log() {
    local ts
    ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    echo "[$ts] $*"
}

err() {
    echo "ERROR: $*" >&2
    exit 1
}

check_cmd() {
    for cmd in "$@"; do
        if ! command -v "${cmd}" >/dev/null 2>&1; then
            err "缺少命令: ${cmd}. 请安装或加载相应环境."
        fi
    done
}

# 下载 EBI bam （按你给的模式）
download_ebi_bam_from_list() {
    local list_file=$1
    local outdir="${My_DATA_DIR}/ebi_bam"
    mkdir -p "${outdir}"
    while IFS= read -r file; do
        local ftpdir="ftp://ftp.sra.ebi.ac.uk/vol1/err/${file:0:6}/00${file:9:9}/${file}.bam"
        log "wget -P ${outdir} ${ftpdir}"
        wget -P "${outdir}" "${ftpdir}"
    done <"${list_file}"
}

# 使用 prefetch 下载 SRR 列表（支持并行）
download_sra_with_prefetch() {
    local list_file=$1
    local threads=${2:-4}
    check_cmd prefetch parallel
    parallel -j "${threads}" --joblog "${My_DATA_DIR}/log_prefetch" --progress "prefetch {} -O ${My_DATA_DIR}" :::: "${list_file}"
}

# 从 GSM 列表查 SRR （依赖 elink/efetch）
gsm_to_srr() {
    local gsm_list=$1
    check_cmd elink efetch grep cut tr
    elink -db gds -input "${gsm_list}" -target sra | efetch > "${My_DATA_DIR}/list_file_SRR_GSM_meta"
    grep -Ff "${gsm_list}" "${My_DATA_DIR}/list_file_SRR_GSM_meta" > "${My_DATA_DIR}/list_file_SRR_GSM_meta_selected" || true
    cut -f 1 -d "," "${My_DATA_DIR}/list_file_SRR_GSM_meta_selected" > "${My_DATA_DIR}/list_file_SRR"
    log "生成 ${My_DATA_DIR}/list_file_SRR"
}

# 单样本的 hisat2 流程（质控/比对/索引/flagstat）
run_hisat2_sample() {
    local sample=$1
    local genome_index=$2
    local threads=${3:-"${My_DEFAULT_THREADS}"}
    local dir="${My_DATA_DIR}/${sample}"
    local logdir="${dir}/logs"
    mkdir -p "${logdir}"
    local report="${dir}/${sample}_report.txt"
    echo "样本: ${sample}" >"${report}"
    echo "基因组: ${genome_index}" >>"${report}"
    echo "线程: ${threads}" >>"${report}"

    local r1="${dir}/${sample}_1.fastq.gz"
    local r2="${dir}/${sample}_2.fastq.gz"
    local data_type=""

    if [[ -f "${r1}" && -f "${r2}" ]]; then
        data_type="paired"
    elif [[ -f "${r1}" ]]; then
        data_type="single"
    else
        err "样本 ${sample} 未找到 fastq 文件 (${r1} 或 ${r2})"
    fi

    # tools checks
    check_cmd fastqc fastp hisat2 samtools

    # fastqc 原始
    if [[ "${data_type}" == "paired" ]]; then
        fastqc -t 2 "${r1}" "${r2}" -o "${dir}" >"${logdir}/fastqc.log" 2>&1 || err "fastqc 失败 ${sample}"
    else
        fastqc -t 2 "${r1}" -o "${dir}" >"${logdir}/fastqc.log" 2>&1 || err "fastqc 失败 ${sample}"
    fi
    echo "fastqc 原始完成" >>"${report}"

    # fastp
    if [[ "${data_type}" == "paired" ]]; then
        fastp --thread 5 -i "${r1}" -I "${r2}" \
              -o "${dir}/${sample}_1.clean.fq.gz" -O "${dir}/${sample}_2.clean.fq.gz" \
              -j "${dir}/${sample}.json" -h "${dir}/${sample}.html" >"${logdir}/fastp.log" 2>&1 || err "fastp 失败"
    else
        fastp --thread 5 -i "${r1}" -o "${dir}/${sample}_1.clean.fq.gz" \
              -j "${dir}/${sample}.json" -h "${dir}/${sample}.html" >"${logdir}/fastp.log" 2>&1 || err "fastp 失败"
    fi
    echo "fastp 完成" >>"${report}"

    # fastqc 清洗后
    if [[ "${data_type}" == "paired" ]]; then
        fastqc -t 2 "${dir}/${sample}_1.clean.fq.gz" "${dir}/${sample}_2.clean.fq.gz" -o "${dir}" >"${logdir}/fastqc_clean.log" 2>&1 || err "fastqc(clean) 失败"
    else
        fastqc -t 2 "${dir}/${sample}_1.clean.fq.gz" -o "${dir}" >"${logdir}/fastqc_clean.log" 2>&1 || err "fastqc(clean) 失败"
    fi
    echo "fastqc(clean) 完成" >>"${report}"

    # hisat2 -> samtools sort
    if [[ "${data_type}" == "paired" ]]; then
        hisat2 -p "${threads}" -x "${genome_index}" \
            -1 "${dir}/${sample}_1.clean.fq.gz" -2 "${dir}/${sample}_2.clean.fq.gz" 2>"${logdir}/hisat2.log" | \
            samtools sort -@ 10 -o "${dir}/${sample}.sorted.bam" - || err "hisat2/samtools 失败"
    else
        hisat2 -p "${threads}" -x "${genome_index}" \
            -U "${dir}/${sample}_1.clean.fq.gz" 2>"${logdir}/hisat2.log" | \
            samtools sort -@ 10 -o "${dir}/${sample}.sorted.bam" - || err "hisat2/samtools 失败"
    fi
    echo "hisat2 比对完成" >>"${report}"

    # samtools index & flagstat
    samtools index -@ 5 "${dir}/${sample}.sorted.bam" >"${logdir}/samtools_index.log" 2>&1 || err "samtools index 失败"
    samtools flagstat -@ 5 "${dir}/${sample}.sorted.bam" >"${dir}/${sample}.sorted.bam.flagstat" 2>"${logdir}/samtools_flagstat.log" || err "samtools flagstat 失败"
    echo "samtools 统计完成" >>"${report}"
    echo "输出: ${dir}/${sample}.sorted.bam" >>"${report}"
}

# featureCounts 批量计数，默认 -t exon -g gene，可通过参数覆盖
# run_featurecounts <gtf> <list_file> [threads] [feature_type] [feature_attr]
run_featurecounts() {
    local gtf=$1
    local list_file=$2
    local threads=${3:-12}
    local feature_type=${4:-exon}
    local feature_attr=${5:-gene}

    check_cmd featureCounts
    mkdir -p "${My_QUANT_DIR}"
    # 组装 bam 列表
    local bams=()
    while IFS= read -r s; do
        local bam="${My_DATA_DIR}/${s}/${s}.sorted.bam"
        if [[ -f "${bam}" ]]; then
            bams+=("${bam}")
        else
            log "警告: 未找到 ${bam}"
        fi
    done <"${list_file}"

    if [[ ${#bams[@]} -eq 0 ]]; then
        err "没有找到任何 BAM 文件可以计数"
    fi

    log "运行 featureCounts: -T ${threads} -p -t ${feature_type} -g ${feature_attr} -a ${gtf}"
    featureCounts -T "${threads}" -p -t "${feature_type}" -g "${feature_attr}" \
        -a "${gtf}" -o "${My_QUANT_DIR}/result_count.tsv" "${bams[@]}" || err "featureCounts 失败"

    # 格式化输出（保留原脚本的处理）
    sed -n '2p' "${My_QUANT_DIR}/result_count.tsv" | cut -f 1,7- | perl -pe 's|\t\S+/(\S+)/[^/]+\.bam|\t\1|g;s|-|_|g' >"${My_QUANT_DIR}/result_count.without_Length.tsv"
    cut -f 1,7- "${My_QUANT_DIR}/result_count.tsv" | sed '1,2d' | sed 's/ //g' >>"${My_QUANT_DIR}/result_count.without_Length.tsv"

    sed -n '2p' "${My_QUANT_DIR}/result_count.tsv" | cut -f 1,6,7- | perl -pe 's|\t\S+/(\S+)/[^/]+\.bam|\t\1|g;s|-|_|g' >"${My_QUANT_DIR}/result_count.with_Length.tsv"
    cut -f 1,6,7- "${My_QUANT_DIR}/result_count.tsv" | sed '1,2d' | sed 's/ //g' >>"${My_QUANT_DIR}/result_count.with_Length.tsv"
    log "featureCounts 完成，输出目录: ${My_QUANT_DIR}"
}

# 小样本差异分析（Trinity 的 run_DE_analysis.pl）
run_de_small_samples() {
    local count_matrix=$1
    local samples_file=$2
    local contrasts_file=$3
    check_cmd run_DE_analysis.pl analyze_diff_expr.pl
    run_DE_analysis.pl --matrix "${count_matrix}" --method edgeR --samples_file "${samples_file}" --contrasts "${contrasts_file}" --output "${My_QUANT_DIR}/edge.result.dir" --dispersion 0.1
    analyze_diff_expr.pl --matrix "${count_matrix}" -P 1e-3 -C 2
    log "小样本差异分析完成"
}

# Trinity de novo 组装（合并多个样本的 left/right）
run_trinity_de_novo() {
    local list_file=$1
    local lefts=()
    local rights=()
    while IFS= read -r sample; do
        lefts+=("$(find "${My_DATA_DIR}" -name "${sample}_1.fastq.gz" -print -quit || true)")
        rights+=("$(find "${My_DATA_DIR}" -name "${sample}_2.fastq.gz" -print -quit || true)")
    done <"${list_file}"

    local LEFTS
    local RIGHTS
    LEFTS=$(IFS=,; echo "${lefts[*]}")
    RIGHTS=$(IFS=,; echo "${rights[*]}")

    check_cmd Trinity
    Trinity --seqType fq --max_memory 500G --left "${LEFTS}" --right "${RIGHTS}" --CPU 80 --output "${My_DATA_DIR}/trinity_out_dir"
    log "Trinity 拼接完成"
}