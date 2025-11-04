#!/usr/bin/env bash
# 主流水线脚本：支持子命令 download|preprocess|align|count|de|trinity|all
set -euo pipefail

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 加载环境与库
source "${SCRIPTDIR}/Script_Transcript_Upstream_0_source_environment.sh"
source "${SCRIPTDIR}/lib_transcriptome.sh"

usage() {
    cat <<EOF
用法: $0 <command> [options]
命令:
  download   -l <list_file> [-t threads]    下载 SRR/EBI 数据 (prefetch 或 wget)
  preprocess -l <list_file> [-t threads]    运行 fastqc/fastp（按样本目录查找 fastq）
  align      -l <list_file> -g <hisat_index> [-t threads]  运行 hisat2 并生成 sorted.bam
  count      -l <list_file> -a <genome.gtf> [-t threads] [-F feature_type] [-G feature_attr]   运行 featureCounts
               默认 feature_type=exon, feature_attr=gene
  de         -c <count_matrix> -s <samples_file> -x <contrasts_file>  小样本差异分析
  trinity    -l <list_file>                     de novo 拼接
  all        -l <list_file> -o <organism> [-t threads]   下载/质控/比对/计数 一键运行
示例:
  $0 count -l list_file -a genome.gtf -t 12 -F exon -G gene
  $0 all -l list_file -o human -t 40
EOF
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

cmd=$1
shift || true

# 全局参数默认值
threads="${My_DEFAULT_THREADS:-40}"
list_file=""
hisat_index=""
genome_gtf=""
organism=""
feature_type="exon"
feature_attr="gene"

# 解析常用参数（子命令内会再次解析必需项）
while getopts ":l:t:g:a:o:c:s:x:i:F:G:" opt; do
    case ${opt} in
        l) list_file=${OPTARG} ;;
        t) threads=${OPTARG} ;;
        g) hisat_index=${OPTARG} ;;  # hisat index for align
        a) genome_gtf=${OPTARG} ;;   # annotation for count
        o) organism=${OPTARG} ;;
        c) count_matrix=${OPTARG} ;;
        s) samples_file=${OPTARG} ;;
        x) contrasts_file=${OPTARG} ;;
        i) hisat_index=${OPTARG} ;;
        F) feature_type=${OPTARG} ;; # featureCounts -t
        G) feature_attr=${OPTARG} ;; # featureCounts -g
        \?) echo "未知参数 -$OPTARG" >&2; usage ;;
    esac
done
shift $((OPTIND-1)) || true

# helper: map organism -> hisat_index & gtf if user used 'all'
set_genome_by_organism() {
    case "${organism}" in
        human)
            hisat_index="/data/shenlab_data/00.Database/14.Human_genome/GRCh38_index"
            genome_gtf="/data/shenlab_data/00.Database/14.Human_genome/GCF_000001405.40_GRCh38.p14_genomic.gtf"
            ;;
        mouse)
            hisat_index="/data/shenlab_data/00.Database/34.Mus_musculus_genome/Mus_musculus"
            genome_gtf="/data/shenlab_data/00.Database/34.Mus_musculus_genome/GCF_000001635.27_GRCm39_genomic.gtf"
            ;;
        rattus|rat)
            hisat_index="/data/shenlab_data/00.Database/36.Rattus_norvegicus_genome/Rattus_norvegicus_genome"
            genome_gtf="/data/shenlab_data/00.Database/36.Rattus_norvegicus_genome/GCF_036323735.1_GRCr8_genomic.gtf"
            ;;
        *)
            if [[ -n "${hisat_index}" && -n "${genome_gtf}" ]]; then
                :
            else
                err "未知 organism 或未提供索引/gtf，请使用 -i/-a 指定."
            fi
            ;;
    esac
}

case "${cmd}" in
    download)
        if [[ -z "${list_file}" ]]; then err "download 需要 -l list_file"; fi
        if grep -q "^GSM" "${list_file}"; then
            log "检测到 GSM 列表，尝试转换为 SRR"
            gsm_to_srr "${list_file}"
            list_file="${My_DATA_DIR}/list_file_SRR"
            log "使用 ${list_file} 进行 prefetch 下载"
        fi
        download_sra_with_prefetch "${list_file}" "${threads}"
        ;;
    preprocess)
        if [[ -z "${list_file}" ]]; then err "preprocess 需要 -l list_file"; fi
        check_cmd parallel
        # 使用用户指定的线程数或默认值，限制并发任务数以避免过载
        local preprocess_jobs=$((threads / 10))
        [[ ${preprocess_jobs} -lt 1 ]] && preprocess_jobs=1
        [[ ${preprocess_jobs} -gt 6 ]] && preprocess_jobs=6
        parallel -j ${preprocess_jobs} --joblog "${My_DATA_DIR}/log_preprocess" --progress "bash -lc 'source ${SCRIPTDIR}/lib_transcriptome.sh; run_hisat2_sample {} dummy_index 2>/dev/null || true'" :::: "${list_file}" || true
        ;;
    align)
        if [[ -z "${list_file}" || -z "${hisat_index}" ]]; then err "align 需要 -l list_file -g hisat_index"; fi
        check_cmd parallel
        parallel -j 6 --joblog "${My_DATA_DIR}/log_hisat2_${list_file}" --progress "bash -lc 'source ${SCRIPTDIR}/lib_transcriptome.sh; run_hisat2_sample {} ${hisat_index} ${threads}'" :::: "${list_file}"
        ;;
    count)
        if [[ -z "${list_file}" || -z "${genome_gtf}" ]]; then err "count 需要 -l list_file -a genome.gtf"; fi
        # 将用户指定的 feature_type/feature_attr 传入 run_featurecounts
        run_featurecounts "${genome_gtf}" "${list_file}" "${threads}" "${feature_type}" "${feature_attr}"
        ;;
    de)
        if [[ -z "${count_matrix}" || -z "${samples_file}" || -z "${contrasts_file}" ]]; then err "de 需要 -c count_matrix -s samples_file -x contrasts_file"; fi
        run_de_small_samples "${count_matrix}" "${samples_file}" "${contrasts_file}"
        ;;
    trinity)
        if [[ -z "${list_file}" ]]; then err "trinity 需要 -l list_file"; fi
        run_trinity_de_novo "${list_file}"
        ;;
    all)
        if [[ -z "${list_file}" || -z "${organism}" ]]; then err "all 需要 -l list_file -o organism"; fi
        set_genome_by_organism
        log "开始 all 流程: 下载 -> 质控/比对 -> 计数"
        if grep -q "^GSM" "${list_file}"; then
            gsm_to_srr "${list_file}"
            list_file="${My_DATA_DIR}/list_file_SRR"
            download_sra_with_prefetch "${list_file}" 6
        else
            download_sra_with_prefetch "${list_file}" 6
        fi
        parallel -j 6 --joblog "${My_DATA_DIR}/log_hisat2_all" --progress "bash -lc 'source ${SCRIPTDIR}/lib_transcriptome.sh; run_hisat2_sample {} ${hisat_index} ${threads}'" :::: "${list_file}"
        # all 流程中使用默认 feature_type/feature_attr，用户可用 count 命令单独运行计数并定制
        run_featurecounts "${genome_gtf}" "${list_file}" 12 "${feature_type}" "${feature_attr}"
        log "all 流程完成"
        ;;
    *)
        usage
        ;;
esac

exit 0