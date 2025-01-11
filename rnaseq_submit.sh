
# rnaseq_submit.sh
#!/bin/bash

#PBS -N rnaseq.pbs
#PBS -l nodes=1:ppn=20
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_pub

source ~/.bashrc
conda activate base-omics

# Load functions from functions.sh
pipeline_dir="/home_data/home/slst/leixy2023/data/project/ref_development_data/GSE162724_RAW/GSE162584/scripts"
source ${pipeline_dir}/rnaseq_functions.sh
export ${pipeline_dir}

# Ensure script immediately exits if any command fails, and unset variables are treated as errors
set -euo pipefail

# Define pipeline configuration parameters
inputdir=/home_data/home/slst/leixy2023/data/project/ref_development_data/GSE162724_RAW/GSE162584/sra      # Input directory
workdir="${inputdir}/.." # Working directory for analysis results
Spe="mouse"             # Species: "human" or "mouse"
dt="PE"              # Data type: "PE" for paired-end, "SE" for single-end
divider="_"                  # 用于提取文件名

# set global variables
export workdir

# Create output directories if they don't exist
mkdir -p "$workdir"/{01_rawdata,02_cleandata,03_qc,04_align,05_counts,log}

# Determine index and GTF paths based on species
if [ "${Spe}" == "human" ]; then
  index=/public/slst/home/zhenghui/guoli/database/index/hisat2_index/hg38/hg38.p14
  gtf=/public/slst/home/zhenghui/guoli/database/gene_annotation/ucsc/hg38.refGene.gtf
else
  index=/public/slst/home/leixy2023/database/mm10/hisat2_index/mm10_hisat2_index
  gtf=/public/slst/home/leixy2023/database/mm10/gencode.vM25.annotation_remove_miRNA_snoRNA_snRNA.gtf
  #index=/public/slst/home/zhenghui/guoli/database/index/hisat2_index/mm9/mm9_hisat2_index
  #gtf=/public/slst/home/zhenghui/guoli/database/gene_annotation/ucsc/mm9.refGene.gtf
fi

# start pipline
echo "RNA-seq Analysis Pipeline start at $(date)"

# sra2fq.gz
# 判断输入目录是否存在，并且是否包含 SRA 文件
if [ -d "$inputdir" ] && find "$inputdir" -type f -name "*.sra" | read; then
    # 如果满足条件，则运行 convert_sra_to_fastq
    #time convert_sra_to_fastq "${inputdir}" "${workdir}/01_rawdata"
    echo -e "转换 SRA 文件为fastq.gz。\n"
else
    # 如果条件不满足，则输出错误信息
    echo -e "输入目录不存在或者不包含 SRA 文件。\n"
fi

# Generate list of raw Fastq files
find "${workdir}/01_rawdata" -name "*.gz" -exec basename {} \; |sort > "${workdir}/01_rawdata/raw_fq_${dt}.txt"

# trim
#time quality_control_and_trim "${dt}" "${divider}"

# QC
#time qc_with_multiqc "01_rawdata"
#time qc_with_multiqc "02_cleandata"

# hisat2
#time align_reads_with_hisat2 "${dt}" "${index}" "${divider}"

# feature count
time count_feature_with_featurecount "${dt}" "${gtf}"

# norm
#time para_my_norm "rnanorm" "${gtf}"
#time para_my_norm "stringtie" "${gtf}"
conda activate py36
time para_my_norm "cufflinks" "${gtf}"

# The pipeline ends here
echo -e "\nRNA-seq Analysis Pipeline Completed at $(date)"
