RNA 数据不像wgbs那样很大需要分开跑，rna-seq指定输入目录后使用单一脚本即可



[TOC]



## Functions

### 1. prefetch下载数据

```sh
prefetch SRR123456
```

```sh
# 下载metadata
geofetch -i GSE162903 -n metadata --just-metadata
# 下载rawdata
nohup geofetch -i GSE162903 > output.log &
```
### 2. SRA 转 fq.gz

从ncbi上下载的原始测序数据一般以SRA（Sequence Read Archive）格式存储，在分析前需要被转换为 fq.gz 格式

参数：
- `indir`: 输入目录
- `outdir`: 输出目录
- `${workdir}/log`: 默认的运行记录输出目录

```sh
# sra2fq.gz
#convert_sra_to_fastq "${work_dir}/sra" "${work_dir}/rawdata" "PE"
convert_sra_to_fastq() {
    local indir="$1"
    local outdir="$2"
    local dt="$3"

    # Convert SRA files to FASTQ files
    find "${indir}" -name "*.sra" | while read -r file; do
        # Extract the base name of the file without the directory and extension
        base_name=$(basename "$file" .sra)

        # Perform the conversion and compression
        fasterq-dump --threads 10 --split-3 --outfile "${outdir}/${base_name}.fastq" "$file"

        if [ "$dt" = "PE" ]; then
            pigz -p 20 "${outdir}/${base_name}_1.fastq"
            pigz -p 20 "${outdir}/${base_name}_2.fastq"
        elif [ "$dt" = "SE" ]; then
            pigz -p 20 "${outdir}/${base_name}.fastq"
        fi
    done
}

```

### 3. 质控

```sh
# trim
# quality_control_and_trim "PE" "/your/work/directory" 
quality_control_and_trim() {
    echo -e "\nStep 1: Trimming"
    echo -e "------------------------------------"
    
    # Parameters
    local dt="$1"
    local file="$2"
    
    # Trimming
    # Function for cut and trimming paired-end reads
    trim_pe_reads() {
        local read1="$1"
        local read2="$2"

        base=$(echo "${read1}" | awk -v re="${RE_basename}" '{gsub(re, ""); print}')

        echo "trim receive paras: read1=$read1, read2=$read2, base=${base}" >> "${clean_dir}/received_paras.log"

        trim_galore --paired "${fq_dir}/${read1}" "${fq_dir}/${read2}" \
            -o "${clean_dir}" --quality 20 --max_n 4 --length 30 --phred33 --cores 8 \
            >> "${log_dir}/trim.log" 2>&1
    }
    
    # Function for trimming single-end reads
    trim_se_read() {
        local file="$1"
        trim_galore "${fq_dir}/${file}" -o "${clean_dir}" --quality 20 --max_n 4 --length 30 --cores 8 \
            >> "${log_dir}/trim.log" 2>&1
    }
    
    # Parallel processing for trimming
    if [ "${dt}" == "PE" ]; then
        export -f trim_pe_reads
        parallel --colsep '\t' --link trim_pe_reads {1} {2} :::: <(awk 'NR%2==1' "${file}") :::: <(awk 'NR%2==0' "${file}")
    
    else
        export -f trim_se_read
        parallel trim_se_read {1} :::: "${file}"
    fi
    
    # Generate list of trimmed files
    find "${clean_dir}" -name "*.gz" -exec basename {} \; | sort > "${clean_dir}/trimgalore_fq_${dt}.txt"
    
    echo "trim finished at $(date)"
}

```

### 4. 质检
对原始测序数据和质控后的数据进行质检
```sh
# QC with multiqc
qc_with_multiqc() {
    echo -e "\n------QC for ${subdir}------"

    # Parameters
    local subdir="$1"

    # Create folder
    mkdir -p "${qc_dir}/${subdir}"

    find "${workdir}/${subdir}" -name "*.gz" | \
        xargs -P 10 -I{} fastqc {} -o "${qc_dir}/${subdir}" >> "${log_dir}/qc.log" 2>&1
    multiqc "${qc_dir}/${subdir}" -o "${qc_dir}/${subdir}" >> "${log_dir}/log/multiqc.log" 2>&1

    echo "qc_with_multiqc finished at $(date)"
}

```
### 5.利用hisat2将测序数据比对到基因组

```sh
# Align by hisat2
align_reads_with_hisat2() {
    echo -e "\nStep 2: Aligning Reads with HISAT2"
    echo -e "------------------------------------"

    # Parameters
    dt="$1"        # Data type: "PE" for paired-end, otherwise single-end
    index="$2"     # Path to HISAT2 index
    divider="$3"   # Character used for splitting, only the first part is used as basename

    file="${clean_dir}/trimgalore_fq_${dt}.txt"

    if [ "${dt}" == "PE" ]; then
        for i in $(seq 1 2 $(cat "${file}" | wc -l)); do
            read1=$(awk -v row=${i} '(NR == row){print $0}' "${file}")
            read2=$(awk -v row=$((i+1)) '(NR == row){print $0}' "${file}")
            base=$(echo "${read1}" | awk -F"${divider}" '{print $1}')

            echo "hisat2 receive parameters: read1=$read1, read2=$read2, divider=$divider, base=${base}" >> "${log_dir}/hisat2.log" "${align_dir}/received_paras.log"

            hisat2 -t -q -p 20 --dta-cufflinks \
                   -x "${index}" -1 "${clean_dir}/${read1}" -2 "${clean_dir}/${read2}" \
                   -S "${align_dir}/${base}.sam" >> "${log_dir}/hisat2.log" 2>&1

            samtools sort -@ 20 -o "${align_dir}/${base}.Hisat_aln.sorted.bam" "${align_dir}/${base}.sam"
            # 生成.bai
            samtools index "${align_dir}/${base}.Hisat_aln.sorted.bam" 
        done
    else
        while read -r file; do
            base=$(echo "${file}" | awk -F"${divider}" '{print $1}')
            
            hisat2 -t -q -p 20 -x "${index}" -U "${clean_dir}/${file}" \
                   --dta-cufflinks -S "${align_dir}/${base}.sam" >> "${log_dir}/hisat2.log" 2>&1

            samtools sort -n -@ 20 -o "${align_dir}/${base}.Hisat_aln.sorted.bam" "${align_dir}/${base}.sam"
            samtools index "${align_dir}/${base}.Hisat_aln.sorted.bam" 
        done < "${file}"
    fi

    find "${align_dir}/" -name "*.bam" -exec basename {} \; | sort > "${align_dir}/align_fq_${dt}.txt"

    echo "align finished at $(date)"
}

```

### 6. featurecount 统计count

```sh
# feature count
count_feature_with_featurecount() {
    echo -e "\nStep 3: Counting Features with FeatureCounts"
    echo -e "------------------------------------"

    # Parameters
    dt="$1"       # Data type: "PE" for paired-end, otherwise single-end
    gtf="$2"      # Path to GTF file
    RE_name="$3"  #.*RNA-|.Hisat.*  Regex pattern for sample name extraction

    if [ "${dt}" == "PE" ]; then
        featureCounts -T 20 -p --countReadPairs -t exon -g gene_name \
                      -a "${gtf}" -o "${count_dir}/counts.txt" "${align_dir}"/*.bam >> "${log_dir}/featureCounts.log" 2>&1

    else
        featureCounts -T 20 -t exon -g gene_id \
                      -a "${gtf}" -o "${count_dir}/counts.txt" "${align_dir}"/*.bam >> "${log_dir}/featureCounts.log" 2>&1
    fi

    # 整理格式
    # 删除第一列的注释 + 删除中间2-6列：Chr   Start   End Strand  Length + 提取RNA-和.Hisat中间的样本名
    tail -n +2 "${count_dir}/counts.txt" | \
    perl -lane 'splice @F,1,5; print join "\t",@F' | \
    awk -F'\t' -v pattern="${RE_name}" 'BEGIN {OFS="\t"} NR==1 {for (i=2; i<=NF; i++) {gsub(pattern, "", $i)}} {print}' > "${count_dir}/extracted_counts.txt"

    echo "featureCounts finished at $(date)"
}


```


### 7. count 标准化（FPKM/TPM）

```sh
# normalizing
my_norm() {
    local type="$1" # cufflinks,stringtie的输入都是align.bam，rnanorm的输入是counts
    local file="$2"
    local gtf="$3"
    local RE_name="$4" #.*RNA-|.Hisat.*

    local base=$(echo "${file}" | awk -v re="${RE_name}" '{gsub(re, ""); print}')
    echo "my_norm receive parameters: type=$type, file=$file, gtf=$gtf, base=${base}" >> "${norm_dir}/received_paras.log"

    if [ "${type}" == "cufflinks" ]; then
        cufflinks -p 20 -G ${gtf} -o ${norm_dir}/cufflinks/${base} ${align_dir}/${file}

    elif [ "${type}" == "stringtie" ]; then
        echo -e "${base}\t${base}.gtf" >> "${norm_dir}/stringtie/sample_list.txt"
        stringtie -e -B -p 20 -G ${gtf} -o ${norm_dir}/stringtie/${base}.gtf ${align_dir}/${file}

    elif [ "${type}" == "rnanorm" ]; then
        # 转置
        python3 -c "import pandas as pd; pd.read_csv('${file}', sep='\t').T.to_csv('${norm_dir}/rnanorm/extracted_counts_t.csv', sep=',', header=False)"
        # 计算
        rnanorm tpm "${norm_dir}/rnanorm/extracted_counts_t.csv" --gtf ${gtf} > ${norm_dir}/rnanorm/tpm.csv
        rnanorm fpkm "${norm_dir}/rnanorm/extracted_counts_t.csv" --gtf ${gtf} > ${norm_dir}/rnanorm/fpkm.csv
    fi
}

# Define a function to perform RNA-seq data normalization
para_my_norm() {
    local type="$1"  # Type of normalization method (e.g., rnanorm, cufflinks, stringtie)
    local gtf="$2"   # Path to the GTF file
    local RE_name="$3"

    # Create the directory for normalized data
    mkdir -p "${norm_dir}/${type}"

    # Print a message indicating the start of normalization
    echo -e "\nStep 4: ${type} normalizing"
    echo -e "------------------------------------"

    # Export the my_norm function for use with parallel
    export -f my_norm

    # Perform normalization based on the type
    if [ "${type}" == "rnanorm" ]; then
        # Perform rnanorm normalization
        my_norm "rnanorm" "${count_dir}/extracted_counts.txt" "${gtf}" "${RE_name}"

    elif [ "${type}" == "cufflinks" ]; then
        # Perform cufflinks normalization in parallel
        parallel --colsep '\t' --link my_norm "${type}" {1} "${gtf}" "${RE_name}" :::: "${align_dir}/align_fq.txt"
        # Extract and merge FPKM data
        python "${scripts_dir}/cufflinks/merged_FPKM.py" "${norm_dir}/cufflinks/"

    elif [ "${type}" == "stringtie" ]; then
        # Perform stringtie normalization in parallel
        parallel --colsep '\t' --link my_norm "${type}" {1} "${gtf}" "${RE_name}" :::: "${align_dir}/align_fq.txt"
        # Count matrix generation
        python3 "${scripts_dir}/stringtie/prepDE.py3" \
            -i "${norm_dir}/stringtie/sample_list.txt" \
            -g gene_count_matrix.csv
        # FPKM matrix generation
        python3 "${scripts_dir}/stringtie/getFPKM.py3" \
            -i "${norm_dir}/stringtie/sample_list.txt" \
            -g gene_fpkm_matrix.csv
        # TPM matrix generation
        python3 "${scripts_dir}/stringtie/getTPM.py3" \
            -i "${norm_dir}/stringtie/sample_list.txt" \
            -g gene_tpm_matrix.csv
    fi
}

```

### 8. 计算repeats的表达量

#### 8.1 bowtie2 法

需要提供序列信息，直接通过bowtie2 map

```sh
# repeat_with_bowtie2.sh
source ~/.bashrc
conda activate base-omics

workdir=/public/slst/home/leixy2023/project/DNMT3C/rna/240519_mainline
bowtie2_index=/public/slst/home/leixy2023/database/mm10/mm10_repeatmasker/mm10_rep_bowtie2
chrom_fai=/public/slst/home/leixy2023/database/mm10_rep/mm_rep.fa.fai
export workdir bowtie2_index chrom_sizes gff

repeat_with_bowtie2(){
  local read1="$1"
  local read2="$2"
  local divider="$3"
  base=$(echo "${read1}" | awk -F"${divider}" '{print $1}')
  echo "repeat_with_bowtie2 receive parameters: read1=$read1, read2=$read2, divider=$divider, base=${base}" >> "${workdir}/08_rep_bowtie2/bowtie2.log" 
  
  # bowtie2
  bowtie2 -p 10 -x ${bowtie2_index} -1 "${workdir}/02_cleandata/${read1}" -2 "${workdir}/02_cleandata/${read2}" -S "${workdir}/08_rep_bowtie2/${base}.sam" >${workdir}/log/bowtie2_${base}.log 2>&1
  
  # sam2bam
  samtools sort -n -@ 20 -o "${workdir}/08_rep_bowtie2/${base}_sorted.bam" "${workdir}/08_rep_bowtie2/${base}.sam" >${workdir}/log/sam2bam.log 2>&1

  # sam2bed
  mkdir -p ${workdir}/08_rep_bowtie2/{bed,count}
  convert2bed -i sam < "${workdir}/08_rep_bowtie2/${base}.sam" > "${workdir}/08_rep_bowtie2/bed/${base}.bed"


}
export -f repeat_with_bowtie2

file="${workdir}/02_cleandata/trimgalore_fq_PE.txt"
divider="_"
paste <(awk 'NR%2==1' "${file}") <(awk 'NR%2==0' "${file}") | awk -v OFS='\t' -v divider="${divider}" '{$3 = $3 divider; print}'  > "${workdir}/08_rep_bowtie2/tmp_parallel.txt"
# parallel
time parallel --colsep '\t' --link repeat_with_bowtie2 {1} {2} {3} :::: "${workdir}/08_rep_bowtie2/tmp_parallel.txt"

# extractcount
rm "${workdir}/08_rep_bowtie2/bed/rep_bed.txt"
find "${workdir}/08_rep_bowtie2/bed/" -type f -name '*.bed' | while read -r file; do
    base=$(echo "$file" | sed 's/.*RNA-\(.*\)\..*/\1/')
    echo -e "$base\t$file" >> "${workdir}/08_rep_bowtie2/bed/rep_bed.txt"
done

python extract.py "${workdir}/08_rep_bowtie2/bed/rep_bed.txt" test.txt
```

#### 8.2 featurecount 法

提供repeatst的gtf，直接从hisat2的结果中利用featurecount提取

```sh
gtf=/home_data/home/slst/leixy2023/data/database/mm10/mm10_repeatmasker/mm10_rep.gtf

# feature count
count_feature_with_featurecount() {
  echo -e "\nStep 3: Counting Features with FeatureCounts"
  echo -e "------------------------------------"
  
  # Parameters
  dt="$1"       # Data type: "PE" for paired-end, otherwise single-end
  gtf="$2"    # Path to HISAT2 index
  count_input="/home_data/home/slst/leixy2023/data/project/DNMT3C/rna/240519_mainline/04_align"
  count_output="/home_data/home/slst/leixy2023/data/project/DNMT3C/rna/240519_mainline/09_count_rep"

  # Counting Features with FeatureCounts  
  mkdir -p ${workdir}/05_counts
  featureCounts -T 20 -p --countReadPairs -t exon -g gene_id \
                  -a ${gtf} -o ${count_output}/counts.txt  ${count_input}/*.bam >> "${count_output}/featureCounts.log" 2>&1


  
  # 整理格式
  # 删除第一列的注释 + 删除中间2-6列：Chr	Start	End	Strand	Length + 提取RNA-和.Hisat中间的样本名
  tail -n +2 ${workdir}/05_counts/counts.txt | perl -lane 'splice @F,1,5; print join "\t",@F' | awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {for (i=2; i<=NF; i++) {gsub(/.*RNA-|.Hisat.*/, "", $i)}} {print}' > ${workdir}/05_counts/extracted_counts.txt
  
  echo "featurecount finished at $(date)"
}
```


## Scripts

在脚本的开头我一般会注明
1. 脚本的名字
2. Shebang信息，指定使用 Bash shell 来执行脚本，用于在本地测试脚本无调用错误
3. pbs标识。pbs 是学校集群的调度系统，通过添加pbs标识将命令上传到云端（节点）运行。 更多关于pbs见：[SHTH-HPC-PBS](../../A.格物致知/04%20软件技术/SHTH-HPC-PBS.md)

```sh
#rnaseq_submit.sh
#!/bin/bash

#PBS -N rnaseq.pbs
#PBS -l nodes=1:ppn=20
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_pub

```

脚本中可设置的参数：
- `inputdir`: 输入数据所在文件夹。可以是sra或fq.gz
- `workdir`: 输出数据所在的主文件夹
- `Spe`: 使用的物种和版本，目前提供hg38, mm10, mm9的基因组和注释文件，如果使用其它物种需自行添加。
- `dt`: 测序数据类型，目前只支持双端测序数据处理
- `divider`: 用于提取文件名

```sh
# Define pipeline configuration parameters
inputdir=/public/slst/home/leixy2023/project/DNMT3C/rna/240519_mainline/01_rawdata      # Input directory
workdir="${inputdir}/.."  # Working directory for analysis results
Spe="mm10"                # Species: "hg38" or "mm10" or "mm9"
dt="PE"                   # Data type: "PE" for paired-end, "SE" for single-end
divider="_"               # Divider for extracting filenames
```

使用的环境配置：
目录结构创建后使用 `export` 导出，以便在函数脚本 `rnaseq_functions.sh` 中直接调用



```sh
#### Configure Environment ####
source ~/.bashrc
conda activate base-omics

# Load functions from functions.sh
source /home_data/home/slst/leixy2023/pipeline/rna/rnaseq_functions.sh

# Ensure script immediately exits if any command fails, and unset variables are treated as errors
set -euo pipefail

# Determine index and GTF paths based on species
case "${Spe}" in
  "hg38")
    index=/public/slst/home/zhenghui/guoli/database/index/hisat2_index/hg38/hg38.p14
    gtf=/public/slst/home/zhenghui/guoli/database/gene_annotation/ucsc/hg38.refGene.gtf
    ;;
  "mm10")
    index=/public/slst/home/leixy2023/database/mm10/hisat2_index/mm10_hisat2_index
    gtf=/public/slst/home/leixy2023/database/mm10/mm10.refGene.gtf
    ;;
  "mm9")
    index=/public/slst/home/zhenghui/guoli/database/index/hisat2_index/mm9/mm9_hisat2_index
    gtf=/public/slst/home/zhenghui/guoli/database/gene_annotation/ucsc/mm9.refGene.gtf
    ;;
  *)
    echo "${Spe} index isn't exist"
    exit 1
    ;;
esac

# Set directory structure
clean_dir="${workdir}/02_cleandata"
qc_dir="${workdir}/03_qc"
align_dir="${workdir}/04_align"
count_dir="${workdir}/05_counts"
log_dir="${workdir}/log"

mkdir -p "${clean_dir}" "${qc_dir}" "${align_dir}" "${count_dir}" "${log_dir}"

expert "${workdir}" "${clean_dir}" "${qc_dir}" "${align_dir}" "${count_dir}" "${log_dir}"
```