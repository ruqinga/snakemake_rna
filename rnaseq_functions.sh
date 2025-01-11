
# rnaseq_functions.sh
#!/bin/bash

# sra2fq.gz
#convert_sra_to_fastq "${work_dir}/sra" "${work_dir}/rawdata" "PE"
convert_sra_to_fastq_abandon() {
    local indir="$1"
    local outdir="$2"
    local dt="$3"

    # 将SRA文件转换为Fastq.gz文件
    find "${indir}" -name "*.sra" | \
    parallel -P 5 'fasterq-dump --threads 10 --split-3 --outfile "${outdir}/{/.}.fastq" {} && \
    pigz -p 20 "${outdir}/{/.}.fastq"'  2>"${workdir}/log/sra2fq_error.log"

}

convert_sra_to_fastq() {
    local indir="$1"
    local outdir="$2"

    # 将SRA文件转换为Fastq.gz文件
    find "${indir}" -name "*.sra" | while read -r file; do
        # Extract the base name of the file without the directory and extension
        base_name=$(basename "$file" .sra)
        
        # Perform the conversion and compression
        fasterq-dump --threads 10 --split-3 --outfile "${outdir}/${base_name}.fastq" "$file" && \
            if [[ -f "${outdir}/${base_name}_1.fastq" ]]; then
            pigz -p 20 "${outdir}/${base_name}_1.fastq"
            pigz -p 20 "${outdir}/${base_name}_2.fastq"
        else
            pigz -p 20 "${outdir}/${base_name}.fastq"
        fi
done
}


# QC with multiqc
qc_with_multiqc() {
  # Parameters
  local subdir="$1"
  echo -e "\n------QC for ${subdir}------"


  # Create the output directory
  mkdir -p ${workdir}/03_qc/${subdir}

  # Find and process files
  find ${workdir}/${subdir} -name "*.gz" | while read -r file; do
    fastqc "$file" -o ${workdir}/03_qc/${subdir} >> ${workdir}/log/qc.log 2>&1
  done

  # Run multiqc
  multiqc ${workdir}/03_qc/${subdir} -o ${workdir}/03_qc/${subdir} >> ${workdir}/log/multiqc.log 2>&1

  echo "qc_with_multiqc finished at $(date)"
}


qc_with_multiqc_abandon() {
  echo -e "\n------QC for ${subdir}------"

  # Parameters
  local subdir="$1"

  # create fold
  mkdir -p ${workdir}/03_qc/${subdir}

  find ${workdir}/${subdir} -name "*.gz" | xargs -P 10 -I{} fastqc {} -o ${workdir}/03_qc/${subdir} >> ${workdir}/log/qc.log 2>&1
  multiqc ${workdir}/03_qc/${subdir} -o ${workdir}/03_qc/${subdir} >> ${workdir}/log/multiqc.log 2>&1

  echo "qc_with_multiqc finished at $(date)"
}

# trim

# quality_control_and_trim "/your/work/directory" "PE"
quality_control_and_trim() {
  echo -e "\nStep 1: Trimming"
  echo -e "------------------------------------"

  # Parameters
  local dt="$1"      # Data type: "PE" for paired-end, otherwise single-end
  file="${workdir}/01_rawdata/raw_fq_${dt}.txt"

  #  trimming
  # Function for cut and trimming paired-end reads
  trim_pe_reads() {
  local read1="$1"
  local read2="$2"
  local divider="$3"
  base=$(echo "${read1}" | awk -F"${divider}" '{print $1}')
  echo "trim receive paras: read1=$read1, read2=$read2, divider=$divider",base=${base} >> "${workdir}/02_cleandata/received_paras.log"
  trim_galore --paired "${workdir}/01_rawdata/${read1}" "${workdir}/01_rawdata/${read2}" \
            -o "${workdir}/02_cleandata" --quality 20 --max_n 4 --length 30 --phred33 --cores 8 \
            >> "${workdir}/log/trim.log" 2>&1
            }

  # Function for trimming single-end reads
  trim_se_read() {
  local readfile="$1"
  trim_galore "${workdir}/01_rawdata/${readfile}" -o "${workdir}/02_cleandata" --quality 20 --max_n 4 \
            --length 30 --cores 8 > "${workdir}/log/trim.log" 2>&1
            }

  # Parallel processing for trimming
   if [ "${dt}" == "PE" ]; then
     export -f trim_pe_reads
     paste <(awk 'NR%2==1' "${file}") <(awk 'NR%2==0' "${file}") | awk -v OFS='\t' -v divider="${divider}" '{$3 = $3 divider; print}'  > "${workdir}/02_cleandata/tmp_parallel.txt"
     parallel --colsep '\t' --link trim_pe_reads {1} {2} {3} :::: "${workdir}/02_cleandata/tmp_parallel.txt"

  else
    export -f trim_se_read
    parallel trim_se_read {} :::: "${file}"
  fi

  # Generate list of trimmed files
  find "${workdir}/02_cleandata" -name "*.gz" -exec basename {} \; |sort > "${workdir}/02_cleandata/trimgalore_fq_${dt}.txt"

  echo "trim finished at $(date)"
}

# Align by hisat2
align_reads_with_hisat2() {
  echo -e "\nStep 2: Aligning Reads with HISAT2"
  echo -e "------------------------------------"

  # Parameters
  dt="$1"       # Data type: "PE" for paired-end, otherwise single-end
  index="$2"    # Path to HISAT2 index
  divider="$3"      # 以这个字符进行分割，只保留第一位作为basename

  file="${workdir}/02_cleandata/trimgalore_fq_${dt}.txt"

  if [ "${dt}" == "PE" ]; then
    for i in $(seq 1 2 $(cat "${file}" | wc -l)); do
      read1=$(awk -v row=${i} '(NR == row){print $0}' "${file}")
      read2=$(awk -v row=$((i+1)) '(NR == row){print $0}' "${file}")
      base=$(echo "${read1}" | awk -F"${divider}" '{print $1}')

      echo "hisat2 receive parameters: read1=$read1, read2=$read2, divider=$divider, base=${base}" >> "${workdir}/log/hisat2.log" "${workdir}/04_align/received_paras.log"

      hisat2 -t -q -p 20 --dta-cufflinks \
             -x "${index}" -1 "${workdir}/02_cleandata/${read1}" -2 "${workdir}/02_cleandata/${read2}" \
             -S "${workdir}/04_align/${base}.sam" >> "${workdir}/log/hisat2.log" 2>&1

       samtools sort -@ 20 -o "${workdir}/04_align/${base}.Hisat_aln.sorted.bam" "${workdir}/04_align/${base}.sam"
      samtools index "${workdir}/04_align/${base}.Hisat_aln.sorted.bam" "${workdir}/04_align/${base}.Hisat_aln.sorted.bam.bai"
    done
  else
    while IFS= read -r read;do
      base=$(echo "${read}" | awk -F"${divider}" '{print $1}')
      #base="${read%.*}"

      echo "hisat2 receive parameters: read=$read, divider=$divider, base=${base}" >> "${workdir}/log/hisat2.log" "${workdir}/04_align/received_paras.log"

      hisat2 -t -q -p 20 -x "${index}" -U "${workdir}/02_cleandata/${read}" \
             --dta-cufflinks -S "${workdir}/04_align/${base}.sam" >> "${workdir}/log/hisat2.log" 2>&1

      samtools sort -@ 20 -o "${workdir}/04_align/${base}.Hisat_aln.sorted.bam" "${workdir}/04_align/${base}.sam"
      samtools index "${workdir}/04_align/${base}.Hisat_aln.sorted.bam" "${workdir}/04_align/${base}.Hisat_aln.sorted.bam.bai"
    done < "${file}"
  fi

  find "${workdir}/04_align/" -name "*.bam" -exec basename {} \; | sort > "${workdir}/04_align/align_fq.txt"

  echo "align finished at $(date)"
}

# feature count
count_feature_with_featurecount() {
  echo -e "\nStep 3: Counting Features with FeatureCounts"
  echo -e "------------------------------------"

  # Parameters
  dt="$1"       # Data type: "PE" for paired-end, otherwise single-end
  gtf="$2"    # Path to HISAT2 index

  # Counting Features with FeatureCounts
  mkdir -p ${workdir}/05_counts

  if [ "${dt}" == "PE" ]; then
    featureCounts -T 20 -p --countReadPairs -t exon -g gene_name \
                  -a ${gtf} -o ${workdir}/05_counts/counts.txt  ${workdir}/04_align/*.bam >> "${workdir}/log/featureCounts.log" 2>&1

  else
    featureCounts -T 20 -t exon -g gene_name \
                  -a ${gtf} -o ${workdir}/05_counts/counts.txt  ${workdir}/04_align/*.bam
  fi

  # 整理格式
  # 删除第一列的注释 + 删除中间2-6列：Chr   Start   End Strand  Length + 提取RNA-和.Hisat中间的样本名
  tail -n +2 ${workdir}/05_counts/counts.txt | perl -lane 'splice @F,1,5; print join "\t",@F' | awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {for (i=2; i<=NF; i++) {gsub(/.*RNA-|.Hisat.*/, "", $i)}} {print}' > ${workdir}/05_counts/extracted_counts.txt

  echo "featurecount finished at $(date)"
}

# normalizing
my_norm() {
  local type="$1" # cufflinks,stringtie的输入都是align.bam，rnanorm的输入是counts
  local file="$2"
  local gtf="$3"
  local base=$(basename "$file" .Hisat_aln.sorted.bam)
  #local base=$(gsub(/.*RNA-|.Hisat.*/, "","$file") '{print }')
  echo "my_norm receive parameters: type=$type, file=$file, gtf=$gtf, base=${base}" >> "${workdir}/06_norm/received_paras.log"

  if [ "${type}" == "cufflinks" ]; then
    cufflinks -p 20 -G ${gtf} -o ${workdir}/06_norm/cufflinks/${base} ${workdir}/04_align/${file}

  elif [ "${type}" == "stringtie" ]; then
    echo -e "${base}\t${base}.gtf" >> "${workdir}/06_norm/stringtie/sample_list.txt"
    stringtie -e -B -p 20 -G ${gtf} -o ${workdir}/06_norm/stringtie/${base}.gtf ${workdir}/04_align/${file}

  elif [ "${type}" == "rnanorm" ]; then
    # 转置
    python3 -c "import pandas as pd; pd.read_csv('${file}', sep='\t').T.to_csv('${workdir}/06_norm/rnanorm/extracted_counts_t.csv', sep=',', header=False)"
    # 计算
    rnanorm tpm "${workdir}/06_norm/rnanorm/extracted_counts_t.csv" --gtf ${gtf} > ${workdir}/06_norm/rnanorm/tpm.csv
    rnanorm fpkm "${workdir}/06_norm/rnanorm/extracted_counts_t.csv" --gtf ${gtf} > ${workdir}/06_norm/rnanorm/fpkm.csv
  fi
}

# Define a function to perform RNA-seq data normalization
para_my_norm() {
  local type="$1"  # Type of normalization method (e.g., rnanorm, cufflinks, stringtie)
  local gtf="$2"   # Path to the GTF file
  #local pipeline_dir="/home_data/home/slst/leixy2023/pipeline/rna"  # Pipeline directory

  # Create the directory for normalized data
  mkdir -p "${workdir}/06_norm/${type}"

  # Print a message indicating the start of normalization
  echo -e "\nStep 4: ${type} normalizing"
  echo -e "------------------------------------"

  # Export the my_norm function for use with parallel
  export -f my_norm

  # Perform normalization based on the type
  if [ "${type}" == "rnanorm" ]; then
    # Perform rnanorm normalization
    my_norm "rnanorm" "${workdir}/05_counts/extracted_counts.txt" "${gtf}"

  elif [ "${type}" == "cufflinks" ]; then
    # Perform cufflinks normalization in parallel
    parallel --colsep '\t' --link my_norm "${type}" {1} "${gtf}" :::: "${workdir}/04_align/align_fq.txt"
    # Extract and merge FPKM data
    python "${pipeline_dir}/merge_FPKM.py" "${workdir}/06_norm/cufflinks/"

  elif [ "${type}" == "stringtie" ]; then
    # Perform stringtie normalization in parallel
    parallel --colsep '\t' --link my_norm "${type}" {1} "${gtf}" :::: "${workdir}/04_align/align_fq.txt"
    # Count matrix generation
    python3 "${pipeline_dir}/stringtie/prepDE.py3" \
      -i "${workdir}/06_norm/stringtie/sample_list.txt" \
      -g gene_count_matrix.csv
    # FPKM matrix generation
    python3 "${pipeline_dir}/stringtie/getFPKM.py3" \
      -i "${workdir}/06_norm/stringtie/sample_list.txt" \
      -g gene_fpkm_matrix.csv
    # TPM matrix generation
    python3 "${pipeline_dir}/stringtie/getTPM.py3" \
      -i "${workdir}/06_norm/stringtie/sample_list.txt" \
      -g gene_tpm_matrix.csv
  fi
}




