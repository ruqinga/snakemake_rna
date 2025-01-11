#rnaseq_submit.sh
#!/bin/bash

#PBS -N sra2fq.pbs
#PBS -l nodes=1:ppn=6
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_pub

source ~/.bashrc
conda activate base-omics


convert_sra_to_fastq() {
    local indir="$1"
    local outdir="$2"

    # 将SRA文件转换为Fastq.gz文件
    find "${indir}" -name "*.sra" | while read -r file; do
        # Extract the base name of the file without the directory and extension
        base_name=$(basename "$file" .sra)
        
        # Perform the conversion and compression
        fasterq-dump --threads 10 --split-3 --outfile "${outdir}/${base_name}.fastq" "$file" && \
        pigz -p 20 "${outdir}/${base_name}_1.fastq"
        pigz -p 20 "${outdir}/${base_name}_2.fastq"
    done

}

convert_sra_to_fastq ~/data/project/ref_development_data/GSE162724_RAW/GSE162584/sra ~/data/project/ref_development_data/GSE162724_RAW/GSE162584/01_rawdata
