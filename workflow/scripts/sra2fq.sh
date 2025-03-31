#PBS -N sra2fq.pbs
#PBS -l nodes=1:ppn=20
#PBS -S /bin/bash
#PBS -j oe
#PBS -q slst_pub

source ~/.bashrc
conda activate base-omics

input_dir="/home_data/home/slst/leixy2023/data/project/ref_development_data/process/E135_185_bs_241105/sra"
out_dir="/home_data/home/slst/leixy2023/data/project/ref_development_data/process/E135_185_bs_241105/01_rawdata"

mkdir -p $out_dir

find "${input_dir}" -wholename "*/*.sra" | while read file; do
    fasterq-dump --threads 10 --split-3 --outdir "${out_dir}" "$file" && \
    pigz -p 10 "${out_dir}/$(basename "$file" .sra)"*.fastq
done

