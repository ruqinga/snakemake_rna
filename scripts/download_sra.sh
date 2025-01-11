#!/bin/bash
# 用法: bash download_sra.sh <srr_list_file>
# 示例内容:
# sample_title   run_accession   sra_md5 PRJ_ID  strategy
# Wildtype_WGBS_E18-5   SRR13202262  5511afdd30dc80ed76844b2eca7601a6 PRJNA682750 WGBS

# 检查是否提供了srr_list_file
if [ -z "$1" ]; then
    echo "错误: 需要提供SRR列表文件 <srr_list_file>"
    echo "用法: bash download_sra.sh <srr_list_file> 数据会自动下载到<srr_list_file>所在文件夹"
    exit 1
fi

# 有报错立马退出
set -euo pipefail

srr_list_file="$1"
output_dir="$(dirname "$srr_list_file")"

# 确保切换到文件所在目录
#cd "$output_dir" || { echo "无法进入目录: $output_dir"; exit 1; }

# 提取 run_accession 列并保存为 sralist.txt
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="run_accession") col=i} NR>1 {print $col}' "$srr_list_file" > "$output_dir/sralist.txt"

# 下载 SRA 文件
echo "准备下载以下数据:"
cat "$output_dir"/sralist.txt
prefetch --max-size 100G -O "$output_dir" $(<"$output_dir/sralist.txt")

# 校验 MD5
echo "开始 MD5 校验..."
awk -F'\t' 'NR>1 {print $2, $4}' "$srr_list_file" | while read -r run_accession sra_md5; do
    sra_file="${output_dir}/${run_accession}/${run_accession}.sra"
    if [ -f "$sra_file" ]; then
        calculated_md5=$(md5sum "$sra_file" | awk '{print $1}')
        if [ "$calculated_md5" == "$sra_md5" ]; then
            echo "$run_accession: MD5 校验通过"
        else
            echo "$run_accession: MD5 校验失败 (预期: $sra_md5, 计算: $calculated_md5)"
        fi
    else
        echo "$run_accession: 文件不存在"
    fi
done

# 文件重命名
echo "开始重命名文件..."
awk -F'\t' 'NR>1 {print $2,$5}' "$srr_list_file" | while read -r run_accession sample_title; do
    sra_file="${output_dir}/${run_accession}/${run_accession}.sra"
    if [ -f "$sra_file" ]; then
        mv "$sra_file" "${output_dir}/${run_accession}/${sample_title}.sra"
        echo "$sra_file 已重命名为 ${sample_title}.sra"
    else
        echo "$run_accession: 文件不存在，跳过重命名"
    fi
done


echo "流程完成！"


