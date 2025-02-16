#!/bin/bash
# 用于处理单端和双端数据并调用 Snakemake 运行流程
# 用法: bash run.sh <work_dir> <fq_dir> -y

# 默认工作目录为上一级目录
work_dir="${1:-../}"
echo "工作目录设置为: $work_dir"

# 默认fq目录为工作目录下的01_rawdata目录
fq_dir="${2:-${work_dir}/01_rawdata}"
echo "原始数据目录: $fq_dir"

# 激活 Snakemake 的 Conda 环境
source ~/.bashrc
conda activate snakemake_env

# 调用脚本生成单端和双端数据列表

# 初始化数组
json_array_pe=()
json_array_se=()

# 遍历所有 .fastq.gz 文件
for file in $(find "$fq_dir" -name "*.fastq.gz" -maxdepth 1 | sort); do
    # 检查是否是带 _1.fastq.gz 或 _2.fastq.gz 的 PE 文件
    if [[ "$file" =~ _1.fastq.gz$ ]]; then
        # 获取对应的 _2 文件
        file2="${file/_1.fastq.gz/_2.fastq.gz}"

        if [ -f "$file2" ]; then
            # 如果 _2 文件存在，则添加到 PE 数组
            json_entry_pe="{\"read1\": \"$file\", \"read2\": \"$file2\"}"
            json_array_pe+=("$json_entry_pe")
        fi
    elif [[ ! "$file" =~ _1.fastq.gz$ && ! "$file" =~ _2.fastq.gz$ ]]; then
        # 如果不是 _1 或 _2，视为 SE 文件
        json_entry_se="{\"read1\": \"$file\"}"
        json_array_se+=("$json_entry_se")
    fi
done

# 将数组转换为 JSON 字符串
json_output_pe=$(IFS=,; echo "[${json_array_pe[*]}]")
json_output_se=$(IFS=,; echo "[${json_array_se[*]}]")

# 检查 JSON 数据是否为空并输出相应的提示信息
if [[ -z "$json_output_se" || "$json_output_se" == "[]" ]]; then
    echo "没有单端数据"
    json_output_se=""
else
    echo "即将处理如下单端数据:"
    echo "$json_output_se" | jq .
fi

if [[ -z "$json_output_pe" || "$json_output_pe" == "[]" ]]; then
    echo "没有双端数据"
    json_output_pe=""
else
    echo "即将处理如下双端数据:"
    echo "$json_output_pe" | jq .
fi



# 运行 Snakemake 工作流（预览模式）
if [[ -n "$json_output_pe" ]]; then
    echo "运行 Snakemake 处理双端数据（仅预览）..."
    snakemake \
        -np \
        --executor cluster-generic \
        --cluster-generic-submit-cmd 'qsub -q slst_pub -N rna_pe.pbs -l nodes=1:ppn=20' \
        --latency-wait 60 \
        --jobs 4 \
        --use-conda \
        --group-components processing=8 \
        --config fq_dir="$fq_dir" work_dir="$work_dir" dt="PE" reads="$json_output_pe"
fi

if [[ -n "$json_output_se" ]]; then
    echo "运行 Snakemake 处理单端数据（仅预览）..."
    snakemake \
        -np \
        --executor cluster-generic \
        --cluster-generic-submit-cmd 'qsub -q slst_pub -N rna_se.pbs -l nodes=1:ppn=20' \
        --latency-wait 60 \
        --jobs 4 \
        --use-conda \
        --group-components processing=8 \
        --config fq_dir="$fq_dir" work_dir="$work_dir" dt="SE" reads="$json_output_se"
fi

# 提示是否确认实际执行任务
# 检查是否传递了 -y 参数
if [[ "$3" == "-y" ]]; then
    confirm_run="y"
else
    read -p "是否确认执行任务（实际提交作业）？(y/n): " confirm_run
fi

if [[ "$confirm_run" != "y" && "$confirm_run" != "Y" ]]; then
    echo "任务已取消！"
    exit 0
fi

# 实际运行 Snakemake 工作流
if [[ -n "$json_output_pe" ]]; then
    echo "开始处理双端数据..."
    snakemake \
        --executor cluster-generic \
        --cluster-generic-submit-cmd 'qsub -q slst_fat -N rna_pe.pbs -l nodes=1:ppn=10 -j oe' \
        --latency-wait 60 \
        --jobs 8 \
        --use-conda \
        --group-components processing=8 \
        --config fq_dir="$fq_dir" work_dir="$work_dir" dt="PE" reads="$json_output_pe"
fi

if [[ -n "$json_output_se" ]]; then
    echo "开始处理单端数据..."
    snakemake \
        --executor cluster-generic \
        --cluster-generic-submit-cmd 'qsub -q slst_fat -N rna_se.pbs -l nodes=1:ppn=10 -j oe' \
        --latency-wait 60 \
        --jobs 8 \
        --use-conda \
        --group-components processing=8 \
        --config fq_dir="$fq_dir" work_dir="$work_dir" dt="SE" reads="$json_output_se"
fi

echo "任务已完成！"
