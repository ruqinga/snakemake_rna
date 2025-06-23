#!/bin/bash
# 用于处理单端和双端数据并调用 Snakemake 运行流程
# 用法: bash run.sh <fq_dir> -y

# 检查输入目录是否存在
if [ -z "$1" ]; then
    echo "错误: 未指定输入目录！"
    exit 1
elif [ ! -d "$1" ]; then
    echo "错误：目录 $1 不存在！"
    exit 1
fi

# 输入原始数据所在文件夹
fq_dir="${1}"
echo "原始数据目录: $fq_dir"

# 激活 Snakemake 的 Conda 环境
source ~/.bashrc
conda activate snakemake_env

# 调用脚本生成单端和双端数据列表

# 初始化数组
json_array_pe=()
json_array_se=()
json_array=()

# 遍历所有 .fastq.gz 文件
for file in $(find "$fq_dir" -maxdepth 1 -name "*.fastq.gz" | sort); do

    file=$(realpath "$file")  # 处理路径中的特殊字符，比如空格

    # 检查是否是带 _1.fastq.gz 或 _2.fastq.gz 的 PE 文件
    if [[ "$file" =~ _1.fastq.gz$ ]]; then
        # 获取对应的 _2 文件
        file2="${file/_1.fastq.gz/_2.fastq.gz}"

        if [ -f "$file2" ]; then
            # 如果 _2 文件存在，则添加到 PE 数组
            json_entry_pe="{\"read1\": \"$file\", \"read2\": \"$file2\"}"
            json_array_pe+=("$json_entry_pe")
            json_array+=("$json_entry_pe")
        fi
    elif [[ ! "$file" =~ _1.fastq.gz$ && ! "$file" =~ _2.fastq.gz$ ]]; then
        # 如果不是 _1 或 _2，视为 SE 文件
        json_entry_se="{\"read1\": \"$file\"}"
        json_array_se+=("$json_entry_se")
        json_array+=("$json_entry_se")
    fi
done

# 将数组转换为 JSON 字符串
json_output_pe=$(IFS=,; echo "[${json_array_pe[*]}]")
json_output_se=$(IFS=,; echo "[${json_array_se[*]}]")
json_output=$(IFS=,; echo "[${json_array[*]}]")

#echo "$json_output"

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
# snakemake --until extract_methylation_pe --touch # 可以更新文件状态，避免重复生成
echo "运行 Snakemake （仅预览）..."
snakemake \
    -np \
    --use-conda \
    --config fq_dir="$fq_dir" reads="$json_output"

# 提示是否确认实际执行任务
# 检查是否传递了 -y 参数
if [[ "$2" == "-y" ]]; then
    confirm_run="y"
else
    read -p "是否确认执行任务（实际提交作业）？(y/n): " confirm_run
fi

if [[ "$confirm_run" != "y" && "$confirm_run" != "Y" ]]; then
    echo "任务已取消！"
    exit 0
fi

# 实际运行 Snakemake 工作流
echo "运行 Snakemake ..."
snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --sample {wildcards} --rule {rule}" \
    --latency-wait 60 \
    --jobs 5 \
    --use-conda \
    --config fq_dir="$fq_dir" reads="$json_output"
#--sample {wildcards.sample}

echo "任务已完成！"


    #--cluster-generic-submit-cmd 'qsub -q slst_pub -N rna_pe.pbs -l nodes=1:ppn=20 -l walltime=100:00:00' \
    #--cluster-generic-submit-cmd "python workflow/scripts/submit_job.py --config config/cluster_config.yaml --sample {wildcards.sample} --rule "test" " \