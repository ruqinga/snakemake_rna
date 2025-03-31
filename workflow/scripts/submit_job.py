#!/usr/bin/env python3
"""
Snakemake-PBS 提交脚本
"""

import os
import re
import argparse
import yaml
from datetime import datetime
import subprocess


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="cluster_config.yaml 路径")
    parser.add_argument("--sample", nargs='?', default="None", help="样本名称")
    parser.add_argument("--rule", help="规则名称")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="snakemake command")
    return parser.parse_args()


def generate_job_id():
    return datetime.now().strftime("%Y%m%d%H%M%S")


def load_config(config_file):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

# 由于snakemake的输入里不能直接进行有无sample的判断，因此直接导入整个wildcards，在这里解析出sample

def extract_sample_value(input_string):
    if input_string is None:
        return None
    match = re.search(r'sample=(\w+)', input_string)
    if match:
        return match.group(1)
    else:
        return None

def generate_qsub_command(config, sample, rule, command):
    job_id = generate_job_id()

    # 读取 config 中的参数
    queue = config.get('queue', 'slst_fat')
    nodes = config.get('nodes', 1)
    ppn = config.get('ppn', 4)
    walltime = config.get('walltime', '25:00:00')
    mem = config.get('mem', '8G')

    # 定义输出文件路径，使用当前日期
    current_date = datetime.now().strftime("%Y%m%d")  # 获取当前日期，格式化为 YYYYMMDD

    # 定义pbs名称
    # 提取sample
    sample = extract_sample_value(sample)
    # 文件名基于 sample 和 rule
    if sample == None:
        name = f"rna_{rule}"
    else:
        name = f"rna_{sample}"

    current_dir = os.getcwd()
    output_dir = os.path.join(current_dir, f"pbs/log/{current_date}/{name}_out.log")   # 输出文件路径
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)

    # 将 command 列表转换为字符串
    command = " ".join(command).strip("[]")
    print(command)

    # qsub 命令模板，使用 -j oe 将标准输出和标准错误合并到同一个文件
    qsub_command = (
        f"qsub -q {queue} -N {name}.pbs "
        f"-l nodes={nodes}:ppn={ppn} -l walltime={walltime} -l mem={mem} "
        f"-j oe -o {output_dir} "
        f"-- <<EOF\n#!/bin/bash\n{command}\nEOF"
    )

    return qsub_command


def main():
    # 解析命令行参数
    args = parse_args()
    print(f"Config: {args.config}, Sample: {args.sample}, Rule: {args.rule}, command: {args.command}")

    # 从 config 文件加载配置
    config = load_config(args.config)

    # 生成 qsub 命令
    qsub_command = generate_qsub_command(config, args.sample, args.rule, args.command)

    # 打印 qsub 命令
    print(f"Generated qsub command: {qsub_command}")

    # 执行 qsub 命令
    try:
        result = subprocess.run(qsub_command, shell=True, check=True, text=True, capture_output=True)
        print("Job submitted successfully!")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Failed to submit job:")
        print(e.stderr)


if __name__ == "__main__":
    main()
