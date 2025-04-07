import re
import logging
import pandas as pd
from pathlib import Path

# 配置
trim_logs_dir = Path(snakemake.params.trim_log_dir)  # 输入目录1
align_logs_dir = Path(snakemake.params.align_log_dir)  # 输入目录2
outputfile = Path(snakemake.output.summary)

# 设置日志
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

def extract_aligned_reads(log_file):
    """Extract aligned read counts from a HISAT2 log file."""
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary pattern: "aligned concordantly exactly 1 time" and ">1 times"
    exactly_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly exactly 1 time", content)
    more_than_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly >1 times", content)

    # Secondary pattern: "aligned exactly 1 time" and ">1 times"
    if not exactly_1_match or not more_than_1_match:
        exactly_1_match = re.search(r"(\d+) \(.*?\) aligned exactly 1 time", content)
        more_than_1_match = re.search(r"(\d+) \(.*?\) aligned >1 times", content)

    # Convert matches to integers and sum them
    aligned_reads = sum(int(m.group(1)) for m in [exactly_1_match, more_than_1_match] if m)
    return aligned_reads


def extract_raw_clean_reads(log_file):
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary regex patterns
    input_reads = re.search(r"Total number of sequences analysed:\s*(\d+)", content)
    less_than_cutoff = re.search(r"Number of sequence pairs removed because at least one read was shorter than the length cutoff \(30 bp\):\s*(\d+)", content)
    N_than_4 = re.search(r"Number of sequence pairs removed because at least one read contained more N\(s\) than the specified limit of 4:\s*(\d+)", content)

    # Secondary fallback patterns
    if not input_reads or not less_than_cutoff or not N_than_4:
        input_reads = re.search(r"(\d+) sequences processed in total", content)
        less_than_cutoff = re.search(r"Sequences removed because they became shorter than the length cutoff of 30 bp:\s*(\d+)", content)
        N_than_4 = re.search(r"Sequences removed because they contained more Ns than the cutoff of 4:\s*(\d+)", content)

    # Convert matches to integers (default to 0 if not found)
    input_reads_value = int(input_reads.group(1)) if input_reads else 0
    less_than_cutoff_value = int(less_than_cutoff.group(1)) if less_than_cutoff else 0
    N_than_4_value = int(N_than_4.group(1)) if N_than_4 else 0

    # Logging warnings if values are missing
    if not input_reads:
        logging.warning(f"{log_file} 没有找到 input_reads")
    if not less_than_cutoff:
        logging.warning(f"{log_file} 没有找到 less_than_cutoff")
    if not N_than_4:
        logging.warning(f"{log_file} 没有找到 N_than_4")

    # Calculate trimmed reads
    trimmed_reads = input_reads_value - less_than_cutoff_value - N_than_4_value
    return input_reads_value, trimmed_reads


def process_logs(trim_logs_dir, align_logs_dir):
    """Process all log files in the input directories."""
    raw_reads = {}
    clean_reads = {}
    mapped_reads = {}
    mapping_rates = {}

    # 处理来自 trim_logs_dir 目录的日志文件
    for log_file in trim_logs_dir.glob("*.log"):
        file_name_without_log = log_file.stem  # 获取去掉 .log 的文件名
        raw, clean = extract_raw_clean_reads(log_file)
        raw_reads[file_name_without_log] = raw
        clean_reads[file_name_without_log] = clean

    # 处理来自 align_logs_dir 目录的日志文件
    for log_file in align_logs_dir.glob("*.log"):
        file_name_without_log = log_file.stem
        aligned_reads = extract_aligned_reads(log_file)
        mapped_reads[file_name_without_log] = aligned_reads
        mapping_rates[file_name_without_log] = aligned_reads / clean_reads.get(file_name_without_log, 1)  # 防止除零错误

    # 合并数据到 DataFrame
    df = pd.DataFrame({
        "raw_reads": raw_reads,
        "clean_reads": clean_reads,
        "mapped_reads": mapped_reads,
        "mapping_rate": mapping_rates
    }).reset_index().rename(columns={"index": "filename"})  # filename 是列，而不是索引

    # 格式化数值列
    df["raw_reads"] = df["raw_reads"].astype(int)  # 保留 0 位小数
    df["clean_reads"] = df["clean_reads"].astype(int)  # 保留 0 位小数
    df["mapped_reads"] = df["mapped_reads"].astype(int)  # 保留 0 位小数
    df["mapping_rate"] = df["mapping_rate"].round(2)  # 保留 2 位小数

    # 结果排序和保存
    df.to_csv(outputfile, sep='\t', index=False, encoding='utf-8')


if __name__ == "__main__":
    process_logs(trim_logs_dir, align_logs_dir)
