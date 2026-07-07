import re
import logging
import pandas as pd
from pathlib import Path

# 配置
trim_logs_dir = Path("Results/02_trim_out/logs/")  # 输入目录1
align_logs_dir = Path("Results/04_align_out/logs/")  # 输入目录2
outputfile = Path("Results/summary.csv")


# trim_logs_dir = Path(snakemake.params.trim_log_dir)  # 输入目录1
# align_logs_dir = Path(snakemake.params.align_log_dir)  # 输入目录2
# outputfile = Path(snakemake.output.summary)

# 设置日志
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

def extract_clean_aligned_reads(log_file):
    """Extract aligned read counts from a HISAT2 log file."""
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary pattern: "aligned concordantly exactly 1 time" and ">1 times"
    clean_reads_match = re.search(r"(\d+) reads; of these:", content)
    exactly_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly exactly 1 time", content)
    more_than_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly >1 times", content)

    # Secondary pattern: "aligned exactly 1 time" and ">1 times"
    if not exactly_1_match or not more_than_1_match:
        exactly_1_match = re.search(r"(\d+) \(.*?\) aligned exactly 1 time", content)
        more_than_1_match = re.search(r"(\d+) \(.*?\) aligned >1 times", content)

    # Convert matches to integers and sum them
    clean_reads = int(clean_reads_match.group(1)) if clean_reads_match else 0
    aligned_reads = sum(int(m.group(1)) for m in [exactly_1_match, more_than_1_match] if m)
    return clean_reads, aligned_reads


def extract_raw_reads(log_file):
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary regex patterns
    input_reads = re.search(r"Total number of sequences analysed:\s*(\d+)", content)

    # Secondary fallback patterns
    if not input_reads:
        input_reads = re.search(r"(\d+) sequences processed in total", content)

    # Convert matches to integers (default to 0 if not found)
    input_reads_value = int(input_reads.group(1)) if input_reads else 0

    # Logging warnings if values are missing
    if not input_reads:
        logging.warning(f"{log_file} 没有找到 input_reads")

    return input_reads_value


def process_logs(trim_logs_dir, align_logs_dir):
    """Process all log files in the input directories."""
    raw_reads = {}
    clean_reads = {}
    mapped_reads = {}
    mapping_rates = {}

    # 处理来自 trim_logs_dir 目录的日志文件
    for log_file in trim_logs_dir.glob("*.log"):
        if log_file.stat().st_size > 0:  # 检查文件大小是否不为 0
            file_name_without_log = log_file.stem  # 获取去掉 .log 的文件名
            raw = extract_raw_reads(log_file)
            raw_reads[file_name_without_log] = raw

    # 处理来自 align_logs_dir 目录的日志文件
    for log_file in align_logs_dir.glob("*.log"):
        if log_file.stat().st_size > 0:  # 检查文件大小是否不为 0
            file_name_without_log = log_file.stem
            clean, aligned_reads = extract_clean_aligned_reads(log_file)
            mapped_reads[file_name_without_log] = aligned_reads
            clean_reads[file_name_without_log] = clean
            mapping_rates[file_name_without_log] = aligned_reads / clean if clean > 0 else 0 # 防止除零错误

    # 合并数据到 DataFrame
    df = pd.DataFrame({
        "raw_reads": raw_reads,
        "clean_reads": clean_reads,
        "mapped_reads": mapped_reads,
        "mapping_rate": mapping_rates
    }).reset_index().rename(columns={"index": "filename"})  # filename 是列，而不是索引

    #print(f"{df.head()}")

    # 格式化数值列
    df["raw_reads"] = df["raw_reads"].astype("Int64")
    df["clean_reads"] = df["clean_reads"].astype("Int64")
    df["mapped_reads"] = df["mapped_reads"].astype("Int64")
    df["mapping_rate"] = df["mapping_rate"].round(4)  # 保留 4 位小数

    # 结果排序和保存
    df.to_csv(outputfile, sep='\t', index=False, encoding='utf-8')
    print(f"Saved summary to {outputfile}")


if __name__ == "__main__":
    process_logs(trim_logs_dir, align_logs_dir)
