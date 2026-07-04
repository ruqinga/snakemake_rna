import sys
import pandas as pd
from pathlib import Path
import logging
import argparse

def setup_logging(log_path: str = None) -> None:
    """配置日志格式和级别。"""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_path:
        handlers.append(logging.FileHandler(log_path, encoding='utf-8'))
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers
    )

def parse_args():
    """命令行模式参数解析。"""
    parser = argparse.ArgumentParser(description='合并多个基因计数文件')
    parser.add_argument('-i','--input_dir', type=str, required=True, help='输入文件夹')
    parser.add_argument('--suffix', type=str, required=True, help='目标文件后缀（如 .txt）')
    parser.add_argument('-o','--output', type=str, default="merged_count.txt", help='输出文件路径')
    parser.add_argument('--id_column', type=str, default="Geneid",  help='用于合并的列名 (default: Geneid)')
    parser.add_argument('--log', type=str, default=None, help='日志输出路径')
    return parser.parse_args()

def get_files_from_dir(input_dir, suffix):
    """从目录中获取指定后缀的文件列表，升序排序。"""
    files = sorted([str(f) for f in Path(input_dir).glob(f'*' + suffix)])
    logging.info(f'共找到{len(files)}个{suffix}文件: {files}')
    return files

def extract_prefix(filename, delimiter):
    """从文件名中提取前缀（基于分隔符分割的第一部分）"""
    parts = filename.split(delimiter)
    if len(parts) < 2:
        raise ValueError(f"文件名未包含分隔符 '{delimiter}': {filename}")
    return parts[0]

def merge_count_files(input_files, output_path, merge_column, suffix):
    """主合并逻辑。"""
    merged_df = None

    for fp in sorted(input_files):
        try:
            df = pd.read_csv(fp, sep='\t', comment='#')
        except Exception as e:
            logging.error(f'文件 {fp} 读取失败: {e}')
            continue

        # 提取样本ID
        sample_id = extract_prefix(Path(fp).name, suffix)
        sample_id = sample_id.replace('_pe', '').replace('_se', '')

        required_columns = [merge_column, df.columns[-1]]

        if not all(col in df.columns for col in required_columns):
            logging.error(f"文件 {fp} 缺少必要列: {required_columns}")
            raise ValueError(f"文件 {fp} 缺少必要列: {required_columns}")

        df = df[required_columns].copy()
        df.rename(columns={df.columns[-1]: sample_id}, inplace=True)

        # 合并
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(
                merged_df,
                df,
                on=merge_column,
                how='outer',
                validate='one_to_one'
            )
        logging.info(f'文件 {fp} 合并完成。当前维度: {merged_df.shape}')

    # 排序与输出
    merged_df.sort_values(merge_column).to_csv(
        output_path,
        sep='\t',
        index=False,
        encoding='utf-8'
    )
    logging.info(f'合并结果已保存到 {output_path}')

def main():
    """兼容snakemake和命令行两种模式的主入口。"""
    # 判断是否为snakemake调用
    if 'snakemake' in globals():
        smk = snakemake
        input_dir = smk.params.input_dir if hasattr(smk.params, 'input_dir') else None
        suffix = smk.params.suffix if hasattr(smk.params, 'suffix') else None
        output_path = Path(smk.output[0]) if isinstance(smk.output, list) else Path(smk.output.merged_count)
        merge_column = smk.params.id_column
        log_path = smk.log[0] if hasattr(smk, 'log') and isinstance(smk.log, list) and len(smk.log) > 0 else None
        setup_logging(log_path)
        input_files = get_files_from_dir(input_dir, suffix)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merge_count_files(input_files, output_path, merge_column)
    else:
        args = parse_args()
        setup_logging(args.log)
        input_files = get_files_from_dir(args.input_dir, args.suffix)
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merge_count_files(input_files, output_path, args.id_column, args.suffix)

if __name__ == '__main__':
    main()