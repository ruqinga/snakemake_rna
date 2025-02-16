import pandas as pd
from pathlib import Path

# config
input_files = snakemake.input # 输入文件
output_path = Path(snakemake.output.merged_count)  # 合并结果输出文件路径
merge_column = snakemake.params.id_column  # 合并依据的列名


def merge_count_files():

    merged_df = None

    for fp in sorted(input_files):
        # 安全读取文件，跳过注释行
        df = pd.read_csv(fp, sep='\t', comment='#')

        # 提取文件名作为样本ID
        sample_id = '.'.join(Path(fp).stem.split('.')[:-1]) #移除所有拓展名

        print(f"Processing sample_id: {sample_id}")

        # 确定需要的列：基因ID列和表达量列（最后一列）
        required_columns = [merge_column, df.columns[-1]]

        # 检查文件是否包含必要的列
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"文件 {fp} 缺少必要列: {required_columns}")

        # 提取所需列并重命名表达量列为样本ID
        df = df[required_columns].copy()
        df.rename(columns={df.columns[-1]: sample_id}, inplace=True)

        # 合并数据框
        if merged_df is None:
            # 如果是第一次循环，直接初始化合并结果
            merged_df = df
        else:
            # 否则，基于基因ID列进行外连接合并
            merged_df = pd.merge(
                merged_df,
                df,
                on=merge_column,
                how='outer',
                validate='one_to_one'  # 确保基因ID在每个文件中唯一
            )

    # 结果排序和保存
    merged_df.sort_values(merge_column).to_csv(
        output_path,
        sep='\t',
        index=False,
        encoding='utf-8'
    )


if __name__ == '__main__':
    # 创建输出目录（如果不存在）
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merge_count_files()