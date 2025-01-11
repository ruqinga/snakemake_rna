
import os
import sys
import pandas as pd

# 获取工作目录
workdir = sys.argv[1]

# 获取工作目录下所有子文件夹，并按名称升序排序
subfolders = sorted([f.path for f in os.scandir(workdir) if f.is_dir()])

# 初始化一个空的DataFrame用于存储数据
merged_data = pd.DataFrame()

# 遍历每个子文件夹
for folder in subfolders:
    folder_name = os.path.basename(folder)
    file_path = os.path.join(folder, "genes.fpkm_tracking")

    # 检查文件是否存在
    if os.path.isfile(file_path):
        # 读取文件
        temp_data = pd.read_csv(file_path, sep='\t', header=0)

        # 提取第5和第10列，并指定列名
        temp_data = temp_data[['gene_short_name', 'FPKM']]
        temp_data.columns = ['gene', folder_name]

        # 如果存在重复的基因条目，对其进行平均处理
        temp_data = temp_data.groupby('gene').mean().reset_index()

        # 如果是第一个文件，直接赋值给merged_data
        if merged_data.empty:
            merged_data = temp_data
        else:
            # 将数据合并到merged_data中，按照基因列进行合并
            merged_data = pd.merge(merged_data, temp_data, on='gene', how='outer')

# 将合并后的数据保存到文件中
merged_data.to_csv(os.path.join(workdir, "merged_FPKM.txt"), sep='\t', index=False)

