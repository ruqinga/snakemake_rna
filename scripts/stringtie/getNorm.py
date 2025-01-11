# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

def process_file(file_path, file_type, name):
    temp_data = pd.DataFrame(columns=['gene', name])

    print(name)
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith("#") and "transcript" in line:
                fields = line.strip().split("\t")
                if fields[2] == "transcript":
                    attributes = dict(item.strip().split(" ") for item in fields[8].split(";") if item.strip())
                    gene_id = attributes.get("gene_id", "").strip('"')
                    value = float(attributes.get(file_type, 0).strip('"'))
                    temp_data = pd.concat([temp_data, pd.DataFrame([[gene_id, value]], columns=['gene', name])], ignore_index=True)
    return temp_data

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <sample_list_file> <file_type> <output_file>")
        sys.exit(1)

    sample_list_file = sys.argv[1]
    file_type = sys.argv[2]
    output_file = sys.argv[3]

    sample_list = pd.read_csv(sample_list_file, sep='\t', header=None, names=['name', 'file_path'])
    merged_data = pd.DataFrame()

    for _, row in sample_list.iterrows():
        name, file_path = row['name'], row['file_path']
        temp_data = process_file(file_path, file_type, name)
        temp_data.to_csv('test1.csv', sep='\t', index=False)
        # 如果存在重复的基因条目，对其进行平均处理
        temp_data = temp_data.groupby('gene').mean().reset_index()
        temp_data.to_csv('test2.csv', sep='\t', index=False)

         # Merge with existing data
        if merged_data.empty:
            merged_data = temp_data
            print("firt one")
        else:
            merged_data = pd.merge(merged_data, temp_data, on='gene', how='outer')

    merged_data.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

