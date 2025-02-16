# 生成目录
def get_directories(config):
    """
    根据配置文件中的 work_dir 和目录结构返回所有目录的路径。
    """
    work_dir = config["work_dir"]
    directories = {
        key: value.format(work_dir=work_dir) if "{work_dir}" in value else value
        for key, value in config["directories"].items()
    }
    # 打印调试信息
    print(f"Resolved directories: {directories}")
    return directories

# 获取输入文件的basename（不带路径和扩展名）
def get_sample_name(filepath, dt):
    """
    从文件路径中提取样本名，基于数据类型（单端或双端）。
    """
    try:
        if dt == "PE":
            return filepath.split("/")[-1].replace("_1.fastq.gz", "").replace("_2.fastq.gz", "")
        elif dt == "SE":
            return filepath.split("/")[-1].replace(".fastq.gz", "")
        else:
            raise ValueError(f"Unsupported 'dt': {dt}")
    except Exception as e:
        raise ValueError(f"Error in get_sample_name for {filepath}: {e}")

def get_sample_list(config):
    """
    从配置文件中提取样本列表。
    """
    dt = config.get("dt")
    reads = config.get("reads", [])
    samples = {get_sample_name(read["read1"], dt) for read in reads}
    samples_list = list(samples)
    # 打印调试信息
    print("Samples:", samples_list)
    return samples_list

def get_all(directories, samples):
    """
    根据配置生成所有需要的目标文件路径列表。
    """
    dt = config.get("dt")
    is_pe = dt == "PE"

    # 定义不同模式下的文件扩展名
    trim_ext = "_1_val_1.fq.gz" if is_pe else "_trimmed.fq.gz"

    # # 通用文件路径
    # targets = [
    #     "{trim_out}/{sample}" + trim_ext,
    #     "{align_out}/{sample}.Hisat_aln.sorted.bam",
    #     "{align_out}/{sample}.Hisat_aln.sorted.bw",
    #     "{count_out}/{sample}.txt",
    #     "{rep_out}/{sample}_rep.txt",
    #     "{norm_out}/{sample}/genes.fpkm_tracking",
    #     "{norm_out}/merged_fpkm.txt",
    #     "{count_out}/merged_count.txt"
    #     #"{rep_out}/merged_count_rep.txt"
    # ]
    #
    # # 根据模式展开所有目标文件路径
    # all_targets = [
    #     expand(path, **directories, sample=samples) for path in targets
    # ]
    #
    # # 展平列表并打印调试信息
    # flattened_targets = [item for sublist in all_targets for item in sublist]
    # print(flattened_targets)  # 打印调试信息

    flattened_targets = "RNA_results/01_cleandata/{sample}_1_val_1.fq.gz"

    return flattened_targets



def get_fq_list(wildcards):
    if config["dt"] == "SE":
        return f"{config['fq_dir']}/{wildcards.sample}.fastq.gz"
    elif config["dt"] == "PE":
        return [
            f"{config['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{config['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_trimmed_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['trim_out']}/{wildcards.sample}_trimmed.fq.gz"
    elif config["dt"] == "PE":
        return [
            f"{directories['trim_out']}/{wildcards.sample}_1_val_1.fq.gz",
            f"{directories['trim_out']}/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_alined_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['align_out']}/{wildcards.sample}.Hisat_aln.sorted.bam"
    elif config["dt"] == "PE":
        return f"{directories['align_out']}/{wildcards.sample}.Hisat_aln.sorted.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")