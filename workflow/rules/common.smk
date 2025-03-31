# 解析input，生成output
class SampleProcessor:
    def __init__(self, config):
        self.reads = config.get("reads", [])
        self.sample_names = []
        self.sample_info = {}
        self.process_reads()

    # 生成sample和对应sample_info
    def process_reads(self):
        for read in self.reads:
            sample_name = read["read1"].split("/")[-1].replace("_1.fastq.gz", "").replace(".fastq.gz", "")
            if sample_name in self.sample_names:
                print(f"Error: Duplicate sample name '{sample_name}' found.")
                exit(1)
            self.sample_names.append(sample_name)

            # 判断是否包含read2
            if "read2" in read and read["read2"]:
                self.sample_info[sample_name] = "PE"
            else:
                self.sample_info[sample_name] = "SE"

    # 根据样品类型（SE/PE）确定trim之后得到的后缀
    def get_trim_ext(self, sample):
        return "1_val_1" if self.sample_info[sample] == "PE" else "trimmed"

    # 所有要生成的文件
    def generate_targets(self, sample):
        dt = 'pe' if self.sample_info[sample] == 'PE' else 'se'
        return [
            f"Results/02_trim_out/{sample}_{self.get_trim_ext(sample)}.fq.gz",
            f"Results/03_qc/rawdata/{sample}_tmp.txt",
            f"Results/03_qc/cleandata/{sample}_tmp.txt",
            f"Results/03_qc/rawdata/multiqc_report.html",
            f"Results/03_qc/cleandata/multiqc_report.html",
            f"Results/04_align_out/{sample}_{dt}.Hisat_aln.sorted.bam",
            f"Results/05_count/{sample}_{dt}.txt",
            f"Results/06_norm/{sample}_{dt}/genes.fpkm_tracking",
            f"Results/07_repeats_count/{sample}_{dt}.txt",
            f"Results/05_count/merged_count.txt",
            f"Results/06_norm/merged_fpkm.txt",
            f"Results/07_repeats_count/merged_count.rep.txt"
        ]

    # 获取所有目标路径
    def get_all_targets(self):
        all_paths = []
        for sample in self.sample_names:
            paths = self.generate_targets(sample)
            all_paths.extend(paths)
        return all_paths


def get_fq_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"{config['fq_dir']}/{wildcards.sample}.fastq.gz"
    else:
        return [
            f"{config['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{config['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]


def get_trimmed_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/02_trim_out/{wildcards.sample}_trimmed.fq.gz"
    else:
        return [
            f"Results/02_trim_out/{wildcards.sample}_1_val_1.fq.gz",
            f"Results/02_trim_out/{wildcards.sample}_2_val_2.fq.gz"
        ]


def get_tidy_all(output, samples, sample_info):
    return [
        f"{output}/{sample}_{ext}_fastqc.zip"
        for sample in samples
        for ext in (["1"] if sample_info[sample] == "PE" else [""])
    ]

def get_qc_trim_all(samples, sample_info):
    return [
        f"Results/03_qc/rawdata/{sample}_{trim}_fastqc.zip"
        for sample in samples
        for trim in (["1_val_1"] if sample_info[sample] == "PE" else ["trimmed"])
    ]