configfile: "config.yaml"

include: "rules/common.smk"

# 获取全局参数
print(f"Config: {config}")

directories = get_directories(config)
samples = get_sample_list(config)


rule all:
    input: get_all(directories, samples)


#load rules
include: "rules/trim.smk"
include: "rules/hisat2.smk"
include: "rules/extract_count.smk"
include: "rules/norm.smk"
include: "rules/tidy.smk"
