# 基于Snakemake的RNA-Seq数据处理

修改日期：2025.03.31

使用方法：

```sh
bash run.sh <fq_dir> -y
# `-y`:确定提交到pbs
```

下面是测试数据的运行情况：

![image-20250331170152131](https://raw.githubusercontent.com/ruqinga/picture/main/2024/image-20250331170152131.png)

确定没问题则通过nohup提交

```
nohup bash run.sh <fq_dir> -y &
```

**pipline: fq.gz → trim → qc → hisat2 → featurecount + cufflink**

目录结构：

```txt
.
├── config # <--- ！！！你需要检查的
│   ├── cluster_config.yaml 
│   └── config.yaml
├── pbs # <--- pbs的输出log
│   └── log
├── Results # <--- 输出的结果文件
│   ├── 02_trim_out
│   ├── 03_qc
│   ├── 04_align_out
│   ├── 05_count
│   ├── 06_norm
│   └── 07_repeats_count
├── run.sh # <--- 提交命令的脚本
└── workflow
    ├── env # <--- 对应conda环境（待完成）
    ├── rules # <--- 定义数据处理规则
    ├── scripts # <--- 存放辅助代码
    └── Snakefile # <--- 流程控制

```

通过简单三步即可部署成功：

## Step1: install snakemake

使用 conda 或 [mamba](https://github.com/mamba-org/micromamba-releases/releases)（推荐） 安装 snakemake，使 Snakemake 能够处理工作流程的软件依赖 。

​    

1）创建一个新的conda环境并安装snakemake和snakedeploy

```
conda create -n snakemake_env -c conda-forge -c bioconda snakemake snakedeploy
# mamba的命令与conda一致，把开头的conda改为mamba就可以了
mamba create -n snakemake_env -c conda-forge -c bioconda snakemake snakedeploy
```

这个命令指定了从 conda-forge 和 bioconda 安装 snakemake

​    

2）激活新创建的环境

```
conda activate snakemake_env
```

​    

3）验证 snakemake 的安装。安装完成后你可以通过运行以下命令来检查 snakemake 是否安装正确以及对应的版本号

```
snakemake --version
```

​    

4）在8.0版本之后都需要安装 `snakemake-executor-plugin-cluster-generic` 插件用于提交任务到 cluster。 👉[发布的声明](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html#id6)

```
pip install snakemake-executor-plugin-cluster-generic
```

​    

## Step2: 下载workflow

​    

从GitHub下载workflow

```
git clone https://github.com/ruqinga/snakemake_rna
```

​    

## Step3: 配置workflow

### 3.1 配置config.yaml

​    

修改conda环境为你的conda环境名称

```
# Environment
conda_env: "base-omics"
conda_cufflink: "py36"
```

​     

修改index路径以及物种信息

```
Spe: "mm10" # Species: "human" or "mm10" or "mm9"
```

​     

其它参数可视需求调整

​    

### 配置cluster_config.yaml

默认配置：

```sh
__default__:
  queue: "slst_pub"      # 默认队列 slst_pub slst_fat  pub_fast
  nodes: 1               # 默认节点数
  ppn: 20                # 默认每节点核心数
  walltime: "30:00:00"   # 默认运行时间
  mem: "100G"              # 默认内存
```

如要添加更多参数，需要修改`workflow/scripts/submit_job.py`