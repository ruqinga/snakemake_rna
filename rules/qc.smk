rule fastqc_rawdata:
    input:
        rawdata=get_fq_list
    output:
        directory("{qc_out}/rawdata/fastqc_done")
    group: "processing_group"
    params:
        qc_out=lambda wildcards: f"{directories['qc_out']}/rawdata"
    log:
        "{qc_out}/rawdata/logs/fastqc_raw.log"
    shell:
        """
        mkdir -p {params.qc_out}/logs
        fastqc {input.rawdata} -o {params.qc_out} > {log} 2>&1
        touch {output}
        """

rule fastqc_trim:
    input:
        trimmed_fq=get_trimmed_list
    output:
        directory("{qc_out}/trimmed/fastqc_done")
    group: "processing_group"
    params:
        qc_out=lambda wildcards: f"{directories['qc_out']}/trimmed"
    log:
        "{qc_out}/trimmed/logs/fastqc_trim.log"
    shell:
        """
        mkdir -p {params.qc_out}/logs
        fastqc {input.trimmed_fq} -o {params.qc_out} > {log} 2>&1
        touch {output}
        """

rule multiqc_rawdata:
    input:
        fastqc_done="{qc_out}/rawdata/fastqc_done"
    output:
        report="{qc_out}/rawdata/multiqc_report.html"
    log:
        "{qc_out}/rawdata/logs/multiqc_raw.log"
    params:
        qc_out=lambda wildcards: f"{directories['qc_out']}/rawdata"
    shell:
        """
        multiqc {params.qc_out} -o {params.qc_out} > {log} 2>&1
        """

rule multiqc_trim:
    input:
        fastqc_done="{qc_out}/trimmed/fastqc_done"
    output:
        report="{qc_out}/trimmed/multiqc_report.html"
    log:
        "{qc_out}/trimmed/logs/multiqc_trim.log"
    params:
        qc_out=lambda wildcards: f"{directories['qc_out']}/trimmed"
    shell:
        """
        multiqc {params.qc_out} -o {params.qc_out} > {log} 2>&1
        """
