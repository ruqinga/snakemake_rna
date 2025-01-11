rule align_reads_with_hisat2:
    input:
        read = get_trimmed_list
    output:
        sam = temp("{align_out}/{sample}.sam"),
        sorted_bam = "{align_out}/{sample}.Hisat_aln.sorted.bam",
        sorted_bam_bai = "{align_out}/{sample}.Hisat_aln.sorted.bam.bai",
        sorted_bw = "{align_out}/{sample}.Hisat_aln.sorted.bw"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["hisat2"]["params"],
        index= config["hisat2"]["index"][config["Spe"]],
        align_out = directories["align_out"]
    log:
        log = "{align_out}/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
           hisat2 {params.option} "{params.index}" -U "{input.read}" -S "{output.sam}" > {log.log} 2>&1
        else
           hisat2 {params.option} "{params.index}" -1 {input.read[0]} -2 {input.read[1]} -S "{output.sam}" > {log.log} 2>&1
        fi
 
        samtools sort -@ 20 -o "{output.sorted_bam}" "{output.sam}" > {log.log} 2>&1
        samtools index "{output.sorted_bam}" "{output.sorted_bam_bai}" > {log.log} 2>&1
        
        bamCoverage --bam "{output.sorted_bam}" --outFileName "{output.sorted_bw}" --outFileFormat bigwig --binSize 100 --normalizeUsing RPKM --numberOfProcessors 20 > {log.log} 2>&1
        echo "alignment finished at $(date)" > {log.log} 2>&1
        """
