rule make_minimap2_index:
    input:
        reference
    output:
        reference.replace("fa", "mmi")
    params:
        "-ax map-ont --MD -Y"
    shell:
        "minimap2 {params} {input} -d {output}"

rule minimap2_alignment:
    input:
        reference.replace("fa", "mmi"),
        fastq
    output:
        "{workdir}/{sample}/2.alignment/aln.bam"
    log:
        "{workdir}/{sample}/log/minimap2_alignment.log"
    benchmark:
        "{workdir}/{sample}/log/minimap2_alignment.benchmark"
    threads:
        threads
    params:
        minimap2_align = "-ax map-ont --MD -Y",
        samtools_sort = "-O BAM"
    shell:
        "minimap2 -t {threads} {params.minimap2_align} {input} | \
        samtools sort -@ {threads} {params.samtools_sort} -o {output} - \
        && samtools index -@ {threads} {output} 2> {log}"
        
rule bam2bed:
    input:
        "{workdir}/{sample}/2.alignment/aln.bam"
    output:
        "{workdir}/{sample}/2.alignment/aln.bed"
    log:
        "{workdir}/{sample}/log/bam2bed.log"
    benchmark:
        "{workdir}/{sample}/log/bam2bed.benchmark"
    run:
        if target_interval == "":
            shell("bedtools bamtobed -i {input} | \
                bedtools merge -i - | \
                awk '$1!~/_random|chrUn_/' > {output} 2> {log}")
        else:
            shell("cp {TARGET_INTERVAL} {output}")