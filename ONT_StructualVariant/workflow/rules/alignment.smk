rule make_lra_index:
    input:
        reference
    output:
        "{}.gli".format(reference)
    shell:
        "lra index -ONT {input}"
            
rule lra_alignment:
    input:
        fq = fastq,
        ref = reference,
        index = "{}.gli".format(reference)
    output:
        "{workdir}/{sample}/2.alignment/aln.bam"
    threads:
        threads
    params:
        lra_align = "-ONT -p s"
    log:
        "{workdir}/{sample}/log/lra_alignment.log"
    benchmark:
        "{workdir}/{sample}/log/lra_alignment.benchmark"
    shell:
        "seqtk seq -A {input.fq} | \
        lra align {params.lra_align} -t {threads} {input.ref} - | \
        samtools sort -@ 4 -o {output} - && \
        samtools index -@ 4 {output}"
        
rule target_bed:
    input:
        "{workdir}/{sample}/2.alignment/aln.bam"
    output:
        "{workdir}/{sample}/2.alignment/aln.bed"
    params:
        interval = target_interval
    log:
        "{workdir}/{sample}/log/target_bed.log"
    benchmark:
        "{workdir}/{sample}/log/target_bed.benchmark"
    run:
        if target_interval == "":
            shell("bedtools bamtobed -i {input} | \
                bedtools merge -i - | \
                awk '$1!~/_random|chrUn_/' > {output} 2> {log}")
        else:
            shell("cp {params.interval} {output}")