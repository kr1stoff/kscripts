rule run_cuteSV:
    input:
        bam = "{workdir}/{sample}/2.alignment/aln.bam",
        ref = reference
    output:
        vcf = "{workdir}/{sample}/3.variant/cutesv.vcf",
        dir = directory("{workdir}/{sample}/3.variant/cutesv_calls")
    threads:
        threads
    log:
        "{workdir}/{sample}/log/run_cuteSV.log"
    benchmark:
        "{workdir}/{sample}/log/run_cuteSV.benchmark"
    params:
        "--min_size 50 --max_size 1000000 --retain_work_dir \
        --report_readid --min_support 1 --genotype \
        --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3"
    shell:
        "mkdir {output.dir} && \
        cuteSV -t {threads} {params} {input.bam} {input.ref} {output.vcf} {output.dir} 2> {log}"

rule SURVIVOR_lowMQ:
    input:
        "{workdir}/{sample}/2.alignment/aln.bam"
    output:
        lowmq_sam = temp("{workdir}/{sample}/2.alignment/lowMQ.sam"),
        lowmq_bam = temp("{workdir}/{sample}/2.alignment/lowMQ.bam"),
        lowmq_cov = temp("{workdir}/{sample}/2.alignment/lowMQ.cov"),
        lowmq_bed = "{workdir}/{sample}/2.alignment/lowMQ.bed"
    log:
        "{workdir}/{sample}/log/SURVIVOR_lowMQ.log"
    benchmark:
        "{workdir}/{sample}/log/SURVIVOR_lowMQ.benchmark"
    threads:
        threads
    params:
        mapq_threshold = "5",
        bincov = "10 1"
    shell:
        "samtools view -H {input} > {output.lowmq_sam} && \
        samtools view -@ {threads} {input} | awk '$5<{params.mapq_threshold}' >> {output.lowmq_sam} && \
        samtools view -@ {threads} -Sb -h {output.lowmq_sam} > {output.lowmq_bam} && \
        samtools depth {output.lowmq_bam} > {output.lowmq_cov} && \
        SURVIVOR bincov {output.lowmq_cov} {params.bincov} | \
        bedtools sort -i - | bedtools merge -i - > {output.lowmq_bed}"

rule filter_by_region:
    input:
        vcf = "{workdir}/{sample}/3.variant/cutesv.vcf",
        interval = "{workdir}/{sample}/2.alignment/aln.bed",
        lowmq_bed = "{workdir}/{sample}/2.alignment/lowMQ.bed"
    output:
        "{workdir}/{sample}/3.variant/cutesv_filter.vcf"
    log:
        "{workdir}/{sample}/log/filter_by_region.log"
    benchmark:
        "{workdir}/{sample}/log/filter_by_region.benchmark"
    ## [20210621] bedtools [sort -> merge] eliminate redundancy
    params:
        gap = "-header -f 0.5 -v -b " + ucsc + "/UCSC_gap_non_redundancy.hg38.bed",
        centro = "-header -f 0.5 -v -b " + ucsc + "/UCSC_centromeres_non_redundancy.hg38.bed",
        telo = "-header -f 0.5 -v -b " + ucsc + "/UCSC_telomeres_non_redundancy.hg38.bed",
        seg_dup = "-header -f 0.5 -v -b " + ucsc + "/UCSC_segmental_dups_non_redundancy.hg38.bed",
        trf = "-header -f 0.5 -v -b " + ucsc + "/UCSC_single_repeats_non_redundancy.hg38.bed"
    ## filter by region
    ## [1] retaion target region
    ## [2] discard genome special region opverlap 50%: 
    ## lowMQ, gap, centromeres, telomeres, segmental duplications, tanden repeats
    shell:
        "grep -Ev '_alt|_random|chrUn_' {input.vcf} | \
        bedtools intersect -u -header -b {input.interval} -a - | \
        bedtools intersect -v -header -b {input.lowmq_bed} -a - | \
        bedtools intersect {params.gap} -a - | \
        bedtools intersect {params.centro} -a - | \
        bedtools intersect {params.telo} -a - | \
        bedtools intersect {params.seg_dup} -a - | \
        bedtools intersect {params.trf} -a - > {output}"
