rule remake_vcf:
    input:
        "{workdir}/{sample}/3.variant/cutesv_filter.vcf"
    output:
        "{workdir}/{sample}/3.variant/cutesv_remake.vcf"
    log:
        "{workdir}/{sample}/log/remake_vcf.log"
    benchmark:
        "{workdir}/{sample}/log/remake_vcf.benchmark"
    params:
        "{sample}"
    script:
        "../scripts/remake_SVs_vcf.py"
    
rule run_annotsv:
    input:
        "{workdir}/{sample}/3.variant/cutesv_remake.vcf"
    output:
        full = "{workdir}/{sample}/4.annotation/hg38_annotsv_full.tsv",
        split = "{workdir}/{sample}/4.annotation/hg38_annotsv_split.tsv"
    log:
        "{workdir}/{sample}/log/run_annotsv.log"
    benchmark:
        "{workdir}/{sample}/log/run_annotsv.benchmark"
    params:
        full = "-genomeBuild GRCh38 -SVinputInfo 1 -annotationMode full",
        split = "-genomeBuild GRCh38 -SVinputInfo 1 -annotationMode split"
    ## AnnotSV path need export to environment variants
    ## see: https://github.com/lgmgeo/AnnotSV
    shell:
        "$ANNOTSV/bin/AnnotSV {params.full} -SVinputFile {input} -outputFile {output.full} 2> {log} && \
        $ANNOTSV/bin/AnnotSV {params.split} -SVinputFile {input} -outputFile {output.split} 2>> {log}"