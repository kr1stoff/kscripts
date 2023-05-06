rule run_nanoplot:
    input:
        fastq
    output:
        directory("{workdir}/{sample}/1.qc_report/rawdata_nanoplot")
    params:
        "--maxlength 30000"
    log:
        "{workdir}/{sample}/log/run_nanoplot.log"
    benchmark:
        "{workdir}/{sample}/log/run_nanoplot.benchmark"
    shell:
        "NanoPlot -t 4 {params} --fastq {input} -o {output} 2> {log}"

rule run_seqkit_stats:
    input:
        fastq
    output:
        "{workdir}/{sample}/1.qc_report/seqkit_stats.txt"
    params:
        "-T -a -b"
    log:
        "{workdir}/{sample}/log/run_seqkit_stats.log"
    benchmark:
        "{workdir}/{sample}/log/run_seqkit_stats.benchmark"
    shell:
        "seqkit stats {params} -j 4 {input} > {output}"

rule run_fastcat:
    input:
        fastq
    output:
        "{workdir}/{sample}/1.qc_report/per_read_stats.txt"
    log:
        "{workdir}/{sample}/log/per_read_stats.log"
    benchmark:
        "{workdir}/{sample}/log/per_read_stats.benchmark"
    shell:
        "fastcat --read {output} {input} > /dev/null"
