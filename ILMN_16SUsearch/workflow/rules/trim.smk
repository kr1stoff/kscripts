rule run_trimmomatic:
    input:
        config["fastq_dir"] + "/{sample_id}.R1.fq",
        config["fastq_dir"] + "/{sample_id}.R2.fq"
    output:
        "{result_dir}/1.trim/{sample_id}.1P.fastq",
        "{result_dir}/1.trim/{sample_id}.1U.fastq",
        "{result_dir}/1.trim/{sample_id}.2P.fastq",
        "{result_dir}/1.trim/{sample_id}.2U.fastq"
    params:
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:200"
    threads:
        config["threads"]
    log:
        "{result_dir}/log/1.trim/{sample_id}_trimmomatic.log"
    benchmark:
        "{result_dir}/log/1.trim/{sample_id}_trimmomatic.benchmark"
    shell:
        "trimmomatic PE -threads {threads} -phred33 -summary {log} \
        {input} {output} {params}"