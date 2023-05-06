rule merge_pairs:
    input:
        r1="{result_dir}/1.trim/{sample_id}.1P.fastq",
        r2="{result_dir}/1.trim/{sample_id}.2P.fastq"
    output:
        merge="{result_dir}/2.usearch_pipe/{sample_id}_merged.fasta",
        nmerge_f="{result_dir}/2.usearch_pipe/{sample_id}.notmerged_fwd.fastq",
        nmerge_r="{result_dir}/2.usearch_pipe/{sample_id}.notmerged_rev.fastq"
    params:
        "-fastq_maxdiffs 5 -fastq_pctid 90 -fastq_minovlen 16"
    log:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_merge_pairs.log"
    benchmark:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_merge_pairs.benchmark"
    shell:
        "usearch10 {params} -fastq_mergepairs {input.r1} -reverse {input.r2} \
        -fastaout {output.merge} -report {log} \
        -fastqout_notmerged_fwd {output.nmerge_f} -fastqout_notmerged_rev {output.nmerge_r}"

rule dereplication:
    input:
        "{result_dir}/2.usearch_pipe/{sample_id}_merged.fasta"
    output:
        "{result_dir}/2.usearch_pipe/{sample_id}_uniques.fasta"
    threads:
        config["threads"]
    log:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_dereplication.log"
    benchmark:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_dereplication.benchmark"
    shell:
        "usearch10 -threads {threads} -fastx_uniques {input} -fastaout {output} -sizeout &> {log}"

rule unoise3:
    input:
        "{result_dir}/2.usearch_pipe/{sample_id}_uniques.fasta"
    output:
        zotus="{result_dir}/2.usearch_pipe/{sample_id}_zotus.fa",
        tabbedout="{result_dir}/2.usearch_pipe/{sample_id}_unoise3.txt"
    threads:
        config["threads"]
    log:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_unoise3.log"
    benchmark:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_unoise3.benchmark"
    shell:
        "usearch10 -threads {threads} -unoise3 {input} \
        -zotus {output.zotus} -tabbedout {output.tabbedout} &> {log}"

rule usearch_global:
    input:
        fa="{result_dir}/2.usearch_pipe/{sample_id}_merged.fasta",
        db="{result_dir}/2.usearch_pipe/{sample_id}_zotus.fa"
    output:
        "{result_dir}/2.usearch_pipe/{sample_id}_hits.uc"
    threads:
        config["threads"]
    params:
        "-strand both  -id 0.97 "
    log:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_usearch_global.log"
    benchmark:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_usearch_global.benchmark"
    shell:
        "usearch10 -threads {threads} -usearch_global {input.fa} \
        {params} -db {input.db} -uc {output} &> {log}"

rule sintax:
    input:
        fa="{result_dir}/2.usearch_pipe/{sample_id}_zotus.fa",
        udb=config["usearch_udb"]
    output:
        "{result_dir}/2.usearch_pipe/{sample_id}_reads.sintax"
    threads:
        config["threads"]
    params:
        "-strand both -sintax_cutoff 0.8"
    log:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_sintax.log"
    benchmark:
        "{result_dir}/log/2.usearch_pipe/{sample_id}_sintax.benchmark"
    shell:
        "usearch10 -threads {threads} {params} -sintax {input.fa} -db {input.udb} \
        -tabbedout {output}"