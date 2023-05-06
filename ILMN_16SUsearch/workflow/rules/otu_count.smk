rule otu_reads_count:
    input:
        sintax="{result_dir}/2.usearch_pipe/{sample_id}_reads.sintax",
        hitsuc="{result_dir}/2.usearch_pipe/{sample_id}_hits.uc"
    output:
        otu_lineage_reads="{result_dir}/3.otu_count/{sample_id}_otu_lineage_reads.tsv",
        genus_reads="{result_dir}/3.otu_count/{sample_id}_genus_reads.tsv"
    benchmark:
        "{result_dir}/log/3.otu_count/{sample_id}_otu_reads_count.benchmark"
    script:
        "../scripts/single_otu_reads_table.py"