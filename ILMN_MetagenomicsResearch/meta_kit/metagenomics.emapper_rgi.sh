#!/usr/bin/env bash

source /nfs-test/mxf/Software/miniconda3/bin/activate base  ## base
#~ emapper
conda activate microbe_analysis_py27
#~ emapper diamond : 64 threads > 4 parallel * 16 threads > 16 threads 
emapper.py -m diamond --data_dir ${eggnog_db} --override --no_annot --no_file_comments --cpu ${Max_Threads}  -i ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_protein_merge.fa -o ${Outpath}/5.function_annotation/merge_out
emapper.py --annotate_hits_table ${Outpath}/5.function_annotation/merge_out.emapper.seed_orthologs --no_file_comments --cpu ${Threads} --data_dir ${local_eggnog_datadir} --override -o ${Outpath}/5.function_annotation/out
conda deactivate  ## microbe_analysis_py27
#~ RGI CARD
conda activate microbe_analysis_py36
#~ cp localDB from Database/CARD
cp -r /nfs-test/mxf/Database/CARD/localDB/ ./
#~ main program
rgi main --input_sequence ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_protein_merge.fa \
  --output_file ${Outpath}/5.function_annotation/gene_catalogue.uniq_protein_merge.rgi \
  --num_threads ${Max_Threads} --input_type protein --local --alignment_tool DIAMOND --clean # --include_loose
#~ clean local loads
rgi clean --local
rm -r ./localDB
conda deactivate  ## microbe_analysis_py36
conda deactivate  ## base