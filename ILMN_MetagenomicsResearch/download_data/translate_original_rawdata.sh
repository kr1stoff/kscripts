#!/usr/bin/env bash

source /nfs-test/mxf/Software/miniconda3/bin/activate base  ## base
conda activate microbe_analysis  ## microbe_analysis
parallel=/nfs-test/mxf/Software/parallel/bin/parallel

cd rawdata_original/
$parallel -j 0 --xapply echo '{1} ../rawdata/{1}' ::: *.fastq.gz
$parallel -j 0 --xapply "zcat {1} | sed 's/ //g' | bgzip > ../rawdata/{1}" ::: *.fastq.gz
cd -

conda deactivate ## microbe_analysis
conda deactivate  ## base
