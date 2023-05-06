#!/usr/bin/env bash


function print_usage {
  echo 'PROG <work_dir> <group_file>'
  echo -e '\twork_dir  WORK_DIR\t\twork directory, as "./results/5.function_annotation". Include KEGG, eggNOG, CAZy, GO directory, lefse in file'
  echo -e '\tgroup_file  GROUP_FILE\t\tgroup file, as "./config/group.txt"'
  exit 1
}

if test $# -ne 2 || test ! -f $1/KEGG/kegg.ko.merge_table.xls || test ! -f $2;then
  print_usage
fi
#~ global group file
group_file=$2

## activate conda env 
source /nfs-test/mxf/Software/miniconda3/bin/activate base
conda activate lefse
## function run lefse pipeline
function run_lefse {
  merged_table_file=$1
  lefse_in_file=$(echo ${merged_table_file} | python3 -c "a=input();print('.'.join(a.strip().split('.')[:-2] + ['lefse.txt']))")
  lefse_workdir=${lefse_in_file%/*}/lefse
  #~ prepare lefse_in file
  python3 /nfs-test/mxf/Project/2.metagenomics/202010_metapipe/scripts/meta_kit/prepare_lefse_in.py ${merged_table_file} ${group_file}
  ### 创建 lefse 工作目录
  mkdir -p ${lefse_workdir}
  ### 1. 格式化输入文件
  lefse-format_input.py ${lefse_in_file} ${lefse_workdir}/A_lefse.in -c 1 -u 2 -o 1000000
  ### 2. lefse 运算  [default: -l 2.0] Novogene -l 4.0
  run_lefse.py ${lefse_workdir}/A_lefse.in ${lefse_workdir}/B_lefse.res -l 2.0
  # awk -F '\t' '$4>3 && $5!="-"' ${lefse_workdir}/B_lefse.res > ${lefse_workdir}/B_lefse.clean.res
  awk -F '\t' '$5!="-" && $3!=""' ${lefse_workdir}/B_lefse.res > ${lefse_workdir}/B_lefse.clean.res
  ### 3. 绘制 LEfSe 得分值
  lefse-plot_res.py ${lefse_workdir}/B_lefse.clean.res ${lefse_workdir}/C_lefse.lda.pdf --format pdf --dpi 150 --width 16
  lefse-plot_res.py ${lefse_workdir}/B_lefse.clean.res ${lefse_workdir}/C_lefse.lda.png --format png --dpi 150 --width 16
  ### 4. 绘制进化分支图
  # lefse-plot_cladogram.py ${lefse_workdir}/B_lefse.clean.res ${lefse_workdir}/D_lefse.cladogram.pdf --format pdf --dpi 150
  # lefse-plot_cladogram.py ${lefse_workdir}/B_lefse.clean.res ${lefse_workdir}/D_lefse.cladogram.png --format png --dpi 150
}

## run function lefse pipeline
#~ kegg
run_lefse $1/KEGG/kegg.ko.merge_table.xls
#~ eggnog
run_lefse $1/eggNOG/eggnog.og.merge_table.xls
#~ cazy
run_lefse $1/CAZy/cazy_func_lv2.merge_table.xls
#~ go
run_lefse $1/GO/go.merge_table.xls

## deactivate conda env
conda deactivate  ## lefse
conda deactivate  ## base