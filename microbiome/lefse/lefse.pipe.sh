#!/usr/bin/env bash

## Add command-line options
PrintUsage(){
  echo -e "\nUasge:\n\tPROGRAM -i <filename> -o <filepath> [args...]"
  echo -e "\nExample:\n/nfs-test/mxf/Software/Kscripts/lefse.pipe.sh -i otu_table.txt -t tax_assignment.txt -o lefse\n"
  exit 1
}
## 0 args
if [ $# -eq 0 ];then
  PrintUsage
fi
## get option
while getopts "i:t:o:h" args
do
  case $args in
     i)
        echo "Input file name: $OPTARG"
        infile=$OPTARG
        ;;
     t)
        echo "Tax assignment file name: $OPTARG"
        tax_assignment=$OPTARG
        ;;
     o)
        echo "Output directory path:$OPTARG"
        outpath=$OPTARG
        ;;
     h)
        PrintUsage
        ;;
     ?)
        PrintUsage
        ;;
  esac
done

## 需要进入环境 lefse
source /nfs-test/mxf/Software/miniconda3/bin/activate base
## prepare input file
conda activate qiime
lefse_raw=${infile%/*}/lefse
mkdir -p ${lefse_raw}
sed '1i\# Constructed from biom file' ${infile} > ${infile%.*}.tsv
### otu_table.tsv 转换为 biom 格式
biom convert -i ${infile%.*}.tsv -o ${lefse_raw}/lefse.otu_table.biom --table-type="OTU table" --to-json
### 添加 otu 注释信息至 biom 格式的 otu 表 (qiime)
biom add-metadata -i ${lefse_raw}/lefse.otu_table.biom --observation-metadata-fp ${tax_assignment} -o ${lefse_raw}/lefse.otu_table.16S.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
### 对 otu 表的 "界门纲目科属种" 统计 (qiime)
summarize_taxa.py -L 1,2,3,4,5,6,7 -i ${lefse_raw}/lefse.otu_table.16S.biom -o ${lefse_raw}/taxonomy
### 转换为 LEfSe 的输入格式, 合并"界门纲目科属"统计的类群，并排序
python2 /nfs-test/mxf/Software/Kscripts/lefse/1-summarize_trans.py -i ${lefse_raw}/taxonomy --prefix lefse.otu_table.16S
python2 /nfs-test/mxf/Software/Kscripts/lefse/2-lefse_trans.py ${lefse_raw}/taxonomy/lefse.otu_table.16S_all.txt raw_data/lefse.group.txt ${lefse_raw}/lefse_in.txt
conda deactivate  ## qiime

## run lefse
conda activate lefse
### 创建 lefse 工作目录
mkdir -p ${outpath}
### 1. 格式化输入文件
lefse-format_input.py ${lefse_raw}/lefse_in.txt ${outpath}/A_lefse.in -c 1 -u 2 -o 1000000
### 2. lefse 运算  [default: -l 2.0]
run_lefse.py ${outpath}/A_lefse.in ${outpath}/B_lefse.res -l 3.0
### 3. 绘制 LEfSe 得分值
awk -F '\t' '$4>3 && $5!="-"' ${outpath}/B_lefse.res > ${outpath}/B_lefse.clean.res
lefse-plot_res.py ${outpath}/B_lefse.clean.res ${outpath}/C_lefse.lda.pdf --format pdf --dpi 150 --width 16
### 4. 绘制进化分支图
lefse-plot_cladogram.py ${outpath}/B_lefse.clean.res ${outpath}/D_lefse.cladogram.pdf --format pdf --dpi 150
conda deactivate  ## lefse
conda deactivate  ## base
