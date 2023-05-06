# Sars_Cov2_Amplicon

## TODO


## 更新内容
- [x] [221101] 流程优化  
    选项加入"选择引物软件"  
- [221031] 流程优化  
    段落icon优化, 过滤前后fastqc  
- [221009] 流程优化  
    - haplotypecaller 改回 freebayes
    - 二聚体长度, 引物去除，1)fastp参数；2)可选去引物软件; 3)质控数据补充; 4)软剪切硬剪切参数;
                            5)pangolin加线程; 6)upload文件选择; 7)原始FQ质控表格统一  
- [220928] 报告结果优化  
    1.fastqc,bam文件目录给链接; 2.变异表格,变异深度/测序深度/变异频率  
- [220801] 发现一分钟的力量  
    fastqc和genome_coverage.R挂后台运行加速  
- [220628] 速度提升  
    1)删除进化树; 2)freebayes换成gatk HaplotypeCaller, gatk过滤使用python加快速度  
- [220613] python 重写 Sars_Cov2_Amplicon 流程  
    1)新冠全球库; 2)支持单端FASTQ输入   
- [2022.03.01 更新]  
    1. 加入 snpEff 变异注释  
    2. 沈博建议流程在 temp 文件中  
- [2022.02.23 重大更新]  
    1. 加入 ivar trim 流程，重要参数 -e  
    2. 开发测试沈杨博士的去 primer 方法  
    3. mask reference 细节
        a. bedtools 统计低深度区域（length>20bp），使用低深度区域 mask reference，基于 masked reference  
        b. 使用 freebayes 进行 call variant，使用上面得到的 vcf 构建 consensus  
- [2022.02.17 重大更新]   
    1. 提高覆盖度阈值  
        a.低覆盖区域屏蔽(<10x)  
        b.freebayes min alt count > 10x min alt fraction > 0.2  
    2. fastp 参数，增加 Q 值及滑窗剪切参数  
- [2022.02.11 重大更新]  
    1. .tsv 后缀改为 .xls  
    2. 进化树图优化，单图长度  
    3. *删除单样本流程中的 InDelRealign 步骤（多样本出现 JAVA 内存不足报错，考虑 gatk3 版本较老，放弃改步骤）  
    4. cov lineage 结果表格 .xlsx  
- [2022.02.07 更新]  
    1. 北京疾控潘阳老师要求 concensus fastq  
    2. 进化树图调整高度，样本多太挤了  
    3. pangolin trans 表格 tsv 后缀改成 xls  
- [2022.01.27 更新]  
    1. 报告细节优化，表头替换为中文，并解释;  
    2. 删除不必要的说明描述?  
[2022.01.21 更新]  
    1. 报告系统升级为公卫生信组统一的 perl 报告框架  
    2. sars_cov2_pipe.sh: trim_galore 替换为 fastp, 需要过滤前后质控信息；fastp json 结果转成 table 格式  
    3. genome_coverage.R 加入新冠基因组层数  
    4. 重新整理 VCF 输出  
    5. pangolin 结果优化整理脚本，加入中文配置文件  
