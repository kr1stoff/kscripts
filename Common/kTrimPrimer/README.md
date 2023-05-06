# kTrimPrimer
## 介绍  
自编去引物脚本  
- `FQTrimPrimer.py` 支持单端双端FQ数据去引物,双端数据合并pair-end为单端,然后继续分析  
- `SAMTrimPrimer.py` SAM数据去引物  

## 安装
推荐`Conda`安装以下软件  
```
python=3.7.12
pyyaml=6.0
bwa=0.7.17
vsearch=2.22.1
```
配置 `kTrimPrimer/conf/software.yaml` 文件, 把 `bwa,vsearch` 的绝对路径替换上去  

## 更新
- [230421] 参数优化  
    - 降低参数 `--fastq_minovlen` 10 > 5  
    - `--fastqout_notmerged_fwd` | `--fastqout_notmerged_rev` 直接合并到 `fastqout`  
- [230314] 逻辑优化  
    单端测序如果read没有两边都识别到同一对引物, 默认引物都是从5'端开始切 (`FLAG` 0|16 正反向 左或右端切)  
- [221010] 增加部分注释  
- [221011] 逻辑优化   
    - 考虑软剪切问题, seq,qual按照 `CIGAR S` 标识截断  
    - 短read过滤,原始数据和剪切的数据都要在这个阈值内  
