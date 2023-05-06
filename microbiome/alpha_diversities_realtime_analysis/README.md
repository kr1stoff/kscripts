## realtime_single_16Sanalysis.r
Data: 2020/10/20
#### 用法： 
```
Rscript --vanilla realtime_single_16Sanalysis.r [-[-input|i] <character>] [-[-outpath|o] [<character>]] [-[-help|h]]
```
#### I/O
**rawdata 目录下为 example**  
输入问件：每个循环产生的 species abundance table，第一列 science names / tax ids，第二列 abundance(depth)  
输出文件：alpha diversities 文件（每个 cycle 逐行）和 alpha diversities indexes plot 多图
