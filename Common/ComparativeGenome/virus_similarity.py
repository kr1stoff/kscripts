#!/usr/bin/env python
# @CreateTime       : 2023/02/24
# @Author           : mengxf
# @version          : v1.0.0
# @LastModified     : 2023/02/24
# @Describtion      : 仅限病毒. 检索目录下所有FASTA, 并与参考进行比对, 查看相似度

import click
from pathlib import Path
from Bio import Align
from Bio import SeqIO

@click.command(context_settings={'help_option_names': ['-h','--help']})
@click.option('-w', '--workdir', required=True, type=click.Path(exists=True), help='基因组集文件夹, .fa/.fna后缀.')
@click.option('-r', '--reference', required=True, type=click.Path(exists=True), help='参考基因组FA序列.')
@click.option('-o', '--outfile', default='similarity.txt', show_default=True, help='输出相似性数值文件.')
def main(workdir, reference, outfile):
    """病毒基因组相似度计算工具."""
    workdir = Path(workdir)
    aligner = Align.PairwiseAligner() # 比对器
    ref_rec = list(SeqIO.parse(reference , 'fasta'))[0]
    ref_seq = ref_rec.seq
    #out
    fastas = list(workdir.glob('*.fa')) + list(workdir.glob('*.fna'))
    with open(outfile, 'w') as gh:
        for fa in fastas:
            rec = list(SeqIO.parse(fa, 'fasta'))[0]
            sid, seq = rec.id, rec.seq
            similarity = aligner.score(ref_seq, seq) / len(ref_seq)
            gh.write(f'{sid}\t{similarity}\n')

if __name__ == '__main__':
    main()
