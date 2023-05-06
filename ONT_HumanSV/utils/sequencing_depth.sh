#!/usr/bin/env bash


### Add command-line options
PrintUsage(){
	echo -e "BAM stats program"
  echo -e "Uasge:\n\tPROG <INPUT_BAM> <PREFIX>"
	echo -e "Arguments:\n\t<INPUT_BAM> : Input sorted BAM file (reference: hg38)\n\t\
<PREFIX>    : Sample name, as test.\n"
	echo -e "Notice: Work directory must be current directory, because of BAMStats."
	echo -e "requires:\n\t- samtools\n\t- jdk6"
  exit 1
}
#! 0 args
if [ $# -ne	 2 ];then
  PrintUsage
fi
#~ get arguments
INPUT_BAM=$1
CLEANED_BAM=${INPUT_BAM/.bam/.cleaned_chrom.bam}
SAMPLE_NAME=$2
#~ check output, if exists, delete! (BAMStats needs)
if [ -f ${SAMPLE_NAME}.html ];then
	rm -r ${SAMPLE_NAME}*
fi

### BAMStats
# ~ jdk path
java=/home/daruibio/test77/Software/jdk6/jdk1.6.0_45/bin/java
# ~ BAMStats jar
BAMStats_jar=/home/daruibio/test77/Software/BAMStats-1.25/BAMStats-1.25.jar

### run
# ~ clean reference include "alt, random, chrUn"
samtools view -@ 8 -h $INPUT_BAM | grep -Ev 'alt|random|chrUn' | samtools view -@ 8 -h -Sb -o $CLEANED_BAM -
# ~ samtools coverage
samtools coverage -o ${SAMPLE_NAME}.cov $CLEANED_BAM
# ~ samtools depth
# samtools depth --reference ~/test77/Database/references/homo_sapiens/hg38.fa -o sequencing_depth/test.depth test/minimap2.cleaned_chrom.bam
# ~  work directory must be current directory
$java -jar -Xmx8g $BAMStats_jar -d -l -m -q -s --view html -i $CLEANED_BAM -o $SAMPLE_NAME
mv $SAMPLE_NAME ${SAMPLE_NAME}.html
