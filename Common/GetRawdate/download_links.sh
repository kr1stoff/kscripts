#!/usr/bin/env bash


printUsage(){
	echo -e "PROG <DownloadLinksFile> <OutputDirectory>\n"
	exit 1
}

if [ $# != 2 ];then
	printUsage
elif [ ! -f $1 ];then
	printUsage
fi

# example/PRJNA649101_download_links.txt
# fasp.sra.ebi.ac.uk:/vol1/fastq/SRR123/041/SRR12336741/SRR12336741_1.fastq.gz
# fasp.sra.ebi.ac.uk:/vol1/fastq/SRR123/041/SRR12336741/SRR12336741_2.fastq.gz
DownloadLinks=$1
OutputDirectory=$2
mkdir -p $OutputDirectory

# for line in $(head -n 40 ${DownloadLinks});
# do
#   echo ${line}
#   /home/gwadmin/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i /home/gwadmin/asperaweb_id_dsa.openssh era-fasp@${line} ${OutputDirectory}
# done

parallel -j 8 --xapply 'echo {1} && \
	/home/gwadmin/.aspera/connect/bin/ascp -QT -l 300m -P33001 \
	-i /home/gwadmin/asperaweb_id_dsa.openssh era-fasp@{1} {2}' \
	:::: $DownloadLinks ::: $OutputDirectory
