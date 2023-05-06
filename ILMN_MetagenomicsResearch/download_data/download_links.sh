#!/usr/bin/env bash
echo -e "Usage: download_links.sh <filename> \n"
echo -e "Example: download_links.sh /nfs-test/mxf/Project/2.metagenomics/202010_metapipe/config/download_links.txt\n"

cat $1 | while read id
do
  echo "Download link: "${id}
  /home/mxf/.aspera/connect/bin/ascp -QT -l 300m -P33001 \
  -i /home/mxf/.aspera/connect/etc/asperaweb_id_dsa.openssh \
  era-fasp@${id} ./rawdata
done
