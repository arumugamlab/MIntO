#!/bin/bash

#profile_tpm_genes=$1
profile_tpm_list=$1
ident=$2
DB=$3
#working_dir=$4
tmp_dir=$4
#output_file=$6

#sh scripts/msamtools_stats.sh '/emc/cbmr/users/rzv923/SHIME_Cdiff//metaG/6-mapping-profiles/BWA_reads-MAGs/SN17/SN17.8-1-binning_unique_best.p95.profile.abund.all.txt.gz,/emc/cbmr/users/rzv923/SHIME_Cdiff//metaG/6-mapping-profiles/BWA_reads-MAGs/SN41/SN41.8-1-binning_unique_best.p95.profile.abund.all.txt.gz' 8-1-binning_unique_best.p95.profile.abund.all /data/rzv923/
# (zgrep -v '#' ${profile_tpm_genes} |cut -f1) > ${tmp_dir}/${DB}_names.txt # row names
# var="Genes" 
# sed -i "1s/.*/$var/" ${tmp_dir}/${DB}_names.txt #change column name

echo ${profile_tpm_list} | tr ',' "\n" > ${tmp_dir}/${DB}_lists #split list of files by /n

for line in $(cat ${tmp_dir}/${DB}_lists); do
  file=$(basename ${line}); sample=${file%%.p${ident}*}
  echo -n "$sample "
  echo -n $(zcat $line | head | grep "Mapped inserts" | cut -f2 -d'(' | sed "s/%.*//")
  echo -e -n '\n'
done >> ${tmp_dir}/${DB}.maprate.txt;

for line in $(cat ${tmp_dir}/${DB}_lists); do
  file=$(basename ${line}); sample=${file%%.p${ident}*}
  echo -n "$sample "
  echo -n $(zcat $line | head | grep "Mapped inserts" | cut -f2 -d':' | sed "s/^ //")
  echo -e -n '\n'
done >> ${tmp_dir}/${DB}.mapstats.txt;

for line in $(cat ${tmp_dir}/${DB}_lists); do
  file=$(basename ${line}); sample=${file%%.p${ident}*}
  echo -n "$sample "
  echo -n $(zcat $line | head | grep "Multiple mapped" | cut -f2 -d'(' | sed "s/%.*//")
  echo -e -n '\n'
done >> ${tmp_dir}/${DB}.multimap.txt

rm ${tmp_dir}/${DB}_lists

