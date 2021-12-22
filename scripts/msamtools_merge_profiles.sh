#!/bin/bash

profile_tpm_genes=$1
profile_tpm_list=$2
DB=$3
tmp_dir=$4
output_file=$5


(zgrep -v '#' ${profile_tpm_genes} |cut -f1) > ${tmp_dir}/${DB}_names.txt # row names
var="Genes" 
sed -i "1s/.*/$var/" ${tmp_dir}/${DB}_names.txt #change column name

echo ${profile_tpm_list} | tr ',' "\n" > ${tmp_dir}/${DB}_lists #split list of files by /n
#ls ${working_dir}/*.profile_TPM.txt.gz > ${tmp_dir}/${DB}_lists

echo "paste \\" > ${tmp_dir}/${DB}_profiles
cat ${tmp_dir}/${DB}_lists| while read line; do echo "<(zgrep -v '#' $line |cut -f2) \\" >> ${tmp_dir}/${DB}_profiles; done  # get 2nd column values
echo "> ${tmp_dir}/${DB}_merged.txt" >> ${tmp_dir}/${DB}_profiles

bash ${tmp_dir}/${DB}_profiles

paste ${tmp_dir}/${DB}_names.txt ${tmp_dir}/${DB}_merged.txt > ${tmp_dir}/${output_file}

rm ${tmp_dir}/${DB}_lists
rm ${tmp_dir}/${DB}_merged.txt
rm ${tmp_dir}/${DB}_profiles
rm ${tmp_dir}/${DB}_names.txt
 
