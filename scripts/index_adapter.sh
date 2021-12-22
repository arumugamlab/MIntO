#!/bin/bash

# '''
# Retrieve the adapters by selecting the most abundant index in the first 10,000 headers of the raw fastq files (trimmomatic_adaptors=False)
# An example of the expected header is:
# @XXXXXXXXXXXX X:X:X:TGAACGCTTG+NTCGCTAAGC
# TGAACGCTTG+GTCGCTAAGC are the sequence adapters used to generate the adapters.fa file
# Trimmomatic uses the adapters.fa file in the first step of QC_1.smk - Filter reads by quality step
#
# Authors: Carmen Saenz, Mani Arumugam
# '''

read_fw=$1
sample=$2
outdir=$3
working_dir=$4
omics=$5

zcat ${read_fw} | grep "^@" | head -10000 | cut -f2 -d' ' | cut -f4 -d':' | sort | uniq -c | sed 's/^ \+//' | tr -s ' ' ' ' | sort -k1,1nr -t' ' | head -1 | awk '{printf("'${sample}'\t%s\t%d\n", $2, $1)}' > ${outdir}/index_barcodes.txt
# Get reverse complement of the 1st part of the adapter pair
barcode1=$(grep "^"$sample"[[:space:]]" ${outdir}/index_barcodes.txt | cut -f2 | sed "s/;/+/" | cut -f1 -d'+' | tr 'ATGC' 'TACG' | awk '{for (i=length; i!=0; i--) x=x substr($0,i,1);} END{print x}')
# Get the 2nd part of adapter pair
barcode2=$(grep "^"$sample"[[:space:]]" ${outdir}/index_barcodes.txt | cut -f2 | sed "s/;/+/" | cut -f2 -d'+')
cat ${working_dir}/${omics}.adapters.fa.in | mseqtools subset --input - --list ${working_dir}/${omics}.palindrome.list --uncompressed --output - | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" | sed "s^Adapter.*/^PrefixPE-Ad/^" > ${outdir}/adapters.fa
cat ${working_dir}/${omics}.adapters.fa.in | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" >> ${outdir}/adapters.fa
cat ${working_dir}/${omics}.adapters.fa.in | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" | seqtk seq -r - | tr -s '12' '21' | sed "s^/^_rc/^" >> ${outdir}/adapters.fa

# zcat ${read_fw} | grep "^@" | head -10000 | cut -f2 -d' ' | cut -f4 -d':' | sort | uniq -c | sed 's/^ \+//' | tr -s ' ' ' ' | sort -k1,1nr -t' ' | head -1 | awk '{printf("'${sample}'\t%s\t%d\n", $2, $1)}' > ${outdir}/index_barcodes.txt
# barcode1=$(grep "^"$sample"[[:space:]]" ${outdir}/index_barcodes.txt | cut -f2 | sed "s/;/+/" | cut -f1 -d'+' | perl -ne 'chomp(); my $x=$_; $x =~ tr/ATGC/TACG/; $x=reverse($x); print $x;')
# barcode2=$(grep "^"$sample"[[:space:]]" ${outdir}/index_barcodes.txt | cut -f2 | sed "s/;/+/" | curzv923t -f2 -d'+')
# cat ${working_dir}/adapters.fa.in | ${mseqtools_dir}/mseqtools -i - -l ${working_dir}/palindrome.list -o - | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" | sed "s^Adapter.*/^PrefixPE-Ad/^" > ${outdir}/adapters.fa
# cat ${working_dir}/adapters.fa.in | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" >> ${outdir}/adapters.fa
# cat ${working_dir}/adapters.fa.in | sed "s/N\{6,10\}/${barcode1}/;s/X\{6,10\}/${barcode2}/" | ${mseqtools_dir}/seqUtils -t fasta -m revcomp -w 100 - | sed "s^/1_rc^_rc/2^;s^/2_rc^_rc/1^" >> ${outdir}/adapters.fa
