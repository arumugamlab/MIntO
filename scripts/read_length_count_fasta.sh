#!/bin/bash

# '''
# Count length of fasta sequences
#
# Authors: Carmen Saenz
# '''

fasta_file=$1
read_length_out=$2

gunzip -c s{fasta_file}| awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }'|awk -v OFS="\t" '/^>/ {getline seq}; {gsub(/>/,"",$0); print $0, length(seq)}' > {read_length_out}