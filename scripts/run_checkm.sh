#!/bin/bash
checkm_threads=$1
pplacer_threads=$2


for folder in $(ls .);
do
checkm lineage_wf ${folder} checkm${folder} -x fna --threads ${checkm_threads} --pplacer_threads ${pplacer_threads} --tab_table -f checkm${folder}.tsv
done
