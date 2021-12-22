#!/bin/bash

# parameters
run_taxonomy=$1
taxonomy_cpu=$2
unique_genomes_folder=$3
output_phylophlan=$4
database_folder=$5
taxonomy_ended=$6
taxonomy_database=$7



echo "Running Taxonomy: ${run_taxonomy}"

if [ $run_taxonomy == "yes" ]; then
    
    echo "CPUs: ${taxonomy_cpu}"
    echo "Genomes folder: ${unique_genomes_folder}"
    echo "Taxonomy Database ${taxonomy_database}" 

    phylophlan_metagenomic -i ${unique_genomes_folder} --nproc $taxonomy_cpu -d ${taxonomy_database}  -o  ${output_phylophlan} --database_folder ${database_folder} && echo "Taxonomy Generated using ${taxonomy_database}" >> ${taxonomy_ended}
else
    echo "You did not want to launch the taxonomy, right? :)" >> ${taxonomy_ended}
fi
