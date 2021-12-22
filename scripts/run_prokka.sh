#!/bin/bash

# parameters
run_prokka=$1
prokka_cpu=$2
prokka_folder=$3
unique_genomes_folder=$4
prokka_ended=$5


if [ $run_prokka == "yes" ]; then
    mkdir -p ${prokka_folder}/prokka && cd ${unique_genomes_folder} && for i in $(ls .); do folder=${i%.*}; echo ${folder}  ; prokka --outdir ${prokka_folder}/prokka/${folder} --prefix ${folder} --addgenes --cpus 10 ${i} --centre X --compliant ; done && echo 'prokka ended!' >> ${prokka_ended}
else
    echo "You did not want to launch prokka, right? :)" >> ${prokka_ended}
fi
