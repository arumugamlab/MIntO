#!/bin/bash

gpu=$1
binner=$2
contigs_file=$3
depth_file=$4
threads=$5
output=$6


echo "gpu: ${gpu}"
echo "binner: ${binner}"
echo "contigs file: ${contigs_file}"
echo "depth file:  ${depth_file}"
echo "threads: ${threads}"
echo "output: ${output}" 


gpu_option=""
if [ $gpu == "yes" ]; then
    gpu_option="--cuda"
fi
if [ $binner == "vamb_256" ]; then
    vamb_options="-l 16 -n 256 256"
elif [ $binner == "vamb_384" ]; then
    vamb_options="-l 24 -n 384 384"
elif [ $binner == "vamb_512" ]; then
    vamb_options="-l 32 -n 512 512"
elif [ $binner == "vamb_768" ]; then
    vamb_options="-l 40 -n 768 768"
else
    echo "Something went wrong"
fi
vamb $vamb_options $gpu_option --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --outdir $output/tmp
