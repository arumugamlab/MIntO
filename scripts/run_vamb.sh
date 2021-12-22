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


if [ $gpu == "yes" ]; then
    if [ $binner == "vamb_256" ]; then
        echo "Launching vamb 256"
        vamb -l 16 -n 256 256  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp
    elif [ $binner == "vamb_384" ]; then
        vamb -l 24 -n 384 384  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp
    elif [ $binner == "vamb_512" ]; then
        vamb -l 32 -n 512 512 --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp
    elif [ $binner == "vamb_768" ]; then
        vamb -l 40 -n 768 786 --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp
    else
        echo "Something went wrong"
    fi 
else
    if [ $binner == "vamb_256" ]; then
    	vamb -l 16 -n 256 256  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp
   	elif [ $binner == "vamb_384" ]; then
        vamb -l 24 -n 384 384  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp
    elif [ $binner == "vamb_512" ]; then
        vamb -l 32 -n 512 512  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp
    elif [ $binner == "vamb_768" ]; then
        vamb -l 40 -n 768 768  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp
	else
       	echo "Something went wrong"
    fi
fi
