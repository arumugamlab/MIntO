#!/usr/bin/env python

'''
MAGs recovery and annotation

1) Run the binning program  (avamb in VAE and AAE modes)
2) Run CheckM2 on all the data
3) Copy the HQ genomes in a folder
4) Run Coverm on HQ (why coverm, because it is easier to add a new binner in the case)
5) Retrieve the score for the genomes
6) Retrieve the best and unique set of genomes (with old scored formula)
7) Run prokka on the genomes (prokka) [separate environment] [ moved to annotation.smk ]
8) Run taxonomic label on the genomes (PhyloPhlAn Metagenomic) [separate environment]

Authors: Eleonora Nigro, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

localrules: aae_tsv, vae_tsv, collect_genomes_from_all_binners, copy_best_genomes, prepare_bins_for_checkm, collect_HQ_genomes

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

# Variables from configuration yaml file

# some variables

if 'BINNERS' in config:
    if config['BINNERS'] is None:
        print('ERROR in ', config_path, ': BINNERS list is empty. "BINNERS" variable should be combinations of vae256, vae384, vae512, vae768, aaey and aaez. Please, complete ', config_path)
    else:
        try:
            if 'BINNERS' in config:
                #print("Samples:")
                for bin in config["BINNERS"]:
                    if bin in ('vae256', 'vae384', 'vae512', 'vae768', 'aaey', 'aaez'):
                        pass
                    else:
                        raise TypeError('BINNERS variable is not correct. "BINNERS" variable should be combinations of vae256, vae384, vae512, vae768, aaey and aaez. Please, complete ', config_path)
        except:
            print('ERROR in ', config_path, ': BINNERS variable is not correct. "BINNERS" variable should be combinations of vae256, vae384, vae512, vae768, aaey and aaez.')
else:
    print('ERROR in ', config_path, ': BINNERS list is empty. "BINNERS" variable should be combinations of vae256, vae384, vae512, vae768, aaey and aaez. Please, complete', config_path)


if config['VAMB_THREADS'] is None:
    print('ERROR in ', config_path, ': VAMB_THREADS variable is empty. Please, complete ', config_path)
elif type(config['VAMB_THREADS']) != int:
    print('ERROR in ', config_path, ': VAMB_THREADS variable is not an integer. Please, complete ', config_path)

if config['VAMB_memory'] is None:
    print('ERROR in ', config_path, ': VAMB_memory variable is empty. Please, complete ', config_path)
elif type(config['VAMB_memory']) != int:
    print('ERROR in ', config_path, ': VAMB_memory variable is not an integer. Please, complete ', config_path)

if config['VAMB_GPU'] is None:
    print('ERROR in ', config_path, ': VAMB_GPU variable is empty. "VAMB_GPU" variable should be yes or no')
elif config['VAMB_GPU'] == True:
    vamb_gpu = "yes"
    print('NOTE: MIntO is using the GPU')
elif config['VAMB_GPU'] == False:
    vamb_gpu = "no"
    print('NOTE: MIntO is not using the GPU')
else:
    print('ERROR in ', config_path, ': VAMB_GPU variable is empty. "VAMB_GPU" variable should be yes or no')

if config['MIN_FASTA_LENGTH'] is None:
    print('ERROR in ', config_path, ': MIN_FASTA_LENGTH variable is empty. Please, complete ', config_path)
elif type(config['MIN_FASTA_LENGTH']) != int:
    print('ERROR in ', config_path, ': MIN_FASTA_LENGTH variable is not an integer. Please, complete ', config_path)

if config['MIN_MAG_LENGTH'] is None:
    print('ERROR in ', config_path, ': MIN_MAG_LENGTH variable is empty. Please, complete ', config_path)
elif type(config['MIN_MAG_LENGTH']) != int:
    print('ERROR in ', config_path, ': MIN_MAG_LENGTH variable is not an integer. Please, complete ', config_path)

if config['CHECKM_COMPLETENESS'] is None:
    print('ERROR in ', config_path, ': CHECKM_COMPLETENESS variable is empty. Please, complete ', config_path)
elif type(config['CHECKM_COMPLETENESS']) != int:
    print('ERROR in ', config_path, ': CHECKM_COMPLETENESS variable is not an integer. Please, complete ', config_path)

if config['CHECKM_CONTAMINATION'] is None:
    print('ERROR in ', config_path, ': CHECKM_CONTAMINATION variable is empty. Please, complete ', config_path)
elif type(config['CHECKM_CONTAMINATION']) != int:
    print('ERROR in ', config_path, ': CHECKM_CONTAMINATION variable is not an integer. Please, complete ', config_path)

checkm_batch_size = 50
if config['CHECKM_BATCH_SIZE'] is None:
    print('WARNING in ', config_path, ': CHECKM_BATCH_SIZE variable is empty. Using 50', config_path)
elif type(config['CHECKM_CONTAMINATION']) != int:
    print('ERROR in ', config_path, ': CHECKM_CONTAMINATION variable is not an integer. Please, complete ', config_path)
else:
    checkm_batch_size = config['CHECKM_BATCH_SIZE']

checkm_db = None
if config['CHECKM_DATABASE'] is None:
   print('ERROR in', config_path, ': CHECKM_DATABASE variable is empty. Please, complete', config_path)
elif path.exists(config['CHECKM_DATABASE']) is False:
   print('ERROR in', config_path, ': CHECKM_DATABASE variable path does not exit. Please, complete', config_path)
elif path.exists(config['CHECKM_DATABASE']) is True:
   checkm_db = config["CHECKM_DATABASE"]

if checkm_db is None:
    raise Exception("CheckM2 database needs to be provided via CHECKM_DATABASE")

if config['COVERM_THREADS'] is None:
    print('ERROR in ', config_path, ': COVERM_THREADS variable is empty. Please, complete ', config_path)
elif type(config['COVERM_THREADS']) != int:
    print('ERROR in ', config_path, ': COVERM_THREADS variable is not an integer. Please, complete ', config_path)

if config['COVERM_memory'] is None:
    print('ERROR in ', config_path, ': COVERM_memory variable is empty. Please, complete ', config_path)
elif type(config['COVERM_memory']) != int:
    print('ERROR in ', config_path, ': COVERM_memory variable is not an integer. Please, complete ', config_path)

if config['SCORE_METHOD'] == 'checkm':
    pass
else:
    print('ERROR in ', config_path, ': SCORE_METHOD variable can only be checkm at the moment!')

run_taxonomy="no"
if config['RUN_TAXONOMY'] is None:
    print('WARNING in ', config_path, ': RUN_TAXONOMY variable is empty. Setting "RUN_TAXONOMY=no"')
elif config['RUN_TAXONOMY'] == True:
    run_taxonomy = "yes"
elif config['RUN_TAXONOMY'] == False:
    run_taxonomy = "no"
    print('NOTE: MIntO is not running taxonomy labelling of the unique MAGs.')
else:
    print('ERROR in ', config_path, ': RUN_TAXONOMY variable is empty. "RUN_TAXONOMY" variable should be yes or no')

if run_taxonomy == "yes":
    if config['TAXONOMY_CPUS'] is None:
        print('ERROR in ', config_path, ': TAXONOMY_CPUS variable is empty. Please, complete ', config_path)
    elif type(config['TAXONOMY_CPUS']) != int:
        print('ERROR in ', config_path, ': TAXONOMY_CPUS variable is not an integer. Please, complete ', config_path)

    if config['TAXONOMY_memory'] is None:
        print('ERROR in ', config_path, ': TAXONOMY_memory variable is empty. Please, complete ', config_path)
    elif type(config['TAXONOMY_memory']) != int:
        print('ERROR in ', config_path, ': TAXONOMY_memory variable is not an integer. Please, complete ', config_path)

    if config['TAXONOMY_DATABASE_FOLDER'] is None:
       print('ERROR in ', config_path, ': TAXONOMY_DATABASE_FOLDER variable is empty. Please, complete ', config_path)
    elif path.exists(config['TAXONOMY_DATABASE_FOLDER']) is False:
       print('ERROR in ', config_path, ': TAXONOMY_DATABASE_FOLDER variable path does not exit. Please, complete ', config_path)
    elif path.exists(config['TAXONOMY_DATABASE_FOLDER']) is True:
       taxonomy_db_folder = config["TAXONOMY_DATABASE_FOLDER"]

    allowed = ('phylophlan', 'gtdb')
    flags = [0 if x in allowed else 1 for x in config['TAXONOMY_NAME'].split(",")]
    if sum(flags) == 0:
        taxonomy=config["TAXONOMY_NAME"]
    else:
        print('ERROR in ', config_path, ': TAXONOMY_NAME variable is not correct. "TAXONOMY_NAME" variable should be phylophlan, gtdb, or combinations thereof.')

    taxonomies_versioned = list()
    taxonomies = taxonomy.split(",")
    for t in taxonomies:
        version="unknown"
        if t == "phylophlan":
            version=config["PHYLOPHLAN_TAXONOMY_VERSION"]
        elif t == "gtdb":
            version=config["GTDB_TAXONOMY_VERSION"]
        taxonomies_versioned.append(t+"."+version)
    print('NOTE: MIntO is running taxonomy labelling of the unique MAGs using [{}].'.format(", ".join(taxonomies_versioned)))


def mags_recovery():
    result = expand("{wd}/{omics}/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt", wd = working_dir, omics = config['omics'])
    if (run_taxonomy == "yes"):
        result.append(expand("{wd}/{omics}/8-1-binning/mags_generation_pipeline/taxonomy.{taxonomy}.tsv", wd = working_dir, omics = config['omics'], taxonomy = taxonomies_versioned))
    return(result)

rule all:
    input:
        mags_recovery()

##############################

### Run Vamb
rule run_vamb_vae:
    input:
        contigs_file = lambda wildcards: expand("{wd}/{omics}/8-1-binning/scaffolds.{min_fasta_length}.fasta",
                                wd = wildcards.wd,
                                omics = wildcards.omics,
                                min_fasta_length = config['MIN_FASTA_LENGTH']),
        rpkm_file    = lambda wildcards: expand("{wd}/{omics}/8-1-binning/scaffolds.{min_fasta_length}.abundance.npz",
                                wd = wildcards.wd,
                                omics = wildcards.omics,
                                min_fasta_length = config['MIN_FASTA_LENGTH'])
    output:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/vae{vbinner}/vae_clusters.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if vamb_gpu == "yes" else ""),
        latent=lambda wildcards: int(int(wildcards.vbinner)/16)
    log:
        "{wd}/logs/{omics}/mags_generation/run_vamb_vae{vbinner}.log"
    resources:
        mem=config['VAMB_memory'],
        gpu=1 if vamb_gpu == "yes" else 0
    threads:
        4 if vamb_gpu == "yes" else config["VAMB_THREADS"]
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        vamb --fasta {input.contigs_file} \
                --rpkm {input.rpkm_file} \
                --seed 1234 \
                -p {threads} \
                {params.cuda} \
                --outdir out \
                --model vae \
                -l {params.latent} \
                -n {wildcards.vbinner} {wildcards.vbinner}
        rsync -a out/* $(dirname {output.tsv})
        """

# TODO: adjust mem based on #samples or better yet #contigs. With 22M contigs, vamb uses 150G.
rule run_vamb_aae:
    input:
        contigs_file = lambda wildcards: expand("{wd}/{omics}/8-1-binning/scaffolds.{min_fasta_length}.fasta",
                                wd = wildcards.wd,
                                omics = wildcards.omics,
                                min_fasta_length = config['MIN_FASTA_LENGTH']),
        rpkm_file    = lambda wildcards: expand("{wd}/{omics}/8-1-binning/scaffolds.{min_fasta_length}.abundance.npz",
                                wd = wildcards.wd,
                                omics = wildcards.omics,
                                min_fasta_length = config['MIN_FASTA_LENGTH'])
    output:
        tsv_y="{wd}/{omics}/8-1-binning/mags_generation_pipeline/aae/aae_y_clusters.tsv",
        tsv_z="{wd}/{omics}/8-1-binning/mags_generation_pipeline/aae/aae_z_clusters.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if vamb_gpu == "yes" else "")
    log:
        "{wd}/logs/{omics}/mags_generation/run_vamb_aae.log"
    resources:
        mem=config['VAMB_memory'],
        gpu=1 if vamb_gpu == "yes" else 0
    threads:
        4 if vamb_gpu == "yes" else config["VAMB_THREADS"]
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        vamb --fasta {input.contigs_file} \
                --rpkm {input.rpkm_file} \
                --seed 1234 \
                -p {threads} \
                {params.cuda} \
                --outdir out \
                --model aae
        rsync -a out/* $(dirname {output.tsv_y})
        """

# Rename bins by stripping 'vae_', 'aae_y_', 'aae_z_' in the bin id. We will add our own prefix in the MAG outputs.
# clusters.tsv file orders entries by bins, but the contigs are not sorted.
# This means that the fna file has contigs in different order each time.
# Sort the clusters.tsv file by (bin, sample, contig_len), so that final fna is reproducible.

rule aae_tsv:
    input:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/aae/aae_{latent_type}_clusters.tsv",
    output:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/aae{latent_type}_clusters.tsv",
    shell:
        """
        cat {input} | sed "s/^aae_{wildcards.latent_type}_//" | sort -k1,1 -k5,5nr -t '_' > {output}
        """

rule vae_tsv:
    input:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/vae{vbinner}/vae_clusters.tsv",
    output:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/vae{vbinner}_clusters.tsv",
    shell:
        """
        cat {input} | sed "s/^vae_//" | sort -k1,1 -k5,5nr -t '_' > {output}
        """

### Select MAGs that satisfy min_fasta_length criterion
# this is on vamb, if there are other binners, depending on the output, the bins should be processed differently
rule make_avamb_mags:
    input:
        tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}_clusters.tsv",
        contigs_file = lambda wildcards: expand("{wd}/{omics}/8-1-binning/scaffolds.{min_fasta_length}.fasta",
                                wd = wildcards.wd,
                                omics = wildcards.omics,
                                min_fasta_length = config['MIN_FASTA_LENGTH']),
    output:
        discarded_genomes = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}_discarded_genomes.txt",
        bin_folder = directory("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/bins"),
    params:
        min_mag_length = config["MIN_MAG_LENGTH"],
        binsplit_char = config["BINSPLIT_CHAR"]
    log:
        "{wd}/logs/{omics}/mags_generation/avamb_{binner}.take_all_genomes_for_each_run.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (\
                mkdir -p {output.bin_folder}
                python {script_dir}/take_all_genomes.py \
                    --vamb_cluster_tsv {input.tsv} \
                    --binsplit_char {params.binsplit_char} \
                    --contigs_file {input.contigs_file} \
                    --assembly_method_name {wildcards.binner} \
                    --min_fasta_length {params.min_mag_length} \
                    --output_folder {output.bin_folder} \
                    --discarded_genomes_info {output.discarded_genomes}
            ) &> {log}
        """

###############################
# Prepare batches for checkM
###############################

checkpoint prepare_bins_for_checkm:
    input:
        bin_folder = rules.make_avamb_mags.output.bin_folder
    output:
        checkm_groups = directory("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/checkm")
    params:
        batch_size = checkm_batch_size
    shell:
        """
        mkdir -p {output.checkm_groups}
        cd {output.checkm_groups}
        rm -rf batch.*
        files=$(find {input.bin_folder}/ -name "*.fna"  | wc -l)
        batches=$(( $files/{params.batch_size} + 1 ))
        suffix_len=$(echo -n $batches | wc -c)
        find {input.bin_folder}/ -name "*.fna" | split - batch. --lines {params.batch_size} --suffix-length=$suffix_len --numeric-suffixes=1
        """

###############################
# Get the batches per binner
###############################

def get_checkm_output_for_batches(wildcards):
    #Collect the genome bins from previous step
    checkpoint_output = checkpoints.prepare_bins_for_checkm.get(**wildcards).output[0]
    result = expand("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/checkm/{batch}.out/quality_report.tsv",
                    wd=wildcards.wd,
                    binner=wildcards.binner,
                    omics=wildcards.omics,
                    batch=glob_wildcards(os.path.join(checkpoint_output, 'batch.{batch}')).batch)
    return(result)

########################
# CheckM on a batch
########################

rule checkm_batch:
    input:
        '{somewhere}/batch.{something}'
    output:
        '{somewhere}/{something}.out/quality_report.tsv'
    log:
        '{somewhere}/{something}.checkM.log'
    shadow:
        "minimal"
    params:
        checkm_db = checkm_db
    conda:
        config["minto_dir"]+"/envs/checkm2.yml"
    threads: 16
    resources:
        mem = 32
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir tmp
        time ( \
        checkm2 predict --quiet --database_path {params.checkm_db} -x fna --remove_intermediates --threads {threads} --input $(cat {input}) --tmpdir tmp -o out
        ) >& {log}
        rsync -a out/quality_report.tsv {output}
        """

########################
# Merge batches of checkm for a single binner
########################

rule merge_checkm_batches:
    input:
        get_checkm_output_for_batches
    output:
        binner_combined = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}.checkM.txt"
    log:
        "{wd}/logs/{omics}/mags_generation/{binner}.checkM.merge.log"
    resources:
        mem=10
    threads:
        2
    run:
        import pandas as pd
        # concatenate all the .tsv file in the folder in order to create a comphresenive file
        li = []
        for filename in input:
            df = pd.read_csv(filename, index_col=None, header=0, sep = "\t")
            li.append(df)
        all_checkm_output = pd.concat(li, axis=0, ignore_index=True)
        # save the file with all the checkm in the same file
        all_checkm_output.to_csv("{}".format(output.binner_combined), sep = "\t", index = False)

########################
# Collect genomes from multiple binners into one location.
# Copying using 'cp location/*.fna' might fail if string consisting list of fna files goes over bash limit for argument length (ARG_MAX).
# ARG_MAX on our system was 2^21, and with our projects we hit it around 12,000 files.
# Therefore, we make the list of files using 'find' and process in batches of 5000 so that ARG_MAX is not reached.
# If this rule fails in your hands with 'cp: Argument list too long', that means that your file paths are longer than 2000 characters.
# It is strange, but not wrong. In that case, please reduce params.batch_size=5000 below and you should be fine.
# Use of readarray was inspired by: https://stackoverflow.com/a/41268405
########################

rule collect_genomes_from_all_binners:
    input:
        checkm_out = lambda wildcards: expand("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}.checkM.txt",
                                wd = wildcards.wd,
                                omics=wildcards.omics,
                                binner = config['BINNERS'])
    output:
        all_genomes = directory("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/all"),
        collected = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/collect_genomes.done"
    log:
        collected = "{wd}/logs/{omics}/mags_generation/collect_genomes.log"
    params:
        batch_size=5000
    shell:
        """
        rm -rf {output.all_genomes}
        mkdir {output.all_genomes}
        time (
        for i in {input.checkm_out}; do
          location=$(dirname $i)
          echo "Collecting $(basename $location) in batches of {params.batch_size}"
          find $location/bins/ -name "*.fna" |
                  while readarray -t -n {params.batch_size} FILES && ((${{#FILES[@]}})); do
                      ln --symbolic --relative ${{FILES[@]}} {output.all_genomes}/
                  done
        done
        ) >& {log}
        touch {output.collected}
        """

########################
# Create a comprehensive table with checkm from all binners
########################

rule make_comprehensive_table:
    input:
        lambda wildcards: expand("{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}.checkM.txt",
                                wd = wildcards.wd,
                                omics=wildcards.omics,
                                binner = config['BINNERS'])
    output:
        checkm_total = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/checkm/checkm-comprehensive.tsv"
    log:
        "{wd}/logs/{omics}/mags_generation/make_comprehensive_table.log"
    resources:
        mem=10
    threads:
        2
    run:
        import pandas as pd
        # concatenate all the .tsv file in the folder in order to create a comphresenive file
        li = []
        for filename in input:
            df = pd.read_csv(filename, index_col=None, header=0, sep = "\t")
            li.append(df)
        all_checkm_output = pd.concat(li, axis=0, ignore_index=True)
        # save the file with all the checkm in the same file
        all_checkm_output.to_csv("{}".format(output.checkm_total), sep = "\t", index = False)

## Copy HQ genomes inside HQ_genomes folder
rule collect_HQ_genomes:
    input:
        checkm_total = rules.make_comprehensive_table.output,
        collected = rules.collect_genomes_from_all_binners.output.collected
    output:
        HQ_table="{wd}/{omics}/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv",
        HQ_folder=directory("{wd}/{omics}/8-1-binning/mags_generation_pipeline/HQ_genomes")
    params:
        all_genomes_folder = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/avamb/all/",
        completeness = config["CHECKM_COMPLETENESS"],
        contamination = config["CHECKM_CONTAMINATION"]
    log:
        "{wd}/logs/{omics}/mags_generation/collect_HQ_genomes.log"
    resources:
        mem=10
    threads:
        2
    run:
        import subprocess
        import pandas as pd
        # open the checkm_comprhrensive table
        checkm_results=pd.read_csv(str(input.checkm_total), sep = "\t")
        # take and save the HQ table
        HQ_checkm_results = checkm_results[(checkm_results["Completeness"] >= params.completeness) & (checkm_results["Contamination"] <= params.contamination)]
        HQ_checkm_results.to_csv(output.HQ_table, sep = "\t", index = False)
        # create the path for copying the genomes
        try:
            os.mkdir(output.HQ_folder)
        except OSError as error:
            print(error)
        # take the bins
        hq_bins = list(HQ_checkm_results["Name"])
        with open(str(log), 'w') as f:
            for bin_id in hq_bins:
                source_file = params.all_genomes_folder +"/{}.fna".format(bin_id)
                destination_file = output.HQ_folder  + "/{}.fna".format(bin_id)
                print("[rule collect_HQ_genomes] Copying {} to {}".format(source_file, destination_file), file=f)
                subprocess.call(["ln", "-s", source_file, destination_file])

## Run coverm on HQ genomes to create the .tsv file
rule run_coverm:
    input:
        HQ_table=rules.collect_HQ_genomes.output.HQ_table,
        checkm_total = rules.make_comprehensive_table.output 
    output:
        cluster_tsv="{wd}/{omics}/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv"
    params:
        HQ_folder="{wd}/{omics}/8-1-binning/mags_generation_pipeline/HQ_genomes"
    log:
        "{wd}/logs/{omics}/mags_generation/run_coverm.log"
    resources:
        mem=config["COVERM_memory"]
    threads:
        config["COVERM_THREADS"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" # coverm
    shell:
        """
        time (coverm cluster --genome-fasta-directory {params.HQ_folder} --checkm2-quality-report {input.checkm_total} -x fna --cluster-method fastani --ani 99 --fragment-length 2500 --min-aligned-fraction 30 --output-cluster-definition {output.cluster_tsv} --threads {threads} --precluster-method finch --precluster-ani 93) &> {log}
        """

## Run retrieving scored
rule calculate_score_genomes:
    input:
        cluster_tsv = rules.run_coverm.output.cluster_tsv,
        HQ_table = rules.collect_HQ_genomes.output.HQ_table
    output:
        scored_genomes = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm_scored.tsv"
    params:
        HQ_folder="{wd}/{omics}/8-1-binning/mags_generation_pipeline/HQ_genomes",
        #calculate_genomes_score="{script_dir}/calculate_genomes_score.py"
        score_method = config["SCORE_METHOD"]
    log:
        "{wd}/logs/{omics}/mags_generation/calculate_score_genomes.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (python {script_dir}/calculate_genomes_score.py --checkm_output {input.HQ_table} --fasta_folder {params.HQ_folder} --output_file {output.scored_genomes} --score_method {params.score_method}) &> {log}
        """


## Run retrieved the best unique genomes
rule find_unique_and_best_genomes:
    input:
        scored_genomes = rules.calculate_score_genomes.output,
        coverm = rules.run_coverm.output
    output:
        scored = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/coverm_unique_cluster_scored.tsv",
        best_unique_genomes = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"
    log:
        "{wd}/logs/{omics}/mags_generation/find_unique_and_best_genomes.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    run:
        import pandas as pd

        # read the table for the score
        score_table = pd.read_csv(str(input.scored_genomes), sep = "\t", index_col = "Bin_id", comment = "#") # we skip the first line with the --score_method

        # read coverm table
        coverm_table = pd.read_csv(str(input.coverm), sep = "\t", names = ["ref_cluster", "cluster_members"])


        # list of best genomes that should be written in the output
        best_genomes_list = []

        # create a dictionary of cluster
        d_cluster = {}

        for i in range(len(coverm_table)):

            ref_cluster = coverm_table["ref_cluster"][i]
            cluster_members = coverm_table["cluster_members"][i].split("/")[-1].replace(".fna", "") # it will append also the name of the genome withput the path

            if ref_cluster not in d_cluster:
                d_cluster[ref_cluster] = [cluster_members]

            else:
                d_cluster[ref_cluster].append(cluster_members)

        # now we take the best genome based on the score
        for cluster in d_cluster:
            genomes_in_the_cluster = d_cluster[cluster]
            dataframe_score = score_table.loc[genomes_in_the_cluster, ["Score"]].sort_values(by=["Score"], ascending=False) # we take the genomes from the score table
            best_genome = dataframe_score.index[0]
            #best_genomes = dataframe_score["Bin_id"][0] # we take the best genome
            best_genomes_list.append(best_genome)

        # Subset score table in order to have the best genomes only
        best_genomes_scored = score_table.loc[best_genomes_list]

        # Save the file
        best_genomes_scored = best_genomes_scored.sort_values(by=["Score"], ascending=False)
        best_genomes_scored.to_csv(output.scored, sep = "\t", index = True)

        # Create the file
        best_unique_genomes_list= list(best_genomes_scored.index)

        with open(output.best_unique_genomes, "w") as fh:
            for genome in best_unique_genomes_list:
                fh.write("{}\n".format(genome))

## Run copy the best genomes
checkpoint copy_best_genomes:
    input:
        best_unique_genomes = "{wd}/{omics}/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"
    output:
        genome_dir = directory("{wd}/{omics}/8-1-binning/mags_generation_pipeline/unique_genomes")
    log:
        "{wd}/logs/{omics}/mags_generation/copy_best_genomes.log"
    resources:
        mem=10
    threads:
        1 # Decide number of threads
    shell:
        """
        time (mkdir -p {output.genome_dir}
        while read line; do
          cp --dereference {wildcards.wd}/{wildcards.omics}/8-1-binning/mags_generation_pipeline/HQ_genomes/${{line}}.fna {output.genome_dir}/ ;
        done < {input.best_unique_genomes}
        )&> {log}
        """

########################
# PhyloPhlAn on fna files
########################

rule phylophlan_taxonomy_for_genome_collection:
    input:
        genomes="{wd}/{omics}/8-1-binning/mags_generation_pipeline/unique_genomes",
        db_folder=lambda wildcards: "{location}/{db_version}".format(location=taxonomy_db_folder, db_version=wildcards.db_version)
    output:
        "{wd}/{omics}/8-1-binning/mags_generation_pipeline/taxonomy.phylophlan.{db_version}.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/mags_generation_pipeline/taxonomy.phylophlan.{db_version}.log"
    resources:
        mem=config["TAXONOMY_memory"]
    threads:
        config["TAXONOMY_CPUS"]
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (
        phylophlan_metagenomic -i {input.genomes} --nproc {threads} -d {wildcards.db_version} -o taxonomy --database_folder {input.db_folder}
        ) >& {log}
        rsync -a taxonomy.tsv {output}
        """

########################
# PhyloPhlAn on a fna file
########################

rule gtdb_taxonomy_for_genome_collection:
    input:
        genomes="{wd}/{omics}/8-1-binning/mags_generation_pipeline/unique_genomes",
        db_folder=lambda wildcards: "{location}/GTDB.{db_version}/taxonomy/gtdb_taxonomy.tsv".format(location=taxonomy_db_folder, db_version=wildcards.db_version)
    output:
        "{wd}/{omics}/8-1-binning/mags_generation_pipeline/taxonomy.gtdb.{db_version}.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/mags_generation_pipeline/taxonomy.gtdb.{db_version}.log"
    resources:
        mem=70
    threads:
        config["TAXONOMY_CPUS"]
    conda:
        config["minto_dir"]+"/envs/gtdb.yml"
    shell:
        """
        time (
        export GTDBTK_DATA_PATH={input.db_folder}
        gtdbtk classify_wf --genome_dir {input.genomes} -x fna --out_dir tmp --cpus {threads} --skip_ani_screen
        cat tmp/gtdbtk.*.summary.tsv | grep "user_genome" | head -1 > taxonomy.tsv
        cat tmp/gtdbtk.*.summary.tsv | grep -v "user_genome" >> taxonomy.tsv
        ) >& {log}
        rsync -a taxonomy.tsv {output}
        """
