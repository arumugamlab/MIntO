#!/usr/bin/env python

'''
MAGs recovery and annotation

1) Run the binning program  (vamb in different option)
2) Run Checkm on all the data
3) Copy the HQ genomes in a folder
4) Run Coverm on HQ (why coverm, becasue it is easier to add a new binner in the case)
5) Retrieving the score for the genomes
6) Retrieving the best and unique set of genomes (with old scored formula)
7) Run prokka on the genomes (prokka) [separate environment]
8) Run taxonomic label on the genomes (PhylopHlan Metagenomic) [separate environment]

Authors: Eleonora Nigro, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

localrules: copy_genomes_in_all, copy_best_genomes, check_prokka_output

# args = sys.argv
# config_path = args[args.index("--configfile") + 1]
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file

# some variables

if config['working_dir'] is None:
    print('ERROR in ', config_path, ': working_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['working_dir']) is False:
    print('ERROR in ', config_path, ': working_dir variable path does not exit. Please, complete ', config_path)
else:
    working_dir = config['working_dir']
    wdir = config['working_dir']

if config['omics'] in ('metaG'):
    omics = config['omics']
else:
    print('ERROR in ', config_path, ': omics variable is not correct. "omics" variable should be metaG.')

if config['minto_dir'] is None:
    print('ERROR in ', config_path, ': minto_dir variable in configuration yaml file is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]
    script_dir=config["minto_dir"]+"/scripts/"

if 'BINNERS' in config:
    if config['BINNERS'] is None:
        print('ERROR in ', config_path, ': BINNERS list is empty. "BINNERS" variable should be vamb_256, vamb_384, vamb_512 and/or vamb_768. Please, complete ', config_path)
    else:
        try:
            if 'BINNERS' in config:
                #print("Samples:")
                for bin in config["BINNERS"]:
                    if bin in ('vae256', 'vae384', 'vae512', 'vae768', 'aaey', 'aaez'):
                        pass
                    else:
                        raise TypeError('BINNERS variable is not correct. "BINNERS" variable should be vamb_256, vamb_384, vamb_512 or vamb_768. Please, complete ', config_path)
        except:
            print('ERROR in ', config_path, ': BINNERS variable is not correct. "BINNERS" variable should be vamb_256, vamb_384, vamb_512 or vamb_768.')
else:
    print('ERROR in ', config_path, ': BINNERS list is empty. "BINNERS" variable should be vamb_256, vamb_384, vamb_512 and/or vamb_768. Please, complete', config_path)


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
    print('WARNING in ', config_path, ': MIntO is using the GPU')
elif config['VAMB_GPU'] == False:
    vamb_gpu = "no"
    print('WARNING in ', config_path, ': MIntO is not using the GPU')
else:
    print('ERROR in ', config_path, ': VAMB_GPU variable is empty. "VAMB_GPU" variable should be yes or no')

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

if config['CHECKM_DATABASE'] is None:
   print('ERROR in ', config_path, ': CHECKM_DATABASE variable is empty. Please, complete ', config_path)
elif path.exists(config['CHECKM_DATABASE']) is False:
   print('ERROR in ', config_path, ': CHECKM_DATABASE variable path does not exit. Please, complete ', config_path)
elif path.exists(config['CHECKM_DATABASE']) is True:
   checkm_db = config["CHECKM_DATABASE"]

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

if config['RUN_PROKKA'] is None:
    print('ERROR in ', config_path, ': RUN_PROKKA variable is empty. "RUN_PROKKA" variable should be yes or no')
elif config['RUN_PROKKA'] == True:
    run_prokka = "yes"
    print('WARNING in ', config_path, ': MIntO is running Prokka on the unique genomes retrieved.')
elif config['RUN_PROKKA'] == False:
    run_prokka = "no"
    print('WARNING in ', config_path, ': MIntO is not running Prokka on the unique genomes retrieved.')
else:
    print('ERROR in ', config_path, ': RUN_PROKKA variable is empty. "RUN_PROKKA" variable should be yes or no')

if config['PROKKA_CPUS'] is None:
    print('ERROR in ', config_path, ': PROKKA_CPUS variable is empty. Please, complete ', config_path)
elif type(config['PROKKA_CPUS']) != int:
    print('ERROR in ', config_path, ': PROKKA_CPUS variable is not an integer. Please, complete ', config_path)

if config['PROKKA_memory'] is None:
    print('ERROR in ', config_path, ': PROKKA_memory variable is empty. Please, complete ', config_path)
elif type(config['PROKKA_memory']) != int:
    print('ERROR in ', config_path, ': PROKKA_memory variable is not an integer. Please, complete ', config_path)

if config['RUN_TAXONOMY'] is None:
    print('ERROR in ', config_path, ': RUN_TAXONOMY variable is empty. "RUN_TAXONOMY" variable should be yes or no')
elif config['RUN_TAXONOMY'] == True:
    print('WARNING in ', config_path, ': MIntO is running taxonomy labelling of the unique set of genomes using PhyloPhlAn3.')
    run_taxonomy = "yes"
elif config['RUN_TAXONOMY'] == False:
    run_taxonomy = "no"
    print('WARNING in ', config_path, ': MIntO is not running taxonomy labelling of the unique set of genomes using PhyloPhlAn3.')
else:
    print('ERROR in ', config_path, ': RUN_TAXONOMY variable is empty. "RUN_TAXONOMY" variable should be yes or no')

if config['TAXONOMY_DATABASE'] is None:
    print('ERROR in ', config_path, ': TAXONOMY_DATABASE variable is empty. Please, complete ', config_path)

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
   #print(taxonomy_db_folder)


def mags_recovery():
    result = expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt", wd = working_dir)
    if (run_prokka == "yes"):
        result.append(expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka.done", wd = working_dir))
    if (run_taxonomy == "yes"):
        result.append(expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy.tsv", wd = working_dir))
    return(result)

rule all:
    input:
        mags_recovery()

##############################

### Run Vamb
rule run_vamb_vae:
    input:
        contigs_file = "{wd}/metaG/8-1-binning/scaffolds.2500.fasta",
        rpkm_file = "{wd}/metaG/8-1-binning/scaffolds.2500.abundance.npz"
    output:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/vae{vbinner}/vae_clusters.tsv"
    params:
        cuda="{}".format("--cuda" if vamb_gpu == "yes" else ""),
        latent=lambda wildcards: int(int(wildcards.vbinner)/16)
    log:
        "{wd}/logs/metaG/mags_generation/run_vamb_vae{vbinner}.log"
    resources:
        mem=config['VAMB_memory'],
        gpu=1 if vamb_gpu == "yes" else 0
    threads:
        4 if vamb_gpu == "yes" else config["VAMB_THREADS"]
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        rmdir $(dirname {output.tsv})
        vamb --fasta {input.contigs_file} --rpkm {input.rpkm_file} --seed 1234 -p {threads} {params.cuda} --outdir $(dirname {output.tsv}) --model vae -l {params.latent} -n {wildcards.vbinner} {wildcards.vbinner}
        """

rule run_vamb_aae:
    input:
        contigs_file = "{wd}/metaG/8-1-binning/scaffolds.2500.fasta",
        rpkm_file = "{wd}/metaG/8-1-binning/scaffolds.2500.abundance.npz"
    output:
        tsv_y="{wd}/metaG/8-1-binning/mags_generation_pipeline/aae/aae_y_clusters.tsv",
        tsv_z="{wd}/metaG/8-1-binning/mags_generation_pipeline/aae/aae_z_clusters.tsv"
    params:
        cuda="{}".format("--cuda" if vamb_gpu == "yes" else "")
    log:
        "{wd}/logs/metaG/mags_generation/run_vamb_aae.log"
    resources:
        mem=config['VAMB_memory'],
        gpu=1 if vamb_gpu == "yes" else 0
    threads:
        4 if vamb_gpu == "yes" else config["VAMB_THREADS"]
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        rmdir $(dirname {output.tsv_y})
        vamb --fasta {input.contigs_file} --rpkm {input.rpkm_file} --seed 1234 -p {threads} {params.cuda} --outdir $(dirname {output.tsv_y}) --model aae
        """

rule aae_tsv:
    input:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/aae/aae_{latent_type}_clusters.tsv",
    output:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/aae{latent_type}_clusters.tsv",
    shell:
        """
        cat {input} | sed "s/^aae_{wildcards.latent_type}_//" > {output}
        """

rule vae_tsv:
    input:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/vae{vbinner}/vae_clusters.tsv",
    output:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/vae{vbinner}_clusters.tsv",
    shell:
        """
        cat {input} | sed "s/^vae_//" > {output}
        """

### Select MAGs that satisfy min_mag_length criterion
# this is on vamb, if there are other binners, depending on the output, the bins should be processed differently
rule make_avamb_mags:
    input:
        tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}_clusters.tsv",
        contigs_file = "{wd}/metaG/8-1-binning/scaffolds.2500.fasta"
    output:
        discarded_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}_discarded_genomes.txt",
        bin_folder = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/bins"),
    params:
        min_mag_length = config["MIN_MAG_LENGTH"],
        binsplit_char = config["BINSPLIT_CHAR"]
    log:
        "{wd}/logs/metaG/mags_generation/avamb{binner}.take_all_genomes_for_each_run.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (\
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
        checkm_groups = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/checkm")
    params:
        batch_size = checkm_batch_size
    shell:
        """
        cd {input.bin_folder}
        files=$(ls *.fna | wc -l)
        batches=$(($files/{params.batch_size}))
        batches=$(expr $batches + 1)
        for i in $(seq -w 1 $batches); do
          mkdir -p {output.checkm_groups}/$i
          touch {output.checkm_groups}/$i.batch
        done
        for i in $(seq -w 1 $batches); do
          list=$(ls *.fna | head -{params.batch_size})
          mv $list {output.checkm_groups}/$i/
        done
        """

###############################
# Get the batches per binner
###############################

def get_checkm_output_for_batches(wildcards):
    #Collect the genome bins from previous step
    checkpoint_output = checkpoints.prepare_bins_for_checkm.get(**wildcards).output[0]
    result = expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/checkm/{batch}.out/quality_report.tsv",
                    wd=wildcards.wd,
                    binner=wildcards.binner,
                    batch=glob_wildcards(os.path.join(checkpoint_output, '{batch}.batch')).batch)
    return(result)

########################
# CheckM on a batch
########################

rule checkm_batch:
    input:
        '{somewhere}/{something}.batch'
    output:
        '{somewhere}/{something}.out/quality_report.tsv'
    log:
        '{somewhere}/{something}.checkM.log'
    params:
        checkm_db = checkm_db
    conda:
        config["minto_dir"]+"/envs/checkm2.yml"
    threads: 16
    resources:
        mem = 32
    shell:
        """
        tmp=$(mktemp -d)
        time ( \
        checkm2 predict --quiet --database_path {params.checkm_db} -x fna --remove_intermediates --threads {threads} --input {wildcards.somewhere}/{wildcards.something} --tmpdir $tmp -o $(dirname {output}) 
        ) >& {log}
        rm -rf $tmp
        """

########################
# Merge batches of checkm for a single binner
########################

rule merge_checkm_batches:
    input:
        get_checkm_output_for_batches
    output:
        binner_combined = "{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}.checkM.txt"
    log:
        "{wd}/logs/metaG/mags_generation/{binner}.checkM.merge.log"
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

rule move_bins_after_checkm:
    input:
        checkm_report=rules.merge_checkm_batches.output.binner_combined,
        fna_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/checkm/"
    output:
        moved="{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/fna.moved",
    params:
        batch_size = checkm_batch_size
    shell:
        """
        cd {input.fna_folder}
        for i in *.batch; do
          batch=${{i%%.batch}}
          mv $batch/*.fna ../bins/
          rmdir $batch
          rm $batch.checkM.log
          rm -rf $batch.out
          rm $i
        done
        touch {output.moved}
        """

rule copy_genomes_in_all:
    input:
        bin_folder = lambda wildcards: expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/bins",
                                wd = wildcards.wd,
                                binner = config['BINNERS']),
        moved = lambda wildcards: expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/fna.moved",
                                wd = wildcards.wd,
                                binner = config['BINNERS'])
    output:
        all_genomes = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/all"),
        copied = "{wd}/metaG/8-1-binning/mags_generation_pipeline/copy_genomes_all_finished.txt"
    shell:
        """
        mkdir -p {output.all_genomes}
        for i in {input.bin_folder}; do
          cp $i/*.fna {output.all_genomes}/
        done
        touch {output.copied}
        """

########################
# Create a comprehensive table with checkm from all binners
########################

rule make_comprehensive_table:
    input:
        lambda wildcards: expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/{binner}/{binner}.checkM.txt",
                                wd = wildcards.wd,
                                binner = config['BINNERS'])
    output:
        checkm_total = "{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm/checkm-comprehensive.tsv"
    log:
        "{wd}/logs/metaG/mags_generation/make_comprehensive_table.log"#.format(wdir = config['working_dir'])
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
rule copy_HQ_genomes:
    input:
        checkm_total = rules.make_comprehensive_table.output,
        copied = rules.copy_genomes_in_all.output.copied
    output:
        HQ_table="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv",
        HQ_folder=directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes")
    params:
        all_genomes_folder = "{wd}/metaG/8-1-binning/mags_generation_pipeline/avamb/all/",
        completeness = config["CHECKM_COMPLETENESS"],
        contamination = config["CHECKM_CONTAMINATION"]
    log:
        "{wd}/logs/metaG/mags_generation/copy_HQ_genomes.log"
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
                print("[rule copy_HQ_genomes] Copying {} to {}".format(source_file, destination_file), file=f)
                subprocess.call(["cp", source_file, destination_file])

## Run coverm on HQ genomes to create the .tsv file
rule run_coverm:
    input:
        HQ_table=rules.copy_HQ_genomes.output.HQ_table
    output:
        cluster_tsv="{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv"
    params:
        HQ_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes"
    log:
        "{wd}/logs/metaG/mags_generation/run_coverm.log"
    resources:
        mem=config["COVERM_memory"]
    threads:
        config["COVERM_THREADS"]
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (coverm cluster --genome-fasta-directory {params.HQ_folder} -x fna --ani 99 --output-cluster-definition {output.cluster_tsv} --threads {threads} --precluster-method finch) &> {log}
        """

## Run retrieving scored
rule calculate_score_genomes:
    input:
        cluster_tsv = rules.run_coverm.output.cluster_tsv,
        HQ_table = rules.copy_HQ_genomes.output.HQ_table
    output:
        scored_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm_scored.tsv"
    params:
        HQ_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes",
        #calculate_genomes_score="{script_dir}/calculate_genomes_score.py"
        score_method = config["SCORE_METHOD"]
    log:
        "{wd}/logs/metaG/mags_generation/calculate_score_genomes.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (python {script_dir}calculate_genomes_score.py --checkm_output {input.HQ_table} --fasta_folder {params.HQ_folder} --output_file {output.scored_genomes} --score_method {params.score_method}) &> {log}
        """


## Run retrieved the best unique genomes
rule find_unique_and_best_genomes:
    input:
        scored_genomes = rules.calculate_score_genomes.output,
        coverm = rules.run_coverm.output
    output:
        scored = "{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster_scored.tsv",
        best_unique_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"
    log:
        "{wd}/logs/metaG/mags_generation/find_unique_and_best_genomes.log"
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
        best_unique_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"
    output:
        genome_dir = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes")
    log:
        "{wd}/logs/metaG/mags_generation/copy_best_genomes.log"
    resources:
        mem=10
    threads:
        1 # Decide number of threads
    shell:
        """
        time (mkdir -p {output.genome_dir}
        while read line; do
          cp {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes/${{line}}.fna {output.genome_dir}/ ;
        done < {input.best_unique_genomes}
        )&> {log}
        """

########################
# Prokka on a fna file
########################

rule prokka_for_genome:
    input:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes/{mag}.fna"
    output:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka/{mag}/{mag}.gbk"
    log:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka/{mag}.log"
    resources:
        mem=config["PROKKA_memory"]
    threads:
        config["PROKKA_CPUS"]
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        rm -rf $(dirname {output})
        prokka --outdir $(dirname {output}) --prefix {wildcards.mag} --addgenes --cdsrnaolap --cpus {threads} --centre X --compliant {input} >& {log}
        """

###############################
# Get prokka output list for all MAGs in checkpoint above
###############################

def get_hq_mag_prokka_files(wildcards):
    #Collect the genome bins from previous step
    checkpoint_output = checkpoints.copy_best_genomes.get(**wildcards).output[0]
    result = expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka/{mag}/{mag}.gbk",
                    wd=wildcards.wd,
                    mag=glob_wildcards(os.path.join(checkpoint_output, '{mag}.fna')).mag)
    return(result)

rule check_prokka_output:
    input:
        get_hq_mag_prokka_files
    output:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka.done"
    shell:
        """
        touch {output}
        """

########################
# PhyloPhlAn on a fna file
########################

rule taxonomy_for_genome_collection:
    input:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes"
    output:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy.tsv"
    log:
        "{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy.log"
    params:
        run_taxonomy = "{run_taxonomy}".format(run_taxonomy = run_taxonomy),
        taxonomy_database_folder = "{taxonomy_db_folder}".format(taxonomy_db_folder = taxonomy_db_folder),
        taxonomy_database = config["TAXONOMY_DATABASE"],
    resources:
        mem=config["TAXONOMY_memory"]
    threads:
        config["TAXONOMY_CPUS"]
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        cd $(dirname {output})
        phylophlan_metagenomic -i {input} --nproc {threads} -d {params.taxonomy_database} -o taxonomy --taxonomy_database_folder {params.taxonomy_database_folder}
        """
