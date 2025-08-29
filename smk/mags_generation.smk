#!/usr/bin/env python

'''
MAGs recovery and annotation

1) Run the binning program  (vamb in VAE, AAE and VAEVAE modes: default, avamb and taxvamb)
2) Run CheckM2 on all the data
3) Copy the HQ genomes in a folder
4) Run Coverm on HQ (why coverm, because it is easier to add a new binner in the case)
5) Retrieve the score for the genomes
6) Retrieve the best and unique set of genomes (with old scored formula)
7) Run prokka on the genomes (prokka) [separate environment] [ moved to annotation.smk ]
8) Run taxonomic label on the genomes (PhyloPhlAn Metagenomic) [separate environment]

Authors: Eleonora Nigro, Mani Arumugam
'''

import os.path
import math

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/resources.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule mags_base, mags_rpkg, mags_checkm2 from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# Variables from configuration yaml file

# Binners
BINNERS              = validate_required_key(config, 'BINNERS')
for x in BINNERS:
    check_allowed_values('BINNERS', x, ['vae256', 'vae384', 'vae512', 'vae768',
                                        'vaevae256', 'vaevae384', 'vaevae512', 'vaevae768',
                                        'aaey', 'aaez']
                        )

# COVERM params
COVERM_THREADS       = validate_required_key(config, 'COVERM_THREADS')
COVERM_memory        = validate_required_key(config, 'COVERM_memory')

# binning params
MIN_FASTA_LENGTH     = validate_required_key(config, 'MIN_FASTA_LENGTH')
MIN_MAG_LENGTH       = validate_required_key(config, 'MIN_MAG_LENGTH')

# CHECKM params
CHECKM_COMPLETENESS  = validate_required_key(config, 'CHECKM_COMPLETENESS')
CHECKM_CONTAMINATION = validate_required_key(config, 'CHECKM_CONTAMINATION')
CHECKM_BATCH_SIZE    = 50
if (x := validate_optional_key(config, 'CHECKM_BATCH_SIZE')):
    CHECKM_BATCH_SIZE = x

# VAMB params
VAMB_THREADS         = validate_required_key(config, 'VAMB_THREADS')
VAMB_GPU             = validate_required_key(config, 'VAMB_GPU')
if VAMB_GPU:
    print('NOTE: MIntO is using the GPU')
else:
    print('NOTE: MIntO is not using the GPU')

# TAXVAMB annotator
TAXVAMB_ANNOTATOR    = 'metabuli'
if (x := validate_optional_key(config, 'TAXVAMB_ANNOTATOR')):
    TAXVAMB_ANNOTATOR = x
check_allowed_values('TAXVAMB_ANNOTATOR', TAXVAMB_ANNOTATOR, ['mmseqs', 'metabuli'])

# METABULI memory in GB
METABULI_MEM_GB      = 128
if (x := validate_optional_key(config, 'METABULI_MEM_GB')):
    METABULI_MEM_GB = x

# Scoring MAGs
SCORE_METHOD         = validate_required_key(config, 'SCORE_METHOD')
check_allowed_values('SCORE_METHOD', SCORE_METHOD, ['checkm'])


rule all:
    input:
        f"{working_dir}/{omics}/8-1-binning/mags/unique_genomes",
        print_versions.get_version_output(snakefile_name)
    default_target: True

##############################
### Run mmseqs taxonomy
### It needs 800 GB memory!
##############################

rule mmseqs_taxonomy:
    input:
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
        mmseqs_db    = f"{minto_dir}/data/mmseqs/{{mmseqs_tax_db}}/{{mmseqs_tax_db}}"
    output:
        mmseqs_c = f"{{wd}}/{{omics}}/8-1-binning/scaffold_taxonomy/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{mmseqs_tax_db}}.mmseqs.classification.tsv",
        taxonomy = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{mmseqs_tax_db}}.mmseqs.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if VAMB_GPU else ""),
    log:
        "{wd}/logs/{omics}/8-1-binning/mmseqs.taxonomy.{mmseqs_tax_db}.log"
    threads:
        VAMB_THREADS
    resources:
        mem = 800
    conda:
        minto_dir + "/envs/contig_taxonomy.yml"
    shell:
        """
        time (
            mmseqs createdb {input.contigs_file} seqDB
            mmseqs taxonomy seqDB {input.mmseqs_db} results_dir tmp_dir --threads {threads} --tax-lineage 1
            mmseqs createtsv seqDB results_dir taxonomy_raw.tsv --threads {threads}
            taxconverter mmseqs2 -i taxonomy_raw.tsv -o taxonomy_processed.tsv
            rsync -a taxonomy_raw.tsv       {output.mmseqs_c}
            rsync -a taxonomy_processed.tsv {output.taxonomy}
        ) >& {log}
        """

##############################
### Run metabuli taxonomy
##############################

rule metabuli_taxonomy:
    input:
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
        metabuli_db  = f"{minto_dir}/data/metabuli/{{metabuli_tax_db}}"
    output:
        metabuli_c = f"{{wd}}/{{omics}}/8-1-binning/scaffold_taxonomy/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{metabuli_tax_db}}.metabuli.classification.tsv",
        metabuli_r = f"{{wd}}/{{omics}}/8-1-binning/scaffold_taxonomy/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{metabuli_tax_db}}.metabuli.report.tsv",
        metabuli_k = f"{{wd}}/{{omics}}/8-1-binning/scaffold_taxonomy/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{metabuli_tax_db}}.metabuli.krona.html",
        taxonomy   = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.{{metabuli_tax_db}}.metabuli.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/8-1-binning/metabuli.taxonomy.{metabuli_tax_db}.log"
    threads:
        VAMB_THREADS
    resources:
        mem = METABULI_MEM_GB
    conda:
        minto_dir + "/envs/contig_taxonomy.yml"
    shell:
        """
        time (
            metabuli classify {input.contigs_file} {input.metabuli_db}/ results_dir outprefix --threads {threads} --max-ram {resources.mem} --seq-mode 1 --lineage 1
            taxconverter metabuli \
                    -c results_dir/outprefix_classifications.tsv \
                    -r results_dir/outprefix_report.tsv \
                    -o taxonomy_processed.tsv
            rsync -a results_dir/outprefix_classifications.tsv {output.metabuli_c}
            rsync -a results_dir/outprefix_report.tsv          {output.metabuli_r}
            rsync -a results_dir/outprefix_krona.html          {output.metabuli_k}
            rsync -a taxonomy_processed.tsv                    {output.taxonomy}
        ) >& {log}
        """

##############################
### Run Vamb
##############################

# taxvamb mode
# We ask specifically for GTDB taxonomy from either mmseqs or metabuli
# Memory estimates:
# GPU: '2.345e+06 + 1.635e-03*npzsize + 2.051e-04*fastasize' in memKB
# CPU: '2.448e+05 + 1.688e-03*npzsize + 2.457e-04*fastasize' in memKB
# Simplified to 'ceil(2.35 + 1.69e-09*npzsize + 2.46e-10*fastasize)' in memGB
# And add 10GB with each new attempt
rule run_taxvamb:
    input:
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
        rpkm_file    = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.abundance.npz",
        tax_file     = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.taxonomy.GTDB.{TAXVAMB_ANNOTATOR}.tsv",
    output:
        tsv="{wd}/{omics}/8-1-binning/mags/vaevae{vbinner}/vaevae_clusters_unsplit.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if VAMB_GPU else ""),
        latent=lambda wildcards: int(int(wildcards.vbinner)/16)
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/run_taxvamb_vaevae{vbinner}.log"
    resources:
        mem = lambda wildcards, input, attempt: 10*(attempt-1) + math.ceil(2.35 + 1.69e-9*get_file_size(input.rpkm_file) + 2.46e-10*get_file_size(input.contigs_file)),
        gpu=1 if VAMB_GPU else 0
    threads:
        4 if VAMB_GPU else VAMB_THREADS
    conda:
        minto_dir + "/envs/vamb.yml"
    shell:
        """
        vamb bin taxvamb \
                --fasta {input.contigs_file} \
                --abundance {input.rpkm_file} \
                --taxonomy {input.tax_file} \
                -o '' \
                --seed 1234 \
                -p {threads} \
                {params.cuda} \
                -l {params.latent} \
                -n {wildcards.vbinner} {wildcards.vbinner} \
                --outdir out >& {log}
        rsync -a out/* $(dirname {output.tsv})
        """

# VAE mode
# Memory estimates:
# GPU: '2.345e+06 + 1.635e-03*npzsize + 2.051e-04*fastasize' in memKB
# CPU: '2.448e+05 + 1.688e-03*npzsize + 2.457e-04*fastasize' in memKB
# Simplified to 'ceil(2.35 + 1.69e-09*npzsize + 2.46e-10*fastasize)' in memGB
# And add 10GB with each new attempt
rule run_vamb_vae:
    input:
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
        rpkm_file    = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.abundance.npz",
    output:
        tsv="{wd}/{omics}/8-1-binning/mags/vae{vbinner}/vae_clusters_unsplit.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if VAMB_GPU else ""),
        latent=lambda wildcards: int(int(wildcards.vbinner)/16)
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/run_vamb_vae{vbinner}.log"
    resources:
        mem = lambda wildcards, input, attempt: 10*(attempt-1) + math.ceil(2.35 + 1.69e-9*get_file_size(input.rpkm_file) + 2.46e-10*get_file_size(input.contigs_file)),
        gpu=1 if VAMB_GPU else 0
    threads:
        4 if VAMB_GPU else VAMB_THREADS
    conda:
        minto_dir + "/envs/vamb.yml"
    shell:
        """
        vamb bin default \
                --fasta {input.contigs_file} \
                --abundance {input.rpkm_file} \
                -o '' \
                --seed 1234 \
                -p {threads} \
                {params.cuda} \
                -l {params.latent} \
                -n {wildcards.vbinner} {wildcards.vbinner} \
                --outdir out >& {log}
        rsync -a out/* $(dirname {output.tsv})
        """

# AAE mode
# Memory estimates:
# GPU: '2.176e+06 + 1.905e-03*npzsize + 4.859e-04*fastasize' in memKB
# CPU: '1.232e+05 + 1.920e-03*npzsize + 4.662e-04*fastasize' in memKB
# Simplified to 'ceil(2.18 + 1.92e-09*npzsize + 4.86e-10*fastasize)' in memGB
# And add 10GB with each new attempt
rule run_vamb_aae:
    input:
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
        rpkm_file    = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.abundance.npz",
    output:
        tsv_y="{wd}/{omics}/8-1-binning/mags/aae/aae_y_clusters_unsplit.tsv",
        tsv_z="{wd}/{omics}/8-1-binning/mags/aae/aae_z_clusters_unsplit.tsv"
    shadow:
        "minimal"
    params:
        cuda="{}".format("--cuda" if VAMB_GPU else "")
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/run_vamb_aae.log"
    resources:
        mem = lambda wildcards, input, attempt: 10*(attempt-1) + math.ceil(2.18 + 1.85e-9*get_file_size(input.rpkm_file) + 1.73e-10*get_file_size(input.contigs_file)),
        gpu=1 if VAMB_GPU else 0
    threads:
        4 if VAMB_GPU else VAMB_THREADS
    conda:
        minto_dir + "/envs/vamb.yml"
    shell:
        """
        vamb bin avamb \
                --fasta {input.contigs_file} \
                --abundance {input.rpkm_file} \
                -o '' \
                --seed 1234 \
                -p {threads} \
                {params.cuda} \
                --outdir out >& {log}
        rsync -a out/* $(dirname {output.tsv_y})
        """

# Rename bins by stripping 'vae_', 'aae_y_', 'aae_z_' in the bin id. We will add our own prefix in the MAG outputs.
# clusters.tsv file orders entries by bins, but the contigs are not sorted.
# This means that the fna file has contigs in different order each time.
# Sort the clusters.tsv file by (bin, sample, contig_len), so that final fna is reproducible.

rule aae_tsv:
    localrule: True
    input:
        tsv="{wd}/{omics}/8-1-binning/mags/aae/aae_{latent_type}_clusters_unsplit.tsv",
    output:
        tsv="{wd}/{omics}/8-1-binning/mags/vamb/aae{latent_type}_clusters.tsv",
    shell:
        """
        tail -n +2 {input} | sed "s/^aae_{wildcards.latent_type}_//" | sort -k1,1n -k5,5nr -t '_' > {output}
        """

rule vae_tsv:
    localrule: True
    input:
        tsv=rules.run_vamb_vae.output.tsv
    output:
        tsv="{wd}/{omics}/8-1-binning/mags/vamb/vae{vbinner}_clusters.tsv",
    wildcard_constraints:
        vbinner='256|384|512|768'
    shell:
        """
        tail -n +2 {input} | sort -k1,1n -k5,5nr -t '_' > {output}
        """

rule taxvamb_tsv:
    input:
        tsv=rules.run_taxvamb.output.tsv
    output:
        tsv="{wd}/{omics}/8-1-binning/mags/vamb/vaevae{vbinner}_clusters.tsv",
    wildcard_constraints:
        vbinner='256|384|512|768'
    shell:
        """
        tail -n +2 {input} | sort -k1,1n -k5,5nr -t '_' > {output}
        """

# TODO: estimate memory requirement
### Select MAGs that satisfy min_fasta_length criterion
# this is on vamb, if there are other binners, depending on the output, the bins should be processed differently
rule make_vamb_mags:
    input:
        tsv          = "{wd}/{omics}/8-1-binning/mags/vamb/{binner}_clusters.tsv",
        contigs_file = f"{{wd}}/{{omics}}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.fasta.gz",
    output:
        discarded_genomes = "{wd}/{omics}/8-1-binning/mags/vamb/{binner}/{binner}_discarded_genomes.txt",
        bin_folder = directory("{wd}/{omics}/8-1-binning/mags/vamb/{binner}/bins"),
    params:
        min_mag_length = MIN_MAG_LENGTH
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/vamb/{binner}.take_all_genomes_for_each_run.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        minto_dir + "/envs/mags.yml"
    shell:
        """
        time (
            mkdir -p {output.bin_folder}
            python {script_dir}/take_all_genomes.py \
                    --vamb_cluster_tsv {input.tsv} \
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
    localrule: True
    input:
        bin_folder = rules.make_vamb_mags.output.bin_folder
    output:
        checkm_groups = directory("{wd}/{omics}/8-1-binning/mags/vamb/{binner}/checkm")
    params:
        batch_size = CHECKM_BATCH_SIZE
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
    result = expand("{wd}/{omics}/8-1-binning/mags/vamb/{binner}/checkm/{batch}.out/quality_report.tsv",
                    wd=wildcards.wd,
                    omics=wildcards.omics,
                    binner=wildcards.binner,
                    batch=glob_wildcards(os.path.join(checkpoint_output, 'batch.{batch}')).batch)
    return(result)

########################
# CheckM on a batch
########################

rule checkm_batch:
    input:
        fna_list  = "{somewhere}/batch.{something}",
        checkm_db = "{minto_dir}/data/CheckM2_database/uniref100.KO.1.dmnd".format(minto_dir = minto_dir)
    output:
        '{somewhere}/{something}.out/quality_report.tsv'
    log:
        '{somewhere}/{something}.checkM.log'
    shadow:
        "minimal"
    conda:
        minto_dir + "/envs/checkm2.yml"
    threads: 16
    resources:
        mem = 32
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir tmp
        time (
            checkm2 predict --quiet --database_path {input.checkm_db} -x fna --remove_intermediates --threads {threads} --input $(cat {input.fna_list}) --tmpdir tmp -o out
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
        binner_combined = "{wd}/{omics}/8-1-binning/mags/vamb/{binner}/{binner}.checkM.txt"
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/vamb/{binner}.checkM.merge.log"
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
# Create a comprehensive table with checkm from all binners
########################

rule make_comprehensive_table:
    input:
        lambda wildcards: expand("{wd}/{omics}/8-1-binning/mags/vamb/{binner}/{binner}.checkM.txt",
                                wd     = wildcards.wd,
                                omics  = wildcards.omics,
                                binner = BINNERS)
    output:
        checkm_all = "{wd}/{omics}/8-1-binning/mags/checkm/checkm-comprehensive.tsv"
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/make_comprehensive_table.log"
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
        all_checkm_output.to_csv("{}".format(output.checkm_all), sep = "\t", index = False)

## Copy HQ genomes inside HQ_genomes folder
rule collect_HQ_genomes:
    localrule: True
    input:
        checkm_all = rules.make_comprehensive_table.output,
    output:
        checkm_HQ = "{wd}/{omics}/8-1-binning/mags/HQ_genomes_checkm.tsv",
        HQ_folder = directory("{wd}/{omics}/8-1-binning/mags/HQ_genomes")
    params:
        bin_folder    = lambda wildcards: f"{wildcards.wd}/{wildcards.omics}/8-1-binning/mags/vamb/",
        completeness  = CHECKM_COMPLETENESS,
        contamination = CHECKM_CONTAMINATION
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/collect_HQ_genomes.log"
    resources:
        mem=10
    threads:
        2
    run:
        import subprocess
        import pandas as pd

        # open the checkm_comprehensive table
        checkm_results=pd.read_csv(str(input.checkm_all), sep = "\t")

        # take and save the HQ table
        HQ_checkm_results = checkm_results[(checkm_results["Completeness"] >= params.completeness) & (checkm_results["Contamination"] <= params.contamination)]
        HQ_checkm_results.to_csv(output.checkm_HQ, sep = "\t", index = False)

        # create the path for copying the genomes
        try:
            os.mkdir(output.HQ_folder)
        except OSError as error:
            print(error)

        # take the bins
        hq_bins = list(HQ_checkm_results["Name"])
        with open(str(log), 'w') as f:
            for bin_id in hq_bins:
                if '.' not in bin_id:
                    raise Exception(f"Invalid bin_id {bin_id}")
                binner_name      = bin_id.split('.')[0]
                source_file      = params.bin_folder + "/{}/bins/{}.fna".format(binner_name, bin_id)
                destination_file = output.HQ_folder  + "/{}.fna".format(bin_id)
                print("[rule collect_HQ_genomes] Hardlinking {} to {}".format(source_file, destination_file), file=f)
                subprocess.run(args=["ln", source_file, destination_file], stdout=f, stderr=f)

## Run coverm on HQ genomes to create the .tsv file
rule run_coverm:
    input:
        checkm_HQ = rules.collect_HQ_genomes.output.checkm_HQ,
        HQ_folder = rules.collect_HQ_genomes.output.HQ_folder
    output:
        cluster_tsv = "{wd}/{omics}/8-1-binning/mags/coverm_unique_cluster.tsv"
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/run_coverm.log"
    resources:
        mem=COVERM_memory
    threads:
        COVERM_THREADS
    conda:
        minto_dir + "/envs/MIntO_base.yml" # coverm
    shell:
        """
        time (
            export OMP_NUM_THREADS=1
            coverm cluster --genome-fasta-directory {input.HQ_folder} --checkm2-quality-report {input.checkm_HQ} -x fna --cluster-method fastani --ani 99 --fragment-length 2500 --min-aligned-fraction 30 --output-cluster-definition {output.cluster_tsv} --threads {threads} --precluster-method finch --precluster-ani 93
        ) &> {log}
        """

## Run retrieving scored
rule calculate_score_genomes:
    input:
        cluster_tsv = rules.run_coverm.output.cluster_tsv,
        checkm_HQ = rules.collect_HQ_genomes.output.checkm_HQ,
        HQ_folder = rules.collect_HQ_genomes.output.HQ_folder
    output:
        genome_scores = "{wd}/{omics}/8-1-binning/mags/HQ_genomes_checkm_scored.tsv"
    params:
        score_method = SCORE_METHOD
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/calculate_score_genomes.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    conda:
        minto_dir + "/envs/mags.yml"
    shell:
        """
        time (
            python {script_dir}/calculate_genomes_score.py --checkm_output {input.checkm_HQ} --fasta_folder {input.HQ_folder} --output_file {output.genome_scores} --score_method {params.score_method}
        ) &> {log}
        """


## Run retrieved the best unique genomes
rule find_unique_and_best_genomes:
    input:
        genome_scores = rules.calculate_score_genomes.output.genome_scores,
        coverm_tsv    = rules.run_coverm.output.cluster_tsv
    output:
        genome_report = "{wd}/{omics}/8-1-binning/mags/best_unique_genomes.report.tsv",
        genome_list  = "{wd}/{omics}/8-1-binning/mags/best_unique_genomes.list"
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/find_unique_and_best_genomes.log"
    resources:
        mem=10
    threads:
        2 # Decide number of threads
    run:
        import pandas as pd

        # read the table for the score
        score_table = pd.read_csv(str(input.genome_scores), sep = "\t", index_col = "Bin_id", comment = "#") # we skip the first line with the --score_method

        # read coverm table
        coverm_table = pd.read_csv(str(input.coverm_tsv), sep = "\t", names = ["ref_cluster", "cluster_members"])


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
        best_genomes_scored.to_csv(output.genome_report, sep = "\t", index = True)

        # Create the file
        best_unique_genomes_list= list(best_genomes_scored.index)

        with open(output.genome_list, "w") as fh:
            for genome in best_unique_genomes_list:
                fh.write("{}\n".format(genome))

## Run copy the best genomes
checkpoint copy_best_genomes:
    localrule: True
    input:
        genome_list  = rules.find_unique_and_best_genomes.output.genome_list,
        genome_scores = rules.calculate_score_genomes.output.genome_scores,
    output:
        genome_dir = directory("{wd}/{omics}/8-1-binning/mags/unique_genomes"),
    log:
        "{wd}/logs/{omics}/8-1-binning/mags/copy_best_genomes.log"
    resources:
        mem=10
    threads:
        1 # Decide number of threads
    shell:
        """
        time (
            mkdir -p {output.genome_dir}
            while read line; do
                ln {wildcards.wd}/{wildcards.omics}/8-1-binning/mags/HQ_genomes/${{line}}.fna {output.genome_dir}/ ;
            done < {input.genome_list}
        ) &> {log}
        """
