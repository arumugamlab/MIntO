#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step - quality reads filtering

Authors: Carmen Saenz, Mani Arumugam, Judit Szarvas
'''

#
# configuration yaml file
# import sys
import re
from os import path, scandir
import pandas as pd
from glob import glob

localrules: qc0_fake_move_for_multiqc, qc0_create_multiqc, \
            qc1_check_read_length_merge, qc1_cumulative_read_len_plot, \
            qc2_config_yml_file

# Get common config variables
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule QC_0_base, QC_0_rpkg from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

if config['raw_reads_dir'] is None:
    print('ERROR in ', config_path, ': raw_reads_dir variable in configuration yaml file is empty. Please, complete ', config_path)
else:
    raw_dir = config['raw_reads_dir']

if config['FASTP_threads'] is None:
    print('ERROR in ', config_path, ': FASTP_threads variable is empty. Please, complete ', config_path)
elif type(config['FASTP_threads']) != int:
    print('ERROR in ', config_path, ': FASTP_threads variable is not an integer. Please, complete ', config_path)

if config['FASTP_memory'] is None:
    print('ERROR in ', config_path, ': FASTP_memory variable is empty. Please, complete ', config_path)
elif type(config['FASTP_memory']) != int:
    print('ERROR in ', config_path, ': FASTP_memory variable is not an integer. Please, complete ', config_path)

if config['FASTP_adapters'] in ('Skip', 'Quality', 'Overlap', 'Detect'):
    pass
elif config['FASTP_adapters'] is None:
    print('ERROR in ', config_path, ': FASTP_adapters variable is empty. Please, complete ', config_path)
elif path.exists(config['FASTP_adapters']) is False:
    print('ERROR in ', config_path, ': FASTP_adapters variable path does not exit. Please, complete ', config_path)

# file suffixes
ilmn_suffix = ["1.fq.gz", "2.fq.gz"]
if 'ILLUMINA_suffix' in config and config["ILLUMINA_suffix"] is not None:
    ilmn_suffix = config["ILLUMINA_suffix"]

# Make list of illumina samples, if ILLUMINA in config
ilmn_samples = list()
ilmn_samples_organisation = "folder"
ilmn_runs_df = None
if 'ILLUMINA' in config:
    # column_name specified option
    if isinstance(config["ILLUMINA"], str):
        ilmn_samples_organisation = "bulk"
        if path.exists(config["ILLUMINA"]) and path.isfile(config["ILLUMINA"]):
            # extra runs sheet
            print('MIntO uses', config["ILLUMINA"], 'as sample list')
            col_name = "sample"
            ilmn_runs_df = pd.read_table(config["ILLUMINA"])
            for sampleid in ilmn_runs_df['sample'].unique():
                for runid in ilmn_runs_df.loc[ilmn_runs_df['sample'] == sampleid]['run'].to_list():
                    #print(sampleid, runid)
                    sample_pattern = "{}/{}[.-_]{}".format(raw_dir, runid, ilmn_suffix[0])
                    if glob(sample_pattern):
                        if sampleid not in ilmn_samples:
                            ilmn_samples.append(sampleid)
                    else:
                        print('WARNING: sample ', sample_pattern, 'not in bulk data folder ', raw_dir)
        else:
            # column name in metadata sheet
            col_name = config["ILLUMINA"]
            md_df = pd.read_table(metadata)
            if not col_name in md_df.columns:
                raise Exception(f"ERROR in {config_path}: column name specified for ILLUMINA does not exist in metadata or runs sheet. Please, complete {config_path}")
            for sampleid in md_df[col_name].to_list():
                sample_pattern = "{}/{}[.-_]{}".format(raw_dir, sampleid, ilmn_suffix[0])
                if glob(sample_pattern) and sampleid not in ilmn_samples:
                    ilmn_samples.append(sampleid)
                else:
                    print('WARNING in ', metadata, ': sample ', sampleid, 'not in bulk data folder ', raw_dir)
    # listed samples
    else:
        for ilmn in config["ILLUMINA"]:
            location = "{}/{}".format(raw_dir, ilmn)
            if path.exists(location):
                ilmn_samples.append(ilmn)
            else:
                raise Exception(f"ERROR: {location} in raw_dir does not exist. Please, locate {ilmn}")

# trimming options
adapter_trimming_args = ""
if config['FASTP_adapters'] == 'Quality':
    adapter_trimming_args = "--disable_adapter_trimming"
elif config['FASTP_adapters'] == 'Overlap':
    adapter_trimming_args = "--overlap_diff_percent_limit 10"
elif config['FASTP_adapters'] == 'Detect':
    adapter_trimming_args = "--detect_adapter_for_pe"

# Define all the outputs needed by target 'all'

def qc1_read_length_cutoff_output():
    result = expand("{wd}/output/1-trimmed/{omics}_cumulative_read_length_cutoff.pdf",
                    wd = working_dir,
                    omics = omics)
    return(result)

def qc1_config_yml_output():
    result = expand("{wd}/{omics}/QC_2.yaml",
                    wd = working_dir,
                    omics = omics)
    return(result)

def qc0_multiqc_summary_output():
    result = expand("{wd}/output/1-trimmed/{omics}_multiqc.html",
                    wd = working_dir,
                    omics = omics)
    return(result)

rule all:
    input:
        qc1_read_length_cutoff_output(),
        qc1_config_yml_output(),
        qc0_multiqc_summary_output(),
        print_versions.get_version_output(snakefile_name)
    default_target: True


###############################################################################################
# Pre-processing of metaG and metaT data
# Filter reads by quality
###############################################################################################

# Get a sorted list of runs for a sample

def get_runs_for_sample(sample):
    runs = [sample]
    if ilmn_samples_organisation == "folder":
        sample_dir = '{raw_dir}/{sample}'.format(raw_dir=raw_dir, sample=sample)
        runs = [ f.name.split(ilmn_suffix[0])[0][:-1] for f in scandir(sample_dir) if f.is_file() and f.name.endswith(ilmn_suffix[0]) ]
    elif ilmn_runs_df is not None:
        runs = ilmn_runs_df.loc[ilmn_runs_df['sample'] == sample]['run'].to_list()
    return(sorted(runs))

# Get files for a run

def get_raw_reads_for_sample_run(wildcards):
    prefix = '{raw_dir}/{run}'.format(raw_dir=raw_dir, run=wildcards.run)
    if ilmn_samples_organisation == "folder":
        prefix = '{raw_dir}/{sample}/{run}'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run)
    raw_sample_run = {}
    for i, k in enumerate(['read_fw', 'read_rv']):
        raw_sample_run[k] = glob("{}[.-_]{}".format(prefix, ilmn_suffix[i]))[0]
    return raw_sample_run


##########
# FastQC of the initial raw data
##########

rule initial_fastqc:
    input:
        unpack(get_raw_reads_for_sample_run)
    output:
        flag="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastqc.done",
    params:
        outdir="{wd}/{omics}/1-0-qc/{sample}"
    log:
        "{wd}/logs/{omics}/1-0-qc/{sample}_{run}_qc0_fastqc.log"
    resources:
        mem=4
    threads:
        2
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        time (
            fastqc --noextract -f fastq -t {threads} -q -o {params.outdir} {input.read_fw} {input.read_rv} && \
        echo "Fastqc done" > {output.flag} 
        ) >& {log}
        """

##########
# Adaptor-trimming + Quality-trimming or Quality-trimming only?
##########

if config['FASTP_adapters'] == 'Skip':

    # Fake a trim
    rule qc0_trim_quality:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
            json="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastp.json"
        localrule: True
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc0_trim_quality.log"
        shell:
            """
            (ln -s --force {input.read_fw} {output.pairead1}
            ln -s --force {input.read_rv} {output.pairead2}
            touch {output.json}
            echo "Skipped trimming and linked raw files to trimmed files.") >& {log}
            """
elif adapter_trimming_args:
    rule qc0_trim_quality:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
            json="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastp.json",
            html="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastp.html"
        shadow:
            "minimal"
        params:
            mq_5=config['FASTP_front_mean_qual'],
            mq_3=config['FASTP_tail_mean_qual'],
            ml=config['FASTP_min_length'],
            adapter_args=adapter_trimming_args
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc0_trim_quality.log"
        resources:
            mem=config['FASTP_memory']
        threads:
            config['FASTP_threads']
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml"
        shell:
            """
            time ( \
                fastp -i {input.read_fw} --in2 {input.read_rv} \
                -o {output.pairead1} --out2 {output.pairead2} \
                --cut_window_size 4 --cut_front --cut_front_mean_quality {params.mq_5} --cut_tail --cut_tail_mean_quality {params.mq_3} --length_required {params.ml} \
                {params.adapter_args} \
                --dont_eval_duplication \
                --thread {threads} --json {output.json} --html {output.html}
            ) >& {log}
            """
else:
    rule qc0_trim_quality_and_custom_adapter:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
            json="{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastp.json"
        shadow:
            "minimal"
        params:
            mq_5=config['FASTP_front_mean_qual'],
            mq_3=config['FASTP_tail_mean_qual'],
            ml=config['FASTP_min_length'],
            adapter=config['FASTP_adapters']
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc0_trim_customadapter.log"
        resources:
            mem=config['FASTP_memory']
        threads:
            config['FASTP_threads']
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml"
        shell:
            """
            time ( \
                fastp -i {input.read_fw} --in2 {input.read_rv} \
                -o {output.pairead1} --out2 {output.pairead2} \
                --cut_window_size 4 --cut_front --cut_front_mean_quality {params.mq_5} --cut_tail --cut_tail_mean_quality {params.mq_3} --length_required {params.ml} \
                --dont_eval_duplication \
                --adapter_fasta {params.adapter} \
                --thread {threads} --json {output.json}
            ) >& {log}
            """

########
# Multiqc aggregation and read lengths
########

rule qc0_fake_move_for_multiqc:
    input:
        fastqc_html=lambda wildcards: expand("{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastqc.done",
                                            wd=wildcards.wd,
                                            omics=wildcards.omics,
                                            sample=wildcards.sample,
                                            run=get_runs_for_sample(wildcards.sample)
                                            ),
        fastp_json=lambda wildcards: expand("{wd}/{omics}/1-0-qc/{sample}/{sample}-{run}_fastp.json",
                                            wd=wildcards.wd,
                                            omics=wildcards.omics,
                                            sample=wildcards.sample,
                                            run=get_runs_for_sample(wildcards.sample)
                                            )
    output:
        temp("{wd}/{omics}/1-0-qc/{sample}_fake.moved")
    threads:
        2
    shell:
        """
        echo "{wildcards.sample} moved" > {output}
        """

rule qc0_create_multiqc:
    input:
        mqc_flag=lambda wildcards: expand("{wd}/{omics}/1-0-qc/{sample}_fake.moved",
                                            wd=wildcards.wd,
                                            omics=wildcards.omics,
                                            sample=ilmn_samples
                                            )
    output:
        "{wd}/output/1-trimmed/{omics}_multiqc.html"
    params:
        indir="{wd}/{omics}/1-0-qc",
        outdir="{wd}/output/1-trimmed"
    log:
        "{wd}/logs/{omics}/1-trimmed/{omics}_multiqc.log"
    threads:
        2
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        time (
            multiqc --filename {wildcards.omics}_multiqc.html --outdir {params.outdir} -d --zip-data-dir -q --no-ansi --interactive {params.indir}
        ) >& {log}
        """

rule qc1_check_read_length:
    input:
        pairead=lambda wildcards: expand("{wd}/{omics}/1-trimmed/{sample}/{run}.{group}.paired.fq.gz",
                                            wd=wildcards.wd,
                                            omics=wildcards.omics,
                                            sample=wildcards.sample,
                                            run=get_runs_for_sample(wildcards.sample),
                                            group=wildcards.group
                                            )
    output:
        length="{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt"
    log:
        "{wd}/logs/{omics}/1-trimmed/{sample}.{group}.check_read_length.log"
    resources:
        mem=4
    threads:
        2
    shell:
        """
        time (
            sh -c 'gzip -cd {input.pairead} | awk -v f="{wildcards.sample}_{wildcards.group}" "{{if(NR%4==2) print length(\$1),f}}" | sort -n | uniq -c > {output.length}'
        ) >& {log}
        """

rule qc1_check_read_length_merge:
    input:
        length=lambda wildcards: expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt", wd=wildcards.wd, omics=wildcards.omics, sample=ilmn_samples, group=['1', '2'])
    output:
        readlen_dist="{wd}/{omics}/1-trimmed/samples_read_length.txt",
    log:
        "{wd}/logs/{omics}/1-trimmed/qc1_check_read_length_merge.log"
    threads: 1
    shell:
        """
        time (
            cat {input.length} > {output.readlen_dist}
        ) >& {log}
        """

rule qc1_cumulative_read_len_plot:
    input:
        readlen_dist=rules.qc1_check_read_length_merge.output.readlen_dist
    output:
        plot="{wd}/output/1-trimmed/{omics}_cumulative_read_length_cutoff.pdf",
        cutoff_file="{wd}/{omics}/1-trimmed/QC_1_min_len_read_cutoff.txt"
    params:
        cutoff=config["perc_remaining_reads"],
    log:
        "{wd}/logs/{omics}/1-trimmed/plot_cumulative_read_len.log"
    threads: 1
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (
            Rscript {script_dir}/QC_cumulative_read_length_plot.R --input {input.readlen_dist} --frac {params.cutoff} --out_plot {output.plot} --out_cutoff {output.cutoff_file}
        ) >& {log}
        """

##########################################################################################################
# Generate configuration yml file for the next step in pre-processing of metaG and metaT data
##########################################################################################################

rule qc2_config_yml_file:
    input:
        cutoff_file=rules.qc1_cumulative_read_len_plot.output.cutoff_file
    output:
        config_file="{wd}/{omics}/QC_2.yaml"
    params:
        trim_threads=config['FASTP_threads'],
        trim_memory=config['FASTP_memory']
    log:
        "{wd}/logs/{omics}/qc2_config_yml_file.log"
    shell:
        """
        cat > {output.config_file} <<___EOF___
######################
# General settings
######################

PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
minto_dir: {minto_dir}
METADATA: {metadata}

########################################
# Program settings
########################################

#########################
# Read length filtering
#########################

READ_minlen: $(cat {input.cutoff_file})

#########################
# Host genome filtering
#########################

# bwa-mem2 index files will be stored at: <PATH_host_genome>/BWA_index/<NAME_host_genome>.*
# If it already exists, then it will be used directly.
# If not, a fasta file should exist as: <PATH_host_genome>/<NAME_host_genome>
# This will be build into index files using:
#    bwa-mem2 index -p <PATH_host_genome>/BWA_index/<NAME_host_genome> <PATH_host_genome>/<NAME_host_genome>

PATH_host_genome:
NAME_host_genome:
BWA_index_host_memory: 100
BWA_host_threads: 8
BWA_host_memory: 40
___EOF___

        if [ "{omics}" == "metaT" ]; then
            cat >> {output.config_file} <<___EOF___

#########################
# Ribosomal RNA depletion
#########################

sortmeRNA_threads: 8
sortmeRNA_memory: 10
sortmeRNA_db: {minto_dir}/data/rRNA_databases
sortmeRNA_db_idx: {minto_dir}/data/rRNA_databases/idx
___EOF___

        fi

        cat >> {output.config_file} <<___EOF___

##################################
# Assembly-free taxonomy profiling
##################################

# Following values for 'TAXA_profiler' are supported:
#    1. metaphlan - relative abundance using MetaPhlAn
#    2. motus_raw - read counts using mOTUs
#    3. motus_rel - relative abundance using mOTUs
# Comma-delimited combination of multiple options also supported
# Eg:
#    TAXA_profiler: metaphlan,motus_rel
TAXA_threads: 8
TAXA_memory: 10
TAXA_profiler: motus_rel,metaphlan
metaphlan_version: 4.0.6
motus_version: 3.0.3

#########################
# K-mer based comparison
#########################

# FracMinHash comparisons by sourmash
# SOURMASH_min_abund - Minimum count of each k-mer for filtering the sketch (integer)
# SOURMASH_max_abund - Maximum count of each k-mer for filtering the sketch (integer)
# SOURMASH_cutoff    - Dissimilarity cutoff for subclusters via hierarchical clustering

SOURMASH_min_abund: 2
SOURMASH_max_abund: 1000
SOURMASH_cutoff: 0.40

#####################
# Analysis parameters
#####################

# MAIN_factor  - the main factor in the metadata file to differentiate in visualization (using color)
# PLOT_factor2 - the second factor in the metadata file to differentiate in visualization (using shape)
# PLOT_time    - name of the factor in the metadata file denoting time (e.g. hour, day)

MAIN_factor:
PLOT_factor2:
PLOT_time:

######################
# Optionally, do you want to merge replicates or make pseudo samples
# E.g:
# MERGE_ILLUMINA_SAMPLES:
#  - merged=rep1+rep2+rep3
#
# The above directive will combine 3 samples (rep1, rep2 and rep3)
# after the last step into a new sample called 'merged'. Now you can remove
# rep1, rep2 and rep3 from assembly, MAG generation and profiling steps.
# Please note that METADATA file must have an entry for 'merged' as well,
# otherwise QC_2 step will fail.
# Having extra entries in METADATA file does not affect you in any way.
######################

#MERGE_ILLUMINA_SAMPLES:


######################
# Input data
######################

# ILLUMINA section:
# -----------------
# List of illumina samples that will be filtered by read length.
#
# E.g.:
# - I1
# - I2
#
ILLUMINA:
$(for i in {ilmn_samples}; do echo "- '$i'"; done)
___EOF___

        echo {ilmn_samples} >& {log}
        """
