#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step - quality reads filtering

Authors: Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

# args = sys.argv
# config_path = args[args.index("--configfile") + 1]
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ")#, config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file
if config['PROJECT'] is None:
    print('ERROR in ', config_path, ': PROJECT variable is empty. Please, complete ', config_path)
else:
    project_id = config['PROJECT']

if config['working_dir'] is None:
    print('ERROR in ', config_path, ': working_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['working_dir']) is False:
    print('ERROR in ', config_path, ': working_dir variable path does not exit. Please, complete ', config_path)
else:
    working_dir = config['working_dir']

if config['omics'] in ('metaG', 'metaT'):
    omics = config['omics']
else:
    print('ERROR in ', config_path, ': omics variable is not correct. "omics" variable should be metaG or metaT.')

if config['local_dir'] is None:
    print('ERROR in ', config_path, ': local_dir variable is empty. Please, complete ', config_path)
else:
    local_dir = config['local_dir']

if config['raw_reads_dir'] is None:
    print('ERROR in ', config_path, ': raw_reads_dir variable in configuration yaml file is empty. Please, complete ', config_path)
else:
    raw_dir = config['raw_reads_dir']

if config['minto_dir'] is None:
    print('ERROR in ', config_path, ': minto_dir variable in configuration yaml file is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]
    script_dir=config["minto_dir"]+"/scripts"

if config['METADATA'] is None:
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata=config["METADATA"]
elif config['METADATA'] == "None":
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata=config["METADATA"]
elif path.exists(config['METADATA']) is False:
    print('ERROR in ', config_path, ': METADATA variable path does not exit. Please, complete ', config_path)
else:
    metadata=config["METADATA"]

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
    else:
        try:
            # Make list of illumina samples, if ILLUMINA in config
            ilmn_samples = list()
            if 'ILLUMINA' in config:
                #print("Samples:")
                for ilmn in config["ILLUMINA"]:
                    if path.exists(raw_dir+'/'+ilmn+'/'+ilmn+'.1.fq.gz') is True:
                        #print(ilmn)
                        ilmn_samples.append(ilmn)
                    else:
                        raise TypeError('ERROR in ', config_path, ': ILLUMINA list of samples does not exist. Please, complete ', config_path)
        except: 
            print('ERROR in ', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete ', config_path)
else:
    print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)

if config['trimmomatic_adaptors'] in ('False', 'Skip'):
    pass
elif config['trimmomatic_adaptors'] is None:
    print('ERROR in ', config_path, ': trimmomatic_adaptors variable is empty. Please, complete ', config_path)
elif path.exists(config['trimmomatic_adaptors']) is False:
    print('ERROR in ', config_path, ': trimmomatic_adaptors variable path does not exit. Please, complete ', config_path)

if config['perc_remaining_reads'] is None:
    print('ERROR in ', config_path, ': perc_remaining_reads variable is empty. Please, complete ', config_path)

if config['TRIMMOMATIC_threads'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_threads']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is not an integer. Please, complete ', config_path)

if config['TRIMMOMATIC_memory'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_memory']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is not an integer. Please, complete ', config_path)


# Define all the outputs needed by target 'all'
def qc1_trim_quality_output():
    result = expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.paired.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group=['1', '2']),\
    expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.single.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group=['1', '2']),\
    expand("{wd}/{omics}/1-trimmed/{sample}/{sample}_trimlog_quality_adapter", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [])
    return(result)

def qc1_check_read_length_output():
    result = expand("{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group = ['1', '2']),\
    expand("{wd}/{omics}/1-trimmed/samples_read_length.txt",
                    wd = working_dir,
                    omics = omics)
    return(result)

def qc1_read_length_cutoff_output():
    result = expand("{wd}/output/1-trimmed/{omics}_cumulative_read_lenght_cutoff.pdf", 
                    wd = working_dir,
                    omics = omics)
    return(result)

def qc1_config_yml_output():
    result = expand("{wd}/{omics}/QC_2.yaml", 
                    wd = working_dir,
                    omics = omics)
    return(result)

rule all:
    input:
        qc1_trim_quality_output(),
        qc1_check_read_length_output(),
        qc1_read_length_cutoff_output(),
        qc1_config_yml_output()


###############################################################################################
# Pre-processing of metaG and metaT data 
# Filter reads by quality
###############################################################################################    
rule qc1_trim_quality:
    input:
        read_fw=lambda wildcards: '{raw_dir}/{sample}/{sample}.1.fq.gz'.format(raw_dir=raw_dir, sample = wildcards.sample),
        read_rv=lambda wildcards: '{raw_dir}/{sample}/{sample}.2.fq.gz'.format(raw_dir=raw_dir, sample = wildcards.sample),
    output: 
        pairead1="{wd}/{omics}/1-trimmed/{sample}/{sample}.1.paired.fq.gz", 
        singleread1="{wd}/{omics}/1-trimmed/{sample}/{sample}.1.single.fq.gz", 
        pairead2="{wd}/{omics}/1-trimmed/{sample}/{sample}.2.paired.fq.gz", 
        singleread2="{wd}/{omics}/1-trimmed/{sample}/{sample}.2.single.fq.gz", 
        log_adapter="{wd}/{omics}/1-trimmed/{sample}/{sample}_trimlog_quality_adapter"
    params:
        tmp_quality=lambda wildcards: "{local_dir}/{omics}_{sample}_qc1_trim_quality/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
        adapters=config['trimmomatic_adaptors'],
    log: 
        "{wd}/logs/{omics}/1-trimmed/{sample}_qc1_trim_quality.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads']      
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml"#trimmomatic
    shell: 
        """mkdir -p {params.tmp_quality}1-trimmed/
        remote_dir=$(dirname {output.pairead1})
        adapter_file="{params.adapters}"
time (if [ "{params.adapters}" == "Skip" ]
    then echo ${{adapter_file}}; \
cat {input.read_fw} > {output.pairead1}; \
> {output.singleread1}; \
cat {input.read_rv} > {output.pairead2}; \
> {output.singleread2}; \
> {output.log_adapter}
elif [ "{params.adapters}" == "False" ]
    then echo 'Generate adapter sequences file'; \
sh {script_dir}/index_adapter.sh {input.read_fw} {wildcards.sample} {params.tmp_quality}1-trimmed/ {wildcards.wd} {wildcards.omics}; \
trimmomatic PE -threads {threads} -trimlog {params.tmp_quality}1-trimmed/{wildcards.sample}_trimlog_quality_adapter -phred33 {input.read_fw} {input.read_rv} \
{params.tmp_quality}1-trimmed/{wildcards.sample}.1.paired.fq.gz {params.tmp_quality}1-trimmed/{wildcards.sample}.1.single.fq.gz \
{params.tmp_quality}1-trimmed/{wildcards.sample}.2.paired.fq.gz {params.tmp_quality}1-trimmed/{wildcards.sample}.2.single.fq.gz \
TRAILING:5 LEADING:5 SLIDINGWINDOW:4:20 ILLUMINACLIP:{params.tmp_quality}1-trimmed/adapters.fa:2:30:10:2:keepBothReads; \
rsync {params.tmp_quality}1-trimmed/* $remote_dir
else echo 'Filter by quality'
    trimmomatic PE -threads {threads} -trimlog {params.tmp_quality}1-trimmed/{wildcards.sample}_trimlog_quality_adapter -phred33 {input.read_fw} {input.read_rv} \
{params.tmp_quality}1-trimmed/{wildcards.sample}.1.paired.fq.gz {params.tmp_quality}1-trimmed/{wildcards.sample}.1.single.fq.gz \
{params.tmp_quality}1-trimmed/{wildcards.sample}.2.paired.fq.gz {params.tmp_quality}1-trimmed/{wildcards.sample}.2.single.fq.gz \
TRAILING:5 LEADING:5 SLIDINGWINDOW:4:20 ILLUMINACLIP:${{adapter_file}}:2:30:10:2:keepBothReads; \
rsync {params.tmp_quality}1-trimmed/* $remote_dir
        fi ) >& {log}
        rm -rf {params.tmp_quality}"""
 
rule qc1_check_read_length:
    input: 
        pairead="{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.paired.fq.gz"
    output: 
        length="{wd}/{omics}/1-trimmed/{sample}/{sample}.{group}.read_length.txt"
    params:
        tmp_read_length=lambda wildcards: "{local_dir}/{omics}_{sample}.{group}_qc1_check_read_length/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample, group = wildcards.group),
    log: 
        "{wd}/logs/{omics}/1-trimmed/{sample}.{group}.check_read_length.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads'] 
    shell: 
        """ mkdir -p {params.tmp_read_length}
        time (sh {script_dir}/QC_trim_read_length.sh {input.pairead} {params.tmp_read_length}/{wildcards.sample}.{wildcards.group}.read_length.txt {wildcards.sample}_{wildcards.group}
        rsync {params.tmp_read_length}/{wildcards.sample}.{wildcards.group}.read_length.txt {output.length}) >& {log}
        rm -rf {params.tmp_read_length} """

rule qc1_check_read_length_merge:
    input: 
        length=expand("{{wd}}/{{omics}}/1-trimmed/{sample}/{sample}.{group}.read_length.txt", sample=ilmn_samples, group=['1', '2'])
    output: 
        list_samples="{wd}/{omics}/1-trimmed/samples_read_length.txt", 
    log: 
        "{wd}/logs/{omics}/1-trimmed/qc1_check_read_length_merge.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads'] 
    shell: 
        """time (cat {input.length} >> {output.list_samples}) >& {log}"""  
    
rule qc1_cumulative_read_len_plot:
    input:
        list_samples="{wd}/{omics}/1-trimmed/samples_read_length.txt", 
    output:
        plot="{wd}/output/1-trimmed/{omics}_cumulative_read_lenght_cutoff.pdf", 
        cutoff_file="{wd}/{omics}/1-trimmed/QC_1_min_len_read_cutoff.txt"
    params: 
        cutoff=config["perc_remaining_reads"], 
    log: 
        "{wd}/logs/{omics}/1-trimmed/plot_cumulative_read_len.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads'] 
    conda: 
        config["minto_dir"]+"/envs/taxa_env.yml"
    shell:
        """ time (Rscript {script_dir}/QC_cumulative_read_length_plot.R {input.list_samples} {params.cutoff} {output.plot} {output.cutoff_file}) >& {log}""" #; rm {input.list_length}
            
##########################################################################################################
# Generate configuration yml file for the next step in pre-processing of metaG and metaT data 
##########################################################################################################

rule qc1_config_yml_file:
    input: 
        cutoff_file="{wd}/{omics}/1-trimmed/QC_1_min_len_read_cutoff.txt",
    output: 
        config_file="{wd}/{omics}/QC_2.yaml"
    params: 
        tmp_qc1_yaml=lambda wildcards: "{local_dir}/{omics}_qc1_config_yml_file/".format(local_dir=local_dir, omics = omics),
        trim_threads=config['TRIMMOMATIC_memory'],
        trim_memory=config['TRIMMOMATIC_threads'] 
    resources:
        mem=2
    threads: 
        2
    log: 
        "{wd}/logs/{omics}/qc1_config_yml_file.log"
    shell: 
        """
        mkdir -p {params.tmp_qc1_yaml}
        time (files='{ilmn_samples}'; echo $files; echo $files| tr ' ' '\\n'| sed 's/^/- /' - > {params.tmp_qc1_yaml}/samples_illumina.txt
echo "######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
local_dir: {local_dir}
minto_dir: {minto_dir}
METADATA: {metadata}

######################
# Program settings
######################
# Pre-processing - length reads filtering
TRIMMOMATIC_threads: {params.trim_threads}
TRIMMOMATIC_memory: {params.trim_memory}" > {params.tmp_qc1_yaml}QC_2.yaml
cat {input.cutoff_file} >> {params.tmp_qc1_yaml}QC_2.yaml

echo "
# Pre-processing - host genome filtering
PATH_host_genome:
NAME_host_genome:
BWA_index_host_threads:
BWA_index_host_memory:
BWA_host_threads:
BWA_host_memory: " >> {params.tmp_qc1_yaml}QC_2.yaml

if [ "{omics}" == "metaG" ];\
    then echo "
# Assembly-free taxonomy profiling
TAXA_threads:
TAXA_memory:
taxa_profile: metaphlan">> {params.tmp_qc1_yaml}QC_2.yaml
elif [ "{omics}" == "metaT" ];\
    then echo "
# ribosomal RNA depletion
sortmeRNA_threads:
sortmeRNA_memory:
sortmeRNA_db:
sortmeRNA_db_idx: ">> {params.tmp_qc1_yaml}QC_2.yaml
fi
echo "
# Input data
# ILLUMINA section:
# -----------------
# List of illumina samples that will be prerpocessed.
#
# E.g.:
# - I1
# - I2
#
ILLUMINA:" >> {params.tmp_qc1_yaml}QC_2.yaml
cat {params.tmp_qc1_yaml}/samples_illumina.txt >> {params.tmp_qc1_yaml}QC_2.yaml
rsync {params.tmp_qc1_yaml}QC_2.yaml {output.config_file}) >& {log}
rm -rf {params.tmp_qc1_yaml}
        """ 
        
        