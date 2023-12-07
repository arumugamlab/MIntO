#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step - quality reads filtering

Authors: Carmen Saenz, Mani Arumugam
'''

#
# configuration yaml file
# import sys
import re
from os import path

localrules: qc1_check_read_length_merge, qc1_cumulative_read_len_plot, qc2_config_yml_file

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

if config['raw_reads_dir'] is None:
    print('ERROR in ', config_path, ': raw_reads_dir variable in configuration yaml file is empty. Please, complete ', config_path)
else:
    raw_dir = config['raw_reads_dir']

multiplex_tech = 'TruSeq'
if config['MULTIPLEX_TECH'] is None:
    print('ERROR in ', config_path, ': MULTIPLEX_TECH variable in configuration yaml file is empty. Please, complete ', config_path)
else:
    multiplex_tech = config['MULTIPLEX_TECH']

if config['TRIMMOMATIC_threads'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_threads']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is not an integer. Please, complete ', config_path)

if config['TRIMMOMATIC_memory'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_memory']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is not an integer. Please, complete ', config_path)

if config['TRIMMOMATIC_adaptors'] in ('Skip', 'Quality'):
    pass
elif config['TRIMMOMATIC_adaptors'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_adaptors variable is empty. Please, complete ', config_path)
elif path.exists(config['TRIMMOMATIC_adaptors']) is False:
    print('ERROR in ', config_path, ': TRIMMOMATIC_adaptors variable path does not exit. Please, complete ', config_path)

if 'TRIMMOMATIC_palindrome' in config and config['TRIMMOMATIC_palindrome'] is not None:
    if path.exists(config['TRIMMOMATIC_palindrome']) is False:
        print('ERROR in ', config_path, ': TRIMMOMATIC_palindrome variable path does not exit. Please, complete ', config_path)

trimmomatic_simple_clip_threshold = 10
if 'TRIMMOMATIC_simple_clip_threshold' in config:
    trimmomatic_simple_clip_threshold = config['TRIMMOMATIC_simple_clip_threshold']

if config['perc_remaining_reads'] is None:
    print('ERROR in ', config_path, ': perc_remaining_reads variable is empty. Please, complete ', config_path)

# Make list of illumina samples, if ILLUMINA in config
try:
    # Make list of illumina samples, if ILLUMINA in config
    ilmn_samples = list()
    if 'ILLUMINA' in config:
        #print("Samples:")
        for ilmn in config["ILLUMINA"]:
            location = "{}/{}".format(raw_dir, ilmn)
            if path.exists(location) is True:
                #print(fwd)
                ilmn_samples.append(ilmn)
            else:
                raise TypeError('ERROR in', config_path, ':', 'sample', ilmn, 'in ILLUMINA list does not exist. Please, complete', config_path)
except:
    print('ERROR in', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete', config_path)


# Define all the outputs needed by target 'all'

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
    result = expand("{wd}/output/1-trimmed/{omics}_cumulative_read_length_cutoff.pdf",
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
        #qc1_trim_quality_output(),
        #qc1_check_read_length_output(),
        qc1_read_length_cutoff_output(),
        qc1_config_yml_output()


###############################################################################################
# Pre-processing of metaG and metaT data
# Filter reads by quality
###############################################################################################

# Get a sorted list of runs for a sample

def get_runs_for_sample(sample):
    sample_dir = '{raw_dir}/{sample}'.format(raw_dir=raw_dir, sample=sample)
    runs = [ re.sub(".1\.(fq|fastq)\.gz", "", path.basename(f)) for f in os.scandir(sample_dir) if f.is_file() and (f.name.endswith('1.fq.gz') or f.name.endswith('1.fastq.gz')) ]
    #print(runs)
    return(sorted(runs))

# Get files for a run

def get_raw_reads_for_sample_run(wildcards):
    prefix = '{raw_dir}/{sample}/{run}'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run)
    for ext in ['fq', 'fastq']:
        if path.exists(prefix+'.1.'+ext+'.gz'):
            return {
                    'read_fw': '{raw_dir}/{sample}/{run}.1.{ext}.gz'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run, ext=ext),
                    'read_rv': '{raw_dir}/{sample}/{run}.2.{ext}.gz'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run, ext=ext)
                    }
        elif path.exists(prefix+'_1.'+ext+'.gz'):
            return {
                    'read_fw': '{raw_dir}/{sample}/{run}_1.{ext}.gz'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run, ext=ext),
                    'read_rv': '{raw_dir}/{sample}/{run}_2.{ext}.gz'.format(raw_dir=raw_dir, sample=wildcards.sample, run=wildcards.run, ext=ext)
                    }

# Get file for infering index

def get_example_to_infer_index(sample):
    sample_dir = '{raw_dir}/{sample}'.format(raw_dir=raw_dir, sample=sample)
    runs = [ os.path.normpath(f) for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith('1.fq.gz') ]
    #print(runs)
    return(runs[0])

##########
# Index barcodes should be used for adaptor trimming
##########

if 'TRIMMOMATIC_index_barcodes' in config and config['TRIMMOMATIC_index_barcodes'] != None:

    # Create a file with fwd/rev index sequences in the format
    #   GTTACGGA+CTTGGCTA
    # If there is a global index file for all samples given in config['TRIMMOMATIC_index_barcodes'], use it
    # If it does not exist, then infer from fastq file

    localrules: qc1_get_index_barcode, qc1_make_custom_adapter_file, qc1_make_custom_adapter_file_with_palindrome
    ruleorder: qc1_make_custom_adapter_file_with_palindrome > qc1_make_custom_adapter_file

    if config['TRIMMOMATIC_index_barcodes'].lower() == "infer":
        # Infer from fastq file
        rule qc1_get_index_barcode:
            input:
                read_fw=lambda wildcards: get_example_to_infer_index(wildcards.sample),
            output:
                barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
            shell:
                """
                outdir=$(dirname {output.barcodes})
                mkdir -p $outdir
                # Get the most frequent barcode-pair
                # I am invoking a shell within, otherwise snakemake's strict mode catches SIG_PIPE from head to gzip and fails
                barcodes=$(sh -c 'gzip -cd {input.read_fw} | grep "^@" | head -10000 | cut -f2 -d" " | cut -f4 -d":" | sort | uniq -c | sed "s/^ \+//" | tr -s " " " " | sort -k1,1nr -t" " | head -1 | cut -f2 -d" "')
                # Write barcodes for this sample
                echo $barcodes > {output.barcodes}
                """

    else:
        # From global index file
        rule qc1_get_index_barcode:
            input:
                barcodes=config['TRIMMOMATIC_index_barcodes'],
            output:
                barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
            shell:
                """
                outdir=$(dirname {output.barcodes})
                mkdir -p $outdir
                # Get barcodes for this sample
                barcodes=$(grep "^"{wildcards.sample}"[[:space:]]" {input.barcodes} | cut -f2)
                # Write barcodes for this sample
                echo $barcodes > {output.barcodes}
                """

    # Create a custom adapter file to be used in Trimmomatic, given the custom index file above
    rule qc1_make_custom_adapter_file_with_palindrome:
        input:
            barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
            palindrome=config['TRIMMOMATIC_palindrome'],
            template=config['TRIMMOMATIC_adaptors']
        output:
            adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
        params:
            multiplex=multiplex_tech
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml" # seqtk
        shell:
            """
            # Get the 1st part of the adapter pair
            barcode1=$(cat {input.barcodes} | sed "s/;/+/" | cut -f1 -d'+')
            # Get the 2nd part of adapter pair
            barcode2=$(cat {input.barcodes} | sed "s/;/+/" | cut -f2 -d'+')

            # For MGI, revcomp barcode2
            if [ "{params.multiplex}" == "MGIEasy" ]; then
                barcode2=$(echo $barcode2 | tr 'ATGC' 'TACG' | rev)
            fi

            # Make custom adapters for this sample using its index sequences
            #  1. Make palindromes
            cat {input.template} | seqkit grep --quiet -n -f {input.palindrome} -w 1000 -o - | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | sed "s^>.*/^>PrefixPE-Ad/^" | tr -s '12' '21' > {output.adapter}
            #  2. Make 5-prime adapters
            cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" >> {output.adapter}
            #  3. Make 3-prime adapters
            cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | tr -s '12' '21' | sed "s^/^_rc/^" >> {output.adapter}
            """

    # Create a custom adapter file to be used in Trimmomatic, given the custom index file above
    rule qc1_make_custom_adapter_file:
        input:
            barcodes="{wd}/{omics}/1-trimmed/{sample}/index_barcodes.txt",
            template=config['TRIMMOMATIC_adaptors']
        output:
            adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
        params:
            multiplex=multiplex_tech
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml" # seqtk
        shell:
            """
            # Get the 1st part of the adapter pair
            barcode1=$(cat {input.barcodes} | sed "s/;/+/" | cut -f1 -d'+')
            # Get the 2nd part of adapter pair
            barcode2=$(cat {input.barcodes} | sed "s/;/+/" | cut -f2 -d'+')

            # For MGI, revcomp barcode2
            if [ "{params.multiplex}" == "MGIEasy" ]; then
                barcode2=$(echo $barcode2 | tr 'ATGC' 'TACG' | rev)
            fi

            # Make custom adapters for this sample using its index sequences
            #  2. Make 5-prime adapters
            cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" > {output.adapter}
            #  3. Make 3-prime adapters
            cat {input.template} | sed "s/>\{{6,10\}}/${{barcode1}}/;s/<\{{6,10\}}/${{barcode2}}/" | seqtk seq -r - | tr -s '12' '21' | sed "s^/^_rc/^" >> {output.adapter}
            """

##########
# No index barcodes but adaptor trimming needed
##########

elif config['TRIMMOMATIC_adaptors'] != 'Skip':

    localrules: qc1_copy_fixed_adapter_file

    # Symlink standard file if there is no index_barcode file
    rule qc1_copy_fixed_adapter_file:
        input:
            template=config['TRIMMOMATIC_adaptors']
        output:
            adapter="{wd}/{omics}/1-trimmed/{sample}/adapters.fa"
        shell:
            """
            ln -s --force {input.template} {output.adapter}
            """

##########
# Adaptor-trimming + Quality-trimming or Quality-trimming only?
##########

# If there is adapter file, then priority is to use it
ruleorder: qc1_trim_quality_and_adapter > qc1_trim_quality

if config['TRIMMOMATIC_adaptors'] == 'Skip':

    localrules: qc1_trim_quality

    # Fake a trim
    rule qc1_trim_quality:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
        shell:
            """
            (ln -s --force {input.read_fw} {output.pairead1}
            ln -s --force {input.read_rv} {output.pairead2}
            echo "Skipped trimming and linked raw files to trimmed files.") >& {log}
            """

else:
    rule qc1_trim_quality:
        input:
            unpack(get_raw_reads_for_sample_run)
        output:
            pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
            singleread1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.single.fq.gz",
            pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
            singleread2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.single.fq.gz",
            summary="{wd}/{omics}/1-trimmed/{sample}/{run}.trim.summary"
        shadow:
            "minimal"
        log:
            "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
        resources:
            mem=config['TRIMMOMATIC_memory']
        threads:
            config['TRIMMOMATIC_threads']
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml"#trimmomatic
        shell:
            """
            remote_dir=$(dirname {output.pairead1})
            time ( \
                trimmomatic PE -threads {threads} \
                    -summary {output.summary} \
                    -phred33 \
                    {input.read_fw} {input.read_rv} \
                    {wildcards.run}.1.paired.fq.gz {wildcards.run}.1.single.fq.gz \
                    {wildcards.run}.2.paired.fq.gz {wildcards.run}.2.single.fq.gz \
                    TRAILING:20 LEADING:5 SLIDINGWINDOW:4:20 \
                && rsync -a * $remote_dir/
            ) >& {log}
            """

# Trim with Trimmomatic
rule qc1_trim_quality_and_adapter:
    input:
        unpack(get_raw_reads_for_sample_run),
        adapter='{wd}/{omics}/1-trimmed/{sample}/adapters.fa'
    output:
        pairead1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.paired.fq.gz",
        singleread1="{wd}/{omics}/1-trimmed/{sample}/{run}.1.single.fq.gz",
        pairead2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.paired.fq.gz",
        singleread2="{wd}/{omics}/1-trimmed/{sample}/{run}.2.single.fq.gz",
        summary="{wd}/{omics}/1-trimmed/{sample}/{run}.trim.summary"
    shadow:
        "minimal"
    params:
        simple_clip_threshold=trimmomatic_simple_clip_threshold
    log:
        "{wd}/logs/{omics}/1-trimmed/{sample}_{run}_qc1_trim_quality.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"#trimmomatic
    shell:
        """
        remote_dir=$(dirname {output.pairead1})
        time ( \
            trimmomatic PE -threads {threads} \
                -summary {output.summary} \
                -phred33 \
                {input.read_fw} {input.read_rv} \
                {wildcards.run}.1.paired.fq.gz {wildcards.run}.1.single.fq.gz \
                {wildcards.run}.2.paired.fq.gz {wildcards.run}.2.single.fq.gz \
                TRAILING:20 LEADING:5 SLIDINGWINDOW:4:20 \
                ILLUMINACLIP:{input.adapter}:2:30:{params.simple_clip_threshold}:1:TRUE \
            && rsync -a * $remote_dir/
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
        time (sh -c 'gzip -cd {input.pairead} | awk -v f="{wildcards.sample}_{wildcards.group}" "{{if(NR%4==2) print length(\$1),f}}" | sort -n | uniq -c > {output.length}') >& {log}
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
        time (cat {input.length} > {output.readlen_dist}) >& {log}
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
    resources:
        mem=config['TRIMMOMATIC_memory']
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (Rscript {script_dir}/QC_cumulative_read_length_plot.R --input {input.readlen_dist} --frac {params.cutoff} --out_plot {output.plot} --out_cutoff {output.cutoff_file}) >& {log}
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
        trim_threads=config['TRIMMOMATIC_threads'],
        trim_memory=config['TRIMMOMATIC_memory']
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
# Please do not use '.fasta' or '.fa' or '.fna' extension for the fasta file. Name it without extension.

PATH_host_genome:
NAME_host_genome:
BWA_index_host_memory: 40
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
$(for i in {ilmn_samples}; do echo "- $i"; done)
___EOF___

        echo {ilmn_samples} >& {log}
        """
