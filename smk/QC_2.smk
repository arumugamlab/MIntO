#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step 
    - read length filtering
    - host genome filtering
    - rRNA filtering
Assembly-free taxonomy profiling step

Authors: Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, local_dir, minto_dir, script_dir, metadata
include: 'config_parser.smk'

localrules: qc2_filter_config_yml_assembly, qc2_filter_config_yml_mapping, metaphlan_combine_profiles, motus_combine_profiles, plot_taxonomic_profile, motus_calc_motu

##############################################
# Get sample list
##############################################

# Make list of illumina samples, if ILLUMINA in config
try:
    # Make list of illumina samples, if ILLUMINA in config
    ilmn_samples = list()
    if 'ILLUMINA' in config:
        #print("Samples:")
        for ilmn in config["ILLUMINA"]:
            file_found = False
            for loc in ['5-1-sortmerna', '4-hostfree']:
                fwd = "{}/{}/{}/{}/{}.1.fq.gz".format(working_dir, omics, loc, ilmn, ilmn)
                rev = "{}/{}/{}/{}/{}.2.fq.gz".format(working_dir, omics, loc, ilmn, ilmn)
                if (path.exists(fwd) is True) and \
                   (path.exists(rev) is True):
                       file_found = True
            for loc in ['3-minlength', '1-trimmed']:
                fwd = "{}/{}/{}/{}/{}.1.paired.fq.gz".format(working_dir, omics, loc, ilmn, ilmn)
                rev = "{}/{}/{}/{}/{}.2.paired.fq.gz".format(working_dir, omics, loc, ilmn, ilmn)
                if (path.exists(fwd) is True) and \
                   (path.exists(rev) is True):
                       file_found = True
            if file_found == True:
                #print(ilmn)
                ilmn_samples.append(ilmn)
            else:
                raise TypeError('ERROR in ', config_path, ': ILLUMINA list of samples does not exist. Please, complete ', config_path)
    n_samples=len(ilmn_samples)+3
except: 
    print('ERROR in ', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete ', config_path)

##############################################
# Host genome filtering
##############################################

if config['PATH_host_genome'] is None:
    print('ERROR in ', config_path, ': PATH_host_genome variable is empty. Please, complete ', config_path)
elif path.exists(config['PATH_host_genome']) is False:
    print('ERROR in ', config_path, ': PATH_host_genome variable path does not exit. Please, correct ', config_path)
else:
    host_genome_path=config["PATH_host_genome"]

if config['NAME_host_genome'] is None:
    print('ERROR in ', config_path, ': NAME_host_genome variable is empty. Please, complete ', config_path)
elif path.exists(host_genome_path+'/BWA_index/'+config['NAME_host_genome']+'.pac') is True:
    host_genome_name=config["NAME_host_genome"]
    print('Host genome db {dir}/BWA_index/{db} will be used'.format(db=host_genome_name, dir=host_genome_path))
elif path.exists(host_genome_path+'/'+config['NAME_host_genome']) is True:
    host_genome_name=config["NAME_host_genome"]
    print('Host genome sequence {dir}/{db} will be used to create BWA db'.format(db=host_genome_name, dir=host_genome_path))
else:
    print('ERROR in {}: NAME_host_genome={} does not exist as fasta or BWA db in PATH_host_genome={}. Please, correct!'.format(config_path, NAME_host_genome, PATH_host_genome))

##############################################
# Trimmomatic
##############################################

if config['TRIMMOMATIC_threads'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_threads']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_threads variable is not an integer. Please, complete ', config_path)

if config['TRIMMOMATIC_memory'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_memory']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_memory variable is not an integer. Please, complete ', config_path)

if config['TRIMMOMATIC_minlen'] is None:
    print('ERROR in ', config_path, ': TRIMMOMATIC_minlen variable is empty. Please, complete ', config_path)
elif type(config['TRIMMOMATIC_minlen']) != int:
    print('ERROR in ', config_path, ': TRIMMOMATIC_minlen variable is not an integer. Please, complete ', config_path)

##############################################
# BWA for host genome filtering
##############################################

if config['BWA_index_host_memory'] is None:
    print('ERROR in ', config_path, ': BWA_index_host_memory variable is empty. Please, complete ', config_path)
elif type(config['BWA_index_host_memory']) != int:
    print('ERROR in ', config_path, ': BWA_index_host_memory variable is not an integer. Please, complete ', config_path)

if config['BWA_host_threads'] is None:
    print('ERROR in ', config_path, ': BWA_host_threads variable is empty. Please, complete ', config_path)
elif type(config['BWA_host_threads']) != int:
    print('ERROR in ', config_path, ': BWA_host_threads variable is not an integer. Please, complete ', config_path)

if config['BWA_host_memory'] is None:
    print('ERROR in ', config_path, ': BWA_host_memory variable is empty. Please, complete ', config_path)
elif type(config['BWA_host_memory']) != int:
    print('ERROR in ', config_path, ': BWA_host_memory variable is not an integer. Please, complete ', config_path)

##############################################
# taxonomy
##############################################

if config['TAXA_threads'] is None:
    print('ERROR in ', config_path, ': TAXA_threads variable is empty. Please, complete ', config_path)
elif type(config['TAXA_threads']) != int:
    print('ERROR in ', config_path, ': TAXA_threads variable is not an integer. Please, complete ', config_path)
else:
    TAXA_threads=config["TAXA_threads"]

if config['TAXA_memory'] is None:
    print('ERROR in ', config_path, ': TAXA_memory variable is empty. Please, complete ', config_path)
elif type(config['TAXA_memory']) != int:
    print('ERROR in ', config_path, ': TAXA_memory variable is not an integer. Please, complete ', config_path)
else:
    TAXA_memory=config["TAXA_memory"]

allowed = ('metaphlan', 'motus_rel', 'motus_raw')
flags = [0 if x in allowed else 1 for x in config['taxa_profile'].split(",")]
if sum(flags) == 0:
    taxonomy=config["taxa_profile"]
else:
    print('ERROR in ', config_path, ': taxa_profile variable is not correct. "taxa_profile" variable should be metaphlan, motus_rel or motus_raw, or combinations thereof.')

##############################################
# metaG - QC plots
##############################################

plot_args_list = list()
main_factor = None
if config['MAIN_factor'] is not None:
    main_factor = config['MAIN_factor']
    plot_args_list.append('--factor ' + main_factor) 
else:
    print('ERROR in ', config_path, ': "MAIN_factor" variable cannot be empty.')

if config['PLOT_factor2'] is not None:
    plot_args_list.append('--factor2 ' + config['PLOT_factor2']) 
if config['PLOT_time'] is not None:
    plot_args_list.append('--time ' + config['PLOT_time']) 

plot_args_str = ' '.join(plot_args_list)


##############################################
# metaT - rRNA removal
##############################################

sortmeRNA_db = ''
sortmeRNA_db_idx = ''
sortmeRNA_memory = ''
sortmeRNA_threads = 1
if omics == 'metaT':
    if config['sortmeRNA_threads'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_threads variable is empty. Please, complete ', config_path)
    elif type(config['sortmeRNA_threads']) != int:
        print('ERROR in ', config_path, ': sortmeRNA_threads variable is not an integer. Please, complete ', config_path)
    else:
        sortmeRNA_threads=config["sortmeRNA_threads"]

    if config['sortmeRNA_memory'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_memory variable is empty. Please, complete ', config_path)
    elif type(config['sortmeRNA_memory']) != int:
        print('ERROR in ', config_path, ': sortmeRNA_memory variable is not an integer. Please, complete ', config_path)
    else:
        sortmeRNA_memory=config["sortmeRNA_memory"]

    if config['sortmeRNA_db'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_db variable is empty. Please, complete ', config_path)
    elif path.exists(config['sortmeRNA_db']) is False:
        print('ERROR in ', config_path, ': sortmeRNA_db variable path does not exit. Please, complete ', config_path)
    else:
        sortmeRNA_db=config["sortmeRNA_db"]

    if config['sortmeRNA_db_idx'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_db_idx variable is empty. Please, complete ', config_path)
    elif path.exists(config['sortmeRNA_db_idx']) is False:
        print('WARNING in ', config_path, ': sortmeRNA_db_idx variable path does not exit. Indexed SortMeRNA database will be generated in', config['sortmeRNA_db_idx'])
        sortmeRNA_db_idx=config["sortmeRNA_db_idx"]
    else:
        sortmeRNA_db_idx=config["sortmeRNA_db_idx"]


##############################################
# Define all the outputs needed by target 'all'
##############################################

def filter_trim_length_output():
    result = expand("{wd}/{omics}/3-minlength/{sample}/{sample}.{group}.{out}.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    group = ['1', '2'],
                    out = ['single', 'paired'],
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
            expand("{wd}/{omics}/3-minlength/{sample}/{sample}.trim.summary", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [])
    return(result)

def filter_host_genome_output():
    result = expand("{wd}/{omics}/4-hostfree/{sample}/{sample}.{group}.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group = ['1', '2']),\
            expand("{wd}/{omics}/4-hostfree/{sample}/{sample}.host.reads", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [])
    return(result)

def taxonomy_plot_output():
    taxonomies = taxonomy.split(",")
    results = list()
    profiles = expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxonomy}", 
                wd = working_dir,
                omics = omics,
                sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                taxonomy = taxonomies),\
            expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table_species.txt", 
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies)
    plots = expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.PCoA.Bray_Curtis.pdf", 
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies),\
            expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.Top15genera.pdf", 
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies)
    results.append(profiles)
    results.append(plots)
    return(results)

def next_step_config_yml_output():
    result = expand("{wd}/{omics}/assembly.yaml", 
                wd = working_dir,
                omics = omics),\
    expand("{wd}/{omics}/mapping.yaml", 
                wd = working_dir,
                omics = omics)
    return(result)

def extra_output():
    return()

if omics == 'metaT': 
    def extra_output():
        result = expand("{sortmeRNA_db_idx}", 
                    sortmeRNA_db_idx=sortmeRNA_db_idx),\
        expand("{wd}/{omics}/5-1-sortmerna/{sample}/out/aligned.blast.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
        expand("{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.{group}.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group=['1', '2'])
        return(result)
    def next_step_config_yml_output():
        result = expand("{wd}/{omics}/mapping.yaml", 
                    wd = working_dir,
                    omics = omics)
        return(result)


rule all:
    input: 
        #filter_trim_length_output(),
        #filter_host_genome_output(),
        taxonomy_plot_output(),
        extra_output(),
        next_step_config_yml_output()

###############################################################################################
# Pre-processing of metaG and metaT data step 
# Read length filtering using the MINLEN 
###############################################################################################

rule qc2_length_filter:
    input:
        read_fw='{wd}/{omics}/1-trimmed/{sample}/{sample}.1.paired.fq.gz',
        read_rv='{wd}/{omics}/1-trimmed/{sample}/{sample}.2.paired.fq.gz',
    output: 
        paired1="{wd}/{omics}/3-minlength/{sample}/{sample}.1.paired.fq.gz", 
        single1="{wd}/{omics}/3-minlength/{sample}/{sample}.1.single.fq.gz", 
        paired2="{wd}/{omics}/3-minlength/{sample}/{sample}.2.paired.fq.gz", 
        single2="{wd}/{omics}/3-minlength/{sample}/{sample}.2.single.fq.gz", 
        summary="{wd}/{omics}/3-minlength/{sample}/{sample}.trim.summary"
    params:
        tmp_trim_len=lambda wildcards: "{local_dir}/{omics}_{sample}filter_trim_length".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
        read_length_cutoff=config["TRIMMOMATIC_minlen"]
    log:
        "{wd}/logs/{omics}/3-minlength/{sample}.log"
    resources:
        mem=config['TRIMMOMATIC_memory']
    threads:
        config['TRIMMOMATIC_threads'] 
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #trimmomatic
    shell: 
        """
        mkdir -p {params.tmp_trim_len}
        remote_dir=$(dirname {output.paired1})
        echo {resources.mem}
        time ( \
            trimmomatic PE -threads {threads} \
                -summary {output.summary} \
                -phred33 \
                {input.read_fw} {input.read_rv} \
                {params.tmp_trim_len}/$(basename {output.paired1}) {params.tmp_trim_len}/$(basename {output.single1}) \
                {params.tmp_trim_len}/$(basename {output.paired2}) {params.tmp_trim_len}/$(basename {output.single2}) \
                MINLEN:{params.read_length_cutoff} \
            && rsync {params.tmp_trim_len}/* $remote_dir
        ) >& {log}
        rm -rf {params.tmp_trim_len} 
        """

###############################################################################################
# Pre-processing of metaG and metaT data
# Remove host genome sequences 
###############################################################################################

rule bwaindex_host_genome:
    input: 
        host_genome='{somewhere}/{genome}'
    output:
        '{somewhere}/BWA_index/{genome}.0123',
        '{somewhere}/BWA_index/{genome}.amb',
        '{somewhere}/BWA_index/{genome}.ann',
        '{somewhere}/BWA_index/{genome}.bwt.2bit.64',
        '{somewhere}/BWA_index/{genome}.pac',
    params:
        tmp_bwaindex=lambda wildcards: "{local_dir}/{genome}_bwaindex".format(local_dir=local_dir, genome=wildcards.genome),
    log: 
        "{somewhere}/{genome}_BWAindex.log"
    resources:
        mem=config['BWA_index_host_memory']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        mkdir -p {params.tmp_bwaindex}
        time (\
                bwa-mem2 index {input} -p {params.tmp_bwaindex}/{wildcards.genome}
                rsync {params.tmp_bwaindex}/* {wildcards.somewhere}/BWA_index/
                rm -rf {params.tmp_bwaindex}
            ) &> {log} 
        """

rule qc2_host_filter:
    input: 
        pairead_fw=rules.qc2_length_filter.output.paired1,
        pairead_rv=rules.qc2_length_filter.output.paired2,
        bwaindex=lambda wildcards: ancient(expand('{somewhere}/BWA_index/{genome}.{ext}', somewhere=host_genome_path, genome=host_genome_name, ext=['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']))
    output: 
        host_reads="{wd}/{omics}/4-hostfree/{sample}/{sample}.host.reads", 
        host_free_fw="{wd}/{omics}/4-hostfree/{sample}/{sample}.1.fq.gz",
        host_free_rv="{wd}/{omics}/4-hostfree/{sample}/{sample}.2.fq.gz",
    params:
        tmp_bwa=lambda wildcards: "{local_dir}/{omics}_{sample}_filter_host_genome".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
        bwaindex="{host_genome_path}/BWA_index/{host_genome_name}".format(host_genome_path=host_genome_path, host_genome_name=host_genome_name),
    log:
        "{wd}/logs/{omics}/4-hostfree/{sample}_filter_host_genome_BWA.log"
    resources:
        mem=config['BWA_host_memory']
    threads:
        config['BWA_host_threads'] 
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """ 
        mkdir -p {params.tmp_bwa}
        remote_dir=$(dirname {output.host_free_fw})
        host_read_list="{params.tmp_bwa}/$(basename {output.host_reads})"
        time (\
                bwa-mem2 mem -a -t {threads} -v 3 {params.bwaindex} {input.pairead_fw} {input.pairead_rv} \
                | msamtools filter -S -l 30 - \
                | cut -f1 \
                > $host_read_list
                mseqtools subset --exclude --list $host_read_list --paired --input {input.pairead_fw} --output {params.tmp_bwa}/$(basename {output.host_free_fw})
                mseqtools subset --exclude --list $host_read_list --paired --input {input.pairead_rv} --output {params.tmp_bwa}/$(basename {output.host_free_rv})
                rsync {params.tmp_bwa}/* $remote_dir/ 
            ) >& {log}
        rm -rf {params.tmp_bwa}
        """

###############################################################################################
# Pre-processing of metaT data - rRNA filtering - only on metaT data
###############################################################################################

def get_rRNA_db_files(wildcards):
    files = ["rfam-5.8s-database-id98.fasta",
            "rfam-5s-database-id98.fasta",
            "silva-arc-16s-id95.fasta",
            "silva-arc-23s-id98.fasta",
            "silva-bac-16s-id90.fasta",
            "silva-bac-23s-id98.fasta",
            "silva-euk-18s-id95.fasta",
            "silva-euk-28s-id98.fasta"]
    return(expand("{sortmeRNA_db}/{f}",
                    sortmeRNA_db=sortmeRNA_db,
                    f=files))

rule qc2_filter_rRNA_index:
    input:
        rRNA_db=get_rRNA_db_files
    output:
        rRNA_db_index_file = "{sortmeRNA_db_idx}/rRNA_db_index.log".format(sortmeRNA_db_idx=sortmeRNA_db_idx),
        rRNA_db_index = directory(expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx))
    params:
        tmp_working_dir=lambda wildcards: "{local_dir}/{omics}.rRNA_index".format(local_dir=local_dir, omics = omics),
        #sortmeRNA_db_idx = sortmeRNA_db_idx
    resources:
        mem=sortmeRNA_memory
    threads:
        sortmeRNA_threads
    log:
        "{wd}/logs/{omics}/5-1-sortmerna/rRNA_index.log".format(wd=working_dir, omics = omics),
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #sortmerna
    shell: 
        """
        mkdir -p {params.tmp_working_dir}idx/ 
        time (\
            sortmerna --workdir {params.tmp_working_dir}/ --idx-dir {params.tmp_working_dir}/idx/ -index 1 \
                --ref {input.rRNA_db[0]} \
                --ref {input.rRNA_db[1]} \
                --ref {input.rRNA_db[2]} \
                --ref {input.rRNA_db[3]} \
                --ref {input.rRNA_db[4]} \
                --ref {input.rRNA_db[5]} \
                --ref {input.rRNA_db[6]} \
                --ref {input.rRNA_db[7]}
            rsync {params.tmp_working_dir}/idx/* {output.rRNA_db_index}
            echo 'SortMeRNA indexed rRNA_databases done' > {sortmeRNA_db_idx}/rRNA_db_index.log
            ) >& {log}
        rm -rf {params.tmp_working_dir}
        """

rule qc2_filter_rRNA:
    input:
        host_free_fw=rules.qc2_host_filter.output.host_free_fw,
        host_free_rv=rules.qc2_host_filter.output.host_free_rv,
        rRNA_db_index=ancient(expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx))
    output:
        rRNA_out="{wd}/{omics}/5-1-sortmerna/{sample}/out/aligned.blast.gz",
        rRNA_free_fw="{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.1.fq.gz",
        rRNA_free_rv="{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.2.fq.gz"
    params:
        tmp_sortmerna=lambda wildcards: "{local_dir}/{omics}_{sample}.filter_rRNA".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
        db_idx_dir=sortmeRNA_db_idx,
        db_dir=sortmeRNA_db,
    resources:
        mem=sortmeRNA_memory
    threads:
        sortmeRNA_threads
    log:
        "{wd}/logs/{omics}/5-1-sortmerna/{sample}.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #sortmerna
    shell: 
        """
        mkdir -p {params.tmp_sortmerna}
        remote_dir=$(dirname {output.rRNA_free_fw})
        time (sortmerna --workdir {params.tmp_sortmerna} --kvdb {params.tmp_sortmerna}/kvdb/ --idx-dir {params.db_idx_dir}/ --readb {params.tmp_sortmerna}/readb/ --paired_in --fastx false --blast 1 -threads {threads} --num_alignments 1 \
--ref {params.db_dir}/rfam-5.8s-database-id98.fasta \
--ref {params.db_dir}/rfam-5s-database-id98.fasta \
--ref {params.db_dir}/silva-arc-16s-id95.fasta \
--ref {params.db_dir}/silva-arc-23s-id98.fasta \
--ref {params.db_dir}/silva-bac-16s-id90.fasta \
--ref {params.db_dir}/silva-bac-23s-id98.fasta \
--ref {params.db_dir}/silva-euk-18s-id95.fasta \
--ref {params.db_dir}/silva-euk-28s-id98.fasta \
--reads {input.host_free_fw} --reads {input.host_free_rv}
        zcat {params.tmp_sortmerna}/out/aligned.blast.gz | cut -f1 | uniq > {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list
        mseqtools subset --exclude --list {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list --paired --input {input.host_free_fw} --output {params.tmp_sortmerna}/{wildcards.sample}.1.fq.gz
        mseqtools subset --exclude --list {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list --paired --input {input.host_free_rv} --output {params.tmp_sortmerna}/{wildcards.sample}.2.fq.gz
        rsync {params.tmp_sortmerna}/* $remote_dir
        rsync {params.tmp_sortmerna}/out/* $remote_dir/out/ ) >& {log}
        rm -rf {params.tmp_sortmerna}
        """

###############################################################################################
# Assembly-free taxonomy profiling
###############################################################################################

def get_tax_profile_input_files(wildcards):
    if wildcards.omics == "metaT":
        return(expand("{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.{pair}.fq.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    sample = wildcards.sample,
                    pair = [1, 2]))
    else:
        return(expand("{wd}/{omics}/4-hostfree/{sample}/{sample}.{pair}.fq.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    sample = wildcards.sample,
                    pair = [1, 2]))

rule taxonomic_profile_metaphlan_download_db:
    output: 
        metaphlan_db="{minto_dir}/logs/metaphlan_download_db_checkpoint.log".format(minto_dir=minto_dir)
    resources:
        mem=TAXA_memory #lambda wildcards, input: len(input.host_free_fw) + 2
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/metaphlan_download_db.log".format(wd=working_dir, omics = omics)
    conda:
        config["minto_dir"]+"/envs/metaphlan.yml" #metaphlan
    shell:
        """ 
        time (metaphlan --version
        metaphlan --install --bowtie2db {minto_dir}/data/metaphlan/
        echo 'MetaPhlAn database downloaded'
        if [ $? -eq 0 ]; then
        echo OK > {minto_dir}/logs/metaphlan_download_db_checkpoint.log
        else
        echo FAIL
        fi) >& {log}
        """

rule metaphlan_tax_profile:
    input:
        metaphlan_db="{minto_dir}/logs/metaphlan_download_db_checkpoint.log".format(minto_dir=minto_dir), 
        reads=get_tax_profile_input_files
    output:
        ra="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan"
    params:
        tmp_taxa_prof=lambda wildcards: "{local_dir}/{omics}_{sample}.metaphlan_taxonomic_profile".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.metaphlan.log"
    conda:
        config["minto_dir"]+"/envs/metaphlan.yml" #metaphlan
    shell:
        """
        mkdir -p {params.tmp_taxa_prof}
        remote_dir=$(dirname {output.ra})
        time (metaphlan --bowtie2db {minto_dir}/data/metaphlan/ {input.reads[0]},{input.reads[1]} --input_type fastq --bowtie2out {params.tmp_taxa_prof}/{wildcards.sample}.bowtie2.bz2 --nproc {threads} -o {params.tmp_taxa_prof}/{wildcards.sample}.metaphlan -t rel_ab_w_read_stats
        rsync {params.tmp_taxa_prof}/* $remote_dir) >& {log}
        rm -rf {params.tmp_taxa_prof}
        """

rule metaphlan_combine_profiles:
    input: 
        ra=lambda wildcards: expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxonomy}", wd = wildcards.wd, omics = wildcards.omics, taxonomy = wildcards.taxonomy, sample = ilmn_samples),
    output: 
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table.txt",
        species="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table_species.txt",
    wildcard_constraints:
        taxonomy='metaphlan'
    params:
        cut_fields='1,3-{}'.format(len(ilmn_samples)+3)
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}_combine.log"
    conda:
        config["minto_dir"]+"/envs/metaphlan.yml" #metaphlan
    shell:
        """ 
        time (\
                merge_metaphlan_tables.py {input.ra} > {output.merged}
                grep -E "s__|clade_name" {output.merged} | sed 's/^.*s__//' | sed 's/^clade_name/species/' | cut -f{params.cut_fields} > {output.species}
            ) >& {log}
        """

rule motus_map_db:
    input: 
        reads=get_tax_profile_input_files
    output: 
        mgc="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus.mgc",
    params:
        tmp_taxa_prof=lambda wildcards: "{local_dir}/{omics}_{sample}.motus.taxonomic_profile".format(local_dir=local_dir, omics = wildcards.omics, sample = wildcards.sample),
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.motus.mgc.log"
    conda:
        config["minto_dir"]+"/envs/motus_env.yml" #motus3
    shell:
        """
        mkdir -p {params.tmp_taxa_prof}
        time (\
            motus map_tax   -t {threads}          -f {input.reads[0]} -r {input.reads[1]}       -o {params.tmp_taxa_prof}/{wildcards.sample}.motus.bam -b
            motus calc_mgc  -n {wildcards.sample} -i {params.tmp_taxa_prof}/{wildcards.sample}.motus.bam -o {params.tmp_taxa_prof}/{wildcards.sample}.motus.mgc
            rm {params.tmp_taxa_prof}/{wildcards.sample}.motus.bam
            rsync {params.tmp_taxa_prof}/{wildcards.sample}.motus.mgc {output.mgc}
            ) >& {log}
        rm -rf {params.tmp_taxa_prof}
        """

rule motus_calc_motu:
    input: 
        mgc=rules.motus_map_db.output.mgc,
    output: 
        raw="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus_raw",
        rel="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus_rel"
    resources:
        mem=TAXA_memory
    threads: 2
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.motus.log"
    conda:
        config["minto_dir"]+"/envs/motus_env.yml" #motus3
    shell:
        """
        time (\
            motus calc_motu -n {wildcards.sample} -i {input.mgc} -o {output.rel} -p -q
            motus calc_motu -n {wildcards.sample} -i {input.mgc} -o {output.raw} -p -q -c
            ) >& {log}
        """

rule motus_combine_profiles:
    input: 
        profiles=lambda wildcards: expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxonomy}", wd = wildcards.wd, omics = wildcards.omics, taxonomy = wildcards.taxonomy, sample = ilmn_samples),
    output: 
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table.txt",
        species="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table_species.txt",
    wildcard_constraints:
        taxonomy='motus_(raw|rel)'
    params:
        cut_fields='1,3-',
        files=lambda wildcards, input: ",".join(input.profiles)
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}_combine.log"
    conda:
        config["minto_dir"]+"/envs/motus_env.yml" #motus3
    shell:
        """ 
        time (\
                motus merge -i {params.files} | sed 's/^\(\S*\)\s\([^\\t]\+\)\\t/\\2 [\\1]\\t/' | sed -e 's/^consensus_taxonomy \[#mOTU\]/clade_name/' -e 's/NCBI_tax_id/clade_taxid/' > {output.merged}
                grep -E "s__|clade_name" {output.merged} | sed 's/^.*s__//' | sed 's/^clade_name/species/' | cut -f{params.cut_fields} > {output.species}
            ) >& {log}
        """

rule plot_taxonomic_profile:
    input: 
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table.txt",
    output: 
        pcoa="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.PCoA.Bray_Curtis.pdf",
        barplot="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.Top15genera.pdf",
    params:
        plot_args=plot_args_str
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}_plot.log"
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """ 
        time (\
                Rscript {script_dir}/plot_6_taxa_profile.R --table {input.merged} --profiler {wildcards.omics}.{wildcards.taxonomy} --metadata {metadata} --outdir $(dirname {output.pcoa}) {params.plot_args}
            ) >& {log}
        """

##########################################################################################################
# Generate configuration yml file for recovery of MAGs and taxonomic annotation step - assembly/coassembly 
##########################################################################################################

rule qc2_filter_config_yml_assembly:
    input: 
        host_free_fw=expand("{{wd}}/{{omics}}/4-hostfree/{sample}/{sample}.1.fq.gz", sample=ilmn_samples),
    output: 
        config_file="{wd}/{omics}/assembly.yaml"
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/{omics}/config_yml_assembly.log"
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell: 
        """
        cat > {output} <<___EOF___
######################
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
# MetaSPAdes settings
#
METASPADES_qoffset: auto
METASPADES_threads:
METASPADES_memory:
METASPADES_hybrid_max_k:
METASPADES_illumina_max_k:

# MEGAHIT settings
#
# Note on MEGAHIT_memory:
#     MEGAHIT's memory requirement scales with sample count and sample size.
#     Default MEGAHIT_memory below is 10 Gigabytes per sample (plenty for metaG samples of 10M reads).
#     MEGAHIT_memory parameter represents per-sample memory in Gigabytes
#     during the first attempt. If it is not enough and MEGAHIT fails, then
#     each successive attempt will increase per-sample memory by 6 Gibabytes
#     until it succeeds or reaches max_attempts from snakemake.
#     Please make sure that there is enough RAM on the server.
MEGAHIT_memory: 10
MEGAHIT_threads: 32
MEGAHIT_presets:
 - meta-sensitive
 - meta-large

# MetaFlye settings
#
# MetaFlye will be run for each parameter preset listed here.
# By default Flye will be run with these parameters:
#    --meta --genome-size 3.0m
# If you need to add more options, define them here and name them for future reference.
# Notes:
# ------
# 1. Each preset parameter will be applied to each sample. If you only
#    want one parameter to be used, please comment everything else.
# 2. If nothing is listed here, then MetaFlye won't be run.
#    If you just want our default parameters above, then here is a possible option:
#      metaflye-default: ""
# 3. 'tres-o3000-3x' is valid for flye 2.8.3. From 2.9.x, --plasmids and --trestle are
#    not valid. So please use valid options if you are using newer versions of flye.
METAFLYE_presets:
  tres-o3000-3x: --min-overlap 3000 --iterations 3
  #metaflye-default: ""

# BWA settings
# Used when mapping reads back to contigs
#
BWA_threads:

# samtools settings
# Used when sorting bam files
#
SAMTOOLS_sort_threads: 4
SAMTOOLS_sort_memory_gb: 20

# Input data

# HYBRID section:
# ---------------
# MetaSPAdes hybrid assembly will be performed using these definitions.
# Definition format:
#   Each nanopore sample is in the LHS, and corresponding illumina sample(s) are in RHS (delimited by '+').
#   Hybrid assemblies will be performed for each combination of nanopore and illumina samples.
#   E.g.:
#
#   N1: I3+I4 
#
#   The above will result in 2 hybrid assemblies: 'N1-I3' and 'N1-I4'
#
#   This definition was designed based on our practice of pooling multiple samples derived from the same
#   donor for nanopore sequencing to strike a balance between price and sensitivity. If your nanopore samples
#   are paired one-to-one with illumina samples, then the definition gets simpler.
#   E.g.:
#
#   N1: I3
#   N2: I4
#
#   The above will result in 2 hybrid assemblies: 'N1-I3' and 'N2-I4'
#
#HYBRID:

# NANOPORE section:
# -----------------
# List of nanopore samples that will be assembled individually using MetaFlye.
#
#NANOPORE:

# ILLUMINA section:
# -----------------
# List of illumina samples that will be assembled individually using MetaSPAdes.
#
# E.g.:
# - I1
# - I2
#
ILLUMINA:
$(for i in {ilmn_samples}; do echo "- $i"; done)

# COASSEMBLY section:
# -------------------
# MEGAHIT coassembly will be performed using the following definitions.
# Each coassembly is named in the LHS, and corresponding illumina sample(s) are in RHS (delimited by '+').
# One coassembly will be performed for each line. 
# E.g. 'Subject1: I3+I4' will result in 1 coassembly: 'Subject1' using I3 and I4 data.
# Memory per coassembly is calculated to be 10G per sample in the coassembly.
# Please make sure that there is enough RAM on the server.
#
COASSEMBLY:
  Full: $(echo {ilmn_samples} | sed 's/ /+/g')
___EOF___

        R --vanilla --silent --no-echo >> {output} <<___EOF___
library(dplyr)
concat <- function(v) {{
    Reduce(f = paste(sep='+'), x = v)
}}
metadata <- read.table('{metadata}', sep="\\t", header=TRUE) %>% 
    as.data.frame() %>%
    select(sample, {main_factor}) %>% 
    group_by({main_factor}) %>% 
    filter(n() > 1) %>% 
    mutate(co_asm = paste(sample, collapse = "+")) %>% 
    select(-sample) %>%
    slice(1) %>%
    mutate({main_factor}=paste0("  ", {main_factor}))
write.table(metadata, file="", col.names=FALSE, row.names=FALSE, quote=FALSE, sep=": ")
___EOF___


        echo {ilmn_samples} >& {log}
        """

###############################################################################################
# Generate configuration yml file for Alignment, normalization and integration step
###############################################################################################

rule qc2_filter_config_yml_mapping:
    input: 
        host_free_fw=expand("{{wd}}/{{omics}}/4-hostfree/{sample}/{sample}.1.fq.gz", sample=ilmn_samples),
    output: 
        config_file="{wd}/{omics}/mapping.yaml"
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/{omics}/config_yml_mapping.log"
    shell: 
        """
        cat > {output} <<___EOF___
######################
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
# BWA Alignment 
msamtools_filter_length: 50
alignment_identity: 95

# Normalization approach
abundance_normalization: TPM
fetchMGs_dir:

# Map reads to reference
map_reference: MAG
PATH_reference: # path to gene catalog fasta file 
NAME_reference: # file name of gene catalog fasta file (MIntO will generate bwa index with same name)


# List of databases used to performe the genome annotation:
# - dbCAN
# - KEGG
# - eggNOG
ANNOTATION:

BWAindex_threads:
BWAindex_memory:
BWA_threads:
BWA_memory:
MIN_mapped_reads:

# Input data

# ILLUMINA section:
# -----------------
# List of illumina samples.
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
