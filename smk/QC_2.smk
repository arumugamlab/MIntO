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

# args = sys.argv
# config_path = args[args.index("--configfile") + 1]
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file") #: ", config_path)
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
                    if path.exists(working_dir+'/'+omics+'/1-trimmed/'+ilmn+'/'+ilmn+'.1.paired.fq.gz') is True:
                        #print(ilmn)
                        ilmn_samples.append(ilmn)
                    else:
                        raise TypeError('ERROR in ', config_path, ': ILLUMINA list of samples does not exist. Please, complete ', config_path)
            n_samples=len(ilmn_samples)+3
        except: 
            print('ERROR in ', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete ', config_path)
else:
    print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)

if config['PATH_host_genome'] is None:
    print('ERROR in ', config_path, ': PATH_host_genome variable is empty. Please, complete ', config_path)
elif path.exists(config['PATH_host_genome']) is False:
    print('ERROR in ', config_path, ': PATH_host_genome variable path does not exit. Please, complete ', config_path)
else:
    host_genome_db=config["PATH_host_genome"]

if config['NAME_host_genome'] is None:
    print('ERROR in ', config_path, ': NAME_host_genome variable is empty. Please, complete ', config_path)
elif path.exists(host_genome_db+'/'+config['NAME_host_genome']) is False:
    print('ERROR in ', config_path, ': NAME_host_genome variable path does not exit. Please, complete ', config_path)
else:
    host_genome_name=config["NAME_host_genome"]

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

if config['BWA_index_host_threads'] is None:
    print('ERROR in ', config_path, ': BWA_index_host_threads variable is empty. Please, complete ', config_path)
elif type(config['BWA_index_host_threads']) != int:
    print('ERROR in ', config_path, ': BWA_index_host_threads variable is not an integer. Please, complete ', config_path)

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

if omics == 'metaG':
    if config['TAXA_threads'] is None:
        print('ERROR in ', config_path, ': TAXA_threads variable is empty. Please, complete ', config_path)
    elif type(config['TAXA_threads']) != int:
        print('ERROR in ', config_path, ': TAXA_threads variable is not an integer. Please, complete ', config_path)
    else:
        TAXA_threads=config["TAXA_threads"]

if omics == 'metaG':
    if config['TAXA_memory'] is None:
        print('ERROR in ', config_path, ': TAXA_memory variable is empty. Please, complete ', config_path)
    elif type(config['TAXA_memory']) != int:
        print('ERROR in ', config_path, ': TAXA_memory variable is not an integer. Please, complete ', config_path)
    else:
        TAXA_memory=config["TAXA_memory"]

if omics == 'metaG':
    if config['taxa_profile'] in ('metaphlan', 'motus_rel', 'motus_raw'):
        taxa=config["taxa_profile"]
    else:
        print('ERROR in ', config_path, ': taxa_profile variable is not correct. "taxa_profile" variable should be metaphlan, motus_rel or motus_raw.')

if omics == 'metaT':
    if config['sortmeRNA_threads'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_threads variable is empty. Please, complete ', config_path)
    elif type(config['sortmeRNA_threads']) != int:
        print('ERROR in ', config_path, ': sortmeRNA_threads variable is not an integer. Please, complete ', config_path)
    else:
        sortmeRNA_threads=config["sortmeRNA_threads"]

if omics == 'metaT':
    if config['sortmeRNA_memory'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_memory variable is empty. Please, complete ', config_path)
    elif type(config['sortmeRNA_memory']) != int:
        print('ERROR in ', config_path, ': sortmeRNA_memory variable is not an integer. Please, complete ', config_path)
    else:
        sortmeRNA_memory=config["sortmeRNA_memory"]

if omics == 'metaT':
    if config['sortmeRNA_db'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_db variable is empty. Please, complete ', config_path)
    elif path.exists(config['sortmeRNA_db']) is False:
        print('ERROR in ', config_path, ': sortmeRNA_db variable path does not exit. Please, complete ', config_path)
    else:
        sortmeRNA_db=config["sortmeRNA_db"]

if omics == 'metaT':
    if config['sortmeRNA_db_idx'] is None:
        print('ERROR in ', config_path, ': sortmeRNA_db_idx variable is empty. Please, complete ', config_path)
    elif path.exists(config['sortmeRNA_db_idx']) is False:
        print('WARNING in ', config_path, ': sortmeRNA_db_idx variable path does not exit. Indexed SortMeRNA database will be generated in', config['sortmeRNA_db_idx'])
        sortmeRNA_db_idx=config["sortmeRNA_db_idx"]
    else:
        sortmeRNA_db_idx=config["sortmeRNA_db_idx"]


# Define all the outputs needed by target 'all'
def filter_trim_length_output():
    result = expand("{wd}/{omics}/3-minlength/{sample}/{sample}.1.paired.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
    expand("{wd}/{omics}/3-minlength/{sample}/{sample}.1.single.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
    expand("{wd}/{omics}/3-minlength/{sample}/{sample}.2.paired.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
    expand("{wd}/{omics}/3-minlength/{sample}/{sample}.2.single.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
    expand("{wd}/{omics}/3-minlength/{sample}/{sample}_trimlog_length_adapter", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [])
    return(result)

def filter_host_genome_output():
    result = expand("{host_genome_path}/BWA_index/{host_genome_name}.pac", 
                    host_genome_path=host_genome_db, 
                    host_genome_name=host_genome_name),\
    expand("{wd}/{omics}/4-hostfree/{sample}/{sample}.{group}.fq.gz", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    group = ['1', '2']),\
    expand("{wd}/{omics}/4-hostfree/{sample}/{sample}.host.reads", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [])
    return(result)

def filter_config_yml_output():
    result = expand("{wd}/{omics}/assembly.yaml", 
                    wd = working_dir,
                    omics = omics)
    return(result)

if omics == 'metaG': 
    sortmeRNA_db_idx='None'
    sortmeRNA_db='None'
    sortmeRNA_memory=1
    sortmeRNA_threads=1
    #taxa=config["taxa_profile"]
    #TAXA_threads=config["TAXA_threads"]
    #TAXA_memory=config["TAXA_memory"]
    if taxa != 'metaphlan':
        taxa_motus=taxa
        #print({taxa_motus})
        def extra_output():
            result = expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxa_motus}", 
                        wd = working_dir,
                        omics = omics,
                        sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                        taxa_motus = taxa_motus),\
            expand("{wd}/output/6-taxa_profile/{taxa_motus}.PCoA.Bray_Curtis.pdf", 
                        wd = working_dir,
                        taxa_motus = taxa_motus),\
            expand("{wd}/output/6-taxa_profile/{taxa_motus}.Top15genera.pdf", 
                        wd = working_dir,
                        taxa_motus = taxa_motus)
            return(result)
        def filter_config_yml_output():
            result = expand("{wd}/{omics}/assembly.yaml", 
                        wd = working_dir,
                        omics = omics),\
            expand("{wd}/{omics}/mapping.yaml", 
                        wd = working_dir,
                        omics = omics)
            return(result)
    elif taxa == 'metaphlan':
        taxa_motus='None'
        def extra_output():
            result = expand("{minto_dir}/logs/metaphlan_download_db_checkpoint.log",
                        minto_dir=minto_dir),\
            expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan", 
                        wd = working_dir,
                        omics = omics,
                        sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
            expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan.read_stats",
                        wd = working_dir,
                        omics = omics,
                        sample = config["ILLUMINA"] if "ILLUMINA" in config else []),\
            expand("{wd}/output/6-taxa_profile/metaphlan.PCoA.Bray_Curtis.pdf", 
                        wd = working_dir),\
                        #taxa = 'metaphlan'),\
            expand("{wd}/output/6-taxa_profile/metaphlan.Top15genera.pdf", 
                        wd = working_dir)
                        #taxa = 'metaphlan')
            return(result)
        def filter_config_yml_output():
            result = expand("{wd}/{omics}/assembly.yaml", 
                        wd = working_dir,
                        omics = omics),\
            expand("{wd}/{omics}/mapping.yaml", 
                        wd = working_dir,
                        omics = omics)
            return(result)
else:
    #sortmeRNA_db_idx=config["sortmeRNA_db_idx"]
    #sortmeRNA_db=config["sortmeRNA_db"]
    #sortmeRNA_memory=config["sortmeRNA_memory"]
    #sortmeRNA_threads=config["sortmeRNA_threads"]
    taxa='None'
    taxa_motus='None'
    TAXA_threads=1
    TAXA_memory=1
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
    def filter_config_yml_output():
        result = expand("{wd}/{omics}/mapping.yaml", 
                    wd = working_dir,
                    omics = omics)
        return(result)


rule all:
    input: 
        filter_trim_length_output(),
        filter_host_genome_output(),
        extra_output(),
        filter_config_yml_output()

###############################################################################################
# Pre-processing of metaG and metaT data step 
# Read length filtering using the MINLEN 
###############################################################################################

rule qc2_filter_trim_length:
    input:
        read_fw=lambda wildcards: '{wd}/{omics}/1-trimmed/{sample}/{sample}.1.paired.fq.gz',
        read_rv=lambda wildcards: '{wd}/{omics}/1-trimmed/{sample}/{sample}.2.paired.fq.gz',
    output: 
        pairead1="{wd}/{omics}/3-minlength/{sample}/{sample}.1.paired.fq.gz", 
        singleread1="{wd}/{omics}/3-minlength/{sample}/{sample}.1.single.fq.gz", 
        pairead2="{wd}/{omics}/3-minlength/{sample}/{sample}.2.paired.fq.gz", 
        singleread2="{wd}/{omics}/3-minlength/{sample}/{sample}.2.single.fq.gz", 
        log_adapter="{wd}/{omics}/3-minlength/{sample}/{sample}_trimlog_length_adapter"
    params:
        tmp_trim_len=lambda wildcards: "{local_dir}/{omics}_{sample}filter_trim_length/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
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
        remote_dir=$(dirname {output.pairead1})
        echo {resources.mem}
        time (trimmomatic PE -threads {threads} -trimlog {params.tmp_trim_len}{wildcards.sample}_trimlog_length_adapter \
-phred33 {input.read_fw} {input.read_rv} {params.tmp_trim_len}{wildcards.sample}.1.paired.fq.gz {params.tmp_trim_len}{wildcards.sample}.1.single.fq.gz \
{params.tmp_trim_len}{wildcards.sample}.2.paired.fq.gz {params.tmp_trim_len}{wildcards.sample}.2.single.fq.gz MINLEN:{params.read_length_cutoff}
        rsync {params.tmp_trim_len}* $remote_dir) >& {log}
        rm -rf {params.tmp_trim_len} """

###############################################################################################
# Pre-processing of metaG and metaT data
# Remove host genome sequences 
###############################################################################################

rule filter_host_genome_bwaindex:
    input: 
        host_genome=lambda wildcards: "{host_genome_path}/{host_genome_name}".format(host_genome_path=host_genome_db, host_genome_name=host_genome_name),
    output:
        bwaindex="{host_genome_path}/BWA_index/{host_genome_name}.pac".format(host_genome_path=host_genome_db, host_genome_name=host_genome_name),
    params:
        tmp_bwaindex=lambda wildcards: "{local_dir}/{host_genome_name}_bwaindex/".format(local_dir=local_dir, host_genome_name=host_genome_name),
    log: 
        "{wd}/logs/{omics}/4-hostfree/filter_host_genome_BWAindex.log".format(wd=working_dir, omics = omics)
    resources:
        mem=config['BWA_index_host_memory']
    threads:
        config['BWA_index_host_threads'] 
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        mkdir -p {params.tmp_bwaindex}
        time (bwa-mem2 index {host_genome_db}/{host_genome_name} -p {params.tmp_bwaindex}{host_genome_name}
        rsync {params.tmp_bwaindex}* {host_genome_db}/BWA_index/
        rm -rf {params.tmp_bwaindex}) &> {log} """

rule qc2_filter_host_genome:
    input: 
        pairead_fw=lambda wildcards: '{wd}/{omics}/3-minlength/{sample}/{sample}.1.paired.fq.gz',
        pairead_rv=lambda wildcards: '{wd}/{omics}/3-minlength/{sample}/{sample}.2.paired.fq.gz',
        bwaindex="{host_genome_path}/BWA_index/{host_genome_name}.pac".format(host_genome_path=host_genome_db, host_genome_name=host_genome_name),
    output: 
        bwaout_alig="{wd}/{omics}/4-hostfree/{sample}/{sample}.host.reads", 
        host_free_fw="{wd}/{omics}/4-hostfree/{sample}/{sample}.1.fq.gz",
        host_free_rv="{wd}/{omics}/4-hostfree/{sample}/{sample}.2.fq.gz",
    params:
        tmp_bwa=lambda wildcards: "{local_dir}/{omics}_{sample}_filter_host_genome/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
        bwaindex="{host_genome_path}/BWA_index/{host_genome_name}".format(host_genome_path=host_genome_db, host_genome_name=host_genome_name),
    log:
        "{wd}/logs/{omics}/4-hostfree/{sample}_filter_host_genome_BWA.log"
    resources:
        mem=config['BWA_host_memory']
    threads:
        config['BWA_host_threads'] 
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """ mkdir -p {params.tmp_bwa}
        bwaindex_dir=$(dirname {input.bwaindex})
        remote_dir=$(dirname {output.host_free_fw})
        (time bwa-mem2 mem -a -t {threads} -v 3 {params.bwaindex} {input.pairead_fw} {input.pairead_rv}| \
msamtools filter -S -l 30 - | cut -f1 > {params.tmp_bwa}{wildcards.sample}.host.reads
        mseqtools subset --exclude --paired --input {input.pairead_fw} --output {params.tmp_bwa}{wildcards.sample}.1.fq.gz --list {params.tmp_bwa}{wildcards.sample}.host.reads 
        mseqtools subset --exclude --paired --input {input.pairead_rv} --output {params.tmp_bwa}{wildcards.sample}.2.fq.gz --list {params.tmp_bwa}{wildcards.sample}.host.reads  
        rsync {params.tmp_bwa}{wildcards.sample}.host.reads {output.bwaout_alig}
        rsync {params.tmp_bwa}{wildcards.sample}.1.fq.gz {output.host_free_fw}
        rsync {params.tmp_bwa}{wildcards.sample}.2.fq.gz {output.host_free_rv}) >& {log}
        rm -rf {params.tmp_bwa}"""

###############################################################################################
# Assembly-free taxonomy profiling - only on metaG data
###############################################################################################
rule taxonomic_profile_metaphlan_download_db:
    output: 
        metaphlan_db="{minto_dir}/data/metaphlan/mpa_v30_CHOCOPhlAn_201901.fna.bz2".format(minto_dir=minto_dir)
    resources:
        mem=TAXA_memory #lambda wildcards, input: len(input.host_free_fw) + 2
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/metaphlan_download_db.log".format(wd=working_dir, omics = omics)
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml" #metaphlan or motus2
    shell:
        """ 
        time (metaphlan --version
        metaphlan --install --bowtie2db {minto_dir}/data/metaphlan/) >& {log}
        """

rule taxonomic_profile_metaphlan:
    input:
        #metaphlan_db="{minto_dir}/data/metaphlan/mpa_v30_CHOCOPhlAn_201901.fna.bz2".format(minto_dir=minto_dir), 
        host_free_fw= lambda wildcards: "{wd}/{omics}/4-hostfree/{sample}/{sample}.1.fq.gz",
        host_free_rv= lambda wildcards: "{wd}/{omics}/4-hostfree/{sample}/{sample}.2.fq.gz"
    output:
        metaphlan_db="{minto_dir}/logs/metaphlan_download_db_checkpoint.log", 
        ra="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan",
        ra_stats="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan.read_stats"
    params:
        tmp_taxa_prof=lambda wildcards: "{local_dir}/{omics}_{sample}.metaphlan_taxonomic_profile/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
    resources:
        mem=TAXA_memory #lambda wildcards, input: len(input.host_free_fw) + 2
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.metaphlan.log"
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml" #metaphlan or motus2
    shell:
        """ mkdir -p {params.tmp_taxa_prof}
        remote_dir=$(dirname {output.ra})
        time (metaphlan --install --bowtie2db {minto_dir}/data/metaphlan/
        echo 'MetaPhlAn database downloaded'
        if [ $? -eq 0 ]; then
        echo OK > {minto_dir}/logs/metaphlan_download_db_checkpoint.log
        else
        echo FAIL
        fi
        metaphlan --bowtie2db {minto_dir}/data/metaphlan/ {input.host_free_fw},{input.host_free_rv} --input_type fastq --bowtie2out {params.tmp_taxa_prof}{wildcards.sample}.bowtie2.bz2 --nproc {threads} -o {params.tmp_taxa_prof}{wildcards.sample}.metaphlan -t rel_ab_w_read_stats
        metaphlan {params.tmp_taxa_prof}{wildcards.sample}.bowtie2.bz2 --input_type bowtie2out --nproc {threads} -o {params.tmp_taxa_prof}{wildcards.sample}.metaphlan.read_stats -t rel_ab_w_read_stats
        rsync {params.tmp_taxa_prof}* $remote_dir) >& {log}
        rm -rf {params.tmp_taxa_prof}
        """

rule taxonomic_profile_motus:
    input: 
        host_free_fw= lambda wildcards: "{wd}/{omics}/4-hostfree/{sample}/{sample}.1.fq.gz",#.format(wd = working_dir, omics = omics, sample = sample),
        host_free_rv= lambda wildcards: "{wd}/{omics}/4-hostfree/{sample}/{sample}.2.fq.gz",#.format(wd = working_dir, omics = omics, sample = sample),
    output: 
        ra="{{wd}}/{{omics}}/6-taxa_profile/{{sample}}/{{sample}}.{taxa_motus}".format(taxa_motus=taxa_motus),
    params:
        tmp_taxa_prof="{local_dir}/{omics}_{{sample}}.{taxa_motus}.taxonomic_profile/".format(local_dir=local_dir, omics = omics, taxa_motus = taxa_motus),
    resources:
        mem=TAXA_memory #lambda wildcards, input: len(input.host_free_fw) + 2
    threads:
        TAXA_threads
    log:
        "{{wd}}/logs/{{omics}}/6-taxa_profile/{{sample}}.{taxa_motus}.log".format(taxa_motus = taxa_motus),
    conda:
        config["minto_dir"]+"/envs/motus_env.yml"
    shell:
        """ mkdir -p {params.tmp_taxa_prof}
        remote_dir=$(dirname {output.ra})
        time ( if [ "{taxa_motus}" == "motus_rel" ]
            then motus profile -t {threads} -f {input.host_free_fw} -r {input.host_free_rv} -o {params.tmp_taxa_prof}{wildcards.sample}.{taxa_motus} -n {wildcards.sample} -u -q
        fi
        if [ "{taxa_motus}" == "motus_raw" ]
            then motus profile -t {threads} -f {input.host_free_fw} -r {input.host_free_rv} -o {params.tmp_taxa_prof}{wildcards.sample}.{taxa_motus} -n {wildcards.sample} -c -u -q
        fi
        rsync {params.tmp_taxa_prof}* $remote_dir) >& {log}
        rm -rf {params.tmp_taxa_prof}
        """

rule taxonomic_profile_plot:
    input: 
        ra=expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxa}", wd = working_dir, omics = omics, sample=ilmn_samples , taxa= taxa),
    output: 
        pcoa="{wd}/output/6-taxa_profile/{taxa}.PCoA.Bray_Curtis.pdf".format(wd = working_dir, taxa= taxa),
        barplot="{wd}/output/6-taxa_profile/{taxa}.Top15genera.pdf".format(wd = working_dir, taxa= taxa),
    params:
        tmp_taxa_prof_plot=lambda wildcards: "{local_dir}/{omics}.{taxa}_taxonomic_profile_plot/".format(local_dir=local_dir, omics = omics, taxa= taxa),
        sample_number=n_samples
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxa}_plot.log".format(wd = working_dir, omics = omics, taxa= taxa),
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml" #R
    shell:
        """ mkdir -p {params.tmp_taxa_prof_plot}
        remote_dir_in=$(dirname $(dirname {input.ra[0]}))
        remote_dir_out=$(dirname {output.pcoa})
        time ( if [ "{taxa}" == "metaphlan" ]
        then cd ${{remote_dir_in}}/
        merge_metaphlan_tables.py */*.metaphlan.read_stats > {params.tmp_taxa_prof_plot}merged_abundance_table.txt
        grep -E "s__|clade" {params.tmp_taxa_prof_plot}merged_abundance_table.txt | sed 's/^.*s__//g' | cut -f1,3-{params.sample_number} | sed -e 's/clade_name/body_site/g' > {params.tmp_taxa_prof_plot}merged_abundance_table_species.txt
        rsync {params.tmp_taxa_prof_plot}merged_abundance_table.txt $remote_dir_in
        rsync {params.tmp_taxa_prof_plot}merged_abundance_table_species.txt $remote_dir_in
        fi
        Rscript {script_dir}/plot_6_taxa_profile.R $remote_dir_in/ {taxa} $remote_dir_out/ {metadata}) >& {log}
        rm -rf {params.tmp_taxa_prof_plot}"""

###############################################################################################
# Pre-processing of metaG and metaT data - rRNA filtering - only on metaT data
###############################################################################################

rule qc2_filter_rRNA_index:
    input:
        rRNA_db1="{sortmeRNA_db}/rfam-5.8s-database-id98.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db2="{sortmeRNA_db}/rfam-5s-database-id98.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db3="{sortmeRNA_db}/silva-arc-16s-id95.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db4="{sortmeRNA_db}/silva-arc-23s-id98.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db5="{sortmeRNA_db}/silva-bac-16s-id90.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db6="{sortmeRNA_db}/silva-bac-23s-id98.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db7="{sortmeRNA_db}/silva-euk-18s-id95.fasta".format(sortmeRNA_db=sortmeRNA_db),
        rRNA_db8="{sortmeRNA_db}/silva-euk-28s-id98.fasta".format(sortmeRNA_db=sortmeRNA_db)
    output:
        rRNA_db_index_file = "{sortmeRNA_db_idx}/rRNA_db_index.log".format(sortmeRNA_db_idx=sortmeRNA_db_idx),
        rRNA_db_index = directory(expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx))
    params:
        tmp_sortmerna_index=lambda wildcards: "{local_dir}/{omics}.rRNA_index/".format(local_dir=local_dir, omics = omics),
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
        """mkdir -p {params.tmp_sortmerna_index}idx/ 
        time (sortmerna --workdir {params.tmp_sortmerna_index} --idx-dir {params.tmp_sortmerna_index}idx/ -index 1 \
--ref {input.rRNA_db1} \
--ref {input.rRNA_db2} \
--ref {input.rRNA_db3} \
--ref {input.rRNA_db4} \
--ref {input.rRNA_db5} \
--ref {input.rRNA_db6} \
--ref {input.rRNA_db7} \
--ref {input.rRNA_db8}
        rsync {params.tmp_sortmerna_index}idx/* {output.rRNA_db_index}
        echo 'SortMeRNA indexed rRNA_databases done' > {sortmeRNA_db_idx}/rRNA_db_index.log) >& {log}
        rm -rf {params.tmp_sortmerna_index}"""

rule qc2_filter_rRNA:
    input:
        host_free_fw="{wd}/{omics}/4-hostfree/{sample}/{sample}.1.fq.gz",
        host_free_rv="{wd}/{omics}/4-hostfree/{sample}/{sample}.2.fq.gz",
        rRNA_db_index=expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx)
    output:
        rRNA_out="{wd}/{omics}/5-1-sortmerna/{sample}/out/aligned.blast.gz",
        rRNA_free_fw="{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.1.fq.gz",
        rRNA_free_rv="{wd}/{omics}/5-1-sortmerna/{sample}/{sample}.2.fq.gz"
    params:
        tmp_sortmerna=lambda wildcards: "{local_dir}/{omics}_{sample}.filter_rRNA/".format(local_dir=local_dir, omics = omics, sample = wildcards.sample),
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
        """mkdir -p {params.tmp_sortmerna}
        remote_dir=$(dirname {output.rRNA_free_fw})
        time (sortmerna --workdir {params.tmp_sortmerna} --kvdb {params.tmp_sortmerna}kvdb/ --idx-dir {params.db_idx_dir}/ --readb {params.tmp_sortmerna}readb/ --paired_in --fastx false --blast 1 -threads {threads} --num_alignments 1 \
--ref {params.db_dir}/rfam-5.8s-database-id98.fasta \
--ref {params.db_dir}/rfam-5s-database-id98.fasta \
--ref {params.db_dir}/silva-arc-16s-id95.fasta \
--ref {params.db_dir}/silva-arc-23s-id98.fasta \
--ref {params.db_dir}/silva-bac-16s-id90.fasta \
--ref {params.db_dir}/silva-bac-23s-id98.fasta \
--ref {params.db_dir}/silva-euk-18s-id95.fasta \
--ref {params.db_dir}/silva-euk-28s-id98.fasta \
--reads {input.host_free_fw} --reads {input.host_free_rv}
        zcat {params.tmp_sortmerna}/out/aligned.blast.gz| cut -f1 | uniq > {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list
        mseqtools subset --exclude --paired --input {input.host_free_fw} --output {params.tmp_sortmerna}{wildcards.sample}.1.fq.gz --list {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list 
        mseqtools subset --exclude --paired --input {input.host_free_rv} --output {params.tmp_sortmerna}{wildcards.sample}.2.fq.gz --list {params.tmp_sortmerna}/out/{wildcards.sample}.rrna.list
        rsync {params.tmp_sortmerna}* $remote_dir
        rsync {params.tmp_sortmerna}out/* $remote_dir/out/ ) >& {log}
        rm -rf {params.tmp_sortmerna}"""

##########################################################################################################
# Generate configuration yml file for recovery of MAGs and taxonomic annotation step - assembly/coassembly 
##########################################################################################################

rule qc2_filter_config_yml_assembly:
    input: 
        host_free_fw=expand("{{wd}}/{{omics}}/4-hostfree/{sample}/{sample}.1.fq.gz", sample=ilmn_samples),
    output: 
        config_file="{wd}/{omics}/assembly.yaml"
    params: 
        tmp_assembly_yaml=lambda wildcards: "{local_dir}/{omics}_filter_config_yml_assembly/".format(local_dir=local_dir, omics = omics),
        #sample_names=lambda wildcards, input: "\\n-".join(input.host_free_fw),
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/{omics}/config_yml_assembly.log"
    shell: 
        """
        mkdir -p {params.tmp_assembly_yaml}
        time (files='{ilmn_samples}'; echo $files; echo $files| tr ' ' '\\n'| sed 's/^/- /' - > {params.tmp_assembly_yaml}/samples_illumina.txt
        echo $files| tr ' ' '+'| sed 's/^/  Full: /' - > {params.tmp_assembly_yaml}/samples_illumina_coass.txt
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
# MetaSPAdes settings
#
METASPADES_qoffset: auto
METASPADES_threads:
METASPADES_memory:
METASPADES_hybrid_max_k:
METASPADES_illumina_max_k:

# MEGAHIT settings
#
# Total used memory is not given here but calculated to be 10G per sample in the coassembly. 
# Default parameter is 5 for metaG samples with 10M reads.
# Please make sure that there is enough RAM on the server.
MEGAHIT_threads: 32
MEGAHIT_memory: 5
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

# COASSEMBLY section:
# -------------------
# MEGAHIT coassembly will be performed using the following definitions.
# Each coassembly is named in the LHS, and corresponding illumina sample(s) are in RHS (delimited by '+').
# One coassembly will be performed for each line. 
# E.g. 'Subject1: I3+I4' will result in 1 coassembly: 'Subject1' using I3 and I4 data.
# Memory per coassembly is calculated to be 10G per sample in the coassembly.
# Please make sure that there is enough RAM on the server.
#
COASSEMBLY:" > {params.tmp_assembly_yaml}assembly.yaml
cat {params.tmp_assembly_yaml}/samples_illumina_coass.txt >> {params.tmp_assembly_yaml}assembly.yaml
echo "
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
ILLUMINA:" >> {params.tmp_assembly_yaml}assembly.yaml
cat {params.tmp_assembly_yaml}/samples_illumina.txt >> {params.tmp_assembly_yaml}assembly.yaml
rsync {params.tmp_assembly_yaml}assembly.yaml {output.config_file}) >& {log}
rm -rf {params.tmp_assembly_yaml}
        """

###############################################################################################
# Generate configuration yml file for Alignment, normalization and integration step
###############################################################################################

rule qc2_filter_config_yml_mapping:
    input: 
        host_free_fw=expand("{{wd}}/{{omics}}/4-hostfree/{sample}/{sample}.1.fq.gz", sample=ilmn_samples),
    output: 
        config_file="{wd}/{omics}/mapping.yaml"
    params: 
        tmp_map_yaml=lambda wildcards: "{local_dir}/{omics}_filter_config_yml_map/".format(local_dir=local_dir, omics = omics),
        #sample_names=lambda wildcards, input: "\\n-".join(input.host_free_fw),
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/{omics}/config_yml_mapping.log"
    shell: 
        """
        mkdir -p {params.tmp_map_yaml}
        time (files='{ilmn_samples}'; echo $files; echo $files| tr ' ' '\\n'| sed 's/^/- /' - > {params.tmp_map_yaml}/samples_illumina.txt
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
ILLUMINA:" > {params.tmp_map_yaml}assembly.yaml
cat {params.tmp_map_yaml}/samples_illumina.txt >> {params.tmp_map_yaml}assembly.yaml
rsync {params.tmp_map_yaml}assembly.yaml {output.config_file}) >& {log}
rm -rf {params.tmp_map_yaml}
        """ 
