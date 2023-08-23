#!/usr/bin/env python

'''
Download and install dependencies

Authors: Carmen Saenz
'''

# configuration yaml file
#import sys
import os.path
from os import path
import glob

metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'
metaphlan_version = '4.0.6'
motus_version = '3.0.3'

#args = sys.argv
#print(args)
#args_idx = sys.argv.index('--configfile')
#print(args_idx)
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file
if config['minto_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': minto_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]

if config['local_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': local_dir variable is empty. Please, complete ', config_path)
else:
    local_dir = config['local_dir']

if config['download_threads'] is None:
    print('ERROR in ', config_path, ': download_threads variable is empty. Please, complete ', config_path)
elif type(config['download_threads']) != int:
    print('ERROR in ', config_path, ': download_threads variable is not an integer. Please, complete ', config_path)
else:
    download_threads=config["download_threads"]


if config['download_memory'] is None:
    print('ERROR in ', config_path, ': download_memory variable is empty. Please, complete ', config_path)
elif type(config['download_memory']) != int:
    print('ERROR in ', config_path, ': download_memory variable is not an integer. Please, complete ', config_path)
else:
    download_memory=config["download_memory"]

if config['rRNA_index_threads'] is None:
    print('ERROR in ', config_path, ': rRNA_index_threads variable is empty. Please, complete ', config_path)
elif type(config['rRNA_index_threads']) != int:
    print('ERROR in ', config_path, ': rRNA_index_threads variable is not an integer. Please, complete ', config_path)
else:
    index_threads=config["rRNA_index_threads"]

if config['rRNA_index_memory'] is None:
    print('ERROR in ', config_path, ': rRNA_index_memory variable is empty. Please, complete ', config_path)
elif type(config['rRNA_index_memory']) != int:
    print('ERROR in ', config_path, ': rRNA_index_memory variable is not an integer. Please, complete ', config_path)
else:
    index_memory=config["rRNA_index_memory"]


def rRNA_db_out():
    files = ["rfam-5.8s-database-id98.fasta",
                "rfam-5s-database-id98.fasta",
                "silva-arc-16s-id95.fasta",
                "silva-arc-23s-id98.fasta",
                "silva-bac-16s-id90.fasta",
                "silva-bac-23s-id98.fasta",
                "silva-euk-18s-id95.fasta",
                "silva-euk-28s-id98.fasta",
                "idx/rRNA_db_index.log"]
    result = expand("{somewhere}/data/rRNA_databases/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def eggnog_db_out():
    files = [
                "data/eggnog.db",
                "data/eggnog_proteins.dmnd",
                "data/eggnog.taxa.db",
                "data/eggnog.taxa.db.traverse.pkl",
                "data/mmseqs",
                "data/pfam"]
    result = expand("{somewhere}/data/eggnog_data/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def Kofam_db_out():
    files = ["ko_list",
                "profiles",
                "README"]
    result = expand("{somewhere}/data/kofam_db/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def dbCAN_db_out():
    files = ["CAZyDB.09242021.fa",
                "dbCAN.txt",
                "tcdb.fa",
                "tf-1.hmm",
                "tf-2.hmm",
                "stp.hmm"]
    result = expand("{somewhere}/data/dbCAN_db/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def metaphlan_db_out():
    result=expand("{minto_dir}/data/metaphlan/{metaphlan_version}/{metaphlan_index}_VINFO.csv",
        minto_dir=minto_dir,
        metaphlan_version=metaphlan_version,
        metaphlan_index=metaphlan_index)
    return(result)

def motus_db_out():
    result=expand("{minto_dir}/logs/motus_{motus_version}.download_db.log",
        minto_dir=minto_dir,
        motus_version=motus_version)
    return(result)

def checkm2_db_out():
    result=expand("{minto_dir}/data/CheckM2_database/uniref100.KO.1.dmnd",
        minto_dir=minto_dir)
    return(result)

def fetchMGs_out():
    result=expand("{minto_dir}/logs/fetchMGs_download.done",
        minto_dir=minto_dir)
    return(result)

def all_env_out():
    files = ["vamb_env.log",
             "r_pkgs.log",
             "mags_env.log"]
    result = expand("{somewhere}/logs/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

# Define all the outputs needed by target 'all'
rule all:
    input:
        checkm2_db_out(),
        rRNA_db_out(),
        eggnog_db_out(),
        Kofam_db_out(),
        dbCAN_db_out(),
        metaphlan_db_out(),
        motus_db_out(),
        fetchMGs_out(),
        all_env_out()

###############################################################################################
# Download and index rRNA database - SortMeRNA
###############################################################################################
rule rRNA_db_download:
    output:
        "{somewhere}/rRNA_databases/{something}.fasta"
    resources: mem=index_memory
    threads: index_threads
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #sortmerna
    shell:
        """
        mkdir -p {wildcards.somewhere}/rRNA_databases
        cd {wildcards.somewhere}/rRNA_databases
        wget --quiet https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/{wildcards.something}.fasta
        """

def get_rRNA_db_index_input(wildcards):
    files = ["rfam-5.8s-database-id98.fasta",
                "rfam-5s-database-id98.fasta",
                "silva-arc-16s-id95.fasta",
                "silva-arc-23s-id98.fasta",
                "silva-bac-16s-id90.fasta",
                "silva-bac-23s-id98.fasta",
                "silva-euk-18s-id95.fasta",
                "silva-euk-28s-id98.fasta"]
    result = expand("{somewhere}/data/rRNA_databases/{file}",
                somewhere = wildcards.somewhere,
                file = files)
    return(result)


rule rRNA_db_index:
    input:
        get_rRNA_db_index_input
    output:
        rRNA_db_index_file = "{somewhere}/data/rRNA_databases/idx/rRNA_db_index.log",
        rRNA_db_index = directory("{somewhere}/data/rRNA_databases/idx")
    params:
        tmp_sortmerna_index=lambda wildcards: "{local_dir}/MIntO.rRNA_index".format(local_dir=local_dir),
    resources: mem=index_memory
    threads: index_threads
    log:
       "{somewhere}/logs/rRNA_db_index.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #sortmerna
    shell:
        """
        mkdir -p {params.tmp_sortmerna_index}/idx/
        dboption=$(echo {input} | sed "s/ / --ref /g")
        time (sortmerna --workdir {params.tmp_sortmerna_index} --idx-dir {params.tmp_sortmerna_index}/idx/ --index 1 --ref $dboption --threads {threads}
        rsync {params.tmp_sortmerna_index}/idx/* {output.rRNA_db_index}
        echo 'SortMeRNA indexed rRNA_databases done' > {output.rRNA_db_index_file}) >& {log}
        rm -rf {params.tmp_sortmerna_index}
        """

###############################################################################################
# Download Eggnog database - eggNOGmapper
###############################################################################################

rule eggnog_db:
    output:
        eggnog_db1="{minto_dir}/data/eggnog_data/data/eggnog.db",
        eggnog_db2="{minto_dir}/data/eggnog_data/data/eggnog_proteins.dmnd",
        eggnog_db3="{minto_dir}/data/eggnog_data/data/eggnog.taxa.db",
        eggnog_db4="{minto_dir}/data/eggnog_data/data/eggnog.taxa.db.traverse.pkl",
        eggnog_db5=directory("{minto_dir}/data/eggnog_data/data/mmseqs"),
        eggnog_db6=directory("{minto_dir}/data/eggnog_data/data/pfam"),
    params:
        eggnog_db= lambda wildcards: "{minto_dir}/data/eggnog_data/".format(minto_dir = minto_dir)
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/eggnog_db_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/eggnog_data/data
        time (
        download_eggnog_data.py -y --data_dir {minto_dir}/data/eggnog_data/data -P -M -f
        echo 'eggNOG database downloaded') &> {log}
        """

###############################################################################################
# Download KEGG database - KOfamScan
###############################################################################################
rule Kofam_db:
    output:
        kofam_db1="{minto_dir}/data/kofam_db/ko_list",
        #kofam_db2="{minto_dir}/data/kofam_db/profiles.tar",
        kofam_db3=directory("{minto_dir}/data/kofam_db/profiles"), #DIRECTORY xxx
        kofam_db4="{minto_dir}/data/kofam_db/README"
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/kofam_db_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/kofam_db/
        time (cd {minto_dir}/data/kofam_db/
        wget ftp://ftp.genome.jp/pub/db/kofam/*
        gunzip {minto_dir}/data/kofam_db/ko_list.gz
        tar -zxvf {minto_dir}/data/kofam_db/profiles.tar.gz
        echo 'KEGG database downloaded') &> {log}
        """

###############################################################################################
# Download dbCAN database - run_dbcan
###############################################################################################

rule dbCAN_db:
    output:
        dbCAN_db1="{minto_dir}/data/dbCAN_db/CAZyDB.09242021.fa",
        dbCAN_db2="{minto_dir}/data/dbCAN_db/dbCAN.txt",
        dbCAN_db3="{minto_dir}/data/dbCAN_db/tcdb.fa",
        dbCAN_db4="{minto_dir}/data/dbCAN_db/tf-1.hmm",
        dbCAN_db5="{minto_dir}/data/dbCAN_db/tf-2.hmm",
        dbCAN_db6="{minto_dir}/data/dbCAN_db/stp.hmm",
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/dbCAN_db_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/dbCAN_db/
        cd {minto_dir}/data/dbCAN_db/
        time (wget https://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.09242021.fa
        diamond makedb --in CAZyDB.09242021.fa -d CAZy
        wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt
        mv dbCAN-HMMdb-V10.txt dbCAN.txt
        hmmpress dbCAN.txt
        wget https://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa
        diamond makedb --in tcdb.fa -d tcdb
        wget https://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm
        hmmpress tf-1.hmm
        wget https://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm
        hmmpress tf-2.hmm
        wget https://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm
        hmmpress stp.hmm
        echo 'dbCAN database downloaded and installed') &> {log}
        """

## https://github.com/linnabrown/run_dbcan

# Database Installation.
# git clone https://github.com/linnabrown/run_dbcan.git
# cd run_dbcan
# test -d db || mkdir db
# cd db \
#     && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa && diamond makedb --in CAZyDB.09242021.fa -d CAZy \
#     && wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt && mv dbCAN-HMMdb-V10.txt dbCAN.txt && hmmpress dbCAN.txt \
#     && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
#     && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
#     && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
#     && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm


# DATABASES Installation
# https://bcb.unl.edu/dbCAN2/download/Databases/ #Databse -- Database Folder
# https://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.09242021.fa #CAZy.fa--use diamond makedb --in CAZyDB.09242021.fa -d CAZy
# [CAZyme]:included in eCAMI.
# [EC]: included in eCAMI.
# https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt #dbCAN-HMMdb-V10.txt--First use mv dbCAN-HMMdb-V10.txt dbCAN.txt, then use hmmpress dbCAN.txt
# https://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa #tcdb.fa--use diamond makedb --in tcdb.fa -d tcdb
# https://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm #tf-1.hmm--use hmmpress tf-1.hmm
# https://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm #tf-2.hmm--use hmmpress tf-2.hmm
# https://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm #stp.hmm--use hmmpress stp.hmm


###############################################################################################
# Download metaphlan database - MetaPhlAn4
###############################################################################################

rule metaphlan_db:
    output:
        "{minto_dir}/data/metaphlan/{metaphlan_version}/{metaphlan_index}_VINFO.csv"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/metaphlan_{metaphlan_version}_{metaphlan_index}_download_db.log"
    conda:
        config["minto_dir"]+"/envs/metaphlan.yml"
    shell:
        """
        mkdir -p {wildcards.minto_dir}/data/metaphlan/{wildcards.metaphlan_version}
        time (\
                metaphlan --version
                metaphlan --install --index {wildcards.metaphlan_index} --bowtie2db {wildcards.minto_dir}/data/metaphlan/{wildcards.metaphlan_version}/
                if [ $? -eq 0 ]; then
                    echo 'MetaPhlAn database download: OK'
                    echo "{wildcards.metaphlan_index}" > {wildcards.minto_dir}/data/metaphlan/{wildcards.metaphlan_version}/mpa_latest
                else
                    echo 'MetaPhlAn database download: FAIL'
                fi) &> {log}
        """

###############################################################################################
# Download mOTUs database - mOTUs3
###############################################################################################

rule motus_db:
    output:
        "{minto_dir}/data/motus/db.{motus_version}.downloaded"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/motus_{motus_version}.download_db.log"
    conda:
        config["minto_dir"]+"/envs/motus_env.yml"
    shell:
        """
        time (\
            motus downloadDB
            if [ $? -eq 0 ]; then
                echo 'mOTUs3 database download: OK'
                echo OK > {output}
            else
                echo 'mOTUs3 database download: FAIL'
            fi
            ) &> {log}
        """

###############################################################################################
# Download CheckM2 database
###############################################################################################

rule checkm2_db:
    output:
        "{minto_dir}/data/CheckM2_database/uniref100.KO.1.dmnd"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/checkm2_download_db.log"
    conda:
        config["minto_dir"]+"/envs/checkm2.yml"
    shell:
        """
        time (\
            checkm2 database --download --path {minto_dir}/data
            if [ $? -eq 0 ]; then
                echo 'CheckM2 database download: OK'
            else
                echo 'CheckM2 database download: FAIL'
            fi
            ) &> {log}
        """

###############################################################################################
# Download fetchMGs
###############################################################################################

rule download_fetchMGs:
    output:
        done="{minto_dir}/logs/fetchMGs_download.done",
        data=directory("{minto_dir}/data/fetchMGs-1.2")
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/fetchMGs_download.log"
    shell:
        """
        time (\
            cd {minto_dir}/data/
            wget -O fetchMGs-1.2.tar.gz https://github.com/motu-tool/fetchMGs/archive/refs/tags/v1.2.tar.gz
            tar xfz fetchMGs-1.2.tar.gz
            if [ $? -eq 0 ]; then
                echo 'fetchMGs download: OK'
                echo OK > {output.done}
            else
                echo 'fetchMGs download: FAIL'
            fi
            rm fetchMGs-1.2.tar.gz
            ) &> {log}
        """

###############################################################################################
# Generate conda environments
###############################################################################################

rule r_pkgs:
    output:
        r_pkgs="{minto_dir}/logs/r_pkgs.log"
    resources:
        mem=download_memory
    threads:
        download_threads
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (echo 'r_pkgs environment generated') &> {output}
        """

rule mags_gen_vamb:
    output:
        vamb_env="{minto_dir}/logs/vamb_env.log"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/vamb_env.log"
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        time (
        echo 'VAMB environment generated') &> {log}
        """

rule mags_gen:
    output:
        mags_env="{minto_dir}/logs/mags_env.log"
    resources:
        mem=download_memory
    threads:
        download_threads
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (echo 'mags environment generated') &> {output}
        """

rule mags_gen_py36:
    output:
        py36_env="{minto_dir}/logs/py36_env.log"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/py36_env.log"
    conda:
        config["minto_dir"]+"/envs/py36_env.yml"
    shell:
        """
        time (
        echo 'Python 3.6 environment generated') &> {log}
        """
