#!/usr/bin/env python

'''
Download and install dependencies

Authors: Carmen Saenz
'''
#

# configuration yaml file
#import sys
import os.path
from os import path

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



#sortmeRNA_db="{minto_dir}/rRNA_databases/"
#sortmeRNA_db_idx= "{minto_dir}/rRNA_databases/idx"
#eggnog_db="{minto_dir}/data/eggnog_data/"
#kofam_db="{minto_dir}/data/kofam_db/"
#dbCAN_db="{minto_dir}/data/dbCAN_db/"

download_threads=config["download_threads"]
download_memory=config["download_memory"]
index_threads=config["rRNA_index_threads"]
index_memory=config["rRNA_index_memory"]

def rRNA_db_out():
    result = expand("{minto_dir}/data/rRNA_databases/rfam-5.8s-database-id98.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/rfam-5s-database-id98.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-arc-16s-id95.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-arc-23s-id98.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-bac-16s-id90.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-bac-23s-id98.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-euk-18s-id95.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/silva-euk-28s-id98.fasta",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/rRNA_databases/idx/rRNA_db_index.log",
        minto_dir=minto_dir)
    expand("{minto_dir}/data/rRNA_databases/idx",
        minto_dir=minto_dir)
    return(result)

def eggnog_db_out():
    result = expand("{minto_dir}/data/eggnog_data/download_eggnog_data.py",
       minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/eggnog.db",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/eggnog_proteins.dmnd",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/eggnog.taxa.db",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/eggnog.taxa.db.traverse.pkl",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/mmseqs",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/eggnog_data/data/pfam",
        minto_dir=minto_dir)
    return(result)

def Kofam_db_out():
    result = expand("{minto_dir}/data/kofam_db/ko_list",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/kofam_db/profiles",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/kofam_db/README",
        minto_dir=minto_dir)
    return(result)

def dbCAN_db_out():
    result = expand("{minto_dir}/data/dbCAN_db/CAZyDB.09242021.fa",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/dbCAN_db/dbCAN.txt",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/dbCAN_db/tcdb.fa",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/dbCAN_db/tf-1.hmm",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/dbCAN_db/tf-2.hmm",
        minto_dir=minto_dir),\
    expand("{minto_dir}/data/dbCAN_db/stp.hmm",
        minto_dir=minto_dir)
    return(result)

def metaphlan_db_out():
    result=expand("{minto_dir}/data/metaphlan/mpa_v30_CHOCOPhlAn_201901.fna.bz2",
        minto_dir=minto_dir)
    return(result)

# Define all the outputs needed by target 'all'
rule all:
    input: 
        rRNA_db_out(),
        eggnog_db_out(),
        Kofam_db_out(),
        dbCAN_db_out(),
        metaphlan_db_out()

###############################################################################################
# Download and index rRNA database - SortMeRNA
###############################################################################################
rule rRNA_db_download:
    output:
        rRNA_db1="{minto_dir}/data/rRNA_databases/rfam-5.8s-database-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db2="{minto_dir}/data/rRNA_databases/rfam-5s-database-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db3="{minto_dir}/data/rRNA_databases/silva-arc-16s-id95.fasta",#.format(minto_dir=minto_dir),
        rRNA_db4="{minto_dir}/data/rRNA_databases/silva-arc-23s-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db5="{minto_dir}/data/rRNA_databases/silva-bac-16s-id90.fasta",#.format(minto_dir=minto_dir),
        rRNA_db6="{minto_dir}/data/rRNA_databases/silva-bac-23s-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db7="{minto_dir}/data/rRNA_databases/silva-euk-18s-id95.fasta",#.format(minto_dir=minto_dir),
        rRNA_db8="{minto_dir}/data/rRNA_databases/silva-euk-28s-id98.fasta"#.format(minto_dir=minto_dir)
    log:
        "{minto_dir}/logs/rRNA_db_download.log"#.format(minto_dir = config["minto_dir"]),
    resources: mem=download_memory
    threads: download_threads
    shell:
        """ 
        mkdir -p {minto_dir}/data/rRNA_databases
        time (cd {minto_dir}/data/rRNA_databases
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5.8s-database-id98.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/rfam-5s-database-id98.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-16s-id95.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-arc-23s-id98.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-16s-id90.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-bac-23s-id98.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-18s-id95.fasta
        wget https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/silva-euk-28s-id98.fasta
        echo 'SortMeRNA rRNA databases downloaded') &> {log}
        """

rule rRNA_db_index:
    input:
        rRNA_db1="{minto_dir}/data/rRNA_databases/rfam-5.8s-database-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db2="{minto_dir}/data/rRNA_databases/rfam-5s-database-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db3="{minto_dir}/data/rRNA_databases/silva-arc-16s-id95.fasta",#.format(minto_dir=minto_dir),
        rRNA_db4="{minto_dir}/data/rRNA_databases/silva-arc-23s-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db5="{minto_dir}/data/rRNA_databases/silva-bac-16s-id90.fasta",#.format(minto_dir=minto_dir),
        rRNA_db6="{minto_dir}/data/rRNA_databases/silva-bac-23s-id98.fasta",#.format(minto_dir=minto_dir),
        rRNA_db7="{minto_dir}/data/rRNA_databases/silva-euk-18s-id95.fasta",#.format(minto_dir=minto_dir),
        rRNA_db8="{minto_dir}/data/rRNA_databases/silva-euk-28s-id98.fasta"#.format(minto_dir=minto_dir)
    output:
        rRNA_db_index_file = "{minto_dir}/data/rRNA_databases/idx/rRNA_db_index.log",#.format(minto_dir=minto_dir),
        rRNA_db_index = directory("{minto_dir}/data/rRNA_databases/idx")#.format(minto_dir=minto_dir),
    params:
        tmp_sortmerna_index=lambda wildcards: "{local_dir}/MIntO.rRNA_index/".format(local_dir=local_dir),
    resources: mem=index_memory
    threads: index_threads
    log:
       "{minto_dir}/logs/rRNA_db_dindex.log"#.format(minto_dir = config["minto_dir"]),
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
        echo 'SortMeRNA indexed rRNA_databases done' > {output.rRNA_db_index_file}) >& {log}
        rm -rf {params.tmp_sortmerna_index}"""

###############################################################################################
# Download Eggnog database - eggNOGmapper
###############################################################################################

rule eggnog_db:
    output:
        eggnog_py="{minto_dir}/data/eggnog_data/download_eggnog_data.py",
        eggnog_db1="{minto_dir}/data/eggnog_data/data/eggnog.db",
        eggnog_db2="{minto_dir}/data/eggnog_data/data/eggnog_proteins.dmnd",
        eggnog_db3="{minto_dir}/data/eggnog_data/data/eggnog.taxa.db",
        eggnog_db4="{minto_dir}/data/eggnog_data/data/eggnog.taxa.db.traverse.pkl",
        eggnog_db5=directory("{minto_dir}/data/eggnog_data/data/mmseqs"),
        eggnog_db6=directory("{minto_dir}/data/eggnog_data/data/pfam"),
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/eggnog_db_download.log"#.format(minto_dir = config["minto_dir"]),
    conda:
        config["minto_dir"]+"/envs/py38_env.yml"
    shell:
        """ 
        mkdir -p {minto_dir}/data/eggnog_data/data
        time (cd {minto_dir}/data/eggnog_data/
        wget https://raw.githubusercontent.com/eggnogdb/eggnog-mapper/master/download_eggnog_data.py
        printf "y\\ny\\ny\\ny\\ny\\n" |python3 {minto_dir}/data/eggnog_data/download_eggnog_data.py --data_dir {minto_dir}/data/eggnog_data/data -P -M -f
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
        "{minto_dir}/logs/kofam_db_download.log" #.format(minto_dir = config["minto_dir"]),
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
        "{minto_dir}/logs/dbCAN_db_download.log"#.format(minto_dir = config["minto_dir"]),
    conda:
        config["minto_dir"]+"/envs/dbcan_env.yml" #config["dbcan_ironmenv"]
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
# Download metaphlan database - MetaPhlAn3
###############################################################################################

rule metaphlan_db:
    output: 
        metaphlan_db="{minto_dir}/data/metaphlan/mpa_v30_CHOCOPhlAn_201901.fna.bz2"
    resources: 
        mem=download_memory
    threads: 
        download_threads
    log:
        "{minto_dir}/logs/metaphlan_download_db.log"
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml" #metaphlan or motus2
    shell:
        """ 
        mkdir -p {minto_dir}/data/metaphlan/
        time (metaphlan --version
        metaphlan --install --bowtie2db {minto_dir}/data/metaphlan/
        echo 'MetaPhlAn database downloaded') &> {log}
        """








