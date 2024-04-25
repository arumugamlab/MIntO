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

include: 'include/cmdline_validator.smk'

script_dir=workflow.basedir+"/../scripts"

metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'
metaphlan_version = '4.0.6'
phylophlan_db_version = 'Jul20'
motus_version = '3.0.3'
gtdb_release_number = '214'

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
                "KEGG_Pathway2KO.tsv",
                "KEGG_Module2KO.tsv",
                "README"]
    result = expand("{somewhere}/data/kofam_db/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def dbCAN_db_out():
    files = ["CAZyDB.fa",
                "fam-substrate-mapping.tsv",
                "dbCAN.txt",
                "tcdb.fa",
                "tf-1.hmm",
                "tf-2.hmm",
                "stp.hmm"]
    result = expand("{somewhere}/data/dbCAN_db/V12/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def func_db_desc_out():
    files = [
                "KEGG_Pathway.tsv",
                "KEGG_Module.tsv",
                "KEGG_KO.tsv",
                "dbCAN.EC.tsv",
                "eggNOG_OGs.tsv"]
    result = expand("{somewhere}/data/descriptions/{file}",
                somewhere = minto_dir,
                file = files)
    return(result)

def phylophlan_db_out():
    # Hidden feature to turn off phylophlan download.
    # Useful when testing workflows because phylophlan download could take a while and occupies space!
    if ('enable_phylophlan' in config):
        flag = str(config['enable_phylophlan'])
        if (flag.lower() in ('no', 'false', '0')):
            return(list())

    result=expand("{minto_dir}/data/phylophlan/SGB.{version}.txt.bz2",
        minto_dir=minto_dir,
        version=phylophlan_db_version)
    return(result)

def metaphlan_db_out():
    # Hidden feature to turn off metaphlan download.
    # Useful when testing workflows because metaphlan download could take a while and occupies space!
    if ('enable_metaphlan' in config):
        flag = str(config['enable_metaphlan'])
        if (flag.lower() in ('no', 'false', '0')):
            return(list())

    result=expand("{minto_dir}/data/metaphlan/{metaphlan_version}/{metaphlan_index}_VINFO.csv",
        minto_dir=minto_dir,
        metaphlan_version=metaphlan_version,
        metaphlan_index=metaphlan_index)
    return(result)

def motus_db_out():
    # Hidden feature to turn off motus download.
    # Useful when testing workflows because motus download could take a while and occupies space!
    if ('enable_motus' in config):
        flag = str(config['enable_motus'])
        if (flag.lower() in ('no', 'false', '0')):
            return(list())

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

def gtdb_db_out():
    # Hidden feature to turn off GTDB download.
    # Useful when testing workflows because GTDB download could take a full day and occupies space!
    if ('enable_GTDB' in config):
        flag = str(config['enable_GTDB'])
        if (flag.lower() in ('no', 'false', '0')):
            return(list())

    result=expand("{minto_dir}/data/GTDB/r{gtdb_release_number}/taxonomy/gtdb_taxonomy.tsv",
        minto_dir=minto_dir,
        gtdb_release_number=gtdb_release_number)
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
        func_db_desc_out(),
        metaphlan_db_out(),
        phylophlan_db_out(),
        motus_db_out(),
        fetchMGs_out(),
        gtdb_db_out(),
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
    shadow:
        "minimal"
    resources: mem=index_memory
    threads: index_threads
    log:
       "{somewhere}/logs/rRNA_db_index.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #sortmerna
    shell:
        """
        mkdir -p sortmerna/idx/
        dboption=$(echo {input} | sed "s/ / --ref /g")
        time (sortmerna --workdir sortmerna --idx-dir sortmerna/idx/ --index 1 --ref $dboption --threads {threads}
        rsync sortmerna/idx/* {output.rRNA_db_index}
        echo 'SortMeRNA indexed rRNA_databases done' > {output.rRNA_db_index_file}) >& {log}
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
        kofam_db2=directory("{minto_dir}/data/kofam_db/profiles"),
        kofam_db3="{minto_dir}/data/kofam_db/profiles/prokaryote.hal",
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
        time (
            cd {minto_dir}/data/kofam_db/

            # Get kofam databases
            wget ftp://ftp.genome.jp/pub/db/kofam/*
            gunzip ko_list.gz
            tar -zxvf profiles.tar.gz

            echo 'kofam database downloaded'
        ) &> {log}
        """

rule KEGG_maps:
    output:
        kegg_module="{minto_dir}/data/kofam_db/KEGG_Module2KO.tsv",
        kegg_pathway="{minto_dir}/data/kofam_db/KEGG_Pathway2KO.tsv",
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/KEGG_maps_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/kofam_db/
        cd {minto_dir}/data/kofam_db/
        time (

            # Get KEGG Pathway2KO mapping
            curl --silent https://rest.kegg.jp/list/pathway > KEGG_Pathway.tsv
            for i in $(cut -f1 KEGG_Pathway.tsv); do
                curl --silent https://rest.kegg.jp/link/ko/$i
            done | sed 's/^path://;s/ko://' | grep '.' > KEGG_Pathway2KO.tsv

            # Get KEGG Module2KO mapping
            curl --silent https://rest.kegg.jp/list/module > KEGG_Module.tsv
            for i in $(cut -f1 KEGG_Module.tsv); do
                curl --silent https://rest.kegg.jp/link/ko/$i
            done | sed 's/^md://;s/ko://' | grep '.' > KEGG_Module2KO.tsv

            echo 'KEGG mapping downloaded'
        ) &> {log}
        """

rule functional_db_descriptions:
    input:
        kegg_ko="{minto_dir}/data/descriptions/include/KEGG_KO.tsv",
        kegg_module="{minto_dir}/data/descriptions/include/KEGG_Module.tsv",
        kegg_pathway="{minto_dir}/data/descriptions/include/KEGG_Pathway.tsv",
    output:
        kegg_ko="{minto_dir}/data/descriptions/KEGG_KO.tsv",
        kegg_module="{minto_dir}/data/descriptions/KEGG_Module.tsv",
        kegg_pathway="{minto_dir}/data/descriptions/KEGG_Pathway.tsv",
        eggnog_desc="{minto_dir}/data/descriptions/eggNOG_OGs.tsv",
        EC_desc="{minto_dir}/data/descriptions/dbCAN.EC.tsv",
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/func_db_desc_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/descriptions/
        time (

            ## KEGG

            # Get MIntO-provided descriptions that includes descriptions for eggNOG's version-freeze, and update using KEGG
            (cat {input.kegg_pathway}; curl --silent https://rest.kegg.jp/list/pathway) | {script_dir}/merge_function_descriptions.pl > {output.kegg_pathway}
            (cat {input.kegg_ko};      curl --silent https://rest.kegg.jp/list/ko)      | {script_dir}/merge_function_descriptions.pl > {output.kegg_ko}
            (cat {input.kegg_module};  curl --silent https://rest.kegg.jp/list/module)  | {script_dir}/merge_function_descriptions.pl > {output.kegg_module}

            echo 'KEGG descriptions downloaded'

            ## eggNOG

            (echo -e 'Funct\tCategory\tDescription';
             curl --silent http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/1/1_annotations.tsv.gz | gzip -dc | cut -f2-4 | grep -v '^COG';
             curl --silent https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab | cut -f1,2,3) > {output.eggnog_desc}

            echo 'eggNOG descriptions downloaded'

            ## EC numbers for dbCAN

            (echo -e 'Funct\tDescription';
             curl --silent https://rest.kegg.jp/list/enzyme;
             curl --silent https://ftp.expasy.org/databases/enzyme/enzclass.txt | {script_dir}/format_ec_classes.pl) > {output.EC_desc}

            echo 'EC descriptions downloaded'

        ) &> {log}
        """

###############################################################################################
# Download dbCAN database - run_dbcan
###############################################################################################

# Based on instructions at https://github.com/linnabrown/run_dbcan

rule dbCAN_db:
    output:
        dbCAN_db1="{minto_dir}/data/dbCAN_db/V12/CAZyDB.fa",
        dbCAN_db2="{minto_dir}/data/dbCAN_db/V12/dbCAN.txt",
        dbCAN_db3="{minto_dir}/data/dbCAN_db/V12/tcdb.fa",
        dbCAN_db4="{minto_dir}/data/dbCAN_db/V12/tf-1.hmm",
        dbCAN_db5="{minto_dir}/data/dbCAN_db/V12/tf-2.hmm",
        dbCAN_db6="{minto_dir}/data/dbCAN_db/V12/stp.hmm",
        dbCAN_db7="{minto_dir}/data/dbCAN_db/V12/fam-substrate-mapping.tsv"
    resources: mem=download_memory
    threads: download_threads
    log:
        "{minto_dir}/logs/dbCAN_db_download.log"
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {minto_dir}/data/dbCAN_db/V12
        cd {minto_dir}/data/dbCAN_db
        time (
        dbcan_build --cpus 2 --db-dir V12 --clean

        echo 'dbCAN database downloaded and installed') &> {log}
        """

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
# Download phylophlan database
###############################################################################################

rule download_phylophlan_db:
    output:
        "{minto_dir}/data/phylophlan/SGB.{phylophlan_db_version}.txt.bz2"
    shadow:
        "minimal"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/phylophlan.SGB.{phylophlan_db_version}.download.log"
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        time (\
            tar xvfz {minto_dir}/tutorial/genomes.tar.gz
            phylophlan_metagenomic --database_folder {minto_dir}/data/phylophlan -d SGB.{phylophlan_db_version} -i genomes -o tmp
            if [ $? -eq 0 ]; then
                echo 'phylophlan download: OK'
            else
                echo 'phylophlan download: FAIL'
            fi
            rm {minto_dir}/data/phylophlan/SGB.{phylophlan_db_version}.tar
            rm {minto_dir}/data/phylophlan/SGB.{phylophlan_db_version}.md5
            ) &> {log}
        """

###############################################################################################
# Download GTDB database
###############################################################################################

rule download_GTDB_db:
    output:
        "{minto_dir}/data/GTDB/r{gtdb_release_number}/taxonomy/gtdb_taxonomy.tsv"
    resources:
        mem=download_memory
    threads:
        download_threads
    log:
        "{minto_dir}/logs/GTDB.r{gtdb_release_number}.download.log"
    conda:
        config["minto_dir"]+"/envs/gtdb.yml"
    shell:
        """
        time (\
            mkdir -p {minto_dir}/data/GTDB/r{gtdb_release_number}
            cd {minto_dir}/data/GTDB
            wget -O gtdb.tar.gz https://data.gtdb.ecogenomic.org/releases/release{gtdb_release_number}/{gtdb_release_number}.0/auxillary_files/gtdbtk_r{gtdb_release_number}_data.tar.gz
            if [ $? -eq 0 ]; then
                echo 'GTDB download: OK'
            else
                echo 'GTDB download: FAIL'
            fi
            tar xfz gtdb.tar.gz
            mv release{gtdb_release_number}/* r{gtdb_release_number}/
            rmdir release{gtdb_release_number}
            rm gtdb.tar.gz
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
