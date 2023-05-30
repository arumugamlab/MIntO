#!/usr/bin/env python

'''
Gene prediction and functional annotation step

Authors: Vithiagaran Gunalan, Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, local_dir, minto_dir, script_dir, metadata
include: 'config_parser.smk'

if config['map_reference'] in ("MAG", "reference_genome"):
    map_reference=config["map_reference"]
else:
    print('ERROR in ', config_path, ': map_reference variable is not correct. "map_reference" variable should be MAG or reference_genome.')

if map_reference == 'MAG':
    reference_dir="{wd}/{omics}/8-1-binning/mags_generation_pipeline/prokka".format(wd=working_dir, omics=omics)
    print('WARNING in ', config_path, ': MIntO is using "' + reference_dir + '" as PATH_reference variable')
elif map_reference == 'reference_genome':
    reference_dir=config["PATH_reference"]
    print('WARNING in ', config_path, ': MIntO is using "'+ reference_dir+'" as PATH_reference variable')
elif config['PATH_reference'] is None:
    print('ERROR in ', config_path, ': PATH_reference variable is empty. Please, complete ', config_path)
elif path.exists(config['PATH_reference']) is False:
    print('ERROR in ', config_path, ': PATH_reference variable path does not exit. Please, complete ', config_path)

if 'ANNOTATION' in config:
    if config['ANNOTATION'] is None:
        print('ERROR in ', config_path, ': ANNOTATION list is empty. "ANNOTATION" variable should be dbCAN, KEGG and/or eggNOG. Please, complete ', config_path)
    else:
        try:
            # Make list of illumina samples, if ILLUMINA in config
            annot_list = list()
            if 'ANNOTATION' in config:
                #print("Samples:")
                for annot in config["ANNOTATION"]:
                    if annot in ('dbCAN', 'KEGG', 'eggNOG'):
                        annot_list.append(annot)
                        pass
                    else:
                        raise TypeError('ANNOTATION variable is not correct. "ANNOTATION" variable should be dbCAN, KEGG and/or eggNOG. Please, complete ', config_path)
        except: 
            print('ERROR in ', config_path, ': ANNOTATION variable is not correct. "ANNOTATION" variable should be dbCAN, KEGG and/or eggNOG.')
else:
    print('ERROR in ', config_path, ': ANNOTATION list is empty. "ANNOTATION" variable should be dbCAN, KEGG and/or eggNOG. Please, complete', config_path)

# Define all the outputs needed by target 'all'

if map_reference == 'MAG':
    post_analysis_dir="9-MAGs-prokka-post-analysis"
    post_analysis_out="MAGs_genes"
elif map_reference == 'reference_genome':
    post_analysis_dir="9-reference-genes-post-analysis"
    post_analysis_out="reference_genes"
elif map_reference == 'genes_db':
    post_analysis_dir="9-genes-db-post-analysis"
    post_analysis_out="db_genes"

def merge_genes_output():
    result = expand("{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa",
                wd = working_dir,
                post_analysis_dir = post_analysis_dir,
                post_analysis_out = post_analysis_out),\
    expand("{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.bed",
                wd = working_dir,
                post_analysis_dir = post_analysis_dir,
                post_analysis_out = post_analysis_out),\
    expand("{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}_names_modif.bed",
                wd = working_dir,
                post_analysis_dir = post_analysis_dir,
                post_analysis_out = post_analysis_out),\
    expand("{wd}/DB/{post_analysis_dir}/{post_analysis_out}_SUBSET.bed",
                wd = working_dir,
                post_analysis_dir = post_analysis_dir,
                post_analysis_out = post_analysis_out),\
    expand("{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
                wd = working_dir,
                post_analysis_dir = post_analysis_dir,
                post_analysis_out = post_analysis_out)
    return(result)

if map_reference == 'genes_db':
    def merge_genes_output(): # CHECK THIS PART - do not output anything when map_reference == 'genes_db'
        result = expand()
        return(result)

def predicted_genes_annot_out():
        result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_{annot}.tsv",
                        wd = working_dir,
                        post_analysis_dir = post_analysis_dir,
                        post_analysis_out = post_analysis_out,
                        annot = config["ANNOTATION"] if "ANNOTATION" in config else []),
        return(result)

def predicted_genes_collate_out():
    result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.annotations.tsv",
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)
    return(result)

rule all:
    input: 
        merge_genes_output(),
        predicted_genes_annot_out(),
        predicted_genes_collate_out()


###############################################################################################
# Prepare predicted genes/ publicly available genes for annotation
## Combine faa and bed files from the different genomes
###############################################################################################

# Get a sorted list of genomes

def get_genomes_from_refdir():
    genomes = [ f.name for f in os.scandir(reference_dir) if f.is_dir() ]
    return(sorted(genomes))

# Combine FAA files ####

rule make_cd_transl_faa_for_genome:
    input: lambda wildcards: "{reference_dir}/{genome}/{genome}.faa".format(reference_dir=reference_dir, genome=wildcards.genome)
    output: "{wd}/DB/{post_analysis_dir}/CD_transl/{genome}_translated_cds.faa"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.reformat_faa.log"
    params:
        pattern = '^>lcl|' if map_reference=='reference_genome' else '^>'
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #seqtk
    shell:
        """
        time (seqtk seq -l 0 {input} | sed -e "1,$ s/{params.pattern}/>{wildcards.genome}|/g" > {output}) >& {log}
        """

def get_genome_cd_transl_faa(wildcards):
    #Collect the CDS faa files for MAGs
    genomes = get_genomes_from_refdir()
    result = expand("{wd}/DB/{post_analysis_dir}/CD_transl/{genome}_translated_cds.faa",
                    wd=wildcards.wd,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    genome=genomes)
    return(result)

rule make_merged_cds_faa:
    input: get_genome_cd_transl_faa
    output:
        fasta_cd_transl_merge="{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa"
    shell:
        """
        cat {input} > {output}
        """

# Convert GFF file into BED file ####

rule make_bed_for_genome:
    input: lambda wildcards: "{reference_dir}/{genome}/{genome}.gff".format(reference_dir=reference_dir, genome=wildcards.genome)
    output:
        bed=                      "{wd}/DB/{post_analysis_dir}/GFF/{genome}.bed",
        bed_hdr_mod=              "{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.bed",
        bed_hdr_mod_coord_correct="{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.coord_correct.bed"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.make_bed_from_gff.log"
    params:
        pattern = '^' if map_reference=='reference_genome' else '^gnl|X|'
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #gff2bed
    shell:
        """
        time (\
            gff2bed < {input} > {output.bed}
            sed -e "1,$ s/{params.pattern}/{wildcards.genome}|/g" {output.bed} > {output.bed_hdr_mod}
            awk -v FS='\\t' -v OFS='\\t' '{{$2+=1; print $0}}' {output.bed_hdr_mod} > {output.bed_hdr_mod_coord_correct} \
        ) >& {log}
        """

def get_genome_bed(wildcards):
    #Collect the BED files for MAGs
    genomes = get_genomes_from_refdir()
    result = expand("{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.coord_correct.bed",
                    wd=wildcards.wd,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    genome=genomes)
    return(result)

rule make_merged_bed:
    input: get_genome_bed
    output:
        bed_file=       "{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.bed",
        bed_file_header="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.header-modif.coord_correct.bed"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}.merge_bed.log"
    wildcard_constraints:
        post_analysis_out='MAGs_genes|reference_genes'
    shell:
        """
        time (\
            cat {input} > {output.bed_file_header}
            cat {input} \
                    | awk -F'\\t' '{{gsub(/[;]/, "\\t", $10)}} 1' OFS='\\t' \
                    | cut -f 1-10 \
                    > {output.bed_file} \
        ) >& {log}
        """

rule gene_annot_subset:
    input:
        bed_file="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.bed",
        bed_file_header="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.header-modif.coord_correct.bed",
        fasta_cd_transl_merge="{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa",
    output: 
        bed_modif="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}_names_modif.bed",
        bed_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_SUBSET.bed", 
        fasta_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
    params:
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_subset_out.log"
    resources:
        mem=5
    threads: 8
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell: 
        """
        remote_dir={wildcards.wd}/DB/{post_analysis_dir}/
        mkdir -p ${{remote_dir}}/GFF
        mkdir -p ${{remote_dir}}/CD_transl
        time (Rscript {script_dir}/{post_analysis_out}_out_subset.R {threads} {input.bed_file_header} {input.bed_file} {input.fasta_cd_transl_merge} ${{remote_dir}} {post_analysis_out}) &> {log}
        """

################################################################################################
# Genes annotation using 3 different tools to retrieve functions from:
## KEGG
## eggNOG
## dbCAN
## Pfam
################################################################################################

rule gene_annot_kofamscan:
    input: 
        fasta_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
    output: 
        kofam_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_KEGG.tsv",
    params:
        tmp_kofamscan=lambda wildcards: "{local_dir}/{post_analysis_out}_kofamscan/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        kofam_db=lambda wildcards: "{minto_dir}/data/kofam_db/".format(minto_dir = minto_dir) #config["KOFAM_db"]
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_kofamscan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {params.tmp_kofamscan}tmp
        remote_dir=$(dirname {output.kofam_out})
        time (exec_annotation -k {params.kofam_db}/ko_list -p {params.kofam_db}/profiles/prokaryote.hal --tmp-dir {params.tmp_kofamscan}tmp --create-alignment -f mapper --cpu {threads} -o {params.tmp_kofamscan}{params.prefix}_kofam_mapper.txt {input.fasta_subset}
        {script_dir}/kofam_hits.pl {params.tmp_kofamscan}{params.prefix}_kofam_mapper.txt > {params.tmp_kofamscan}{params.prefix}_KEGG.tsv
        sed -i 's/\\;K/\\,K/g' {params.tmp_kofamscan}{params.prefix}_KEGG.tsv
        rm -rf {params.tmp_kofamscan}tmp
        mkdir -p ${{remote_dir}}
        rsync {params.tmp_kofamscan}{params.prefix}_KEGG.tsv {output.kofam_out}) &> {log}
        rm -rf {params.tmp_kofamscan}
        """

rule gene_annot_dbcan:
    input: 
        fasta_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
    output: 
        dbcan_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_dbCAN.tsv",
    params:
        tmp_dbcan=lambda wildcards: "{local_dir}/{post_analysis_out}_dbcan/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        dbcan_db=lambda wildcards: "{minto_dir}/data/dbCAN_db/".format(minto_dir = minto_dir)
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_dbcan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {params.tmp_dbcan}
        time (run_dbcan.py {input.fasta_subset} protein --db_dir {params.dbcan_db} --dia_cpu {threads} --out_pre {params.prefix}_ --out_dir {params.tmp_dbcan}
        {script_dir}/process_dbcan_overview.pl {params.tmp_dbcan}/{params.prefix}_overview.txt > {params.tmp_dbcan}/{params.prefix}_dbCAN.tsv
        rsync {params.tmp_dbcan}* {wildcards.wd}/DB/{post_analysis_dir}/annot/
        rm -rf Hotpep)&> {log}
        rm -rf {params.tmp_dbcan}
        """

rule gene_annot_eggnog:
    input: 
        fasta_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
    output: 
        eggnog_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_eggNOG.tsv",
    params:
        tmp_eggnog=lambda wildcards: "{local_dir}/{post_analysis_out}_eggnog/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        eggnog_db= lambda wildcards: "{minto_dir}/data/eggnog_data/".format(minto_dir = minto_dir) #config["EGGNOG_db"]
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_eggnog.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        mkdir -p {params.tmp_eggnog}
        time (cd {params.eggnog_db}
        [[ -d /dev/shm/eggnog_data ]] || mkdir /dev/shm/eggnog_data
        eggnogdata=(eggnog.db eggnog_proteins.dmnd eggnog.taxa.db eggnog.taxa.db.traverse.pkl)
        for e in "${{eggnogdata[@]}}"
            do [[ -f /dev/shm/eggnog_data/$e ]] || cp data/$e /dev/shm/eggnog_data/
        done
        emapper.py --data_dir /dev/shm/eggnog_data/ -o {params.prefix} \
--no_annot --no_file_comments --report_no_hits --override --output_dir {params.tmp_eggnog} -m diamond -i {input.fasta_subset} --cpu {threads}
        emapper.py --annotate_hits_table {params.tmp_eggnog}/{params.prefix}.emapper.seed_orthologs \
--data_dir /dev/shm/eggnog_data/ -m no_search --no_file_comments --override -o {params.prefix} --output_dir {params.tmp_eggnog} --cpu {threads}
        cut -f 1,5,12,13,14,21 {params.tmp_eggnog}/{params.prefix}.emapper.annotations > {params.tmp_eggnog}/{params.prefix}.eggNOG5
        {script_dir}/process_eggNOG_OGs.pl {params.tmp_eggnog}/{params.prefix}.eggNOG5 > {params.tmp_eggnog}/{params.prefix}_eggNOG.tsv
        rm {params.tmp_eggnog}/{params.prefix}.eggNOG5
        sed -i 's/\#query/ID/' {params.tmp_eggnog}/{params.prefix}_eggNOG.tsv
        sed -i 's/ko\://g' {params.tmp_eggnog}/{params.prefix}_eggNOG.tsv
        rsync {params.tmp_eggnog}* {wildcards.wd}/DB/{post_analysis_dir}/annot/)&> {log}
        rm -rf {params.tmp_eggnog}
        """

################################################################################################
# Combine gene annotation results into one file
################################################################################################

rule predicted_gene_annotation_collate:
    input: 
        #kofam_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_kofam.tsv",
        #dbcan_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_dbCAN.tsv",
        #eggnog_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.eggNOG5.annotations",
        annot_out=expand("{{wd}}/DB/{{post_analysis_dir}}/annot/{{post_analysis_out}}_translated_cds_SUBSET_{annot}.tsv", annot=annot_list)
    output: 
        annot_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.annotations.tsv",
    params:
        #tmp_annot=lambda wildcards: "{local_dir}/{post_analysis_out}_annot/".format(local_dir=local_dir),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_annot.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/mags.yml" # python with pandas
    shell:
        """
        time (python3 {script_dir}/collate.py {wildcards.wd}/DB/{post_analysis_dir}/annot/ {params.prefix})&> {log}
        """
