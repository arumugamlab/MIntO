#!/usr/bin/env python

'''
Gene prediction and functional annotation step

Authors: Vithiagaran Gunalan, Carmen Saenz, Mani Arumugam
'''

# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/gene_annotation.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --use-conda --conda-prefix /data/rzv923/MIntO_snakemake_env/ --configfile mapping.smk.yaml --jobs 5
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/gene_annotation.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A"  --use-conda --conda-prefix /data/rzv923/MIntO_snakemake_env/ --configfile mapping.smk.yaml --jobs 5

# configuration yaml file
# import sys
import os.path
from os import path

# args = sys.argv
# config_path = args[args.index("--configfile") + 1]
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file

# some variables
if config['PROJECT'] is None:
    print('ERROR in ', config_path, ': PROJECT variable is empty. Please, complete ', config_path)
else:
    project_name = config['PROJECT']

if config['working_dir'] is None:
    print('ERROR in ', config_path, ': working_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['working_dir']) is False:
    print('ERROR in ', config_path, ': working_dir variable path does not exit. Please, complete ', config_path)
else:
    working_dir = config['working_dir']

if config['local_dir'] is None:
    prints('ERROR in ', config_path, ': local_dir variable is empty. Please, complete ', config_path)
else:
    local_dir = config['local_dir']

if config['minto_dir'] is None:
    print('ERROR in ', config_path, ': minto_dir variable in configuration yaml file is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]
    script_dir=config["minto_dir"]+"/scripts"

if config['map_reference'] in ("MAG", "reference_genome"):
    map_reference=config["map_reference"]
else:
    print('ERROR in ', config_path, ': map_reference variable is not correct. "map_reference" variable should be MAG or reference_genome.')

if map_reference == 'MAG':
    print('WARNING in ', config_path, ': MIntO is using "'+ working_dir+'/metaG/8-1-binning/mags_generation_pipeline/prokka" as PATH_reference variable')
    reference_dir="{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka".format(wd=working_dir)
elif map_reference == 'reference_genome':
    reference_dir=config["PATH_reference"]
    print('WARNING in ', config_path, ': MIntO is using "'+ reference_dir+'" as PATH_reference variable')
elif config['PATH_reference'] is None:
    print('ERROR in ', config_path, ': PATH_reference variable is empty. Please, complete ', config_path)
elif path.exists(config['PATH_reference']) is False:
    print('ERROR in ', config_path, ': PATH_reference variable path does not exit. Please, complete ', config_path)


# Define all the outputs needed by target 'all'

if map_reference == 'MAG':
    post_analysis_dir="9-MAGs-prokka-post-analysis"
    post_analysis_out="MAGs_genes"
    #reference_dir="{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka".format(wd=working_dir)
    def merge_genes_output():
        result = expand("{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa", 
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out),\
        expand("{wd}/DB/{post_analysis_dir}/genomes_list.txt", 
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir),\
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

if map_reference == 'reference_genome':
    post_analysis_dir="9-reference-genes-post-analysis"
    post_analysis_out="reference_genes"
    #reference_dir=config["reference_dir"]
    def merge_genes_output():
        result = expand("{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa", 
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out),\
        expand("{wd}/DB/{post_analysis_dir}/genomes_list.txt", 
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir),\
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
    post_analysis_dir="9-genes-db-post-analysis"
    post_analysis_out="db_genes"
    def merge_genes_output(): # CHECK THIS PART - do not output anything when map_reference == 'genes_db'
        result = expand()
        return(result)

def predicted_genes_kofam_out():
    result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_kofam.tsv",
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)
    return(result)

def predicted_genes_dbcan_out():
    result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_dbCAN.tsv",
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)
    return(result)

def predicted_genes_eggnog_out():
    result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.eggNOG5.annotations",
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)
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
        predicted_genes_dbcan_out(),
        predicted_genes_kofam_out(),
        predicted_genes_eggnog_out(),
        predicted_genes_collate_out()


###############################################################################################
# Prepare predicted genes/ publicly available genes for annotation
## Combine faa and bed files from the different genomes
###############################################################################################

# Convert GFF file into BED file ####
rule gene_annot_merge:
    input:
        in_dir="{reference_dir}/".format(reference_dir=reference_dir),
        #prova="{wd}"#.format(wd = working_dir)
    output:
        fasta_cd_transl_merge="{{wd}}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa".format(post_analysis_dir = post_analysis_dir, post_analysis_out = post_analysis_out), 
        genomes_list="{{wd}}/DB/{post_analysis_dir}/genomes_list.txt".format(post_analysis_dir = post_analysis_dir), 
        bed_file="{{wd}}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.bed".format(post_analysis_dir = post_analysis_dir, post_analysis_out = post_analysis_out),
        bed_file_header="{{wd}}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.header-modif.coord_correct.bed".format(post_analysis_dir = post_analysis_dir, post_analysis_out = post_analysis_out),
    params:
       tmp_reference=lambda wildcards: "{local_dir}/{post_analysis_dir}_{post_analysis_out}_out/".format(local_dir = local_dir, post_analysis_dir = post_analysis_dir, post_analysis_out = post_analysis_out),
    log:
        "{{wd}}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_out.log".format(post_analysis_dir = post_analysis_dir, post_analysis_out = post_analysis_out)
    resources:
        mem=5
    threads: 4
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #gff2bed
    shell: 
        """
        mkdir -p {params.tmp_reference}CD_transl
        mkdir -p {params.tmp_reference}GFF
        time( cd {reference_dir}
        if [ {map_reference} == 'MAG' ]
            then echo {map_reference}
            for file in */*.faa
                do echo ${{file}}
                genome=$(dirname "${{file}}")
                echo ${{genome}}
                echo ${{genome}} >> {params.tmp_reference}genomes_list.txt
                awk '{{if(NR==1) {{print $0}} else {{if($0 ~ /^>/) {{print "\\n"$0}} else {{printf $0}}}}}}' ${{file}} > {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
                [[ $(tail -c1 {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa) && -f {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa ]]&&echo ''>> {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
                sed -e "1,$ s/^>/>${{genome}}|/g" "{params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa" > {params.tmp_reference}CD_transl/${{genome}}_translated_cds.faa
                rm {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
            done
            for file in {reference_dir}/*/*.gff
                do echo ${{file}}
                genome=$(basename "${{file}}" .gff)
                echo ${{genome}}
                gff2bed < ${{file}} > {params.tmp_reference}GFF/${{genome}}.bed
                file_bed={params.tmp_reference}GFF/${{genome}}.bed
                sed -e "1,$ s/^gnl|X|/${{genome}}|/g" "${{file_bed}}" > {params.tmp_reference}GFF/${{genome}}.header-modif.bed
            done
        elif [ {map_reference} == 'reference_genome' ]
            then echo {map_reference} 
            for file in */*.faa
                do echo ${{file}}
                genome=$(dirname "${{file}}")
                echo ${{genome}}
                echo ${{genome}} >> {params.tmp_reference}genomes_list.txt
                awk '{{if(NR==1) {{print $0}} else {{if($0 ~ /^>/) {{print "\\n"$0}} else {{printf $0}}}}}}' ${{file}} > {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
                [[ $(tail -c1 {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa) && -f {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa ]]&&echo ''>> {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
                sed -e "1,$ s/^>lcl|/>${{genome}}|/g" "{params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa" > {params.tmp_reference}CD_transl/${{genome}}_translated_cds.faa
                rm {params.tmp_reference}CD_transl/${{genome}}_translated_cds_int.faa
            done      
            for file in {reference_dir}/*/*.gff
                do echo ${{file}}
                genome=$(basename "${{file}}" .gff)
                echo ${{genome}}
                gff2bed < ${{file}} > {params.tmp_reference}GFF/${{genome}}.bed
                file_bed={params.tmp_reference}GFF/${{genome}}.bed
                sed -e "1,$ s/^/${{genome}}|/g" "${{file_bed}}" > {params.tmp_reference}GFF/${{genome}}.header-modif.bed
            done
        fi
        for g in $(cat {params.tmp_reference}genomes_list.txt)
            do cat {params.tmp_reference}CD_transl/${{g}}_translated_cds.faa >> {params.tmp_reference}{post_analysis_out}_translated_cds.faa
            cat {params.tmp_reference}GFF/${{g}}.header-modif.bed >> {params.tmp_reference}GFF/{post_analysis_out}.header-modif.bed
        done
        awk -v OFS='\\t' '{{$2+=1; print $0}}' {params.tmp_reference}GFF/{post_analysis_out}.header-modif.bed > {params.tmp_reference}GFF/{post_analysis_out}.header-modif.coord_correct.bed
        awk -F'\\t' '{{gsub(/[;]/, "\\t", $10)}} 1' OFS='\\t' {params.tmp_reference}GFF/{post_analysis_out}.header-modif.coord_correct.bed| cut -f 1-10  > {params.tmp_reference}{post_analysis_out}.bed
        rsync {params.tmp_reference}CD_transl/* {working_dir}/DB/{post_analysis_dir}/CD_transl/
        rsync {params.tmp_reference}{post_analysis_out}_translated_cds.faa {output.fasta_cd_transl_merge}
        rsync {params.tmp_reference}genomes_list.txt {output.genomes_list}
        rsync {params.tmp_reference}GFF/* {working_dir}/DB/{post_analysis_dir}/GFF/
        rsync {params.tmp_reference}{post_analysis_out}.bed {output.bed_file}) &> {log}
        rm -rf {params.tmp_reference}"""

rule gene_annot_subset:
    input:
        bed_file="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.bed",
        bed_file_header="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}.header-modif.coord_correct.bed",
        fasta_cd_transl_merge="{wd}/DB/{post_analysis_dir}/CD_transl/{post_analysis_out}_translated_cds.faa",
        genomes_list="{wd}/DB/{post_analysis_dir}/genomes_list.txt", 
    output: 
        bed_modif="{wd}/DB/{post_analysis_dir}/GFF/{post_analysis_out}_names_modif.bed",
        bed_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_SUBSET.bed", 
        fasta_subset="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_translated_cds_SUBSET.faa",
    params:
        tmp_subset=lambda wildcards: "{local_dir}/{post_analysis_out}_out_subset/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
    log:
        "{wd}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_subset_out.log"
    resources:
        mem=5
    threads: 8
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml"
    shell: 
        """
        mkdir -p {params.tmp_subset}/GFF
        mkdir -p {params.tmp_subset}/CD_transl
        remote_dir={wildcards.wd}/DB/{post_analysis_dir}/
        time (Rscript {script_dir}/{post_analysis_out}_out_subset.R {threads} {input.bed_file_header} {input.bed_file} {input.fasta_cd_transl_merge} {params.tmp_subset} {post_analysis_out} ${{remote_dir}}
        if [ {map_reference} == 'reference_genome' ]
            then rsync {params.tmp_subset}/CD_transl/* ${{remote_dir}}CD_transl/
        fi
        rsync {params.tmp_subset}/GFF/* ${{remote_dir}}GFF/
        rsync {params.tmp_subset}/{post_analysis_out}_SUBSET.bed {output.bed_subset}
        rsync {params.tmp_subset}/{post_analysis_out}_translated_cds_SUBSET.faa {output.fasta_subset}) &> {log}
        rm -rf {params.tmp_subset}
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
        kofam_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_kofam.tsv",
    params:
        tmp_kofamscan=lambda wildcards: "{local_dir}/{post_analysis_out}_kofamscan/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        kofam_db=lambda wildcards: "{minto_dir}/data/kofam_db/".format(minto_dir = minto_dir) #config["KOFAM_db"]
    log:
        "{wd}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_kofamscan.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/kofamscan_env.yml" #config["kofamscan_ironmenv"]
    shell:
        """
        mkdir -p {params.tmp_kofamscan}tmp
        remote_dir=$(dirname {output.kofam_out})
        time (exec_annotation -k {params.kofam_db}/ko_list -p {params.kofam_db}/profiles/prokaryote.hal --tmp-dir {params.tmp_kofamscan}tmp --create-alignment -f mapper --cpu {threads} -o {params.tmp_kofamscan}{params.prefix}_kofam_mapper.txt {input.fasta_subset}
        {script_dir}/kofam_hits.pl {params.tmp_kofamscan}{params.prefix}_kofam_mapper.txt > {params.tmp_kofamscan}{params.prefix}_kofam.tsv
        sed -i 's/\\;K/\\,K/g' {params.tmp_kofamscan}{params.prefix}_kofam.tsv
        rm -rf {params.tmp_kofamscan}tmp
        mkdir -p ${{remote_dir}}
        rsync {params.tmp_kofamscan}{params.prefix}_kofam.tsv {output.kofam_out}) &> {log}
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
        "{wd}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_dbcan.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/dbcan_env.yml" #config["dbcan_ironmenv"]
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
        eggnog_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.eggNOG5.annotations",
    params:
        tmp_eggnog=lambda wildcards: "{local_dir}/{post_analysis_out}_eggnog/".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        eggnog_db= lambda wildcards: "{minto_dir}/data/eggnog_data/".format(minto_dir = minto_dir) #config["EGGNOG_db"]
    log:
        "{wd}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_eggnog.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/py38_env.yml" #config["py38_ironmenv"]
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
        {script_dir}/process_eggNOG_OGs.pl {params.tmp_eggnog}/{params.prefix}.eggNOG5 > {params.tmp_eggnog}/{params.prefix}.eggNOG5.annotations
        rm {params.tmp_eggnog}/{params.prefix}.eggNOG5
        sed -i 's/\#query/ID/' {params.tmp_eggnog}/{params.prefix}.eggNOG5.annotations
        sed -i 's/ko\://g' {params.tmp_eggnog}/{params.prefix}.eggNOG5.annotations
        rm -rf /dev/shm/eggnog_data
        rsync {params.tmp_eggnog}* {wildcards.wd}/DB/{post_analysis_dir}/annot/)&> {log}
        rm -rf {params.tmp_eggnog}
        """

################################################################################################
# Combine gene annotation results into one file
################################################################################################

rule predicted_gene_annotation_collate:
    input: 
        kofam_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_kofam.tsv",
        dbcan_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_dbCAN.tsv",
        eggnog_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.eggNOG5.annotations",
    output: 
        annot_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.annotations.tsv",
    params:
        #tmp_annot=lambda wildcards: "{local_dir}/{post_analysis_out}_annot/".format(local_dir=local_dir),
        prefix="{post_analysis_out}_translated_cds_SUBSET",
    log:
        "{wd}/logs/metaG/{post_analysis_dir}/{post_analysis_out}_annot.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/py36_env.yml" #config["py36_ironmenv"]
    shell:
        """
        time (python3 {script_dir}/collate.py {wildcards.wd}/DB/{post_analysis_dir}/annot/ {params.prefix})&> {log}
        """