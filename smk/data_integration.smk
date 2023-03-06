#!/usr/bin/env python

'''
Gene and function profiling step

Authors: Carmen Saenz
'''

# configuration yaml file
# import sys
import os.path
from os import path

# args = sys.argv
# config_path = args[args.index("--configfile") + 1]
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file")#: ", config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file

# some variables
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

if config['omics'] in ('metaG','metaT', 'metaG_metaT'):
    omics = config['omics']
else:
    print('ERROR in ', config_path, ': omics variable is not correct. "omics" variable should be metaG, metaT or metaG_metaT.')

if config['minto_dir'] is None:
    print('ERROR in ', config_path, ': minto_dir variable in configuration yaml file is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]
    script_dir=config["minto_dir"]+"/scripts"

if config['METADATA'] is None:
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata_file=config["METADATA"]
elif config['METADATA'] == "None":
    print('WARNING in ', config_path, ': METADATA variable is empty. Samples will be analyzed excluding the metadata.')
    metadata_file=config["METADATA"]
elif path.exists(config['METADATA']) is False:
    print('ERROR in ', config_path, ': METADATA variable path does not exit. Please, complete ', config_path)
else:
    metadata_file=config["METADATA"]


if config['map_reference'] in ("MAG", "reference_genome","genes_db"):
    map_reference=config["map_reference"]
else:
    print('ERROR in ', config_path, ': map_reference variable is not correct. "map_reference" variable should be MAG, reference_genome or genes_db.')

if config['abundance_normalization'] in ("MG", "TPM"):
    normalization=config['abundance_normalization']
else:
    print('ERROR in ', config_path, ': abundance_normalization variable is not correct. "abundance_normalization" variable should be MG or TPM.')

if normalization == 'MG' and map_reference in ("genes_db"):
    print('ERROR in ', config_path, ': In "genes_db" mode, TPM nomralization is only allowed.')


if config['alignment_identity'] is None:
    print('ERROR in ', config_path, ': alignment_identity variable is empty. Please, complete ', config_path)
elif type(config['alignment_identity']) != int:
    print('ERROR in ', config_path, ': alignment_identity variable is not an integer. Please, complete ', config_path)
elif type(config['alignment_identity']) == int:
    identity=config['alignment_identity']

if config['MERGE_memory'] is None:
    print('ERROR in ', config_path, ': MERGE_memory variable is empty. Please, complete ', config_path)
elif type(config['MERGE_memory']) != int:
    print('ERROR in ', config_path, ': MERGE_memory variable is not an integer. Please, complete ', config_path)

if config['MERGE_threads'] is None:
    print('ERROR in ', config_path, ': MERGE_threads variable is empty. Please, complete ', config_path)
elif type(config['MERGE_threads']) != int:
    print('ERROR in ', config_path, ': MERGE_threads variable is not an integer. Please, complete ', config_path)

if map_reference == 'genes_db':
    if config['ANNOTATION_file'] is None:
        print('ERROR in ', config_path, ': ANNOTATION_file variable in configuration yaml file is empty. Please, complete ', config_path)
    elif path.exists(config['ANNOTATION_file']) is False:
        print('ERROR in ', config_path, ': ANNOTATION_file variable path does not exit. Please, complete ', config_path)
    elif path.exists(config['ANNOTATION_file']) is True:
        annot_file=config['ANNOTATION_file']
elif map_reference == 'MAG':
    print('WARNING in ', config_path, ': MIntO is using "'+ working_dir+'/DB/9-MAGs-prokka-post-analysis/annot/MAGs_genes_translated_cds_SUBSET.annotations.tsv" as ANNOTATION_file variable.')
elif map_reference == 'reference_genome':
    print('WARNING in ', config_path, ': MIntO is using "'+ working_dir+'/DB/9-reference-genes-post-analysis/annot/reference_genes_translated_cds_SUBSET.annotations.tsv" as ANNOTATION_file variable.')
    "{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv"

if map_reference in ('MAG', 'reference_genome'):
    #print('WARNING in ', config_path, ': MIntO is using eggNOG_OGs, KEGG_Pathway, KEGG_Module, KEGG_KO, PFAMs, dbCAN.mod and dbCAN.enzclass as ANNOTATION_ids variable.')
    if config['ANNOTATION_ids'] is None:
        print('ERROR in ', config_path, ': ANNOTATION_ids variable in configuration yaml file is empty. Please, complete ', config_path)
    else:
        funct_opt=config['ANNOTATION_ids']
elif map_reference in ('genes_db'):
    if config['ANNOTATION_ids'] is None:
        print('ERROR in ', config_path, ': ANNOTATION_ids variable in configuration yaml file is empty. Please, complete ', config_path)
    else:
        funct_opt=config['ANNOTATION_ids']


if omics == 'metaG':
    omics_opt='metaG'
    omics_prof='A'
    #gene_abund="{wd}/{omics_opt}/6-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics_opt = omics_opt, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    #gene_abund_file="metaG/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"
elif omics == 'metaT':
    omics_opt='metaT'
    omics_prof='T'
    #gene_abund_file="metaT/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"
elif omics == 'metaG_metaT':
    omics_opt=['metaG','metaT']
    omics_prof=['A','T','E']
    #gene_abund_file="metaG/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
    #"metaT/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"

if map_reference == 'MAG':
    post_analysis_dir="9-MAGs-prokka-post-analysis"
    post_analysis_out="MAGs_genes"
    post_analysis_genome="MAGs_genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    #funct_opt=('eggNOG_OGs','KEGG_Pathway','KEGG_Module','KEGG_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')
    #funct_opt_list = ','.join(['"' + id + '"' for id in ('eggNOG_OGs','KEGG_Pathway','KEGG_Module','KEGG_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')])
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])
elif map_reference == 'reference_genome':
    post_analysis_dir="9-reference-genes-post-analysis"
    post_analysis_out="reference_genes"
    post_analysis_genome="reference_genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    #funct_opt=('eggNOG_OGs','KEGG_Pathway','KEGG_Module','KEGG_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')
    #funct_opt_list = ','.join(['"' + id + '"' for id in ('eggNOG_OGs','KEGG_Pathway','KEGG_Module','KEGG_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')])
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])
elif map_reference == 'genes_db':
    post_analysis_dir="9-genes-db-post-analysis"
    post_analysis_out="db_genes"
    post_analysis_genome="None"
    #annot_file=config['ANNOTATION_file']
    #funct_opt=config['ANNOTATION_ids']
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])

if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
    post_analysis_TPM=post_analysis_out
    post_analysis_other='None'
else:
    post_analysis_TPM="None"
    post_analysis_other=post_analysis_out

#for id in funct_opt:
#            print(id)
#            funct_opt_list.append(id)
#print(funct_opt_list)
# if normalization == 'MG' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
#     post_analysis_MG=post_analysis_out
# else:
#     post_analysis_MG="None"

def integration_merge_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization)
    return(result)

def integration_gene_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_prof}.csv",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_prof}.PCA.pdf",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof)
    return(result)

def integration_function_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.csv",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCA.pdf",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt)
    return(result)

# Define all the outputs needed by target 'all'
rule all:
    input:
        integration_merge_profiles(),
        integration_gene_profiles(),
        #integration_function_profiles()


###############################################################################################
# Prepare gene profile
## Merge gene abundance and gene transcript profiles
###############################################################################################
rule integration_merge_profiles:
    input:
        gene_abund = lambda wildcards: expand("{wd}/{omics_individual}/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
                                            wd = wildcards.wd,
                                            omics_individual = wildcards.omics.split("_"),
                                            map_reference = wildcards.map_reference,
                                            identity = wildcards.identity,
                                            normalization = wildcards.normalization),
    output:
        gene_abund_merge="{wd}/output/data_integration/{map_reference}/{omics}.genes_abundances.p{identity}.{normalization}.csv"
    log:
        "{wd}/logs/output/data_integration/{map_reference}/{omics}.p{identity}.{normalization}.integration_merge_profiles.log"
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/mags.yml" # python with pandas
    shell:
        """
        remote_dir=$(dirname {output.gene_abund_merge})
        time (if [ {omics} == 'metaG_metaT' ]
        then python3 {script_dir}/gene_abundances_merge_profiles.py {wildcards.map_reference} {output.gene_abund_merge} {input.gene_abund}
        else
        rsync {input.gene_abund} {output.gene_abund_merge}
        fi ) &> {log}
        """

###############################################################################################
# Generate gene expression profile
###############################################################################################
# Change names to R scripts
# what do to with annot file
rule integration_gene_profiles:
    input:
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
        #gene_abund="{input_dir}/".format(input_dir=gene_abund_file),
        #gene_abund="{wd}/{omics}/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
    output:
        gene_abund_prof=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_prof}.csv", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
        gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
        gene_abund_plots=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_prof}.PCA.pdf", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
    params:
        annot_file={annot_file},
        funct_opt= funct_opt_list
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_gene_profiles.{omics}.{normalization}.log".format(wd = working_dir, post_analysis_out = post_analysis_out, normalization = normalization, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time ( echo 'integration of gene profiles'
        if [[ {map_reference} == 'MAG' ]] || [[ {map_reference} == 'reference_genome' ]]
        then Rscript {script_dir}/gene_expression_profile_genome_based.R {threads} $(dirname {output.gene_abund_prof[0]}) {omics} {params.annot_file} {metadata_file} {input.gene_abund_merge} {params.funct_opt}
        elif [[ {map_reference} == 'genes_db' ]]
        then Rscript {script_dir}/gene_expression_profile_gene_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {params.annot_file} {metadata_file} {input.gene_abund_merge} {params.funct_opt}
        fi ) &> {log}
        """

##################################################################################################
# Generate function expression profile for the different databases included in the annotation file
##################################################################################################
rule merge_absolute_counts_TPM:
    input:
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
        absolute_counts = expand("{wd}/{omics_opt}/6-mapping-profiles/BWA_reads-{post_analysis_TPM}/genes_abundances.p{identity}.bed", wd = working_dir, omics_opt = omics_opt, post_analysis_TPM = post_analysis_TPM, identity = identity),
    output:
        absolute_counts_merge="{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.bed".format(wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity),
    log:
        "{wd}/logs/output/data_integration/{post_analysis_TPM}/merge_raw_counts.{omics}.TPM.log".format(wd = working_dir, post_analysis_TPM = post_analysis_TPM, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/mags.yml" # python with pandas
    shell:
        """
        time (python3 {script_dir}/gene_abundances_merge_raw_profiles.py {threads} {resources.mem} {working_dir} {omics} {post_analysis_TPM} {normalization} {identity}) &> {log}
        """

rule integration_function_profiles_TPM:
    input:
        absolute_counts_merge="{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.bed".format(wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity),
        gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, normalization = normalization),
    output:
        func_abund_prof=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/F{omics_prof}.{funct_opt}.csv", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
        func_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
        func_abund_plots=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/plots/F{omics_prof}.{funct_opt}.PCA.pdf", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
    params:
        annot_file={annot_file},
        funct_opt=funct_opt_list,
        mapped_reads_threshold=config["MIN_mapped_reads"]
    log:
        "{wd}/logs/output/data_integration/{post_analysis_TPM}/integration_function_profiles.{omics}.TPM.log".format(wd = working_dir, post_analysis_TPM = post_analysis_TPM, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        time (Rscript {script_dir}/function_expression_profile_genome_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_TPM} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}) &> {log}
        """

rule integration_function_profiles_MG_TPM:
        input:
            gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization),
        output:
            func_abund_prof=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.csv", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
            func_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
            func_abund_plots=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCA.pdf", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
        params:
            annot_file={annot_file},
            funct_opt=funct_opt_list,
            mapped_reads_threshold=config["MIN_mapped_reads"]
        log:
            "{wd}/logs/output/data_integration/{post_analysis_other}/integration_funtion_profiles.{omics}.{normalization}.log".format(wd = working_dir, post_analysis_other = post_analysis_other, normalization = normalization, omics = omics),
        resources:
            mem=config["MERGE_memory"]
        threads: config["MERGE_threads"]
        conda:
            config["minto_dir"]+"/envs/r_pkgs.yml" #R
        shell:
            """
            time ( echo 'integration of function profiles'
            if ( [[ {map_reference} == 'MAG' ]] || [[ {map_reference} == 'reference_genome' ]] ) && [[ {normalization} == 'MG' ]]
            then Rscript {script_dir}/function_expression_profile_genome_MG_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}
            elif [[ {map_reference} == 'genes_db' ]] && [[ {normalization} == 'TPM' ]]
            then Rscript {script_dir}/function_expression_profile_gene_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}
            fi) &> {log}
            """
