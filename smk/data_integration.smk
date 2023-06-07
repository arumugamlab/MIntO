#!/usr/bin/env python

'''
Gene and function profiling step

Authors: Carmen Saenz
'''

# configuration yaml file
# import sys
import os.path
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, local_dir, minto_dir, script_dir, metadata
include: 'config_parser.smk'

# some variables

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
elif omics == 'metaT':
    omics_opt='metaT'
    omics_prof='T'
elif omics == 'metaG_metaT':
    omics_opt=['metaG','metaT']
    omics_prof=['A','T','E']

if map_reference == 'MAG':
    post_analysis_dir="9-MAG-genes-post-analysis"
    post_analysis_out="MAG-genes"
    post_analysis_genome="MAG-genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])
elif map_reference == 'reference_genome':
    post_analysis_dir="9-refgenome-genes-post-analysis"
    post_analysis_out="refgenome_genes"
    post_analysis_genome="refgenome_genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])
elif map_reference == 'genes_db':
    post_analysis_dir="9-db-genes-post-analysis"
    post_analysis_out="db-genes"
    post_analysis_genome="None"
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])

print('WARNING in ', config_path, ': MIntO is using "' + annot_file + '" as ANNOTATION_file variable.')

if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
    post_analysis_TPM=post_analysis_out
    post_analysis_other='None'
else:
    post_analysis_TPM="None"
    post_analysis_other=post_analysis_out

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
        integration_function_profiles()


###############################################################################################
# Prepare gene profile
## Merge gene abundance and gene transcript profiles
###############################################################################################
rule integration_merge_profiles:
    input:
        gene_abund = lambda wildcards: expand("{wd}/{omics_individual}/9-mapping-profiles/{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
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
        time (if [ {wildcards.omics} == 'metaG_metaT' ]
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
        #gene_abund="{wd}/{omics}/9-mapping-profiles/{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
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
        then Rscript {script_dir}/gene_expression_profile_genome_based.R {threads} $(dirname {output.gene_abund_prof[0]}) {omics} {params.annot_file} {metadata} {input.gene_abund_merge} {params.funct_opt}
        elif [[ {map_reference} == 'genes_db' ]]
        then Rscript {script_dir}/gene_expression_profile_gene_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {params.annot_file} {metadata} {input.gene_abund_merge} {params.funct_opt}
        fi ) &> {log}
        """

##################################################################################################
# Generate function expression profile for the different databases included in the annotation file
##################################################################################################
rule merge_absolute_counts_TPM:
    input:
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
        absolute_counts = expand("{wd}/{omics_opt}/9-mapping-profiles/{post_analysis_TPM}/genes_abundances.p{identity}.bed", wd = working_dir, omics_opt = omics_opt, post_analysis_TPM = post_analysis_TPM, identity = identity),
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
        time (Rscript {script_dir}/function_expression_profile_genome_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_TPM} {normalization} {identity} {params.annot_file} {metadata} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}) &> {log}
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
            then Rscript {script_dir}/function_expression_profile_genome_MG_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}
            elif [[ {map_reference} == 'genes_db' ]] && [[ {normalization} == 'TPM' ]]
            then Rscript {script_dir}/function_expression_profile_gene_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold}
            fi) &> {log}
            """
