#!/usr/bin/env python

'''
Gene and function profiling step

Authors: Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

# some variables

main_factor = None
if config['MAIN_factor'] is not None:
    main_factor = config['MAIN_factor']

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
        raise Exception("Gene functional annotation needs to be provided via ANNOTATION_file variable")
    elif path.exists(config['ANNOTATION_file']) is False:
        raise Exception("File specified in ANNOTATION_file variable does not exist")
    elif path.exists(config['ANNOTATION_file']) is True:
        annot_file=config['ANNOTATION_file']

if map_reference in ('MAG', 'reference_genome'):
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
elif map_reference == 'reference_genome':
    post_analysis_dir="9-refgenome-genes-post-analysis"
    post_analysis_out="refgenome-genes"
    post_analysis_genome="refgenome-genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
elif map_reference == 'genes_db':
    post_analysis_dir="9-db-genes-post-analysis"
    post_analysis_out="db-genes"
    post_analysis_genome="None"

print('NOTE: MIntO is using "' + annot_file + '" as ANNOTATION_file variable.')

funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])
print('NOTE: MIntO is using ', funct_opt, ' as ANNOTATION_ids variable.')

if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
    post_analysis_TPM=post_analysis_out
else:
    post_analysis_TPM="None"

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
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.tsv",
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

# Function to column-wise combine profiles, after checking whether they all agree on values in key_columns.
# TODO: This is duplicated between here and gene_abundance.smk. Should be moved to 'include/' directory.

def combine_profiles(input_list, output_file, log_file, key_columns):
    import pandas as pd
    import hashlib
    import pickle
    import datetime

    def logme(stream, msg):
        print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), msg, file=stream)

    df_list = list()
    with open(str(log_file), 'w') as f:
        logme(f, "INFO: reading file 0")
        df = pd.read_csv(input_list[0], comment='#', header=0, sep = "\t", memory_map=True)
        df_list.append(df)
        md5_first = hashlib.md5(pickle.dumps(df[key_columns])).hexdigest() # make hash for sequence feature ID
        for i in range(1, len(input_list)):
            logme(f, "INFO: reading file {}".format(i))
            df = pd.read_csv(input_list[i], comment='#', header=0, sep = "\t", memory_map=True)
            md5_next = hashlib.md5(pickle.dumps(df[key_columns])).hexdigest()
            if md5_next != md5_first:
                raise Exception("combine_profiles: Features don't match between {} and {}, using keys {}".format(input_list[0], input_list[i], ",".join(key_columns)))
            df_list.append(df.drop(columns=key_columns))
        logme(f, "INFO: concatenating {} files".format(len(input_list)))
        df = pd.concat(df_list, axis=1, ignore_index=False, copy=False, sort=False)
        logme(f, "INFO: writing to output file")
        if output_file.endswith('.gz'):
            df.to_csv(output_file, sep = "\t", index = False, compression={'method': 'gzip', 'compresslevel': 1})
        else:
            df.to_csv(output_file, sep = "\t", index = False)
        logme(f, "INFO: done")

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
    threads: 1
    run:
        import shutil

        # Set key columns to merge metaG and metaT profiles
        if (wildcards.map_reference == 'db-genes'):
            key_columns = ['ID']
        else:
            key_columns = ['coord', 'chr', 'start', 'stop', 'name', 'score', 'strand', 'source', 'feature', 'frame', 'info', 'gene_length']

        # Merge metaG and metaT profiles
        if (wildcards.omics == 'metaG_metaT'):
            combine_profiles(input.gene_abund, 'combined.txt', log, key_columns=key_columns)
            shutil.copy2('combined.txt', output.gene_abund_merge)
        else:
            shutil.copy2(input.gene_abund[0], output.gene_abund_merge)

###############################################################################################
# Generate gene expression profile
###############################################################################################
# Change names to R scripts
# what do to with annot file
rule integration_gene_profiles:
    input:
        gene_abund_merge=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv",
                                wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    output:
        gene_abund_prof=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_letter}.csv",
                                wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, omics_letter = omics_prof, normalization = normalization),
        gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_letter}.rds",
                                wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, omics_letter = omics_prof, normalization = normalization),
        gene_abund_plots=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_letter}.PCA.pdf",
                                wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, omics_letter = omics_prof, normalization = normalization),
    params:
        annot_file={annot_file},
        funct_opt= funct_opt_list
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_gene_profiles.{omics}.p{identity}.{normalization}.GX.log".format(
                wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    resources:
        mem=25
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time ( echo 'integration of gene profiles'
            Rscript {script_dir}/gene_expression_profile_genome_based.R {threads} $(dirname {output.gene_abund_prof[0]}) {omics} {params.annot_file} {metadata} {input.gene_abund_merge} {params.funct_opt} {main_factor}
        ) &> {log}
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

def get_function_profile_integration_input(wildcards):
    gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            omics_prof = omics_prof,
            normalization = normalization)
    genome_profile=expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.profile.relabund.prop.genome.txt",
            wd = working_dir,
            omics = 'metaT' if omics=='metaG_metaT' else omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            omics_prof = omics_prof,
            normalization = normalization)
    if (normalization == 'MG'):
        return {
                'gene_abund_phyloseq' : gene_abund_phyloseq,
                'genome_profile'      : genome_profile
            }
    else:
        return {
                'gene_abund_phyloseq' : gene_abund_phyloseq,
            }

# TODO: Call R scripts for each functional category independently, so that all values are passed as wildcards. Also R scripts need to be fixed.
rule integration_function_profiles_MG_TPM:
    input:
        unpack(get_function_profile_integration_input)
    output:
        tsv=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.tsv", 
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                identity = identity,
                omics_prof = omics_prof,
                normalization = normalization,
                funct_opt = funct_opt),
        physeq=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds",
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                identity = identity,
                omics_prof = omics_prof,
                normalization = normalization,
                funct_opt = funct_opt),
        plots=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCA.pdf",
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                identity = identity,
                omics_prof = omics_prof,
                normalization = normalization,
                funct_opt = funct_opt),
    params:
        annot_file={annot_file},
        funct_opt=funct_opt_list,
        mapped_reads_threshold=config["MIN_mapped_reads"],
        genome_weights_arg = lambda wildcards, input: input.genome_profile if (normalization == 'MG') else ""
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_funtion_profiles.{omics}.{normalization}.log".format(wd = working_dir, post_analysis_out = post_analysis_out, normalization = normalization, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        time ( echo 'integration of function profiles'
            Rscript {script_dir}/function_expression_profile_genome_MG_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {params.annot_file} {metadata} {minto_dir} {params.funct_opt} {params.mapped_reads_threshold} {main_factor} {params.genome_weights_arg}
        ) &> {log}
        """
