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
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_prof}_F{omics_prof}_features.pdf",
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),\
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

        # Merge metaG and metaT profiles if the input is of length 2
        # Or just read and write the single input if length is 1
        combine_profiles(input.gene_abund, 'combined.txt', log, key_columns=key_columns)
        shutil.copy2('combined.txt', output.gene_abund_merge)

###############################################################################################
# Generate gene profile
# ~~~~~~~~~~~~~~~~~~~~~
#
# Generate one of {GA, GT, GE}.
# These outputs need metaG, metaT, metaG_metaT, resp. as --omics arg for the R script.
# So, we use param script_omics to figure that out.
# When the main workflow needs GA/GT/GE because the main omics is 'metaG_metaT', this rule will
# be run 3 times, once for each.
###############################################################################################
rule integration_gene_profiles:
    input:
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv"
    output:
        gene_abund_prof="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_alphabet}.csv",
        gene_abund_phyloseq="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_alphabet}.rds",
        gene_abund_plots="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_alphabet}.PCA.pdf",
    params:
        annot_file = annot_file,
        metadata_file = metadata,
        funcat_names = funct_opt_list,
        script_omics = lambda wildcards: 'metaG' if wildcards.omics_alphabet == 'A' \
                                         else ('metaT' if wildcards.omics_alphabet == 'T' \
                                               else 'metaG_metaT')
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_gene_profiles.{omics}.p{identity}.{normalization}.G{omics_alphabet}.log"
    resources:
        mem=25
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time ( echo 'integration of gene profiles'
            Rscript {script_dir}/gene_expression_profile_genome_based.R \
                    --threads {threads} \
                    --outdir $(dirname {output.gene_abund_prof}) \
                    --main-factor {main_factor} \
                    --annotation {params.annot_file} \
                    --metadata {params.metadata_file} \
                    --gene-profile {input.gene_abund_merge} \
                    --funcat-names {params.funcat_names} \
                    --omics {params.script_omics}
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
    gene_abund_phyloseq="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds".format(
            wd = wildcards.wd,
            omics = wildcards.omics,
            post_analysis_out = wildcards.post_analysis_out,
            identity = wildcards.identity,
            omics_prof = wildcards.omics_prof,
            normalization = wildcards.normalization)
    ret_dict = {'gene_abund_phyloseq' : gene_abund_phyloseq}

    # Add genome-weights if it is MG-normalization
    genome_profiles=expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.profile.relabund.prop.genome.txt",
            wd = working_dir,
            omics = ['metaG', 'metaT'],
            post_analysis_out = wildcards.post_analysis_out,
            identity = wildcards.identity,
            omics_prof = wildcards.omics_prof,
            normalization = wildcards.normalization)
    if (normalization == 'MG'):
        if (wildcards.omics.find('metaG') != -1):
            ret_dict['metaG_profile'] = genome_profiles[0]
        if (wildcards.omics.find('metaT') != -1):
            ret_dict['metaT_profile'] = genome_profiles[1]

    return ret_dict

#################################################################
# Quantify functions by functional category of interest.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Handles all different normalizations in one rule using one R script.
#
# Input:
# ------
#  Either TPM normalized or MG normalized gene abundances from metaG and/or metaT
# Output:
# -------
#  Depending on the input, we will generate:
#   FA: quantification of the function from metaG data
#   FT: quantification of the function from metaT data
#   FE: ratio FA/FT
# TPM normalization:
# ------------------
#  Quantification is straightforward: sum up all genes mapped to a given functional unit (e.g. K12345).
#  This also means that functional profiles do not sum to million anymore.
#  Since each gene can map to multiple KEGG KO's for example, total might be over million.
# MG normalization:
# -----------------
#  Quantification is a little bit more involved.
#  Weighted sum of genes belonging to a given functional unit (e.g. K12345), weighted by the relative abundance of the MAG it belongs to.
#  Weighting for metaG/metaT gene profiles uses relative abundance from metaG/metaT space, respectively.
#  Interpretation:
#    If each species in the community carries the functional unit and expresses it at the same level as marker genes, it will be 1.0
#################################################################

rule integration_function_profiles:
    input:
        unpack(get_function_profile_integration_input)
    output:
        abundance="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funcat}.tsv",
        features="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_prof}_F{omics_prof}_features.{funcat}.tsv",
        physeq="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funcat}.rds",
        pca="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funcat}.PCA.pdf",
    wildcard_constraints:
        normalization='MG|TPM',
        omics_prof='A|E|T'
    params:
        funcat_desc_file = lambda wildcards: "{location}/data/descriptions/{name}.tsv".format(
                                                    location=minto_dir,
                                                    name=wildcards.funcat.replace("kofam_", "KEGG_").replace("merged_", "KEGG_")),
        weights_arg = lambda wildcards, input: "" if (wildcards.normalization == 'TPM') \
                                               else (f"--genome-weights-metaG {input.metaG_profile}" if (wildcards.omics == 'metaG') \
                                                     else (f"--genome-weights-metaT {input.metaT_profile}" if (wildcards.omics == 'metaT') \
                                                           else f"--genome-weights-metaG {input.metaG_profile} --genome-weights-metaT {input.metaT_profile}" \
                                                          ) \
                                                    )
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_funtion_profiles.{omics}.p{identity}.{normalization}.F{omics_prof}.{funcat}.log"
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        echo 'integration of function profiles'
        time (
            Rscript {script_dir}/function_abundance_and_expression_profiling.R \
                    --threads {threads} \
                    --outdir $(dirname {output.abundance}) \
                    --main-factor {main_factor} \
                    {params.weights_arg} \
                    --normalization {wildcards.normalization} \
                    --funcat-name {wildcards.funcat} \
                    --funcat-desc {params.funcat_desc_file} \
                    --omics {wildcards.omics}
        ) &> {log}
        """

# Combine individual tsv files listing the number of features within each funcat
# into a single data frame.

rule combine_feature_counts:
    input:
        features=lambda wildcards: expand("{somewhere}/{something}_features.{funcat}.tsv",
                                            somewhere=wildcards.somewhere,
                                            something=wildcards.something,
                                            funcat=funct_opt)
    output:
        combined="{somewhere}/{something}_features.tsv",
    run:
        import pandas as pd

        # Read files and append
        df_list = list()
        for f in input.features:
            df = pd.read_csv(f, comment='#', header=0, sep = "\t", memory_map=True)
            df_list.append(df)
        df = pd.concat(df_list, ignore_index=True, copy=False, sort=False)
        df.to_csv(output.combined, sep = "\t", index = False)

# Make a plot of number of features from each functional category.
# It can make GA_FA, GT_FT, GE_FE, depending on what is being asked.

rule make_feature_count_plot:
    input:
        tsv=rules.combine_feature_counts.output.combined
    output:
        pdf="{somewhere}/plots/{something}_features.pdf"
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        R --vanilla --silent --no-echo <<___EOF___
        library(data.table)
        library(ggplot2)
        library(dplyr)

        count_df <- fread('{input.tsv}', header=T) %>%
                        as.data.frame(stringsAsFactors = F, row.names = T) %>%
                        distinct() %>%
                        mutate(DB = reorder(DB, feature_n))

        pdf('{output.pdf}', width=6, height=5, paper="special" )
        print(ggplot(data=count_df, aes(x=DB, y=feature_n)) +
                geom_bar(stat="identity")+ theme_minimal()+
                facet_wrap(feature~., scales= "free") + labs(y='Number of features', x='')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))+
                geom_text(aes(label=feature_n), position=position_dodge(width=0.9), vjust=-0.25, size = 3)
             )
        dev.off()
___EOF___
        """
