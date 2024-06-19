#!/usr/bin/env python

'''
Gene and function profiling step

Authors: Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import os.path
from os import path
import re

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule integration_rpkg from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# some variables

main_factor = None
if config['MAIN_factor'] is not None:
    main_factor = config['MAIN_factor']

# MIntO mode and database-mapping

# Define the 3 modes
valid_minto_modes = ['MAG', 'refgenome', 'catalog']

# Which database are we mapping reads to?
if 'MINTO_MODE' in config and config['MINTO_MODE'] != None:
    MINTO_MODE=config['MINTO_MODE']
else:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be defined".format(config_path))

# Backward compatibility and common misnomers
if MINTO_MODE in ['db_genes', 'db-genes', 'genes_db', 'gene_catalog', 'gene-catalog']:
    MINTO_MODE = 'catalog'
elif MINTO_MODE in ['reference_genome', 'reference-genome', 'reference', 'refgenomes']:
    MINTO_MODE = 'refgenome'
elif MINTO_MODE in ['MAGs', 'mag', 'mags']:
    MINTO_MODE = 'MAG'

if not MINTO_MODE in valid_minto_modes:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be {}.".format(config_path, " or ".join(valid_minto_modes)))

if config['abundance_normalization'] in ("MG", "TPM"):
    normalization=config['abundance_normalization']
else:
    print('ERROR in ', config_path, ': abundance_normalization variable is not correct. "abundance_normalization" variable should be MG or TPM.')

if normalization == 'MG' and MINTO_MODE in ("catalog"):
    raise Exception("ERROR in {}: In 'catalog' mode, only TPM normalization is allowed.".format(config_path))


# Mapping percent identity
if config['alignment_identity'] is None:
    raise Exception("ERROR in {}: alignment_identity variable is empty. Please, fix.".format(config_path))
elif type(config['alignment_identity']) != int:
    raise Exception("ERROR in {}: alignment_identity variable is not an integer. Please, fix.".format(config_path))
identity=config['alignment_identity']

if config['MERGE_memory'] is None:
    raise Exception("ERROR in {}: MERGE_memory variable is empty. Please, fix.".format(config_path))
elif type(config['MERGE_memory']) != int:
    raise Exception("ERROR in {}: MERGE_memory variable is not an integer. Please, fix.".format(config_path))

if config['MERGE_threads'] is None:
    raise Exception("ERROR in {}: MERGE_threads variable is empty. Please, fix.".format(config_path))
elif type(config['MERGE_threads']) != int:
    raise Exception("ERROR in {}: MERGE_threads variable is not an integer. Please, fix.".format(config_path))

if MINTO_MODE == 'catalog':
    if config['ANNOTATION_file'] is None:
        raise Exception("Gene functional annotation needs to be provided via ANNOTATION_file variable")
    elif path.exists(config['ANNOTATION_file']) is False:
        raise Exception("File specified in ANNOTATION_file variable does not exist")
    elif path.exists(config['ANNOTATION_file']) is True:
        annot_file=config['ANNOTATION_file']

if config['ANNOTATION_ids'] is None:
    raise Exception("ERROR in {}: ANNOTATION_ids variable in configuration yaml file is empty. Please, complete.".format(config_path))
else:
    funct_opt=config['ANNOTATION_ids']

# Define all the outputs needed by target 'all'

if omics == 'metaG':
    omics_prof='A'
elif omics == 'metaT':
    omics_prof='T'
elif omics == 'metaG_metaT':
    omics_prof=['A','T','E']

print('NOTE: MIntO is using ', omics, ' as omics variable.')

for omics_type in omics.split("_"):
    omics_folder = path.join(working_dir, omics_type)
    if not path.exists(omics_folder):
        raise Exception(f"ERROR in {omics} setting, the folder {omics_folder} does not exist.")

GENE_DB_TYPE = MINTO_MODE + '-genes'

if MINTO_MODE in ['MAG', 'refgenome']:
    annot_file="{wd}/DB/{subdir}/4-annotations/combined_annotations.tsv".format(wd = working_dir, subdir = MINTO_MODE)

print('NOTE: MIntO is using ', annot_file ,' as ANNOTATION_file variable.')

print('NOTE: MIntO is using ', funct_opt, ' as ANNOTATION_ids variable.')


def integration_merge_profiles():
    result = expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}.csv",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization)
    return(result)

def integration_gene_profiles():
    result = expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/G{omics_prof}.csv",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),
    expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.qs",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),
    expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/plots/G{omics_prof}.PCA.pdf",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof)
    return(result)

def integration_function_profiles():
    result = expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.tsv",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),
    expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/plots/G{omics_prof}_F{omics_prof}_features.pdf",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),
    expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.qs",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),
    expand("{wd}/output/data_integration/{subdir}/{omics}.gene_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCA.pdf",
            wd = working_dir,
            omics = omics,
            subdir = GENE_DB_TYPE,
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
        integration_function_profiles(),
        print_versions.get_version_output(snakefile_name)
    default_target: True

###############################################################################################
# Prepare gene profile
## Merge gene abundance and gene transcript profiles
###############################################################################################
rule integration_merge_profiles:
    input:
        single=lambda wildcards: expand("{wd}/{omics_individual}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.{normalization}.csv",
                                            wd = wildcards.wd,
                                            omics_individual = wildcards.omics.split("_"),
                                            minto_mode = wildcards.gene_db.replace('-genes', ''),
                                            identity = wildcards.identity,
                                            normalization = wildcards.normalization),
    output:
        merged="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}.csv"
    params:
        files = lambda wildcards, input: ",".join(input.single)
    shadow:
        "minimal"
    log:
        "{wd}/logs/output/data_integration/{gene_db}/{omics}.p{identity}.{normalization}.integration_merge_profiles.log"
    resources:
        mem=config["MERGE_memory"]
    threads: 1
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (
            Rscript {script_dir}/merge_profiles.R --threads {threads} --memory {resources.mem} --input {params.files} --out out.txt --keys ID
            rsync -a out.txt {output.merged}
        ) >& {log}
        """

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
        annot_file = annot_file,
        gene_abund_merge=rules.integration_merge_profiles.output.merged
    output:
        gene_abund_prof=    "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/G{omics_alphabet}.csv",
        gene_abund_phyloseq="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_alphabet}.qs",
        gene_abund_plots=   "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/plots/G{omics_alphabet}.PCA.pdf",
    params:
        metadata_file = metadata,
        funcat_names = ','.join(['"' + id + '"' for id in funct_opt]),
        script_omics = lambda wildcards: 'metaG' if wildcards.omics_alphabet == 'A' \
                                         else ('metaT' if wildcards.omics_alphabet == 'T' \
                                               else 'metaG_metaT')
    log:
        "{wd}/logs/output/data_integration/{gene_db}/integration_gene_profiles.{omics}.p{identity}.{normalization}.G{omics_alphabet}.log"
    resources:
        mem=25
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        echo 'integration of gene profiles'
        time (
            Rscript {script_dir}/gene_abundance_and_expression_profiling.R \
                    --threads {threads} \
                    --outdir $(dirname {output.gene_abund_prof}) \
                    --main-factor {main_factor} \
                    --annotation {input.annot_file} \
                    --metadata {params.metadata_file} \
                    --gene-profile {input.gene_abund_merge} \
                    --funcat-names {params.funcat_names} \
                    --omics {params.script_omics}
        ) &> {log}
        """

##################################################################################################
# Generate function expression profile for the different databases included in the annotation file
##################################################################################################

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

# We replicate it 2 times, so that FA/FT, FA+FT+FE are handled separately
# This is because FE mode generates FA+FT+FE; but FA/FT generate only FA/FT resp.
# And to avoid complications, we define the rules conditionally.

if omics == 'metaG_metaT':
    def get_function_profile_integration_input_FE(wildcards):
        gene_abund_phyloseq="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/GE.qs".format(
                wd = wildcards.wd,
                omics = wildcards.omics,
                gene_db = wildcards.gene_db,
                identity = wildcards.identity,
                normalization = wildcards.normalization)
        ret_dict = {'gene_abund_phyloseq' : gene_abund_phyloseq}

        # Add genome-weights if it is MG-normalization
        genome_profiles=expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/genome_abundances.p{identity}.profile.relabund.prop.genome.txt",
                wd = working_dir,
                omics = ['metaG', 'metaT'],
                minto_mode = wildcards.gene_db.replace('-genes', ''),
                identity = wildcards.identity,
                normalization = wildcards.normalization)
        if (normalization == 'MG'):
            if (wildcards.omics.find('metaG') != -1):
                ret_dict['metaG_profile'] = genome_profiles[0]
            if (wildcards.omics.find('metaT') != -1):
                ret_dict['metaT_profile'] = genome_profiles[1]

        return ret_dict

    rule integration_function_profiles_FE:
        input:
            unpack(get_function_profile_integration_input_FE)
        output:
            abundance=   "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/FE.{funcat}.tsv",
            features=    "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/GE_FE_features.{funcat}.tsv",
            physeq=      "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/FE.{funcat}.qs",
            pca=         "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/plots/FE.{funcat}.PCA.pdf",
            FA_abundance="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/FA.{funcat}.tsv",
            FA_features= "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/GA_FA_features.{funcat}.tsv",
            FA_physeq=   "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/FA.{funcat}.qs",
            FA_pca=      "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/plots/FA.{funcat}.PCA.pdf",
            FT_abundance="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/FT.{funcat}.tsv",
            FT_features= "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/GT_FT_features.{funcat}.tsv",
            FT_physeq=   "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/FT.{funcat}.qs",
            FT_pca=      "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/plots/FT.{funcat}.PCA.pdf",
        wildcard_constraints:
            normalization='MG|TPM',
        params:
            funcat_desc_file = lambda wildcards: "{location}/data/descriptions/{name}.tsv".format(
                                                        location=minto_dir,
                                                        name=re.sub("eggNOG.KEGG_|kofam.KEGG_|merged.KEGG_", "KEGG_", wildcards.funcat)),
            weights_arg = lambda wildcards, input: "" if (wildcards.normalization == 'TPM') else f"--genome-weights-metaG {input.metaG_profile} --genome-weights-metaT {input.metaT_profile}"
        log:
            "{wd}/logs/output/data_integration/{gene_db}/integration_funtion_profiles.{omics}.p{identity}.{normalization}.FE.{funcat}.log"
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
else :
    def get_function_profile_integration_input_FA_FT(wildcards):
        gene_abund_phyloseq="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.qs".format(
                wd = wildcards.wd,
                omics = wildcards.omics,
                gene_db = wildcards.gene_db,
                identity = wildcards.identity,
                omics_prof = wildcards.omics_prof,
                normalization = wildcards.normalization)
        ret_dict = {'gene_abund_phyloseq' : gene_abund_phyloseq}

        # Add genome-weights if it is MG-normalization
        genome_profiles=expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/genome_abundances.p{identity}.profile.relabund.prop.genome.txt",
                wd = working_dir,
                omics = ['metaG', 'metaT'],
                minto_mode = wildcards.gene_db.replace('-genes', ''),
                identity = wildcards.identity,
                normalization = wildcards.normalization)
        if (normalization == 'MG'):
            if (wildcards.omics.find('metaG') != -1):
                ret_dict['metaG_profile'] = genome_profiles[0]
            if (wildcards.omics.find('metaT') != -1):
                ret_dict['metaT_profile'] = genome_profiles[1]

        return ret_dict

    rule integration_function_profiles_FA_FT:
        input:
            unpack(get_function_profile_integration_input_FA_FT)
        output:
            abundance="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/F{omics_prof}.{funcat}.tsv",
            features= "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/G{omics_prof}_F{omics_prof}_features.{funcat}.tsv",
            physeq=   "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funcat}.qs",
            pca=      "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funcat}.PCA.pdf",
        wildcard_constraints:
            normalization='MG|TPM',
            omics_prof='A|T'
        params:
            funcat_desc_file = lambda wildcards: "{location}/data/descriptions/{name}.tsv".format(
                                                        location=minto_dir,
                                                        name=re.sub("eggNOG.KEGG_|kofam.KEGG_|merged_", "KEGG_", wildcards.funcat)),
            weights_arg = lambda wildcards, input: "" if (wildcards.normalization == 'TPM') \
                                                   else (f"--genome-weights-metaG {input.metaG_profile}" if (wildcards.omics == 'metaG') \
                                                         else f"--genome-weights-metaT {input.metaT_profile}" \
                                                        )
        log:
            "{wd}/logs/output/data_integration/{gene_db}/integration_funtion_profiles.{omics}.p{identity}.{normalization}.F{omics_prof}.{funcat}.log"
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
    localrule: True
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
    resources:
        mem=config["MERGE_memory"]
    threads: 1
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
