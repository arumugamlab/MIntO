#!/usr/bin/env python

'''
KEGG module completeness estimation

Authors: Carmen Saenz, Mani Arumugam, Judit Szarvas
'''

# configuration yaml file
# import sys
import os.path
import math
import re

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/resources.smk'

# Define the 3 modes
valid_minto_modes = ['MAG', 'refgenome']

# Which database are we mapping reads to?
if 'MINTO_MODE' in config and config['MINTO_MODE'] != None:
    MINTO_MODE=config['MINTO_MODE']
else:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be defined".format(config_path))

# Backward compatibility and common misnomers
if MINTO_MODE in ['reference_genome', 'reference-genome', 'reference', 'refgenomes']:
    MINTO_MODE = 'refgenome'
elif MINTO_MODE in ['MAGs', 'mag', 'mags']:
    MINTO_MODE = 'MAG'

if not MINTO_MODE in valid_minto_modes:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be {}.".format(config_path, " or ".join(valid_minto_modes)))

if config['abundance_normalization'] in ("MG", "TPM"):
    normalization=config['abundance_normalization']
else:
    print('ERROR in ', config_path, ': abundance_normalization variable is not correct. "abundance_normalization" variable should be MG or TPM.')

# Mapping percent identity
if config['alignment_identity'] is None:
    raise Exception("ERROR in {}: alignment_identity variable is empty. Please, fix.".format(config_path))
elif type(config['alignment_identity']) != int:
    raise Exception("ERROR in {}: alignment_identity variable is not an integer. Please, fix.".format(config_path))
identity=config['alignment_identity']

if config['ANNOTATION'] is None:
    raise Exception("ERROR in {}: ANNOTATION variable in configuration yaml file is empty. Please, complete.".format(config_path))
else:
    if all([x in ('kofam','eggNOG') for x in config['ANNOTATION']]):
        annot=config['ANNOTATION']
    else:
        raise Exception("ERROR in {}: ANNOTATION variable in configuration yaml file should be kofam or eggNOG. Please, correct it.".format(config_path))

# Define all the outputs needed by target 'all'

if omics == 'metaG':
    omics_prof='A'
elif omics == 'metaT':
    omics_prof='T'
elif omics == 'metaG_metaT':
    omics_prof='T'

print('NOTE: MIntO is using ', omics, ' as omics variable.')

for omics_type in omics.split("_"):
    omics_folder = os.path.join(working_dir, omics_type)
    if not os.path.exists(omics_folder):
        raise Exception(f"ERROR in {omics} setting, the folder {omics_folder} does not exist.")

GENE_DB_TYPE = MINTO_MODE + '-genes'


print('NOTE: MIntO is using ', annot ,' as ANNOTATION variable.')

# Define all the outputs needed by target 'all'
rule all:
    input:
        expand("{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/F{omics_alphabet}.per_sample.{annot}.module_comp.tsv",
            wd = working_dir,
            omics = omics,
            gene_db = GENE_DB_TYPE,
            identity = identity,
            normalization = normalization,
            omics_alphabet = omics_prof,
            annot = annot),
        expand("{wd}/DB/{minto_mode}/4-annotations/per_mag.{annot}.module_comp.tsv",
            wd = working_dir,
            minto_mode = MINTO_MODE,
            annot = annot)

###############################################################################################
## Prepare inputs to EBI::kegg-pathways-completeness-tool
###############################################################################################

rule prepare_per_sample_input:
    input:
        pso = "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_alphabet}.{annot}.KEGG_KO.qs"
    output:
        tsv = temp("{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_alphabet}.per_sample.{annot}.tsv")
    shadow:
        "minimal"
    resources:
        mem = 5
    threads: 1
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        R --vanilla --silent --no-echo <<___EOF___
library(phyloseq)
library(qs)
library(data.table)
pso <- qread('{input.pso}')
snames <- sample_names(pso)
kos <- rownames(otu_table(pso))
dt <- data.table(unclass(otu_table(pso)))
dt <- dt[, lapply(.SD, function(x) ifelse(x == 0, NA, x))]
dt <- dt[, lapply(.SD, as.character)]
# kos are now cols
dt <- transpose(dt)
kos_obj <- vector("list", length(kos))
for (i in seq(1,length(kos))){{
    kos_obj[[i]] <- dt[,lapply(.SD, function(x) ifelse(! is.na(x), kos[[i]], x)), .SDcols = paste0("V", i)]
}}
kos_obj_dt <- do.call("cbind", kos_obj)
kos_obj_dt <- cbind(data.table(ID=snames), kos_obj_dt)
fwrite(kos_obj_dt, file = 'per_sample.kos.tsv', sep = "\t", row.names = F,col.names = F)
___EOF___

tr -s '\t' < per_sample.kos.tsv > {output.tsv}
        """

rule prepare_per_mag_input:
    input:
        tsv = "{wd}/DB/{minto_mode}/4-annotations/{annot}.tsv",
        combi = "{wd}/DB/{minto_mode}/4-annotations/combined_annotations.tsv"
    output:
        tsv = temp("{wd}/DB/{minto_mode}/4-annotations/{annot}/per_mag.{annot}.tsv"),
    localrule: True
    run:
        mag_kos = {}
        with open(input.combi, "r") as fp:
            header = fp.readline().strip().split("\t")
            col_name = f"{wildcards.annot}.KEGG_KO"
            col_idx = header.index(col_name)
            for line in fp:
                tmp = line.strip().split("\t")
                mag_id = tmp[0].split("|")[-1].split("_")[0]
                kos = tmp[col_idx].split(",")
                if "-" not in kos:
                    if mag_id in mag_kos:
                        mag_kos[mag_id].extend(kos)
                    else:
                        mag_kos[mag_id] = kos
        with open(output.tsv, "w") as of:
            for m in mag_kos.keys():
                l = list(set(mag_kos.get(m)))
                print(m, "\t".join(l), file=of, sep="\t")

###############################################################################################
# Estimate module completeness
###############################################################################################
rule sample_completeness:
    input:
        tsv="{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_alphabet}.per_sample.{annot}.tsv",
        kpc_pathways="{}/data/kofam_db/all_pathways.txt".format(minto_dir),
        kpc_class="{}/data/kofam_db/all_pathways_class.txt".format(minto_dir),
        kpc_names="{}/data/kofam_db/all_pathways_names.txt".format(minto_dir),
    output:
        "{wd}/output/data_integration/{gene_db}/{omics}.gene_abundances.p{identity}.{normalization}/F{omics_alphabet}.per_sample.{annot}.module_comp.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/output/data_integration/{gene_db}/{omics}.p{identity}.{normalization}.F{omics_alphabet}.per_sample.{annot}.module_comp.log"
    resources:
        mem = 5
    threads: 4
    conda:
        config["minto_dir"]+"/envs/kpc.yml"
    shell:
        """
        time (
            give_pathways -i {input.tsv} --pathways {input.kpc_pathways} --names {input.kpc_names} --classes {input.kpc_class} -o per_sample
            sed -e 's|^contig|sample_alias|' per_sample.summary.kegg_contigs.tsv > {output}
        ) >& {log}
        """

rule mag_completeness:
    input:
        tsv="{wd}/DB/{minto_mode}/4-annotations/{annot}/per_mag.{annot}.tsv",
        kpc_pathways="{}/data/kofam_db/all_pathways.txt".format(minto_dir),
        kpc_class="{}/data/kofam_db/all_pathways_class.txt".format(minto_dir),
        kpc_names="{}/data/kofam_db/all_pathways_names.txt".format(minto_dir),
    output:
        "{wd}/DB/{minto_mode}/4-annotations/per_mag.{annot}.module_comp.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/per_mag.{annot}.module_comp.log"
    resources:
        mem = 5
    threads: 4
    conda:
        config["minto_dir"]+"/envs/kpc.yml"
    shell:
        """
        time (
            give_pathways -i {input.tsv} --pathways {input.kpc_pathways} --names {input.kpc_names} --classes {input.kpc_class} -o per_mag
            sed -e 's|^contig|ID|' per_mag.summary.kegg_contigs.tsv > {output}
        ) >& {log}
        """
