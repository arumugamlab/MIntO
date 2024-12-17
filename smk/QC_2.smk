#!/usr/bin/env python

'''
Pre-processing of metaG and metaT data step
    - read length filtering
    - host genome filtering
    - rRNA filtering
Assembly-free taxonomy profiling step

Authors: Carmen Saenz, Mani Arumugam
'''

import re
import os.path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule QC_2_base, QC_2_rpkg, QC_2_mpl, QC_2_motus from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

localrules: qc2_filter_config_yml_assembly, qc2_filter_config_yml_mapping, \
            metaphlan_combine_profiles, motus_combine_profiles, motus_calc_motu, \
            plot_taxonomic_profile, plot_sourmash_kmers

taxonomies_versioned = list()
ilmn_samples = list()
merged_illumina_samples = dict()

##############################################
# Register composite samples
##############################################

# Make list of illumina coassemblies, if MERGE_ILLUMINA_SAMPLES in config
if (x := validate_optional_key(config, 'MERGE_ILLUMINA_SAMPLES')):
    for m in x:
        #print(" "+m)
        merged_illumina_samples.append(m)

##############################################
# Get sample list
##############################################

# Make list of illumina samples, if ILLUMINA in config
if (x := validate_optional_key(config, 'ILLUMINA')):
    check_fastq_file_locations(x, locations = ['5-1-sortmerna', '4-hostfree', '3-minlength', '1-trimmed'])
    ilmn_samples = x

    # If it's composite sample, then don't need to see them until it gets merged later
    for i in x:
        if i in merged_illumina_samples.keys():
            ilmn_samples.remove(i)

##############################################
# Make list of clean illumina samples after merging reps
##############################################

ilmn_reps_to_delete = [j.strip() for m in merged_illumina_samples.values() for j in m.split('+')]
nonredundant_ilmn_samples = list(merged_illumina_samples.keys())
for i in ilmn_samples:
    if i not in ilmn_reps_to_delete:
        nonredundant_ilmn_samples.append(i)

##############################################
# Host genome filtering
##############################################

host_genome_path = validate_required_key(config, 'PATH_host_genome')
host_genome_name = validate_required_key(config, 'NAME_host_genome')

if os.path.exists(f"{host_genome_path}/BWA_index/{host_genome_name}.pac"):
    print(f"Host genome db {host_genome_path}/BWA_index/{host_genome_name} will be used")
elif os.path.exists(f"{host_genome_path}/{host_genome_name}"):
    print(f"Host genome sequence {host_genome_path}/{host_genome_name} will be used to create BWA db")
else:
    raise Exception(f"NAME_host_genome={host_genome_name} does not exist as fasta or BWA db in PATH_host_genome={host_genome_path}. Please fix {config_path}")

##############################################
# Minimum read-length trimming
##############################################

read_min_len = validate_required_key(config, 'READ_minlen')
if read_min_len < 50:
    read_min_len = 50

##############################################
# BWA for host genome filtering
##############################################

bwa_threads = validate_required_key(config, 'BWA_host_threads')

##############################################
# taxonomy
##############################################

TAXA_threads = validate_required_key(config, 'TAXA_threads')
TAXA_memory  = validate_required_key(config, 'TAXA_memory')

taxonomy = validate_required_key(config, 'TAXA_profiler')
allowed = ('metaphlan', 'motus_rel', 'motus_raw')
for x in taxonomy.split(","):
    check_allowed_values('TAXA_profiler', x, allowed)

metaphlan_version = validate_required_key(config, 'metaphlan_version')
motus_version     = validate_required_key(config, 'motus_version')

taxonomies = taxonomy.split(",")
for t in taxonomies:
    version="unknown"
    if t.startswith("motus"):
        version=motus_version
    elif t.startswith("metaphlan"):
        version=metaphlan_version
    taxonomies_versioned.append(t+"."+version)

##############################################
# metaG - QC plots
##############################################

plot_args_list = list()

main_factor = validate_required_key(config, 'MAIN_factor')
plot_args_list.append('--factor ' + main_factor)

plot_factor2 = validate_optional_key(config, 'PLOT_factor2')
if plot_factor2 is not None:
    plot_args_list.append(f"--factor2 {plot_factor2}")

plot_time = validate_optional_key(config, 'PLOT_time')
if plot_time is not None:
    plot_args_list.append(f"--time {plot_time}")

plot_args_str = ' '.join(plot_args_list)


##############################################
# metaT - rRNA removal
##############################################

if omics == 'metaT':
    sortmeRNA_threads = validate_required_key(config, 'sortmeRNA_threads')
    sortmeRNA_memory = validate_required_key(config, 'sortmeRNA_memory')
    sortmeRNA_db = validate_required_key(config, 'sortmeRNA_db')
    sortmeRNA_db_idx = validate_required_key(config, 'sortmeRNA_db_idx')


##############################################
# Metaphlan DB version
##############################################
with open(minto_dir + "/data/metaphlan/" + metaphlan_version + "/mpa_latest", 'r') as file:
    metaphlan_index = file.read().rstrip()

##############################################
# Co-assembly grouping variable
##############################################

coas_factor = validate_optional_key(config, 'COAS_factor')
if coas_factor is None:
    coas_factor = main_factor

# Site customization for avoiding NFS traffic during I/O heavy steps such as mapping

CLUSTER_NODES            = None
CLUSTER_LOCAL_DIR        = None
CLUSTER_WORKLOAD_MANAGER = None
include: minto_dir + '/site/cluster_def.py'

# Cluster-aware bwa-index rules

include: 'include/bwa_index_wrapper.smk'

##############################################
# Define all the outputs needed by target 'all'
##############################################

def taxonomy_plot_output():

    results = list()
    profiles = [expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.tsv",
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies_versioned),
            expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.merged_abundance_table_species.txt",
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies_versioned)]
    plots = [expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.PCoA.Bray_Curtis.pdf",
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies_versioned),
            expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.Top15genera.pdf",
                wd = working_dir,
                omics = omics,
                taxonomy = taxonomies_versioned)]
    results.append(profiles)
    results.append(plots)
    return(results)

def smash_plot_output():
    if omics == "metaG":
        barplots = "{wd}/output/6-1-smash/{omics}.{taxonomy_versioned}.clusters.pdf".format(
                        wd = working_dir,
                        omics = omics,
                        taxonomy_versioned = taxonomies_versioned[0])
        results = [barplots]
        return(results)
    else:
        return()

def merged_sample_output():
    return()

if merged_illumina_samples:
    def merged_sample_output():
        result = expand("{wd}/{omics}/{location}/{sample}/{sample}.{pair}.fq.gz",
                        wd = working_dir,
                        omics = omics,
                        location = get_qc2_output_location(omics),
                        sample = merged_illumina_samples.keys(),
                        pair = ['1', '2'])
        return(result)

def next_step_config_yml_output():
    result = expand("{wd}/{omics}/{yaml}.yaml",
                wd = working_dir,
                yaml = ['assembly', 'mapping'],
                omics = omics)
    return(result)

def clean_bwa_index_host():
    if CLUSTER_NODES != None and validate_optional_key(config, 'CLEAN_BWA_INDEX'):
        return(f"{host_genome_path}/BWA_index/{host_genome_name}.clustersync/cleaning.done")
    else:
        return([])

rule all:
    input:
        clean_bwa_index_host(),
        merged_sample_output(),
        taxonomy_plot_output(),
        smash_plot_output(),
        next_step_config_yml_output(),
        print_versions.get_version_output(snakefile_name)
    default_target: True

###############################################################################################
# Pre-processing of metaG and metaT data step
# Read length filtering using the MINLEN
###############################################################################################

# We use suffix fq.gz for "seqkit seq" output even though it outputs unzipped fq.
# This is because the next step "seqkit pair" uses extension to decide whether to compress using pigz or not.
rule qc2_length_filter:
    input:
        read_fw='{wd}/{omics}/1-trimmed/{sample}/{run}.1.fq.gz',
        read_rv='{wd}/{omics}/1-trimmed/{sample}/{run}.2.fq.gz',
    output:
        paired1="{wd}/{omics}/3-minlength/{sample}/{run}.1.fq.gz",
        paired2="{wd}/{omics}/3-minlength/{sample}/{run}.2.fq.gz",
        summary="{wd}/{omics}/3-minlength/{sample}/{run}.trim.summary"
    shadow:
        "minimal"
    params:
        read_length_cutoff=read_min_len
    log:
        "{wd}/logs/{omics}/3-minlength/{sample}_{run}.log"
    resources:
        mem=20
    threads:
        16
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        remote_dir=$(dirname {output.paired1})
        time (
            mkfifo {wildcards.run}.1.fq.gz
            mkfifo {wildcards.run}.2.fq.gz
            seqkit seq {input.read_fw} -m {params.read_length_cutoff} | seqkit replace -p '/1' -r '' > {wildcards.run}.1.fq.gz &
            seqkit seq {input.read_rv} -m {params.read_length_cutoff} | seqkit replace -p '/2' -r '' > {wildcards.run}.2.fq.gz &
            seqkit pair -1 {wildcards.run}.1.fq.gz -2 {wildcards.run}.2.fq.gz -O result -u >& {output.summary}
            rsync -a result/{wildcards.run}.* $remote_dir/
        ) >& {log}
        """

###############################################################################################
# Pre-processing of metaG and metaT data
# Remove host genome sequences
###############################################################################################

# If CLUSTER_NODES is defined, then return the file symlink'ed to local drives.
# Otherwise, return the original index files in project work_dir.
def get_host_bwa_index(wildcards):

    # Where are the index files?
    if CLUSTER_NODES != None:
        index_location = 'BWA_index_local'
    else:
        index_location = 'BWA_index'

    # Get all the index files!
    files = expand("{somewhere}/{location}/{genome}.{ext}",
                    somewhere = host_genome_path,
                    location  = index_location,
                    genome    = host_genome_name,
                    ext       = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'])

    # Return them
    return(files)

# BWA mem memory is estimated as 3.1 bytes per base in database (regression: mem = 5.556e+09 + 3.011*input).
rule qc2_host_filter:
    input:
        pairead_fw=rules.qc2_length_filter.output.paired1,
        pairead_rv=rules.qc2_length_filter.output.paired2,
        bwaindex=get_host_bwa_index
    output:
        host_free_fw="{wd}/{omics}/4-hostfree/{sample}/{run}.1.fq.gz",
        host_free_rv="{wd}/{omics}/4-hostfree/{sample}/{run}.2.fq.gz",
    shadow:
        "minimal"
    params:
        bwaindex="{host_genome_path}/BWA_index/{host_genome_name}".format(host_genome_path=host_genome_path, host_genome_name=host_genome_name),
    log:
        "{wd}/logs/{omics}/4-hostfree/{sample}_{run}_filter_host_genome_BWA.log"
    resources:
        mem = lambda wildcards, input, attempt: 10 + int(3.1*os.path.getsize(input.bwaindex[0])/1e9) + 10*(attempt-1)
    threads:
        bwa_threads
    conda:
        minto_dir + "/envs/MIntO_base.yml" #bwa-mem2, msamtools>=1.1.1, samtools
    shell:
        """
        remote_dir=$(dirname {output.host_free_fw})
        time (
                bwa-mem2 mem -t {threads} -v 3 {params.bwaindex} {input.pairead_fw} {input.pairead_rv} \
                  | msamtools filter -S -l 30 --invert --keep_unmapped -bu - \
                  | samtools fastq -1 $(basename {output.host_free_fw}) -2 $(basename {output.host_free_rv}) -s /dev/null -c 6 -N -
                rsync -a * $remote_dir/
        ) >& {log}
        """

###############################################################################################
# Pre-processing of metaT data - rRNA filtering - only on metaT data
###############################################################################################

if omics == 'metaT':

    def get_rRNA_db_files(wildcards):
        files = ["rfam-5.8s-database-id98.fasta",
                "rfam-5s-database-id98.fasta",
                "silva-arc-16s-id95.fasta",
                "silva-arc-23s-id98.fasta",
                "silva-bac-16s-id90.fasta",
                "silva-bac-23s-id98.fasta",
                "silva-euk-18s-id95.fasta",
                "silva-euk-28s-id98.fasta"]
        return(expand("{sortmeRNA_db}/{f}",
                        sortmeRNA_db=sortmeRNA_db,
                        f=files))

    rule qc2_filter_rRNA_index:
        input:
            rRNA_db=get_rRNA_db_files
        output:
            rRNA_db_index_file = "{sortmeRNA_db_idx}/rRNA_db_index.log".format(sortmeRNA_db_idx=sortmeRNA_db_idx),
            rRNA_db_index = directory(expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx))
        shadow:
            "minimal"
        resources:
            mem=sortmeRNA_memory
        threads:
            sortmeRNA_threads
        log:
            "{wd}/logs/{omics}/5-1-sortmerna/rRNA_index.log".format(wd=working_dir, omics = omics),
        conda:
            minto_dir + "/envs/MIntO_base.yml" #sortmerna
        shell:
            """
            time (
                sortmerna --workdir . --idx-dir ./idx/ -index 1 \
                    --ref {input.rRNA_db[0]} \
                    --ref {input.rRNA_db[1]} \
                    --ref {input.rRNA_db[2]} \
                    --ref {input.rRNA_db[3]} \
                    --ref {input.rRNA_db[4]} \
                    --ref {input.rRNA_db[5]} \
                    --ref {input.rRNA_db[6]} \
                    --ref {input.rRNA_db[7]}
                rsync -a ./idx/* {output.rRNA_db_index}
                echo 'SortMeRNA indexed rRNA_databases done' > {sortmeRNA_db_idx}/rRNA_db_index.log
            ) >& {log}
            """

    rule qc2_filter_rRNA:
        input:
            host_free_fw=rules.qc2_host_filter.output.host_free_fw,
            host_free_rv=rules.qc2_host_filter.output.host_free_rv,
            rRNA_db_index=ancient(expand("{sortmeRNA_db_idx}", sortmeRNA_db_idx=sortmeRNA_db_idx))
        output:
            rRNA_out="{wd}/{omics}/5-1-sortmerna/{sample}/out/{run}.aligned.log",
            rRNA_free_fw="{wd}/{omics}/5-1-sortmerna/{sample}/{run}.1.fq.gz",
            rRNA_free_rv="{wd}/{omics}/5-1-sortmerna/{sample}/{run}.2.fq.gz"
        shadow:
            "minimal"
        params:
            db_idx_dir=sortmeRNA_db_idx,
            db_dir=sortmeRNA_db,
        resources:
            mem=sortmeRNA_memory
        threads:
            sortmeRNA_threads
        log:
            "{wd}/logs/{omics}/5-1-sortmerna/{sample}_{run}.log"
        conda:
            minto_dir + "/envs/MIntO_base.yml" #sortmerna
        shell:
            """
            time (
                sortmerna --paired_in --fastx --out2 --other --threads {threads} --no-best --num_alignments 1 --workdir . --idx-dir {params.db_idx_dir}/ \
                            --ref {params.db_dir}/rfam-5.8s-database-id98.fasta \
                            --ref {params.db_dir}/rfam-5s-database-id98.fasta \
                            --ref {params.db_dir}/silva-arc-16s-id95.fasta \
                            --ref {params.db_dir}/silva-arc-23s-id98.fasta \
                            --ref {params.db_dir}/silva-bac-16s-id90.fasta \
                            --ref {params.db_dir}/silva-bac-23s-id98.fasta \
                            --ref {params.db_dir}/silva-euk-18s-id95.fasta \
                            --ref {params.db_dir}/silva-euk-28s-id98.fasta \
                            --reads {input.host_free_fw} --reads {input.host_free_rv}
                parallel --jobs {threads} <<__EOM__
    rsync -a out/other_fwd.fq.gz {output.rRNA_free_fw}
    rsync -a out/other_rev.fq.gz {output.rRNA_free_rv}
    rsync -a out/aligned.log {output.rRNA_out}
__EOM__
            ) >& {log}
            """

###############################################################################################
# Create pseudo-samples that are created by merging multiple samples.
# E.g., for time-series data, we can create a composite sample with all time points
#       and create profiles for this composite sample.
###############################################################################################

if merged_illumina_samples:

    localrules: merge_fastqs_for_composite_samples
    ruleorder: merge_fastqs_for_composite_samples > qc2_host_filter

    if omics == 'metaT':
        ruleorder: merge_fastqs_for_composite_samples > qc2_filter_rRNA

    # Get the individual reps for the sample
    # And concat all the files for each rep into one
    def get_rep_files_for_composite_sample(wildcards):

        files = []
        reps = [x.strip() for x in merged_illumina_samples[wildcards.merged_sample].split('+')]
        for x in reps:
            files.extend(get_qc2_output_files_one_end(wildcards.wd, wildcards.omics, x, wildcards.pair))
        return(files)

    # Merge files for a given sample from all its reps
    # Restrict it to only those appearing in MERGE_ILLUMINA_SAMPLES dict in config file
    rule merge_fastqs_for_composite_samples:
        input:
            fastq=get_rep_files_for_composite_sample
        output:
            fastq="{wd}/{omics}/{location}/{merged_sample}/{merged_sample}.{pair}.fq.gz"
        shadow:
            "minimal"
        wildcard_constraints:
            merged_sample = '|'.join(merged_illumina_samples.keys()),
            location      = r'5-1-sortmerna|4-hostfree'
        shell:
            """
            cat {input.fastq} > combined.fq.gz
            rsync -a combined.fq.gz {output.fastq}
            """

###############################################################################################
# Assembly-free taxonomy profiling
###############################################################################################

# To enable multiple versions of taxonomy profiles for the same project, we include {version} in taxonomy profile output file name.
# But changing '{sample}.{taxonomy}' to '{sample}.{taxonomy}.{version}' leads to trouble as metaphlan's combining script infers the
# sample name by removing the word after the last dot. If we named files as 'D1.metaphlan.4.0.6', then the combined table lists this
# sample as 'D1.metaphlan.4.0'. This disagrees with the sample metadata and difficult to recover. So we now add '.tsv' extension to
# profile output, let metaphlan remove the '.tsv', and then remove '{taxonomy}.{version}' ourselves. This is the history behind
# naming profile outputs as '{sample}.{taxonomy}.{version}.tsv'.

# MetaPhlAn cannot take multiple runs per sample.
# So named pipes will be used to combine runs into one fastq file.

rule metaphlan_tax_profile:
    input:
        metaphlan_db=lambda wildcards: expand("{minto_dir}/data/metaphlan/{version}/{metaphlan_index}_VINFO.csv",
                                                minto_dir=minto_dir,
                                                version=wildcards.version,
                                                metaphlan_index=metaphlan_index),
        fwd=get_qc2_output_files_fwd_only,
        rev=get_qc2_output_files_rev_only,
    output:
        ra="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.metaphlan.{version}.tsv"
    shadow:
        "minimal"
    params:
        input_files = lambda wildcards, input: ','.join(input.fwd + input.rev)
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.metaphlan.{version}.log"
    conda:
        minto_dir + "/envs/metaphlan.yml" #metaphlan
    shell:
        """
        remote_dir=$(dirname {output.ra})

        time (
            metaphlan {params.input_files} \
                    --input_type fastq \
                    -o {wildcards.sample}.metaphlan.{wildcards.version}.tsv \
                    --bowtie2out {wildcards.sample}.metaphlan.{wildcards.version}.bowtie2.bz2 \
                    --bowtie2db {minto_dir}/data/metaphlan/{wildcards.version} \
                    --index {metaphlan_index} \
                    --nproc {threads} \
                    -t rel_ab_w_read_stats \
                    --unclassified_estimation
            rsync -a {wildcards.sample}.metaphlan.* $remote_dir
        ) >& {log}
        """

rule metaphlan_combine_profiles:
    input:
        ra=lambda wildcards: expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxonomy}.{version}.tsv", wd = wildcards.wd, omics = wildcards.omics, taxonomy = wildcards.taxonomy, version = wildcards.version, sample = nonredundant_ilmn_samples),
    output:
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.merged_abundance_table.txt",
        species="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.merged_abundance_table_species.txt",
    wildcard_constraints:
        taxonomy = r'metaphlan'
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}.{version}_combine.log"
    conda:
        minto_dir + "/envs/metaphlan.yml" #metaphlan
    shell:
        """
        time (
            merge_metaphlan_tables.py {input.ra} | sed -e 's/\.{wildcards.taxonomy}\.{wildcards.version}//g' -e 's/^UNCLASSIFIED/Unknown/' > {output.merged}
            grep -E "s__|clade_name|Unknown" {output.merged} | grep -v "t__" | sed 's/^.*s__//' | sed 's/^clade_name/species/' > {output.species}
        ) >& {log}
        """

# motus can take in multiple runs as comma-separated files.
# So we just construct it in {params}.
rule motus_map_db:
    input:
        db = f"{minto_dir}/data/motus/{{version}}/db_mOTU/db_mOTU_versions",
        fwd=get_qc2_output_files_fwd_only,
        rev=get_qc2_output_files_rev_only
    output:
        mgc="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus.{version}.mgc"
    shadow:
        "minimal"
    params:
        fwd_files = lambda wildcards, input: ",".join(input.fwd),
        rev_files = lambda wildcards, input: ",".join(input.rev),
        motus_db  = lambda wildcards, input: os.path.dirname(input.db)
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.motus.{version}.mgc.log"
    conda:
        minto_dir + "/envs/motus_env.yml" #motus3
    shell:
        """
        time (
            motus map_tax   -t {threads} -f {params.fwd_files} -r {params.rev_files} -db {params.motus_db} -o {wildcards.sample}.motus.bam -b
            motus calc_mgc  -n {wildcards.sample}                                    -db {params.motus_db} -i {wildcards.sample}.motus.bam -o {wildcards.sample}.motus.mgc
            rsync -a {wildcards.sample}.motus.mgc {output.mgc}
        ) >& {log}
        """

rule motus_calc_motu:
    input:
        db = f"{minto_dir}/data/motus/{{version}}/db_mOTU/db_mOTU_versions",
        mgc=rules.motus_map_db.output.mgc,
    output:
        raw="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus_raw.{version}.tsv",
        rel="{wd}/{omics}/6-taxa_profile/{sample}/{sample}.motus_rel.{version}.tsv"
    params:
        motus_db = lambda wildcards, input: os.path.dirname(input.db)
    resources:
        mem=TAXA_memory
    threads: 2
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{sample}.motus.{version}.log"
    conda:
        minto_dir + "/envs/motus_env.yml" #motus3
    shell:
        """
        time (
            motus calc_motu -n {wildcards.sample} -db {params.motus_db} -i {input.mgc} -o {output.rel} -p -q
            motus calc_motu -n {wildcards.sample} -db {params.motus_db} -i {input.mgc} -o {output.raw} -p -q -c
        ) >& {log}
        """

rule motus_combine_profiles:
    input:
        db = f"{minto_dir}/data/motus/{{version}}/db_mOTU/db_mOTU_versions",
        profiles=lambda wildcards: expand("{wd}/{omics}/6-taxa_profile/{sample}/{sample}.{taxonomy}.{version}.tsv", wd = wildcards.wd, omics = wildcards.omics, taxonomy = wildcards.taxonomy, version = wildcards.version, sample = nonredundant_ilmn_samples),
    output:
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.merged_abundance_table.txt",
        species="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.merged_abundance_table_species.txt",
    wildcard_constraints:
        taxonomy = r'motus_(raw|rel)'
    params:
        motus_db   = lambda wildcards, input: os.path.dirname(input.db),
        cut_fields = '1,3-',
        files      = lambda wildcards, input: ",".join(input.profiles)
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}.{version}_combine.log"
    conda:
        minto_dir + "/envs/motus_env.yml" #motus3
    shell:
        """
        time (
            motus merge -db {params.motus_db} -i {params.files} | sed 's/^\(\S*\)\s\([^\\t]\+\)\\t/\\2 [\\1]\\t/' | sed -e 's/^consensus_taxonomy \[#mOTU\]/clade_name/' -e 's/NCBI_tax_id/clade_taxid/' -e 's/^unassigned .unassigned./Unknown/' | cut -f{params.cut_fields}  > {output.merged}
            grep -E "s__|clade_name" {output.merged} | sed 's/^.*s__//' | sed 's/^clade_name/species/' > {output.species}
        ) >& {log}
        """

rule plot_taxonomic_profile:
    input:
        merged="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.merged_abundance_table.txt",
    output:
        profile="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.tsv",
        pcoa="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.PCoA.Bray_Curtis.pdf",
        barplot="{wd}/output/6-taxa_profile/{omics}.{taxonomy}.{version}.Top15genera.pdf",
    wildcard_constraints:
        taxonomy = r'motus_(raw|rel)|metaphlan'
    params:
        plot_args=plot_args_str
    threads: 1
    log:
        "{wd}/logs/{omics}/6-taxa_profile/{taxonomy}.{version}_plot.log"
    conda:
        minto_dir + "/envs/r_pkgs.yml" #R
    shell:
        """
        time (
            Rscript {script_dir}/plot_6_taxa_profile.R --table {input.merged} --profiler {wildcards.omics}.{wildcards.taxonomy}.{wildcards.version} --metadata {metadata} --outdir $(dirname {output.pcoa}) {params.plot_args}
        ) >& {log}
        """

###############################################################################################
# Assembly-free sample comparison
###############################################################################################

# We are using sourmash to create weighted sketches from the host-free metaG fastq files.
# If there are multiple runs for a sample, they will be included in the same sketch.

rule sourmash_sketch:
    input:
        fwd=get_qc2_output_files_fwd_only,
        rev=get_qc2_output_files_rev_only,
    output:
        temp("{wd}/{omics}/6-1a-smash/{sample}.sig.gz")
    shadow:
        "minimal"
    params:
        k=21,
        scaled=100
    resources:
        mem=TAXA_memory
    threads:
        2
    log:
        "{wd}/logs/{omics}/6-1-smash/{sample}.sourmash.sketch.log"
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        time (
            sourmash sketch dna -p k={params.k},scaled={params.scaled},abund --name {wildcards.sample} -o {output} {input.fwd} {input.rev}
        ) >& {log}
        """

rule sourmash_filter:
    input:
        "{wd}/{omics}/6-1a-smash/{sample}.sig.gz"
    output:
        "{wd}/{omics}/6-1-smash/{sample}.sig.gz"
    shadow:
        "minimal"
    params:
        k=21,
        m=config['SOURMASH_min_abund'],
        M=config['SOURMASH_max_abund']
    resources:
        mem=TAXA_memory
    threads:
        2
    log:
        "{wd}/logs/{omics}/6-1-smash/{sample}.sourmash.filter.log"
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        time (
            sourmash signature filter -k {params.k} -m {params.m} -M {params.M} -o {output} {input}
        ) >& {log}
        """

rule sourmash_compare:
    input:
        expand("{wd}/{omics}/6-1-smash/{sample}.sig.gz",
                wd = working_dir,
                omics = omics,
                sample = nonredundant_ilmn_samples)
    output:
        csv="{wd}/{omics}/6-1-smash/{omics}.sourmash_cosine_similarity.csv",
        npy="{wd}/{omics}/6-1-smash/{omics}.sourmash_cosine_similarity.npy"
    shadow:
        "minimal"
    params:
        k=21
    resources:
        mem=TAXA_memory
    threads:
        TAXA_threads
    log:
        "{wd}/logs/{omics}/6-1-smash/{omics}.sourmash.compare.log"
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        time (
            sourmash compare -k {params.k} -p {threads} --csv {output.csv} -o {output.npy} {input}
        ) >& {log}
        """

rule plot_sourmash_kmers:
    input:
        csv = expand("{wd}/{omics}/6-1-smash/{omics}.sourmash_cosine_similarity.csv",
                wd = working_dir,
                omics = omics),
        npy = expand("{wd}/{omics}/6-1-smash/{omics}.sourmash_cosine_similarity.npy",
                wd = working_dir,
                omics = omics),
        merged = "{wd}/output/6-taxa_profile/{omics}.{taxonomy_versioned}.merged_abundance_table.txt",
    output:
        barplot="{wd}/output/6-1-smash/{omics}.{taxonomy_versioned}.clusters.pdf",
        tsv="{wd}/output/6-1-smash/{omics}.{taxonomy_versioned}.sourmash_clusters.tsv"
    wildcard_constraints:
        omics = r'metaG|metaT'
    params:
        cutoff=config['SOURMASH_cutoff'],
        plot_args=plot_args_str
    threads:
        1
    log:
        "{wd}/logs/{omics}/6-1-smash/{omics}.{taxonomy_versioned}.sourmash.plot.log"
    conda:
        minto_dir + "/envs/r_pkgs.yml" #R
    shell:
        """
        time (
            Rscript {script_dir}/plot_6-1_sourmash.R --csv {input.csv} --cutoff {params.cutoff} --table {input.merged} --metadata {metadata} --outdir $(dirname {output.barplot}) {params.plot_args}
        ) >& {log}
        """

rule dummy_sourmash_clusters:
    input:
        "{wd}/output/6-1-smash/{omics}.{taxonomy_versioned}.sourmash_clusters.tsv".format(
                        wd = working_dir,
                        omics = omics,
                        taxonomy_versioned = taxonomies_versioned[0])
    output:
        "{wd}/output/6-1-smash/{omics}.sourmash_clusters.tsv"
    threads:
        1
    localrule: True
    shell:
        """
        mv {input} {output}
        """

##########################################################################################################
# Generate configuration yml file for recovery of MAGs and taxonomic annotation step - assembly/coassembly
##########################################################################################################

rule qc2_filter_config_yml_assembly:
    input:
        tax=lambda wildcards: expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.tsv",
                                    wd = wildcards.wd,
                                    omics = wildcards.omics,
                                    taxonomy = taxonomies_versioned),
        metadata=metadata,
        table="{wd}/output/6-1-smash/{omics}.sourmash_clusters.tsv"
    output:
        config_file="{wd}/{omics}/assembly.yaml"
    params:
        sample_list = nonredundant_ilmn_samples,
        sample_string = ', '.join(["'{}'".format(i) for i in nonredundant_ilmn_samples]),
        merge_illumina_samples_directive = '\n'.join([" {} : {}".format(i, merged_illumina_samples[i]) for i in merged_illumina_samples.keys()])
    resources:
        mem=2
    threads: 2
    log:
        "{wd}/logs/{omics}/config_yml_assembly.log"
    conda:
        minto_dir + "/envs/r_pkgs.yml"
    shell:
        """
        cat > {output} <<___EOF___
######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
minto_dir: {minto_dir}
METADATA: {input.metadata}

######################
# Analysis settings
######################

MAIN_factor: {main_factor}

######################
# Program settings
######################
# MetaSPAdes settings
#
METASPADES_qoffset: auto
METASPADES_threads: 16
METASPADES_memory: 10
METASPADES_hybrid_max_k: 99
METASPADES_illumina_max_k: 99

# MEGAHIT settings
#
# Note on MEGAHIT_memory:
#     MEGAHIT's memory requirement scales with sample count and sample size.
#     Default MEGAHIT_memory below is 10 Gigabytes per sample (plenty for metaG samples of 10M reads).
#     MEGAHIT_memory parameter represents per-sample memory in Gigabytes
#     during the first attempt. If it is not enough and MEGAHIT fails, then
#     each successive attempt will increase per-sample memory by 6 Gibabytes
#     until it succeeds or reaches max_attempts from snakemake.
#     Please make sure that there is enough RAM on the server.
# Custom k-list for assembly: k must be odd, in the range 15-255, increment <= 28, fx. 21,29,39,59,79,99,119,141
MEGAHIT_memory: 10
MEGAHIT_threads: 32
MEGAHIT_presets:
 - meta-sensitive
 - meta-large
#MEGAHIT_custom:
# - 21,29,39,59,79,99,119,141

# MetaFlye settings
#
# MetaFlye will be run for each parameter preset listed here.
# By default Flye will be run with these parameters:
#    --meta --genome-size 3.0m
# If you need to add more options, define them here and name them for future reference.
# Notes:
# ------
# 1. Each preset parameter will be applied to each sample. If you only
#    want one parameter to be used, please comment everything else.
# 2. If nothing is listed here, then MetaFlye won't be run.
#    If you just want our default parameters above, then here is a possible option:
#      metaflye-default: ""
# 3. 'tres-o3000-3x' is valid for flye 2.8.3. From 2.9.x, --plasmids and --trestle are
#    not valid. So please use valid options if you are using newer versions of flye.
METAFLYE_presets:
  tres-o3000-3x: --min-overlap 3000 --iterations 3
  #metaflye-default: ""

###############################
# Binning preparation settings
###############################

# Whether to use contigs or scaffolds from SPAdes
SPADES_CONTIGS_OR_SCAFFOLDS: contigs

# minimum contig/scaffold fasta length
MIN_FASTA_LENGTH: 2500

# maximum available run for a job to calculate assembly batch size for mapping reads to combined contig sets
MAX_RAM_GB_PER_JOB: 180

# Should we exclude any assembly type during MAG generation?
# E.g., if you make MAGs from metaT, individual sample assemblies
# will be quite fragmented. If you have several co-assemblies, then
# ignoring 'illumina_single' might improve your MAG quality. This can
# be achieved using:
# ---
# EXCLUDE_ASSEMBLY_TYPES:
#     - illumina_single
# ---
#
# EXCLUDE_ASSEMBLY_TYPES:

# Contig-depth: bwa-mem2 settings
# Used when mapping reads back to contigs
#
# BWA Alignment
BWA_threads: 10

# Contig-depth: samtools sort settings
# Used when sorting bam files using 3 threads
# Memory listed below is PER-THREAD, so please make sure you have enough
SAMTOOLS_sort_perthread_memgb: 10

###############################
# Input data
###############################

# HYBRID section:
# ---------------
# MetaSPAdes hybrid assembly will be performed using these definitions.
# Definition format:
#   Each nanopore sample is in the LHS, and corresponding illumina sample(s) are in RHS (delimited by '+').
#   Hybrid assemblies will be performed for each combination of nanopore and illumina samples.
#   E.g.:
#
#   N1: I3+I4
#
#   The above will result in 2 hybrid assemblies: 'N1-I3' and 'N1-I4'
#
#   This definition was designed based on our practice of pooling multiple samples derived from the same
#   donor for nanopore sequencing to strike a balance between price and sensitivity. If your nanopore samples
#   are paired one-to-one with illumina samples, then the definition gets simpler.
#   E.g.:
#
#   N1: I3
#   N2: I4
#
#   The above will result in 2 hybrid assemblies: 'N1-I3' and 'N2-I4'
#
#HYBRID:

# NANOPORE section:
# -----------------
# List of nanopore samples that will be assembled individually using MetaFlye.
#
#NANOPORE:

######################
# Optionally, do you want to merge replicates or make pseudo samples
# E.g:
# MERGE_ILLUMINA_SAMPLES:
#  sample1: rep1a+rep1b+rep1c
#  sample2: rep2a+rep2b+rep2c
#
# The above directive will make 2 new composite or pseudo samples at the end of QC_2.
# Imagine you had triplicates for sample1 named as rep1a, rep1b and rep1c.
# And likewise for sample2. The directive above will:
#     - combine 3 samples (rep1a, rep1b and rep1c)into a new sample called 'sample1'.
#     - combine 3 samples (rep2a, rep2b and rep2c)into a new sample called 'sample2'.
# For all subsequent steps, namely:
#     - profiling (done within this snakemake script),
#     - assembly (done by assembly.smk),
#     - binning (done by binning_preparation.smk and mags_generation.smk),
# you can just use 'sample2' instead of the replicates rep2a, rep2b and rep2c in the yaml files.
# Please note that METADATA file must have an entry for 'sample2' as well,
# otherwise QC_2 step will fail.
# Having extra entries in METADATA file does not affect you in any way.
# Therefore, it is safe to have metadata recorded for
# rep2a, rep2b, rep2c, sample2 from the beginning.
######################

MERGE_ILLUMINA_SAMPLES:
{params.merge_illumina_samples_directive}

# ILLUMINA section:
# -----------------
# List of illumina samples that will be assembled individually using MetaSPAdes.
#
# E.g.:
# - I1
# - I2
#
ILLUMINA:
$(for i in {params.sample_list}; do echo "- '$i'"; done)

###############################
# COASSEMBLY section:
###############################

# If enable_COASSEMBLY is "yes", MEGAHIT coassembly will be performed using the following definitions.
# Each coassembly is named in the LHS, and corresponding illumina sample(s) are in RHS (delimited by '+').
# One coassembly will be performed for each line.
# E.g. 'Subject1: I3+I4' will result in 1 coassembly: 'Subject1' using I3 and I4 data.
# Memory per coassembly is calculated to be 10G per sample in the coassembly.
# Please make sure that there is enough RAM on the server.
#
enable_COASSEMBLY: no
COASSEMBLY:
  'Full': $(echo {params.sample_list} | sed 's/ /+/g')
___EOF___

        R --vanilla --silent --no-echo >> {output} <<___EOF___
library(dplyr)
metadata <- read.table('{input.metadata}', sep="\\t", header=TRUE) %>%
    as.data.frame() %>%
    select(sample, {coas_factor}) %>%
    filter(sample %in% c({params.sample_string})) %>%
    group_by({coas_factor}) %>%
    filter(n() > 1) %>%
    mutate(co_asm = paste(sample, collapse = "+")) %>%
    select(-sample) %>%
    slice(1) %>%
    mutate({coas_factor}=paste0("  '", {coas_factor}, "'"))
write.table(metadata, file="", col.names=FALSE, row.names=FALSE, quote=FALSE, sep=": ")
___EOF___

if [[ "metaG" == "{wildcards.omics}" ]]; then
        R --vanilla --silent --no-echo >> {output} <<___EOF___
library(dplyr)
metadata <- read.table('{input.table}', sep="\\t", header=TRUE) %>%
    as.data.frame() %>%
    select(sample, clustering) %>%
    group_by(clustering) %>%
    filter(n() > 1) %>%
    mutate(co_asm = paste(sample, collapse = "+")) %>%
    select(-sample) %>%
    slice(1) %>%
    mutate(clustering=paste0("  SCL", clustering))
write.table(metadata, file="", col.names=FALSE, row.names=FALSE, quote=FALSE, sep=": ")
___EOF___
fi

        echo {params.sample_list} >& {log}
        """

###############################################################################################
# Generate configuration yml file for Alignment, normalization and integration step
###############################################################################################

rule qc2_filter_config_yml_mapping:
    input:
        tsv = lambda wildcards: expand("{wd}/output/6-taxa_profile/{omics}.{taxonomy}.tsv",
                                    wd = wildcards.wd,
                                    omics = wildcards.omics,
                                    taxonomy = taxonomies_versioned),
        metadata = metadata
    output:
        config_file="{wd}/{omics}/mapping.yaml"
    params:
        sample_list = nonredundant_ilmn_samples
    resources:
        mem=2
    threads: 2
    log:
        "{wd}/logs/{omics}/config_yml_mapping.log"
    shell:
        """
        cat > {output} <<___EOF___
######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
minto_dir: {minto_dir}
METADATA: {input.metadata}

######################
# Analysis settings
######################

MAIN_factor: {main_factor}
PLOT_factor2: {plot_factor2}
PLOT_time: {plot_time}

######################
# Annotation settings
######################

# Set MIntO mode
# Where should we map reads to? MAG, refgenome, catalog
MINTO_MODE: MAG

# Which omics for MAGs?
MAG_omics: metaG

# path to gene catalog fasta file or refgenome directory
PATH_reference:

# file name of gene catalog fasta file (MIntO will generate bwa index with same name)
NAME_reference:

# List of software used to perform genome function annotation:
# - dbCAN
# - kofam
# - eggNOG
ANNOTATION:
 - dbCAN
 - kofam
 - eggNOG

#########################
# Gene abundance settings
#########################

# BWA Alignment
BWA_threads: 10

# Alignment filtering
# -------------------
# msamtools_filter_length: Default 50 works well for ILLUMINA 2x150bp.
#                          For shorter sequencing, e.g. 2x75bp, please reduce accordingly.
# alignment_identity: Default 97 works well for MAGs that are dereped at 99%.
#                     For refgenomes, modify it according to the derep level of your refgenomes.
#                     For gene catalog, modify it to match the derep level of the gene catalog (usually 95%).
# MIN_mapped_reads:   Default 2 is designed for low-to-medium sequencing depth.
#                     Increase according to depth (e.g., set to 10 when depth is above 5Gb).
msamtools_filter_length: 50
alignment_identity: 97
MIN_mapped_reads: 2

# Normalization approach
# Could be TPM, MG or comma-delimited combinations
abundance_normalization: TPM,MG

##########################
# Genome taxonomy settings
##########################
#
RUN_TAXONOMY: yes
TAXONOMY_NAME: phylophlan,gtdb  # Currently, phylophlan or gtdb or combination
TAXONOMY_CPUS: 8
TAXONOMY_memory: 5

# Taxonomy database versions
#
PHYLOPHLAN_TAXONOMY_VERSION: SGB.Jun23
GTDB_TAXONOMY_VERSION: r220

# Input data

# ILLUMINA section:
# -----------------
# List of illumina samples.
#
# E.g.:
# - I1
# - I2
#
ILLUMINA:
$(for i in {params.sample_list}; do echo "- '$i'"; done)
___EOF___

        echo {params.sample_list} >& {log}
        """
