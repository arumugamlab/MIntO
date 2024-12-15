#!/usr/bin/env python

'''
Alignment, normalization and integration step

Authors: Carmen Saenz, Mani Arumugam
'''

import os.path
import pathlib
import math

localrules: make_merged_genome_fna, make_genome_def, \
            config_yml_integration, read_map_stats

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'
include: 'include/resources.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule abundance_base, abundance_rpkg from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# We prefer original non-error-corrected reads for profiling to preserve strain variation in population
read_dir=get_qc2_output_location(omics)

# If the original reads have been removed for whatever reasons (running out of space, may be?),
# then we look for error-corrected reads
if os.path.exists("{}/{}/{}/".format(working_dir, omics, read_dir)) is False:
    if os.path.exists("{}/{}/{}/".format(working_dir, omics, '6-corrected')) is True:
        read_dir='6-corrected'
    else:
        raise Exception("ERROR in {}: One of {} or 6-corrected must exist. Please correct {}".format(config_path, get_qc2_output_location(omics), working_dir))

print("NOTE: MIntO is using '{}' as read directory".format(read_dir))

main_factor      = validate_required_key(config, 'MAIN_factor')
plot_factor2     = validate_optional_key(config, 'PLOT_factor2')
plot_time        = validate_optional_key(config, 'PLOT_time')
MIN_mapped_reads = validate_required_key(config, 'MIN_mapped_reads')

# Make list of illumina samples, if ILLUMINA in config

ilmn_samples = list()
if (x := validate_required_key(config, 'ILLUMINA')):
    check_fastq_file_locations(x, locations = [read_dir])
    ilmn_samples = x

# MIntO mode and database-mapping

# Which mode are we running?
MINTO_MODE = validate_required_key(config, 'MINTO_MODE')

# Define the 3 modes
valid_minto_modes = ['MAG', 'refgenome', 'catalog']

# Backward compatibility and common misnomers
if MINTO_MODE in ['db_genes', 'db-genes', 'genes_db', 'gene_catalog', 'gene-catalog']:
    MINTO_MODE = 'catalog'
elif MINTO_MODE in ['reference_genome', 'reference-genome', 'reference', 'refgenomes']:
    MINTO_MODE = 'refgenome'
elif MINTO_MODE in ['MAGs', 'mag', 'mags']:
    MINTO_MODE = 'MAG'

check_allowed_values('MINTO_MODE', MINTO_MODE, valid_minto_modes)

# Normalization
normalization = validate_required_key(config, 'abundance_normalization')
normalization_modes = normalization.split(",")
for m in normalization_modes:
    check_allowed_values('abundance_normalization', m, ("MG", "TPM"))

# Alignment filtering
identity = validate_required_key(config, 'alignment_identity')
msamtools_filter_length = validate_required_key(config, 'msamtools_filter_length')

mag_omics = 'metaG'
gene_catalog_path = None
gene_catalog_name = None

if MINTO_MODE == 'MAG':
    mag_omics = 'metaG'
    if (x := validate_optional_key(config, 'MAG_omics')):
        mag_omics = x
    reference_dir = f"{working_dir}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes"
    print('NOTE: MIntO is using "' + reference_dir + '" as PATH_reference variable')
else:
    if (x := validate_required_key(config, 'PATH_reference')):
        if MINTO_MODE == 'refgenome':
            reference_dir = x
            print('NOTE: MIntO is using "'+ reference_dir+'" as PATH_reference variable')
        elif MINTO_MODE == 'catalog':
            gene_catalog_path = x
            if (gene_catalog_name := validate_required_key(config, 'NAME_reference')):
                if os.path.exists(gene_catalog_path + '/' + gene_catalog_name):
                    print(f"NOTE: MIntO is using {gene_catalog_path}/{gene_catalog_name} as gene database.")

BWA_threads = validate_required_key(config, 'BWA_threads')

if 'MG' in normalization_modes and MINTO_MODE == 'catalog':
    raise Exception("ERROR in {}: In 'catalog' mode, only TPM normalization is allowed.".format(config_path))

# Taxonomic profiles from mapping reads to MAGs or refgenomes

taxonomies_versioned = list()

run_taxonomy = False
if (x := validate_optional_key(config, 'RUN_TAXONOMY')):
    run_taxonomy = x

if run_taxonomy:
    TAXONOMY = validate_required_key(config, 'TAXONOMY_NAME')
    allowed = ('phylophlan', 'gtdb')
    for x in TAXONOMY.split(","):
        check_allowed_values('TAXONOMY_NAME', x, allowed)

    taxonomies = TAXONOMY.split(",")
    for t in taxonomies:
        version="unknown"
        if t == "phylophlan":
            version = validate_required_key(config, "PHYLOPHLAN_TAXONOMY_VERSION")
        elif t == "gtdb":
            version = validate_required_key(config, "GTDB_TAXONOMY_VERSION")
        taxonomies_versioned.append(t+"."+version)
    print('NOTE: MIntO is using taxonomy labelling of the unique MAGs from [{}].'.format(", ".join(taxonomies_versioned)))

# Site customization for avoiding NFS traffic during I/O heavy steps such as mapping

CLUSTER_NODES            = None
CLUSTER_LOCAL_DIR        = None
CLUSTER_WORKLOAD_MANAGER = None
include: minto_dir + '/site/cluster_def.py'

# Cluster-aware bwa-index rules

include: 'include/bwa_index_wrapper.smk'

######
# Make a lookup table for sample --> sample_alias using metadata file
######

def get_sample2alias_map(in_file):
    import pandas
    df = pandas.read_csv(in_file, comment='#', header=0, sep = "\t")
    if 'sample' not in df.columns:
        raise Exception("Need 'sample' column in metadata file")
    if 'sample_alias' not in df.columns:
        raise Exception("Need 'sample_alias' column in metadata file")
    df = df.set_index('sample')['sample_alias']
    return(df.to_dict())

sample2alias = get_sample2alias_map(metadata)

# Define all the outputs needed by target 'all'

def combined_genome_profiles_annotated():
    result = expand("{wd}/output/9-mapping-profiles/{mode}/{omics}.{taxonomy}.p{identity}.tsv",
                wd = working_dir,
                omics = omics,
                mode = MINTO_MODE,
                taxonomy = taxonomies_versioned,
                identity = identity)
    return(result)

def combined_genome_profiles():
    folder = "{wd}/{omics}/9-mapping-profiles/{subdir}".format(
                wd = working_dir,
                omics = omics,
                subdir = MINTO_MODE,
                )
    files = expand("{label}.___ID___.{variant}.tsv",
                zip,
                label = ['genome_abundances', 'genome_abundances', 'contig_abundances', 'contig_abundances'],
                variant = ['abund.prop', 'relabund.prop', 'abund.prop', 'abund.all']
                )
    result = [folder + '/' + f.replace('___ID___', f"p{identity}") for f in files]
    return(result)

def combined_gene_abundance_profiles():
    result = expand("{wd}/{omics}/9-mapping-profiles/{subdir}/{label}.p{identity}.{norm}.tsv",
                label = 'gene_abundances',
                wd = working_dir,
                omics = omics,
                subdir = MINTO_MODE,
                norm = normalization_modes,
                identity = identity)
    return(result)

if MINTO_MODE == 'catalog':
    reference_dir = gene_catalog_path
    def combined_genome_profiles():
        return()
    def combined_genome_profiles_annotated():
        return()

def mapping_statistics():
    result = expand("{wd}/{omics}/9-mapping-profiles/{subdir}/mapping.p{identity}.{stats}.txt",
                    wd = working_dir,
                    omics = omics,
                    subdir = MINTO_MODE,
                    identity = identity,
                    stats = ['maprate', 'mapstats', 'multimap'])
    return(result)

def config_yaml():
    result = "{wd}/data_integration.yaml".format(
                wd = working_dir)
    return(result)

# If BWA index clean-up requested, only do cleanup and nothing else
if CLUSTER_NODES != None and validate_optional_key(config, 'CLEAN_BWA_INDEX'):

    print("NOTE: BWA index cleanup mode requested.")

    def clean_bwa_index():
        if MINTO_MODE == 'catalog':
            return(f"{gene_catalog_path}/BWA_index/{gene_catalog_name}.clustersync/cleaning.done")
        else:
            return(f"{working_dir}/DB/{MINTO_MODE}/BWA_index/{MINTO_MODE}.fna.clustersync/cleaning.done")

    rule all:
        input:
            clean_bwa_index()
        default_target: True

else:
    rule all:
        input:
            combined_genome_profiles_annotated(),
            combined_genome_profiles(),
            combined_gene_abundance_profiles(),
            mapping_statistics(),
            config_yaml(),
            print_versions.get_version_output(snakefile_name)
        default_target: True

###############################################################################################
# Functions to get samples, runs and files
###############################################################################################

# Get a sorted list of runs for a sample

def get_runs_for_sample(wildcards):
    if read_dir == '6-corrected':
        runs = [wildcards.sample]
    else:
        sample_dir = '{wd}/{omics}/{location}/{sample}'.format(wd=wildcards.wd, omics=wildcards.omics, location=read_dir, sample=wildcards.sample)
        runs = sorted([ re.sub("\.1\.fq\.gz", "", os.path.basename(f)) for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith('.1.fq.gz') ])
        #print(runs)
    return(runs)

# Get a sorted list of runs for a sample

def get_fwd_files_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = read_dir,
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '1')
    return(files)

def get_rev_files_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = read_dir,
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '2')
    return(files)

###############################################################################################
# MIntO modes: MAG and refgenome
###############################################################################################

#########################
# Generate combined genome
#########################

# Get a sorted list of genomes

def get_genomes_from_refdir(ref_dir):
    genomes = [ pathlib.Path(f).stem for f in os.scandir(ref_dir) if f.is_file() and f.name.endswith('.fna') ]
    return(sorted(genomes))

def get_genome_fna(wildcards):
    #Collect the fna files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    result = expand("{wd}/DB/{minto_mode}/2-postprocessed/{genome}.fna",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(result)

rule make_merged_genome_fna:
    input: get_genome_fna
    output:
        fasta_merge="{wd}/DB/{minto_mode}/{minto_mode}.fna"
    log:
        "{wd}/logs/DB/{minto_mode}/{minto_mode}.merge_genome.log"
    wildcard_constraints:
        minto_mode = r'MAG|refgenome'
    shell:
        """
        time (
            cat {input} > {output}
        ) >& {log}
        """

# For each sequence name '<seqname> = <locustag>_<seqid>', make '<locustag>\t<seqname>'
rule make_genome_def:
    input: get_genome_fna
    output:
        genome_def="{wd}/DB/{minto_mode}/{minto_mode}.genome.def"
    shell:
        """
        cat {input} | grep '^>' | sed -e 's/^>\([^_]\+\)_\(.*\)/\\1\\t\\1_\\2/' > {output}
        """

################################################################################################
# Get BWA index for the filtered contigs
################################################################################################

# If CLUSTER_NODES is defined, then return the file symlink'ed to local drives.
# Otherwise, return the original index files in project work_dir.
def get_genome_bwa_index(wildcards):

    # Where are the index files?
    if CLUSTER_NODES != None:
        index_location = 'BWA_index_local'
    else:
        index_location = 'BWA_index'

    # Get all the index files!
    files = expand("{wd}/DB/{minto_mode}/{location}/{minto_mode}.fna.{ext}",
            wd         = wildcards.wd,
            minto_mode = wildcards.minto_mode,
            location   = index_location,
            ext        = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'])

    # Return them
    return(files)

#########################
# Map reads using BWA2, profile using msamtools, sort bam file, and create gene-level BED file.
# Mapping and computation of read counts using samtools bedcov are done in the same rule.
# Mapping and sorting are in separate processes, so they don't need to share memory.
# Memory estimates:
# 1. bwa-mem2 mem
#   BWA mem memory is estimated as 3.1 bytes per base in database with 5.6GB offset (regression: mem = 5.556e+09 + 3.011*input).
# 2. samtools sort
#   Sorting uses a hard-coded 3 threads, so we divide the total memory by 3 and send it to 'samtools sort'.
#   But we recommend at least 5GB per thread, even if the database is smaller.
#   That's where the 'mem = max(3*5, x)' comes from.
#   And, for whatever reason, 'samtools sort' uses ~10% more memory than what you allocated, so we adjust for this behavior.
#   That's where the 0.9 in 'sort_memory=$((9*{resources.mem}/{resources.sort_threads}/10))' comes from.
#########################

rule genome_mapping_profiling:
    input:
        bwaindex=get_genome_bwa_index,
        genome_def=rules.make_genome_def.output.genome_def,
        fwd=get_fwd_files_only,
        rev=get_rev_files_only,
        bed_mini="{wd}/DB/{minto_mode}/{minto_mode}-genes.bed.mini"
    output:
        bwa_log=        "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.bwa.log",
        raw_all_seq=    "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
        raw_prop_seq=   "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.txt.gz",
        raw_prop_genome="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.genome.txt.gz",
        rel_prop_genome="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.relabund.prop.genome.txt.gz",
        absolute_counts="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.gene_abundances.p{identity}.bed.gz"
    shadow:
        "minimal"
    params:
        length = msamtools_filter_length,
        mapped_reads_threshold = MIN_mapped_reads,
        bedcov_lines = 500000,
        sample_alias = lambda wildcards: sample2alias[wildcards.sample],
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{sample}.p{identity}.map_profile.log"
    wildcard_constraints:
        identity   = r'\d+',
        minto_mode = r'MAG|refgenome'
    threads:
        BWA_threads
    resources:
        sort_threads = 3,
        bedcov_threads = lambda wildcards, threads: min(10, threads),
        mem = lambda wildcards, input, attempt: max(3*5, math.ceil(5.6 + 3.1e-9*get_file_size(input.bwaindex[0]))) + 10*(attempt-1)
    conda:
        minto_dir + "/envs/MIntO_base.yml" # BWA + samtools + msamtools + perl + parallel
    shell:
        """
        # Make named pipes if needed
        if [ "{params.multiple_runs}" == "yes" ]; then
            mkfifo {wildcards.sample}.1.fq.gz
            mkfifo {wildcards.sample}.2.fq.gz
            cat {input.fwd} > {wildcards.sample}.1.fq.gz &
            cat {input.rev} > {wildcards.sample}.2.fq.gz &
            input_files="{wildcards.sample}.1.fq.gz {wildcards.sample}.2.fq.gz"
        else
            input_files="{input.fwd} {input.rev}"
        fi

        # Get index file name
        bwaindex_prefix={input.bwaindex[0]}
        bwaindex_prefix=${{bwaindex_prefix%.0123}}

        # Estimate samtools sort memory
        sort_memory=$((9*{resources.mem}/{resources.sort_threads}/10))
        echo "Using {resources.sort_threads} threads and $sort_memory GB memory per thread for 'samtools sort'"
        echo "Using {resources.bedcov_threads} threads and {params.bedcov_lines} lines per batch-file for 'samtools bedcov'"

        (time (bwa-mem2 mem -a -t {threads} -v 3 ${{bwaindex_prefix}} $input_files | \
                    msamtools filter -S -b -l {params.length} -p {wildcards.identity} -z 80 --besthit - > aligned.bam) >& {output.bwa_log}
            total_reads="$(grep Processed {output.bwa_log} | perl -ne 'm/Processed (\\d+) reads/; $sum+=$1; END{{printf "%d\\n", $sum/2;}}')"
            #echo $total_reads

            # Run multiple msamtools modes + samtools sort in parallel using GNU parallel.
            # We need 4 threads for msamtools and {resources.sort_threads} threads for samtools sort.
            # We cap it by 'threads' limit to parallel.
            common_args="--label {wildcards.omics}.{params.sample_alias} --total=$total_reads --mincount={params.mapped_reads_threshold} --pandas"
            parallel --jobs {threads} <<__EOM__
msamtools profile aligned.bam $common_args -o {output.raw_all_seq}     --multi=all  --unit=abund --nolen
msamtools profile aligned.bam $common_args -o {output.raw_prop_seq}    --multi=prop --unit=abund --nolen
msamtools profile aligned.bam $common_args -o {output.raw_prop_genome} --multi=prop --unit=abund --nolen --genome {input.genome_def}
msamtools profile aligned.bam $common_args -o {output.rel_prop_genome} --multi=prop --unit=rel           --genome {input.genome_def}
samtools sort aligned.bam -o sorted.bam -@ {resources.sort_threads} -m ${{sort_memory}}G --output-fmt=BAM
__EOM__

            # Use samtools bedcov to generate read-counts per gene
            # 1. Index sorted.bam file
            # 2. Enable parallelization by splitting the BED file into smaller files and then giving it to samtools bedcov
            #    a. Create smaller bed files with {params.bedcov_lines} lines each, like x0000000.bedsub, x0000001.bedsub, ...
            #    b. Send 'ls' output to GNU parallel to process these smaller BED files.

            parallel --jobs {threads} <<__EOM__
samtools index sorted.bam sorted.bam.bai -@ {threads}
split -d --suffix-length=7 --additional-suffix=.bedsub --lines={params.bedcov_lines} {input.bed_mini}
__EOM__

            # Pipe it to parallel
            time (ls | grep -E '\.bedsub$' | parallel --jobs {resources.bedcov_threads} --plus "samtools bedcov -c -g SECONDARY {{}} sorted.bam | cut -f4,5,7 | (grep -v $'\\t0$' || true) > {{.}}.bed.part")
            (echo -e 'gene_length\\tID\\t{wildcards.omics}.{params.sample_alias}'; cat *.bed.part) | gzip -c > out.bed.gz

            # Rsync output bed file
            rsync -a out.bed.gz {output.absolute_counts}
        ) >& {log}
        """

###############################################################################################
# MIntO mode: catalog
###############################################################################################

#########################
# Get bwa index for catalog
#########################

# If CLUSTER_NODES is defined, then return the file symlink'ed to local drives.
# Otherwise, return the original index files in project work_dir.
def get_catalog_bwa_index(wildcards):

    # Where are the index files?
    if CLUSTER_NODES != None:
        index_location = 'BWA_index_local'
    else:
        index_location = 'BWA_index'

    # Get all the index files!
    files = expand("{path}/{location}/{name}.{ext}",
            path       = gene_catalog_path,
            name       = gene_catalog_name,
            location   = index_location,
            ext        = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'])

    # Return them
    return(files)

#########################
# Map reads using BWA2
# Mapping, computation of read counts and TPM normalization is done in the same rule
# TPM normalization: sequence depth and genes’ length
# BWA mem memory is estimated as 3.1 bytes per base in database with 5.6GB offset (regression: mem = 5.556e+09 + 3.011*input).
#########################

rule gene_catalog_mapping_profiling:
    input:
        bwaindex=get_catalog_bwa_index,
        fwd=get_fwd_files_only,
        rev=get_rev_files_only,
    output:
        bwa_log=    "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.log",
        profile_tpm="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.TPM.txt.gz",
        map_profile="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz"
    shadow:
        "minimal"
    params:
        sample_alias=lambda wildcards: sample2alias[wildcards.sample],
        length=msamtools_filter_length,
        mapped_reads_threshold=MIN_mapped_reads,
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{sample}.p{identity}_bwa.log"
    wildcard_constraints:
        minto_mode = r'catalog'
    threads:
        BWA_threads
    resources:
        mem = lambda wildcards, input, attempt: math.ceil(5.6 + 3.1e-9*get_file_size(input.bwaindex[0])) + 10*(attempt-1)
    conda:
        minto_dir + "/envs/MIntO_base.yml" #BWA + samtools
    shell:
        """
        # Make named pipes if needed
        if [ "{params.multiple_runs}" == "yes" ]; then
            mkfifo {wildcards.sample}.1.fq.gz
            mkfifo {wildcards.sample}.2.fq.gz
            cat {input.fwd} > {wildcards.sample}.1.fq.gz &
            cat {input.rev} > {wildcards.sample}.2.fq.gz &
            input_files="{wildcards.sample}.1.fq.gz {wildcards.sample}.2.fq.gz"
        else
            input_files="{input.fwd} {input.rev}"
        fi

        bwaindex_prefix={input.bwaindex[0]}
        bwaindex_prefix=${{bwaindex_prefix%.0123}}
        (time (bwa-mem2 mem -a -t {threads} -v 3 ${{bwaindex_prefix}} $input_files | \
                msamtools filter -S -b -l {params.length} -p {wildcards.identity} -z 80 --besthit - > aligned.bam) >& {output.bwa_log}
            total_reads="$(grep Processed {output.bwa_log} | perl -ne 'm/Processed (\\d+) reads/; $sum+=$1; END{{printf "%d\\n", $sum/2;}}')"
            #echo $total_reads

            # Run multiple msamtools modes in parallel using GNU parallel.
            # We cap it by 'threads' limit to parallel.
            common_args="--label {wildcards.omics}.{params.sample_alias} --total=$total_reads --mincount={params.mapped_reads_threshold} --pandas"
            parallel --jobs {threads} <<__EOM__
msamtools profile aligned.bam $common_args -o profile.TPM.txt.gz --multi prop --unit tpm
msamtools profile aligned.bam $common_args -o profile.abund.all.txt.gz --multi all --unit abund --nolen
__EOM__

            # rsync outputs to output dir
            parallel --jobs {threads} <<__EOM__
rsync -a profile.TPM.txt.gz {output.profile_tpm}
rsync -a profile.abund.all.txt.gz {output.map_profile}
__EOM__
        ) >& {log}
        """

###############################################################################################
# Merge gene or genome profiles
###############################################################################################

# Strip the given common upstream path from a list of paths.
# Useful when creating a concatenated string with lots of filenames, where the string could be shortened using this approach.
def strip_topdir(list_of_paths, topdir):
    stripped = [x.replace(topdir, '') if x.startswith(topdir) else x for x in list_of_paths]
    stripped = [re.sub(r'^/+', '', x) for x in stripped]
    return(stripped)

############################
# Merge msamtools profiles
############################

# Get a list of msamtools profiles based on wildcards
def get_msamtools_profiles(wildcards):

    # Ensure right MINTO_MODE and filename combinations
    if     (wildcards.minto_mode in ['MAG', 'refgenome'] and wildcards.filename == 'gene_abundances') \
        or (wildcards.minto_mode == 'catalog'            and wildcards.filename in ['genome_abundances', 'contig_abundances']):
        raise Exception("MIntO mode {} cannot generate {} profiles".format(wildcards.minto_mode, wildcards.filename))

    # Ensure right MINTO_MODE and norm combination
    if (wildcards.minto_mode == 'catalog' and wildcards.type != 'TPM'):
        raise Exception("MIntO mode {} cannot generate {} profiles".format(wildcards.minto_mode, wildcards.type))

    # For genome_abundances, look for file with '<type>.genome' in it
    typespec = wildcards.type
    if (wildcards.filename == 'genome_abundances'):
        typespec = wildcards.type + '.genome'

    profiles = expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.{typespec}.txt.gz",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            minto_mode = wildcards.minto_mode,
                            identity = wildcards.identity,
                            typespec = typespec,
                            sample = ilmn_samples)
    return(profiles)

# Merge individual msamtools profiles from mapping to gene or genome DB.
# Input generated by the function above, which ensures compatibility.
# We set '--zeroes' because msamtools output leaves zero entries in individual profile files.

rule merge_msamtools_profiles:
    input:
        single=get_msamtools_profiles
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{filename}.p{identity}.{type}.tsv"
    params:
        topdir = lambda wildcards, input: os.path.commonpath(input.single),
        files  = lambda wildcards, input: ','.join(strip_topdir(input.single, os.path.commonpath(input.single)))
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/merge_msamtools_profiles.{filename}.p{identity}.{type}.log"
    resources:
        mem=30
    threads: lambda wildcards,input: min(10, len(input.single))
    wildcard_constraints:
        filename = r'gene_abundances|genome_abundances|contig_abundances',
        identity = r'\d+',
    conda:
        minto_dir + "/envs/r_pkgs.yml"
    shell:
        """
        time (
            shadowdir=$(pwd)
            cd {params.topdir}
            Rscript {script_dir}/merge_profiles.R --threads {threads} --memory {resources.mem} --input {params.files} --out $shadowdir/out.txt --keys ID --zeroes
            cd $shadowdir
            rsync -a out.txt {output.combined}
        ) >& {log}
        """

rule add_annotation_to_genome_profiles:
    input:
        profile  = "{wd}/{omics}/9-mapping-profiles/{minto_mode}/genome_abundances.p{identity}.relabund.prop.tsv",
        locusmap = "{wd}/DB/{minto_mode}/1-prokka/locus_id_list.txt",
        taxonomy = "{wd}/DB/{minto_mode}/3-taxonomy/taxonomy.{taxonomy}.{db_version}.tsv"
    output:
        mag_profile = "{wd}/output/9-mapping-profiles/{minto_mode}/{omics}.{taxonomy}.{db_version}.p{identity}.tsv",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{omics}.{taxonomy}.{db_version}.p{identity}.log"
    wildcard_constraints:
        taxonomy = r'gtdb|phylophlan',
        identity = r'\d+',
    resources:
        mem=10
    conda:
        minto_dir + "/envs/r_pkgs.yml"
    shell:
        """
        time (
            R --vanilla --silent --no-echo >> {output} <<___EOF___

library(data.table)

# Read tables
profile  = fread('{input.profile}',  header = TRUE, sep = "\\t")
locusmap = fread('{input.locusmap}', header = TRUE, sep = "\\t")
taxonomy = fread('{input.taxonomy}', header = TRUE, sep = "\\t", fill = TRUE)

# Link locus_tag and mag_id in the profile
# Keep Unmapped via left-join
dt = (
        merge(profile, locusmap, by.x='ID', by.y='locus_id', all.x=TRUE)
        [is.na(mag_id), mag_id := 'Unmapped']
     )

# Add taxonomy using mag_id
# Keep Unmapped via left-join
# MAGs without annotation will show up with empty taxonomy fields
# Annotate Unmapped as 'Unknown'
# Annotate missing taxonomy fields also as 'Unknown'
dt = (
        merge(dt, taxonomy, by='mag_id', all.x=TRUE)
        [, lapply(.SD, function(x) ifelse(x == '' | is.na(x), 'Unknown', x))]
        [, ID := NULL]
     )
setcolorder(dt, 'mag_id')

# Write mag_profile
fwrite(dt, file = '{output.mag_profile}', row.names = F, col.names = T, sep = "\\t", quote = F)

___EOF___
        ) >& {log}
        """

############################
# Merge bedcov profiles
############################

# Merge individual bedcov BED files from genome mapping
# We don't set '--zeroes' because rule 'genome_mapping_profiling' removed all zero entries in individual BED files.
# Memory is estimated based on a regression using 9Gb metagenomes.

rule merge_gene_abund:
    input:
        single=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.gene_abundances.p{identity}.bed.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    minto_mode = wildcards.minto_mode,
                    identity = wildcards.identity,
                    sample=ilmn_samples),
    params:
        topdir = lambda wildcards, input: os.path.commonpath(input.single),
        files  = lambda wildcards, input: ','.join(strip_topdir(input.single, os.path.commonpath(input.single)))
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.bed"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.log"
    threads: lambda wildcards,input: min(10, len(input.single))
    resources:
        mem=lambda wildcards,input: 5 + 0.3*len(input.single)
    conda:
        minto_dir + "/envs/r_pkgs.yml"
    shell:
        """
        time (
            shadowdir=$(pwd)
            cd {params.topdir}
            Rscript {script_dir}/merge_profiles.R --threads {threads} --memory {resources.mem} --input {params.files} --out $shadowdir/out.txt --keys gene_length,ID
            cd $shadowdir
            rsync -a out.txt {output.combined}
        ) >& {log}
        """

###############################################################################################
# Normalize profiles
###############################################################################################

############################
# Normalization of read counts by sequence depth and genes’ length (TPM normalization)
# or marker genes (MG normalization)
# Normalize and add prefix to samples so that metaG and metaT do not clash when combined in future
############################

# Input depends on norm mode, since MG mode needs location of MG files.
# Thus, input will be encapsulated in this function.

def get_gene_abund_normalization_input(wildcards):
    if (wildcards.norm == 'MG'):
        return {
            'absolute_counts' : "{wd}/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.bed".format(
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    minto_mode = wildcards.minto_mode,
                    identity = wildcards.identity),
            'genomes_marker_genes' : "{wd}/DB/{minto_mode}/genomes.marker_genes.table".format(
                    wd = wildcards.wd,
                    minto_mode = wildcards.minto_mode)
            }
    else:
        return {
            'absolute_counts' : "{wd}/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.bed".format(
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    minto_mode = wildcards.minto_mode,
                    identity = wildcards.identity)
            }

# Memory estimates:
# MG:  '1.157e+07 + 1.727e-02*bed - 2.122e+00*mg - 8.937e+04*samples' in memKB
# TPM: '1.990e+06 + 1.541e-02*bed                - 1.783e+04*samples' in memKB
# Simplified to 'ceil(1.73e-08*bed) + MG?12:2' in memGB
# Ignore the reduction of memory per sample
# And add 10GB with each new attempt
rule gene_abund_normalization:
    input:
        unpack(get_gene_abund_normalization_input)
    output:
        norm_counts="{wd}/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.{norm}.tsv"
    params:
        mapped_reads_threshold=MIN_mapped_reads,
        optional_arg_MG = lambda wildcards, input: "" if wildcards.norm == "TPM" else "--MG " + input.genomes_marker_genes
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.{norm}.log"
    wildcard_constraints:
        minto_mode = r'MAG|refgenome'
    threads: 4
    resources:
        mem = lambda wildcards, input, attempt: (12 if wildcards.norm == 'MG' else 2) + math.ceil(1.73e-8*get_file_size(input.absolute_counts)) + 10*(attempt-1)
    conda:
        minto_dir + "/envs/r_pkgs.yml" #R
    shell:
        """
        time (
            Rscript {script_dir}/normalize_profiles.R --normalize {wildcards.norm} --threads {threads} --memory {resources.mem} --bed {input.absolute_counts} {params.optional_arg_MG} --out {output.norm_counts} --min-read-count {params.mapped_reads_threshold}
        ) &> {log}
        """

###############################################################################################
# Summarize mappability rates
###############################################################################################

############################
## Mappability rate
## Mappability stats
## Multimapping read count
############################

rule read_map_stats:
    input:
        map_profile=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            minto_mode = wildcards.minto_mode,
                                            identity = wildcards.identity,
                                            sample=ilmn_samples)
    output:
        maprate= "{wd}/{omics}/9-mapping-profiles/{minto_mode}/mapping.p{identity}.maprate.txt",
        mapstats="{wd}/{omics}/9-mapping-profiles/{minto_mode}/mapping.p{identity}.mapstats.txt",
        multimap="{wd}/{omics}/9-mapping-profiles/{minto_mode}/mapping.p{identity}.multimap.txt"
    params:
        map_profile_list=lambda wildcards, input: ",".join(input.map_profile),
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/mapping.p{identity}.read_map_stats.log"
    wildcard_constraints:
        identity = r'\d+'
    threads: 1
    resources:
        mem=2
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        # Init empty files
        true > {output.maprate}; true > {output.mapstats}; true > {output.multimap}

        # Generate content
        time (
            for file in {input.map_profile}; do
                sample=$(basename $file); sample=${{sample%%.p{wildcards.identity}*}}
                (echo -e -n "$sample\\t"; bash -c "zcat $file | head | grep 'Mapped inserts'  | cut -f2 -d'(' | sed 's/%.*//'" ) >> {output.maprate}
                (echo -e -n "$sample\\t"; bash -c "zcat $file | head | grep 'Mapped inserts'  | cut -f2 -d':' | sed 's/^ //'"  ) >> {output.mapstats}
                (echo -e -n "$sample\\t"; bash -c "zcat $file | head | grep 'Multiple mapped' | cut -f2 -d'(' | sed 's/%.*//'" ) >> {output.multimap}
            done
        ) &> {log}
        """

###############################################################################################
# Generate configuration yml file for data integration - gene and function profiles
###############################################################################################

rule config_yml_integration:
    input: f"{working_dir}/{omics}/9-mapping-profiles/{MINTO_MODE}/genome_abundances.p{identity}.relabund.prop.tsv"
    output:
        config_file="{wd}/data_integration.yaml"
    params:
        mapped_reads_threshold=MIN_mapped_reads,
    resources:
        mem=2
    threads: 1
    log:
        "{wd}/logs/config_yml_integration.log"
    shell:
        """
        time (echo "######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: metaG_metaT
minto_dir: {minto_dir}
METADATA: {metadata}

######################
# Analysis settings
######################

MAIN_factor: {main_factor}
PLOT_factor2: {plot_factor2}
PLOT_time: {plot_time}

######################
# Program settings
######################

alignment_identity: {identity}
abundance_normalization: MG
MINTO_MODE: {MINTO_MODE}

MERGE_threads: 4

ANNOTATION_file:

# List annotation IDs to generate function profiles.
#
# If MINTO_MODE is 'MAG' or 'refgenome', this list could contain elements from:
# 'eggNOG.OGs', 'eggNOG.KEGG_Pathway', 'eggNOG.KEGG_Module', 'eggNOG.KEGG_KO', 'eggNOG.PFAMs',
# 'kofam.KEGG_Pathway', 'kofam.KEGG_Module', 'kofam.KO',
# 'dbCAN.module', 'dbCAN.enzclass', 'dbCAN.subfamily', 'dbCAN.EC', 'dbCAN.eCAMI_subfamily', 'dbCAN.eCAMI_submodule'.
#
#   E.g.:
# - eggNOG.OGs
# - kofam.KEGG_Pathway
#
# If MINTO_MODE is 'catalog', the names should match the ANNOTATION_file column names.

ANNOTATION_ids:
 - eggNOG.OGs
 - eggNOG.PFAMs
 - dbCAN.module
 - dbCAN.enzclass
 - dbCAN.subfamily
 - dbCAN.EC
 - dbCAN.eCAMI_subfamily
 - dbCAN.eCAMI_submodule
 - kofam.KEGG_Pathway
 - kofam.KEGG_Module
 - kofam.KEGG_KO
 - eggNOG.KEGG_Pathway
 - eggNOG.KEGG_Module
 - eggNOG.KEGG_KO" > {output.config_file}
) >& {log}
        """
