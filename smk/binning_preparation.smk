#!/usr/bin/env python

'''
Binning preparation of contigs from assembly/coassembly.

Authors: Carmen Saenz, Mani Arumugam
'''

import os.path
import math

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/fasta_bam_helpers.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule binning_preparation_base, binning_preparation_vamb from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# Variables from configuration yaml file

##############################################
# Get sample list
##############################################

# Make list of illumina samples, if ILLUMINA in config

ilmn_samples = list()
if (x := validate_required_key(config, 'ILLUMINA')):
    check_input_directory(x, locations = ['7-assembly'])
    ilmn_samples = x

###############################
# MetaSPAdes parameters
###############################

METASPADES_illumina_max_k = validate_required_key(config, 'METASPADES_illumina_max_k')
check_number_is_odd('METASPADES_illumina_max_k', METASPADES_illumina_max_k)

METASPADES_hybrid_max_k = validate_optional_key(config, 'METASPADES_hybrid_max_k')
if METASPADES_hybrid_max_k is not None:
    check_number_is_odd('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k)

if (x := validate_optional_key(config, 'METASPADES_custom_build')):
    check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 300)
    if METASPADES_hybrid_max_k is not None:
        check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 300)
else:
    check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 128)
    if METASPADES_hybrid_max_k is not None:
        check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 128)

###############################
# MEGAHIT parameter sets
###############################

MEGAHIT_presets = list()
mega_k_list = list()

# Check for MEGAHIT_presets
if (x := validate_optional_key(config, 'MEGAHIT_presets')):
    MEGAHIT_presets = x

# Also check for MEGAHIT_custom
if (x := validate_optional_key(config, 'MEGAHIT_custom')):
    # if the custom k-s are set, that should be added to the MEGAHIT assembly types
    if isinstance(x, str):
        x = [x]
    for i, k_list in enumerate(x):
        if k_list and not k_list.isspace():
            MEGAHIT_presets.append(f'meta-custom-{i+1}')
            mega_k_list.append(k_list)

###############################
# MetaFlye parameter sets
###############################

METAFLYE_presets = dict()

# Check for METAFLYE_presets
if (x := validate_optional_key(config, 'METAFLYE_presets')):
    METAFLYE_presets = x

max_job_ram_gb = validate_required_key(config, 'MAX_RAM_GB_PER_JOB')
if max_job_ram_gb < 10:
    print('WARNING in ', config_path, ': MAX_RAM_GB_PER_JOB variable has to be minimum 10GB. Value adjusted to 10GB.')
    max_job_ram_gb = 10

###############################
# MAG building parameters
###############################

MIN_FASTA_LENGTH = validate_required_key(config, 'MIN_FASTA_LENGTH')

spades_contigs_or_scaffolds = "scaffolds"
if (x := validate_optional_key(config, 'SPADES_CONTIGS_OR_SCAFFOLDS')):
    check_allowed_values('SPADES_CONTIGS_OR_SCAFFOLDS', x, ('contigs', 'contig', 'scaffolds', 'scaffold'))
    if not x.endswith('s'):
        x += 's'
    spades_contigs_or_scaffolds = x

BWA_threads = validate_required_key(config, 'BWA_threads')
SAMTOOLS_sort_perthread_memgb = validate_required_key(config, 'SAMTOOLS_sort_perthread_memgb')

###############################
# Make a list of assemblies to use
###############################

# Scaffold type
SCAFFOLDS_type = list()

if ilmn_samples:
    print(f"Found ILLUMINA in {config_path}.")
    SCAFFOLDS_type.append('illumina_single')

# Make list of nanopore samples, if NANOPORE in config

nanopore_assemblies = list()
if (x := validate_optional_key(config, 'NANOPORE')):
    print(f"Found NANOPORE in {config_path}.")
    check_input_directory(x, locations = ['7-assembly'])
    nanopore_assemblies = x
    SCAFFOLDS_type.append('nanopore')

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config

hybrid_assemblies = list()
if (x := validate_optional_key(config, 'HYBRID')):
    print(f"Found HYBRID in {config_path}.")
    SCAFFOLDS_type.append('illumina_single_nanopore')
    for nano in x:
        for ilmn in x[nano].split("+"):
            hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if COASSEMBLY in config
co_assemblies = list()
if validate_optional_key(config, 'enable_COASSEMBLY'):
    if (x := validate_optional_key(config, 'COASSEMBLY')):
        print(f"Found COASSEMBLY in {config_path}.")
        SCAFFOLDS_type.append('illumina_coas')
        co_assemblies = x

# Initialize some variables

possible_assembly_types = ['illumina_single_nanopore', 'nanopore', 'illumina_coas', 'illumina_single']

# List of assemblies per type
assemblies = {}
assemblies['illumina_single_nanopore'] = hybrid_assemblies
assemblies['nanopore']                 = nanopore_assemblies
assemblies['illumina_coas']            = co_assemblies
assemblies['illumina_single']          = ilmn_samples

# Remove assembly types if specified

if (x := validate_optional_key(config, 'EXCLUDE_ASSEMBLY_TYPES')):
    for item in x:
        if item in SCAFFOLDS_type:
            print(f"Removing assembly_type={item} on request")
            SCAFFOLDS_type.remove(item)

# Site customization for avoiding NFS traffic during I/O heavy steps such as mapping

CLUSTER_NODES            = None
CLUSTER_LOCAL_DIR        = None
CLUSTER_WORKLOAD_MANAGER = None
include: minto_dir + '/site/cluster_def.py'

# Cluster-aware bwa-index rules

include: 'include/bwa_index_wrapper.smk'

# If BWA index clean-up requested, only do cleanup and nothing else
if CLUSTER_NODES != None and validate_optional_key(config, 'CLEAN_BWA_INDEX'):

    print("NOTE: BWA index cleanup mode requested.")

    def clean_bwa_index():
        # Get the cleanup flag for index files per scaf_type
        files = expand("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/BWA_index.batches.cleaning.done",
                            wd         = working_dir,
                            omics      = omics,
                            scaf_type  = SCAFFOLDS_type,
                            min_length = MIN_FASTA_LENGTH)
        return(files)

    rule all:
        input:
            clean_bwa_index()
        default_target: True

else:

# Define all the outputs needed by target 'all'

    rule all:
        input:
            abundance   = f"{working_dir}/{omics}/8-1-binning/scaffolds.{MIN_FASTA_LENGTH}.abundance.npz",
            config_yaml = f"{working_dir}/{omics}/mags_generation.yaml",
            versions    = print_versions.get_version_output(snakefile_name)
        default_target: True

###############################################################################################
# Filter contigs from
#  1. MetaSPAdes-individual-assembled illumina samples,
#  2. MEGAHIT-co-assembled illumina samples,
#  3. MetaSPAdes-hybrid-assembled illumina+nanopore data,
#  4. MetaFlye-assembled nanopore sequences
###############################################################################################

def get_assemblies_for_scaf_type(wildcards):
    if (wildcards.scaf_type == 'nanopore'):
        asms = expand("{wd}/{omics}/7-assembly/{assembly}/{assembly_preset}/{assembly}.assembly.fasta",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            assembly = assemblies[wildcards.scaf_type],
                            assembly_preset = METAFLYE_presets)
    elif (wildcards.scaf_type == 'illumina_single'):
        asms = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{contig_or_scaffold}.fasta",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            assembly = assemblies[wildcards.scaf_type],
                            kmer_dir = "k21-" + str(METASPADES_illumina_max_k),
                            contig_or_scaffold = spades_contigs_or_scaffolds)
    elif (wildcards.scaf_type == 'illumina_coas'):
        asms = expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            coassembly = assemblies[wildcards.scaf_type],
                            assembly_preset = MEGAHIT_presets)
    elif (wildcards.scaf_type == 'illumina_single_nanopore'):
        asms = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{contig_or_scaffold}.fasta",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            assembly = assemblies[wildcards.scaf_type],
                            kmer_dir = "k21-" + str(hybrid_max_k),
                            contig_or_scaffold = spades_contigs_or_scaffolds)
    else:
        raise Exception(f"MIntO error: scaf_type={wildcards.scaf_type} is not allowed!")

    return(asms)

def calculate_batch_filesize_mb(wildcards):
    filtered_batch_size = (max_job_ram_gb-5)/22*1024
    if wildcards.scaf_type.endswith("nanopore"):
        # assuming that size filtering doesn't drastically reduce file size
        return(round(filtered_batch_size))
    else:
        # 2500bp filtering reduces file size by approx. half,
        # we account for lower min. length and better assemblies here
        return(round((filtered_batch_size-50)/0.6))

# The following checkpoint will make a txt file with the names of assemblies per batch.
# The next rule will make the fasta file for each batch of assemblies.
# This separation is to enable deletion of the batch fasta files when
# a combined file is created, which will save space for large studies.
checkpoint make_assembly_batches:
    input:
        assemblies = get_assemblies_for_scaf_type
    output:
        location = directory("{wd}/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}")
    log:
        "{wd}/logs/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}.batching.log"
    wildcard_constraints:
        min_length = r'\d+',
        scaf_type  = r'illumina_single|illumina_coas|illumina_single_nanopore|nanopore'
    resources:
        mem = 5
    params:
        batch_filesize = calculate_batch_filesize_mb
    run:
        import datetime
        import os

        def logme(stream, msg):
            print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), msg, file=stream)

        # make the checkpoint directory
        os.mkdir(output.location)

        with open(str(log), 'w') as f:
            logme(f, "INFO: assemblies={} batch_size={}Mb".format(len(input.assemblies), params.batch_filesize))
            start = 0
            batch = 1
            current_batchfilesize = 0
            for i in range(0, len(input.assemblies)):
                current_assemblysize = os.path.getsize(input.assemblies[i]) / (1024**2)

                # if adding the current assembly made batch larger than accepted, then write out batch
                if (i and ((current_batchfilesize + current_assemblysize) > params.batch_filesize)):
                    logme(f, "INFO: Making BATCH={} size={:.0f}Mb".format(batch, current_batchfilesize))

                    end = i

                    # batch goes from start to end-1
                    batch_input = input.assemblies[start:end]
                    batch_list  = output.location + f"/batch{batch}.list"
                    logme(f, "INFO:   INPUT ={}".format(batch_input))
                    logme(f, "INFO:   OUTPUT={}".format(batch_list))

                    # write out batch list
                    with open(batch_list, 'w') as listfile:
                        for asm in batch_input:
                            print(asm, file=listfile)

                    # reset variables
                    start = i
                    batch += 1
                    current_batchfilesize = 0
                current_batchfilesize += current_assemblysize

            # last assembly and batch
            batch_input = input.assemblies[start:]
            batch_list  = output.location + f"/batch{batch}.list"
            logme(f, "INFO:   INPUT ={}".format(batch_input))
            logme(f, "INFO:   OUTPUT={}".format(batch_list))
            with open(batch_list, 'w') as listfile:
                for asm in batch_input:
                    print(asm, file=listfile)
            logme(f, "INFO: done")

rule write_assembly_batch_fasta:
    input:
        list  = "{wd}/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}/batch{batch}.list"
    output:
        fasta = temp("{wd}/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}/batch{batch}.fasta.gz")
    log:
        "{wd}/logs/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}/batch{batch}.writing.log"
    wildcard_constraints:
        min_length = r'\d+',
        scaf_type  = r'illumina_single|illumina_coas|illumina_single_nanopore|nanopore'
    resources:
        mem = 5
    params:
        min_length = lambda wildcards: wildcards.min_length,
    run:
        import datetime

        def logme(stream, msg):
            print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), msg, file=stream)

        with open(str(log), 'w') as f:
            with open(input.list) as listfile:
                # get assembly list
                assemblies = listfile.read().splitlines()

                # write out batch
                filter_fasta_list_by_length(assemblies, output.fasta, params.min_length)

                # log progress
                logme(f, "INFO:   INPUT ={}".format(assemblies))
                logme(f, "INFO:   OUTPUT={}".format(output.fasta))

            logme(f, "INFO: done")

################################################################################################
# Create BWA index for the filtered contigs
################################################################################################

# If CLUSTER_NODES is defined, then return the file symlink'ed to local drives.
# Otherwise, return the original index files in project work_dir.
def get_contig_bwa_index(wildcards):

    # Where are the index files?
    if CLUSTER_NODES != None:
        index_location = 'BWA_index_local'
    else:
        index_location = 'BWA_index'

    # Get all the index files!
    files = expand("{wd}/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}/{location}/batch{batch}.fasta.gz.{ext}",
            wd         = wildcards.wd,
            omics      = wildcards.omics,
            scaf_type  = wildcards.scaf_type,
            location   = index_location,
            batch      = wildcards.batch,
            min_length = wildcards.min_length,
            ext        = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'])

    # Return them
    return(files)

# Maps reads to contig-set using bwa2
# Baseline memory usage is:
#   BWA : 3.1 bytes per base in DNA database (regression: mem = 5.556e+09 + 3.011*input)
#   Sort: using config values: SAMTOOLS_sort_threads and SAMTOOLS_sort_perthread_memgb;
#   Total = BWA + Sort
#   E.g. 25GB for BWA, 2 sort threads and 10GB per-thread = 45GB
# If an attempt fails, 30GB extra per sort thread + 30GB extra for BWA.
# Threads for individual tasks: Work backwards from {threads} from snakemake
# Once mapping succeeds, run coverM to get contig depth for this bam

rule map_contigs_BWA_depth_coverM:
    input:
        bwaindex=get_contig_bwa_index,
        fwd='{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz',
        rev='{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz'
    output:
        depth = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batch{batch}/{illumina}.depth.txt.gz")
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/6-mapping/{illumina}/{illumina}.scaffolds_{scaf_type}.{min_length}.batch{batch}.bwa2.log"
    wildcard_constraints:
        batch    = r'\d+',
        illumina = r'[^/]+'
    resources:
        samtools_sort_threads = 3,
        map_threads = lambda wildcards, threads: max(1, threads - 3),
        coverm_threads = lambda wildcards, threads: min(8, threads),
        mem = lambda wildcards, input, attempt: int(10 + 3.1*os.path.getsize(input.bwaindex[0])/1e9 + 1.1*3*(SAMTOOLS_sort_perthread_memgb + 30*(attempt-1))),
        sort_mem = lambda wildcards, attempt: SAMTOOLS_sort_perthread_memgb + 30*(attempt-1)
    threads:
        BWA_threads + 3
    conda:
        minto_dir + "/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        mkdir -p $(dirname {output})
        db_name=$(echo {input.bwaindex[0]} | sed "s/.0123//")
        time (bwa-mem2 mem -P -a -t {resources.map_threads} $db_name {input.fwd} {input.rev} \
                | msamtools filter -buS -p 95 -l 45 - \
                | samtools sort -m {resources.sort_mem}G --threads {resources.samtools_sort_threads} - \
                > {wildcards.illumina}.bam
              coverm contig --methods metabat --trim-min 10 --trim-max 90 --min-read-percent-identity 95 --threads {resources.coverm_threads} --output-file sorted.depth --bam-files {wildcards.illumina}.bam
              gzip -2 sorted.depth
              rsync -a sorted.depth.gz {output.depth}
        ) >& {log}
        """

# Combine coverM depths from individual samples
# Ignore columns ending with '-var'
# Also ignore contigLen, totalAvgDepth
rule colbind_sample_contig_depths_for_batch:
    input:
        depths = lambda wildcards: expand("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batch{batch}/{illumina}.depth.txt.gz",
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            scaf_type = wildcards.scaf_type,
                                            min_length = wildcards.min_length,
                                            batch = wildcards.batch,
                                            illumina=ilmn_samples)
    output:
        depths = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batch{batch}.depth.txt.gz"),
        header = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batch{batch}.header.txt.gz"),
    log:
        "{wd}/logs/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batch{batch}.depth.log"
    wildcard_constraints:
        batch = r'\d+',
    shadow:
        "minimal"
    resources:
        mem = lambda wildcards, input: 5+len(input.depths)*0.1
    threads:
        2
    run:
        import pandas as pd
        import hashlib
        import pickle
        import datetime
        import shutil

        def logme(stream, msg):
            print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), msg, file=stream)

        with open(str(log), 'w') as f:
            df_list = list()
            logme(f, "INFO: reading file 0")

            # Read first file into df
            df = pd.read_csv(input.depths[0], header=0, sep = "\t", memory_map=True)

            # Get md5 of first 2 columns
            md5_first = hashlib.md5(pickle.dumps(df.iloc[:, 0:2])).hexdigest() # make hash for seqname and seqlen

            # Drop columns ending in '-var' as this is not used by vamb, and add to df_list
            df_list.append(df.drop(df.filter(regex='^contigLen$|^totalAvgDepth$|-var$').columns, axis='columns'))

            # Process remaining files
            for i in range(1, len(input.depths)):
                logme(f, "INFO: reading file {}".format(i))

                # Read next file
                df = pd.read_csv(input.depths[i], header=0, sep = "\t", memory_map=True)

                # Make sure md5 matches
                md5_next = hashlib.md5(pickle.dumps(df.iloc[:, 0:2])).hexdigest()
                if md5_next != md5_first:
                    raise Exception("colbind_sample_contig_depths_for_batch: Sequences don't match between {} and {}".format(input.depths[0], input.depths[i]))

                # Drop common columns and the '-var' column
                df_list.append(df.drop(df.filter(regex='^contigName$|^contigLen$|^totalAvgDepth$|-var$').columns, axis='columns'))

            # Concat df_list into single df
            logme(f, "INFO: concatenating {} files".format(len(input.depths)))
            df = pd.concat(df_list, axis=1, ignore_index=False, copy=False, sort=False)

            # Write header
            logme(f, "INFO: writing header file")
            df.head(0).to_csv(output.header, sep = "\t", header = True, index = False, compression={'method': 'gzip', 'compresslevel': 2})

            # Write output
            logme(f, "INFO: writing depths to temporary file")
            df.to_csv('depths.gz', sep = "\t", header = False, index = False, compression={'method': 'gzip', 'compresslevel': 2})
            logme(f, "INFO: copying depth file to final output")
            shutil.copy2('depths.gz', output.depths)
            logme(f, "INFO: done")

##################################################
# Combining across batches within a scaffold_type
##################################################

def get_batches_for_scaf_type_generic(wildcards, filetype):
    chkpnt_output = checkpoints.make_assembly_batches.get(**wildcards).output[0]
    batches       = glob_wildcards(os.path.join(chkpnt_output, "batch{batch,\d+}.list")).batch
    prefix        = 'scaffolds' if (filetype == 'fasta') else 'depth'
    result        = expand("{wd}/{omics}/8-1-binning/{prefix}_{scaf_type}.{min_length}/batch{batch}.{extension}.gz",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            prefix = prefix,
                            scaf_type = wildcards.scaf_type,
                            min_length = wildcards.min_length,
                            extension = filetype,
                            batch = batches)
    return(result)

def get_fasta_batches_for_scaf_type(wildcards):
    return(get_batches_for_scaf_type_generic(wildcards, filetype='fasta'))

def get_depth_batches_for_scaf_type(wildcards):
    return(get_batches_for_scaf_type_generic(wildcards, filetype='depth.txt'))

def get_header_batches_for_scaf_type(wildcards):
    return(get_batches_for_scaf_type_generic(wildcards, filetype='header.txt'))

# Make one header for this scaf_type
# Sanity check to ensure that the order of sample-depths in the per-batch depth files is the same. Otherwise binning will be wrong!
rule combine_contig_depth_header_batches:
    localrule: True
    input:
        header = get_header_batches_for_scaf_type
    output:
        header = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/combined.header.txt.gz")
    log:
        "{wd}/logs/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/batches.check_depth.log"
    shell:
        """
        rm --force {output}
        uniq_headers=$(for file in {input.header}; do zcat $file; done | sort -u | wc -l)
        if [ "$uniq_headers" == "1" ]; then
            rsync {input.header[0]} {output.header}
        else
            # Headers were not unique. So report and set exit code for this rule to non-zero
            echo "Headers for depth files were not unique - please check depth and header files:"
            for file in {input.header}; do zcat $file; done | sort -u
            >&2 echo "Headers for depth files were not unique - please check depth and header files"
            test "1" == "2"
        fi
        """

rule combine_contig_depth_batches:
    localrule: True
    input:
        depths = get_depth_batches_for_scaf_type,
        header = rules.combine_contig_depth_header_batches.output.header
    output:
        depths = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/combined.depth.txt.gz")
    shadow:
        "minimal"
    params:
        num_files = lambda wildcards, input: len(input.depths)
    resources:
       mem = 10
    threads:
        1
    shell:
        """
        if (( {params.num_files} > 1 )); then
            cat {input.depths} > combined.depth.gz
            rsync -a combined.depth.gz {output.depths}
        else
            ln --force {input.depths[0]} {output.depths}
        fi
        """

# Combine multiple fasta files from batches into one per scaf_type
rule combine_fasta_batches:
    localrule: True
    input:
        fasta = get_fasta_batches_for_scaf_type
    output:
        fasta_combined = temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/combined.fasta.gz")
    wildcard_constraints:
        min_length = r'\d+',
        scaf_type  = r'illumina_single|illumina_coas|illumina_single_nanopore|nanopore'
    shadow:
        "minimal"
    params:
        num_files = lambda wildcards, input: len(input.fasta)
    resources:
       mem = 10
    threads:
        1
    shell:
        """
        if (( {params.num_files} > 1 )); then
            cat {input.fasta} > combined.fasta.gz
            rsync -a combined.fasta.gz {output}
        else
            ln --force {input.fasta[0]} {output}
        fi
        """

#####################################
# Combining across scaffold_type
#####################################

def get_files_across_scaffold_types_generic(wildcards, filetype):
    files = expand("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/combined.{extension}.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    min_length = wildcards.min_length,
                    extension = filetype,
                    scaf_type = SCAFFOLDS_type)
    return(files)

def get_fasta_files_across_scaffold_types(wildcards):
    return(get_files_across_scaffold_types_generic(wildcards, 'fasta'))

def get_depth_files_across_scaffold_types(wildcards):
    return(get_files_across_scaffold_types_generic(wildcards, 'depth.txt'))

def get_header_files_across_scaffold_types(wildcards):
    return(get_files_across_scaffold_types_generic(wildcards, 'header.txt'))

# Combine multiple fasta files into one
rule combine_fasta:
    localrule: True
    input:
        fasta = get_fasta_files_across_scaffold_types
    output:
        fasta_combined="{wd}/{omics}/8-1-binning/scaffolds.{min_length}.fasta.gz"
    shadow:
        "minimal"
    params:
        num_files = lambda wildcards, input: len(input.fasta)
    resources:
       mem = 10
    threads:
        1
    shell:
        """
        if (( {params.num_files} > 1 )); then
            cat {input} > combined.fasta.gz
            rsync -a combined.fasta.gz {output}
        else
            ln --force {input[0]} {output}
        fi
        """

# Sanity check to ensure that the order of sample-depths in the depth file is the same. Otherwise binning will be wrong!
rule combine_contig_depth_header:
    localrule: True
    input:
        header = get_header_files_across_scaffold_types
    output:
        header = temp("{wd}/{omics}/8-1-binning/scaffolds.{min_length}.header.txt.gz")
    shell:
        """
        rm --force {output}
        uniq_headers=$(for file in {input.header}; do zcat $file; done | sort -u | wc -l)
        if [ "$uniq_headers" == "1" ]; then
            rsync {input.header[0]} {output.header}
        else
            # Headers were not unique. So report and set exit code for this rule to non-zero
            for file in {input.header}; do zcat $file; done | sort -u
            >&2 echo "Headers in depth files were not unique - please check depth files"
            test "1" == "2"
        fi
        """

### Prepare abundance.npz for vamb v5+
# Memory requirement:
# From regression of 'gzip -2' depth file sizes, number of samples and maxmem from GNU time:
#    memKB = 4.336e+5 + 1.745e-3*filesize - 1.150e+3*samples
#    memGB = 0.43 + 1.745e-9*filesize - 1.2e-3*samples
# Ignored the negative effect of samples to err on cautious side.
# Safely converted to 2 + int(1.75e-9*filesize)
# Using dummy filesize value to allow --dry-run to work when files don't exist yet
rule make_abundance_npz:
    input:
        contigs = rules.combine_fasta.output.fasta_combined,
        header  = rules.combine_contig_depth_header.output.header,
        depths  = get_depth_files_across_scaffold_types,
    output:
        npz="{wd}/{omics}/8-1-binning/scaffolds.{min_length}.abundance.npz"
    log:
        "{wd}/logs/{omics}/8-1-binning/scaffolds.{min_length}.abundance.log"
    shadow:
        "minimal"
    threads:
        1
    resources:
        mem = lambda wildcards, input: 2 + int(1.75e-9*sum((lambda file: os.path.getsize(file) if os.path.exists(file) else 1048576)(file) for file in input.depths))
    conda:
        minto_dir + "/envs/vamb.yml"
    shell:
        """
        time (
            (zcat {input.header} | sed 's/^contigName/contigname/'; zcat {input.depths}) > combined.depth
            python3 {script_dir}/make_vamb_abundance_npz.py --fasta {input.contigs} --abundance-tsv combined.depth --minlength {wildcards.min_length} --output abundance.npz
            rsync -a abundance.npz {output.npz}
        ) >& {log}
        """

##################################################
# Clean up BWA index files in the mirror locations
##################################################

# Get list of clean flags for cleaning up all batches for this scaf_type
def get_flags_to_clean_bwaindex_mirror_for_scaf_type(wildcards):
    chkpnt_output = checkpoints.make_assembly_batches.get(**wildcards).output[0]
    batches       = glob_wildcards(os.path.join(chkpnt_output, "batch{batch,\d+}.list")).batch
    result        = expand("{wd}/{omics}/8-1-binning/scaffolds_{scaf_type}.{min_length}/BWA_index/batch{batch}.fasta.gz.clustersync/cleaning.done",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            scaf_type = wildcards.scaf_type,
                            min_length = wildcards.min_length,
                            batch = batches)
    return(result)

# Clean BWA index for all batches for a given scaf_type
rule clean_bwaindex_mirror_for_scaf_type:
    localrule: True
    input:
        fasta = get_flags_to_clean_bwaindex_mirror_for_scaf_type
    output:
        temp("{wd}/{omics}/8-1-binning/depth_{scaf_type}.{min_length}/BWA_index.batches.cleaning.done")
    wildcard_constraints:
        min_length = r'\d+',
        scaf_type  = r'illumina_single|illumina_coas|illumina_single_nanopore|nanopore'
    shadow:
        "minimal"
    threads:
        1
    shell:
        """
        touch {output}
        """

###############################################################################################
# Generate configuration yml file for recovery of MAGs and taxonomic annotation step - binning
###############################################################################################

rule config_yml_binning:
    localrule: True
    output:
        config_file="{wd}/{omics}/mags_generation.yaml",
    resources:
        mem=2
    threads: 2
    log:
        "{wd}/logs/{omics}/config_yml_mags_generation.log"
    params:
        min_fasta_length=MIN_FASTA_LENGTH
    shell:
        """
        time (echo "######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
minto_dir: {minto_dir}
METADATA: {metadata}

######################
# Program settings
######################
# COMMON PARAMETERS
#
MIN_FASTA_LENGTH: {params.min_fasta_length}
MIN_MAG_LENGTH: 500000

# VAMB settings
#
BINNERS:
- aaey
- aaez
- vae384

VAMB_THREADS: 20

# Use GPU in VAMB:
# could be yes or no
VAMB_GPU: no

# Metabuli for taxvamb
#
TAXVAMB_ANNOTATOR: metabuli
METABULI_MEM_GB: 128

# CHECKM settings
#
CHECKM_COMPLETENESS: 90  # higher than this
CHECKM_CONTAMINATION: 5  # lower than this
CHECKM_BATCH_SIZE: 50    # Process MAGs with this batch size

# COVERM settings
#
COVERM_THREADS: 8
COVERM_memory: 5

# SCORING THE BEST GENOMES settings
#
# this could be checkm or genome
SCORE_METHOD: checkm
" > {output.config_file}

) >& {log}
        """
