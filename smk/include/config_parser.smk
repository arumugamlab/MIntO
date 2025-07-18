#!/usr/bin/env python

'''
Config file parser

Authors: Carmen Saenz, Mani Arumugam
'''

import os.path

##############################################
# Config key dictionary
##############################################

GLOBAL_CONFIG_KEYTYPES = {
    'abundance_normalization' : str,
    'alignment_identity' : int,
    'ANNOTATION_file' : 'file',
    'ANNOTATION_ids' : list,
    'ANNOTATION' : list,
    'BINNERS' : list,
    'BWA_host_threads' : int,
    'BWA_threads' : int,
    'CHECKM_BATCH_SIZE' : int,
    'CHECKM_COMPLETENESS' : int,
    'CHECKM_CONTAMINATION' : int,
    'CLEAN_BWA_INDEX' : bool,
    'COAS_factor' : str,
    'COASSEMBLY' : dict,
    'CONTIG_MAPPING_BATCH_SIZE' : int,
    'COVERM_memory' : int,
    'COVERM_THREADS' : int,
    'download_memory' : int,
    'download_threads' : int,
    'eggNOG_dbmem' : bool,
    'enable_COASSEMBLY' : bool,
    'enable_GTDB' : bool,
    'enable_metaphlan' : bool,
    'enable_motus' : bool,
    'EXCLUDE_ASSEMBLY_TYPES' : list,
    'FASTP_adapters' : ['file', str],
    'FASTP_front_mean_qual' : int,
    'FASTP_memory' : int,
    'FASTP_min_length' : int,
    'FASTP_tail_mean_qual' : int,
    'FASTP_threads' : int,
    'GTDB_TAXONOMY_VERSION' : str,
    'HYBRID' : dict,
    'ILLUMINA' : [list, 'file', str],
    'ILLUMINA_suffix' : list,
    'MAG_omics' : str,
    'MAIN_factor' : str,
    'MAX_RAM_GB_PER_JOB' : int,
    'MEGAHIT_custom' : [str, list],
    'MEGAHIT_memory' : int,
    'MEGAHIT_presets' : list,
    'MEGAHIT_threads' : int,
    'MERGE_ILLUMINA_SAMPLES': dict,
    'MERGE_memory' : int,
    'MERGE_threads' : int,
    'METABULI_MEM_GB' : int,
    'METADATA' : 'file',
    'METAFLYE_presets' : dict,
    'metaphlan_version' : str,
    'METASPADES_custom_build' : 'file',
    'METASPADES_hybrid_max_k' : int,
    'METASPADES_illumina_max_k' : int,
    'METASPADES_memory' : int,
    'METASPADES_qoffset' : str,
    'METASPADES_threads' : int,
    'MIN_KEGG_PATHWAY_COMPLETENESS' : int,
    'MIN_FASTA_LENGTH' : int,
    'MIN_MAG_LENGTH' : int,
    'MIN_mapped_reads' : int,
    'minto_dir' : 'directory',
    'MINTO_MODE' : str,
    'motus_version' : str,
    'msamtools_filter_length' : int,
    'MULTIPLEX_TECH' : str,
    'NAME_host_genome' : str,
    'NAME_reference' : str,
    'NANOPORE' : list,
    'omics' : str,
    'PATH_host_genome' : 'directory',
    'PATH_reference' : 'directory',
    'perc_remaining_reads' : int,
    'PHYLOPHLAN_TAXONOMY_VERSION' : str,
    'PLOT_factor2' : str,
    'PLOT_time' : str,
    'PROJECT' : str,
    'raw_reads_dir' : 'directory',
    'READ_minlen' : int,
    'rRNA_index_memory' : int,
    'rRNA_index_threads' : int,
    'RUN_TAXONOMY' : bool,
    'SAMTOOLS_sort_perthread_memgb' : int,
    'SAMTOOLS_sort_threads' : int,
    'SCORE_METHOD' : str,
    'sortmeRNA_memory' : int,
    'sortmeRNA_threads' : int,
    'sortmeRNA_db' : 'directory',
    'sortmeRNA_db_idx' : 'directory',
    'SOURMASH_cutoff' : float,
    'SOURMASH_max_abund' : int,
    'SOURMASH_min_abund' : int,
    'SPADES_CONTIGS_OR_SCAFFOLDS' : str,
    'TAXA_memory' : int,
    'TAXA_profiler' : str,
    'TAXA_threads' : int,
    'TAXONOMY_CPUS' : int,
    'TAXONOMY_memory' : int,
    'TAXONOMY_NAME' : str,
    'TAXVAMB_ANNOTATOR' : str,
    'TRIMMOMATIC_adaptors' : ['file', str],
    'TRIMMOMATIC_index_barcodes' : str,
    'TRIMMOMATIC_memory' : int,
    'TRIMMOMATIC_palindrome' : 'file',
    'TRIMMOMATIC_simple_clip_threshold' : int,
    'TRIMMOMATIC_threads' : int,
    'VAMB_GPU' : bool,
    'VAMB_memory' : int,
    'VAMB_THREADS' : int,
    'working_dir' : 'directory'
}

# Aliases for keys.
# E.g., This module will always refer to minimum read length with 'READ_minlen'.
#       But if it cannot find 'READ_minlen', it will look for 'TRIMMOMATIC_minlen'.
ALIASES = {
    'READ_minlen' : ['TRIMMOMATIC_minlen'],
    'FASTP_adapters' : ['FASTP_adaptors'],
    'TRIMMOMATIC_adaptors' : ['TRIMMOMATIC_adapters']
}

##############################################
# Functions to parse the keys
##############################################

def validate_keytype(config, type_dict, key):
    # Get all possible aliases
    pos_keys = [key]
    if key in ALIASES:
        pos_keys += ALIASES[key]

    # Check for key in global type registry
    if key not in type_dict:
        raise KeyError(f"Configuration key not found in global registry: {key}")

    # Check for all aliases
    # Get value
    value = "MINTO_KEY_NOT_FOUND_EXCEPTION"
    for x in pos_keys:
        if x in config:
            # return None if the actual value in config file is 'None'
            if config[x] == None:
                return None
            value = config[x]
            break

    # If this key was not found, raise error for missing key
    if value == "MINTO_KEY_NOT_FOUND_EXCEPTION":
        return value

    # Get list of expected types
    expected_types = type_dict[key]
    if not isinstance(expected_types, list):
        expected_types = [expected_types]

    # Set up type discovery
    valid_type = False
    found_type = None

    # Check for each possible valid type

    for etype in expected_types:

        # Handle simple cases: etype is inbuilt python type
        if not isinstance(etype, str):
            if isinstance(value, etype):
                valid_type = True
                found_type = etype
                break

        # Handle special cases: file and directory
        if isinstance(etype, str):
            if etype == 'file':
                found_type = etype
                if os.path.isfile(value):
                    valid_type = True
            if etype == 'directory':
                found_type = etype
                if os.path.isdir(value):
                    valid_type = True

    # Raise error if right type was not found

    if not valid_type:
        if found_type == 'directory':
            raise TypeError(f"Directory mapped to configuration key '{key}' does not exist: {value}")
        elif found_type == 'file':
            raise TypeError(f"File mapped to configuration key '{key}' does not exist: {value}")
        else:
            raise TypeError("Configuration key '{}' should be of type {}, but got {}".format(key, ' or '.join(x.__name__ for x in expected_types), type(value).__name__))

    # return
    return(value)

def validate_required_key(config, key):

    # Check if key is present and the right type
    val = validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key)

    # If key was present but was explicitly 'None' then raise specific error
    if val is None:
        raise KeyError(f"Value for required configuration key '{key}' cannot be 'None'")

    # If key was missing then raise specific error
    if isinstance(val, str) and val == "MINTO_KEY_NOT_FOUND_EXCEPTION":
        raise KeyError(f"Missing required configuration key: {key}")

    # return
    return(val)

def validate_optional_key(config, key):

    # Check if key is present and the right type
    val = validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key)

    # If key was missing return None

    if isinstance(val, str) and val == "MINTO_KEY_NOT_FOUND_EXCEPTION":
        return None

    # return
    return(val)

##############################################
# Functions to verify the keys
##############################################

def check_allowed_values(key, value, allowed):
    if value not in allowed:
        raise ValueError("Invalid variable '{}={}' : must be one of {}".format(key, value, ', '.join(allowed)))

def check_number_is_odd(key, value):
    if value%2 == 0:
        raise ValueError(f"Invalid variable '{key}={value}' : must be an odd integer!")

def check_number_is_between(key, value, n_from, n_to):
    if value < n_from or value > n_to:
        raise ValueError(f"Invalid variable '{key}={value}' : must be in range [{n_from}, {n_to}]")

def check_input_directory(samples, locations):
    for x in samples:
        file_found = False
        for loc in locations:
            folder = "{}/{}/{}/{}".format(working_dir, omics, loc, x)
            if (os.path.exists(folder)):
                   file_found = True
        if not file_found:
            raise Exception(f"ERROR in {config_path}: fastq directory for sample {x} does not exist.")

##############################################
# Special functions to handle MINTO_MODE
##############################################

def get_minto_mode(config):
    MINTO_MODE = validate_required_key(config, 'MINTO_MODE')

    # Backward compatibility and common misnomers
    if MINTO_MODE in ['reference_genome', 'reference-genome', 'reference', 'refgenomes']:
        MINTO_MODE = 'refgenome'
    elif MINTO_MODE in ['MAGs', 'mag', 'mags']:
        MINTO_MODE = 'MAG'
    elif MINTO_MODE in ['db_genes', 'db-genes', 'genes_db', 'gene_catalog', 'gene-catalog']:
        MINTO_MODE = 'catalog'

    return MINTO_MODE

################################
# config file name:
# Getting the configfile's name is simple but just not obvious.
# See: https://stackoverflow.com/a/72096363
################################
config_path = 'configuration yaml file'
if '--configfile' in sys.argv:
    i = sys.argv.index('--configfile')
    config_path = sys.argv[i+1]
elif '--configfiles' in sys.argv:
    i = sys.argv.index('--configfiles')
    config_path = sys.argv[i+1]

print(" *******************************")
print(" Reading configuration yaml file: ", config_path)
print(" *******************************")
print("  ")

# Set script_dir relative to calling snakemake script

script_dir = os.path.join(os.path.dirname(workflow.basedir), 'scripts')

# Variables from configuration yaml file

minto_dir   = validate_required_key(config, 'minto_dir')

if NEED_PROJECT_VARIABLES:
    project_id  = validate_required_key(config, 'PROJECT')
    metadata    = validate_optional_key(config, 'METADATA')
    working_dir = validate_required_key(config, 'working_dir')
    omics       = validate_required_key(config, 'omics')

    check_allowed_values('omics', omics, ('metaG', 'metaT', 'metaG_metaT'))
