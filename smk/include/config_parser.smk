#!/usr/bin/env python

'''
Config file parser

Authors: Carmen Saenz, Mani Arumugam
'''

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
    'enable_COASSEMBLY' : bool,
    'EXCLUDE_ASSEMBLY_TYPES' : list,
    'FASTP_adapters' : str,
    'FASTP_front_mean_qual' : int,
    'FASTP_memory' : int,
    'FASTP_min_length' : int,
    'FASTP_tail_mean_qual' : int,
    'FASTP_threads' : int,
    'GTDB_TAXONOMY_VERSION' : str,
    'HYBRID' : dict,
    'ILLUMINA' : list,
    'ILLUMINA_suffix' : list,
    'MAG_omics' : str,
    'MAIN_factor' : str,
    'MAX_RAM_GB_PER_JOB' : int,
    'MEGAHIT_memory' : int,
    'MEGAHIT_presets' : list,
    'MEGAHIT_threads' : int,
    'MERGE_ILLUMINA_SAMPLES': dict,
    'MERGE_memory' : int,
    'MERGE_threads' : int,
    'METADATA' : 'file',
    'METAFLYE_presets' : dict,
    'metaphlan_version' : str,
    'METASPADES_custom_build' : 'file',
    'METASPADES_hybrid_max_k' : int,
    'METASPADES_illumina_max_k' : int,
    'METASPADES_memory' : int,
    'METASPADES_qoffset' : str,
    'METASPADES_threads' : int,
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
    'TRIMMOMATIC_adaptors' : str,
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

ALIASES = {
    'READ_minlen' : ['TRIMMOMATIC_minlen']
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
    value = None
    for x in pos_keys:
        if x in config:
            value = config[x]
    if not value:
        return None

    expected_type = type_dict[key]

    # Handle simple cases: expected_type is inbuilt python type
    if not isinstance(expected_type, str):
        if not isinstance(value, expected_type):
            raise TypeError(f"Configuration key '{key}' should be of type {expected_type.__name__}, but got {type(value).__name__}")

    # Handle special cases: file and directory
    if isinstance(expected_type, str):
        if expected_type == 'file':
            if not os.path.isfile(value):
                raise TypeError(f"File mapped to configuration key '{key}' does not exist: {value}")
        if expected_type == 'directory':
            if not os.path.isdir(value):
                raise TypeError(f"Directory mapped to configuration key '{key}' does not exist: {value}")

    # return
    return(value)

def validate_required_key(config, key):

    # Check if key is present and the right type
    val = validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key)

    if not val:
        raise KeyError(f"Missing required configuration key: {key}")

    # return
    return(val)

def validate_optional_key(config, key):

    # Check if key is present and the right type
    val = validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key)

    # return
    return(val)

################################
# config file name:
# Getting the configfile's name is simple but just not obvious.
# See: https://stackoverflow.com/a/72096363
################################
config_path = 'configuration yaml file' #args[args_idx+1]
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
project_id  = validate_required_key(config, 'PROJECT')
metadata    = validate_optional_key(config, 'METADATA')
working_dir = validate_required_key(config, 'working_dir')
omics       = validate_required_key(config, 'omics')

check_allowed_values('omics', omics, ('metaG', 'metaT', 'metaG_metaT'))
