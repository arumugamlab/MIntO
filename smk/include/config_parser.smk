#!/usr/bin/env python

'''
Config file parser

Authors: Carmen Saenz, Mani Arumugam
'''

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
    'ILLUMINA' : list,
    'ILLUMINA_suffix' : list,
    'MAG_omics' : str,
    'MAIN_factor' : str,
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
    'TRIMMOMATIC_adaptors' : 'file',
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

def validate_keytype(config, type_dict, key):
    # Check for key in global type registry
    if key not in type_dict or config[key] is None:
        raise KeyError(f"Configuration key not found in global registry: {key}")
    expected_type = type_dict[key]

    # Get value
    value = config[key]

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
    if key not in config:
        raise KeyError(f"Missing required configuration key: {key}")

    # return
    return(validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key))

def validate_optional_key(config, key):

    # Check if key is present and the right type
    if key not in config or config[key] is None:
        return(None)

    # return
    return(validate_keytype(config, GLOBAL_CONFIG_KEYTYPES, key))

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

script_dir=workflow.basedir+"/../scripts"

# Variables from configuration yaml file

project_id  = validate_required_key(config, 'PROJECT')
omics       = validate_required_key(config, 'omics')
working_dir = validate_required_key(config, 'working_dir')
minto_dir   = validate_required_key(config, 'minto_dir')

metadata    = validate_optional_key(config, 'METADATA')

# Sanity check for omics
if omics not in ('metaG', 'metaT', 'metaG_metaT'):
    raise Exception(f"ERROR in {config_path}: omics={omics} is not allows. 'omics' variable should be metaG, metaT or metaG_metaT.")

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
else:
    print('WARNING in ', config_path, ': ILLUMINA list of samples is missing. Proceed with caution.')
