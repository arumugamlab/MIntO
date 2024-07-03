#!/usr/bin/env python

'''
Alignment, normalization and integration step

Authors: Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
from os import path
import pathlib

localrules: make_merged_genome_fna, make_genome_def, merge_MG_tables, fetchMG_genome_cds_faa, \
            config_yml_integration, read_map_stats

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'

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
if path.exists("{}/{}/{}/".format(working_dir, omics, read_dir)) is False:
    if path.exists("{}/{}/{}/".format(working_dir, omics, '6-corrected')) is True:
        read_dir='6-corrected'
    else:
        raise Exception("ERROR in {}: One of {} or 6-corrected must exist. Please correct {}".format(config_path, get_qc2_output_location(omics), working_dir))

print("NOTE: MIntO is using '{}' as read directory".format(read_dir))

main_factor = None
if config['MAIN_factor'] is not None:
    main_factor = config['MAIN_factor']

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
    else:
        # Make list of illumina samples, if ILLUMINA in config
        ilmn_samples = list()
        if 'ILLUMINA' in config:
            #print("Samples:")
            for ilmn in config["ILLUMINA"]:
                if path.exists("{}/{}/{}/{}".format(working_dir, omics, read_dir, ilmn)) is True:
                    #print(ilmn)
                    ilmn_samples.append(ilmn)
                else:
                    raise Exception('ERROR in ', config_path, ': ILLUMINA sample', ilmn, 'does not exist. Please, complete ', config_path)
        n_samples=len(ilmn_samples)+3
else:
    print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)

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
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be {}.".format(config_path, valid_minto_modes))

# Normalization
if config['abundance_normalization'] is None:
    print('ERROR in ', config_path, ': abundance_normalization variable is not set. "abundance_normalization" variable should be MG or TPM.')
else:
    normalization=config['abundance_normalization']
    normalization_modes=normalization.split(",")
    for m in normalization_modes:
        if m not in ("MG", "TPM"):
            print('ERROR in ', config_path, ': abundance_normalization variable is not correct. "abundance_normalization" variable should be MG or TPM.')

if config['alignment_identity'] is None:
    print('ERROR in ', config_path, ': alignment_identity variable is empty. Please, complete ', config_path)
elif type(config['alignment_identity']) != int:
    print('ERROR in ', config_path, ': alignment_identity variable is not an integer. Please, complete ', config_path)
elif type(config['alignment_identity']) == int:
    identity=config['alignment_identity']

if config['msamtools_filter_length'] is None:
    print('ERROR in ', config_path, ': msamtools_filter_length variable is empty. Please, complete ', config_path)
elif type(config['msamtools_filter_length']) != int:
    print('ERROR in ', config_path, ': msamtools_filter_length variable is not an integer. Please, complete ', config_path)

if config['NAME_reference'] is None and MINTO_MODE == 'catalog':
    raise Exception("ERROR in {}: 'NAME_reference' variable must be defined".format(config_path))

mag_omics = 'metaG'
gene_catalog_db="None"
gene_catalog_name="None"
if MINTO_MODE == 'MAG':
    if 'MAG_omics' in config and config['MAG_omics'] != None:
        mag_omics = config['MAG_omics']
    reference_dir="{wd}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes".format(wd = working_dir, mag_omics = mag_omics)
    print('NOTE: MIntO is using "'+ reference_dir+'" as PATH_reference variable')
else:
    if config['PATH_reference'] is None:
        print('ERROR in ', config_path, ': PATH_reference variable is empty. Please, complete ', config_path)
    elif path.exists(config['PATH_reference']) is False:
        print('ERROR in ', config_path, ': PATH_reference variable path does not exit. Please, complete ', config_path)
    else:
        if MINTO_MODE == 'refgenome':
            reference_dir=config["PATH_reference"]
            print('NOTE: MIntO is using "'+ reference_dir+'" as PATH_reference variable')
        elif MINTO_MODE == 'catalog':
            if path.exists(config['PATH_reference']+'/'+config['NAME_reference']) is True:
                print('NOTE: MIntO is using "'+ config['PATH_reference']+'/'+config['NAME_reference']+'" as PATH_reference and NAME_reference variables.')
                gene_catalog_db=config["PATH_reference"]
                gene_catalog_name=config["NAME_reference"]
            else:
                print('ERROR in ', config_path, ': NAME_reference variable does not exit. Please, complete ', config_path)

if config['BWA_threads'] is None:
    print('ERROR in ', config_path, ': BWA_threads variable is empty. Please, complete ', config_path)
elif type(config['BWA_threads']) != int:
    print('ERROR in ', config_path, ': BWA_threads variable is not an integer. Please, complete ', config_path)

if config['BWA_memory'] is None:
    print('ERROR in ', config_path, ': BWA_memory variable is empty. Please, complete ', config_path)
elif type(config['BWA_memory']) != int:
    print('ERROR in ', config_path, ': BWA_memory variable is not an integer. Please, complete ', config_path)

if 'MG' in normalization_modes and MINTO_MODE == 'catalog':
    raise Exception("ERROR in {}: In 'catalog' mode, only TPM normalization is allowed.".format(config_path))

# Define all the outputs needed by target 'all'

GENE_DB_TYPE = MINTO_MODE + '-genes'

bwaindex_db="DB/{subdir}/BWA_index/{analysis_name}".format(subdir=MINTO_MODE, analysis_name=MINTO_MODE)

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

def combined_genome_profiles():
    result = expand("{wd}/{omics}/9-mapping-profiles/{subdir}/{label}.p{identity}.{variant}.tsv",
                label = 'genome_abundances',
                wd = working_dir,
                omics = omics,
                subdir = MINTO_MODE,
                identity = identity,
                variant = ['abund.prop', 'abund.prop.genome', 'relabund.prop.genome']
                )
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
    reference_dir=config["PATH_reference"]
    def combined_genome_profiles():
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

rule all:
    input:
        combined_genome_profiles(),
        combined_gene_abundance_profiles(),
        mapping_statistics(),
        config_yaml(),
        print_versions.get_version_output(snakefile_name)
    default_target: True

###############################################################################################
# Useful python functions
###############################################################################################

# Get a sorted list of runs for a sample

def get_runs_for_sample(wildcards):
    if read_dir == '6-corrected':
        runs = [wildcards.sample]
    else:
        sample_dir = '{wd}/{omics}/{location}/{sample}'.format(wd=wildcards.wd, omics=wildcards.omics, location=read_dir, sample=wildcards.sample)
        runs = sorted([ re.sub("\.1\.fq\.gz", "", path.basename(f)) for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith('.1.fq.gz') ])
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
# Prepare genes for mapping to MAGs or publicly available genomes
## Generate MAGs or publicly available genomes index
###############################################################################################

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
        minto_mode='MAG|refgenome'
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

# Memory is based on number of MAGs - 50 MB per genome; increase by 50 MB each new attempt.
rule genome_bwaindex:
    input:
        fna=get_genome_fna,
        fasta_merge=rules.make_merged_genome_fna.output.fasta_merge
    output:
        "{wd}/DB/{minto_mode}/BWA_index/{minto_mode}.0123",
        "{wd}/DB/{minto_mode}/BWA_index/{minto_mode}.amb",
        "{wd}/DB/{minto_mode}/BWA_index/{minto_mode}.ann",
        "{wd}/DB/{minto_mode}/BWA_index/{minto_mode}.bwt.2bit.64",
        "{wd}/DB/{minto_mode}/BWA_index/{minto_mode}.pac",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{minto_mode}_bwaindex.log"
    threads: 2
    resources:
        mem = lambda wildcards, input, attempt: 0.05*len(input.fna)*attempt
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        time (
                bwa-mem2 index {input.fasta_merge} -p {wildcards.minto_mode}
                ls {wildcards.minto_mode}.*
                rsync -a {wildcards.minto_mode}.* $(dirname {output[0]})/
            ) >& {log}
        """

rule genome_mapping_profiling:
    input:
        bwaindex=rules.genome_bwaindex.output,
        genome_def=rules.make_genome_def.output.genome_def,
        fwd=get_fwd_files_only,
        rev=get_rev_files_only,
        bed_mini="{wd}/DB/{minto_mode}/{minto_mode}-genes.bed.mini"
    output:
        sorted=         "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.sorted.bam",
        index=          "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.sorted.bam.bai",
        bwa_log=        "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.log",
        raw_all_seq=    "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
        raw_prop_seq=   "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.txt.gz",
        raw_prop_genome="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.genome.txt.gz",
        rel_prop_genome="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.relabund.prop.genome.txt.gz",
        absolute_counts="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.gene_abundances.p{identity}.bed.gz"
    shadow:
        "minimal"
    params:
        sample_alias=lambda wildcards: sample2alias[wildcards.sample],
        length=config["msamtools_filter_length"],
        sort_threads=lambda wildcards, threads, resources: int(1+threads/4),
        sort_memory=lambda wildcards, threads, resources: int(resources.mem/int(1+threads/4)),
        mapped_reads_threshold=config["MIN_mapped_reads"],
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{sample}.p{identity}.map_profile.log"
    wildcard_constraints:
        identity=r'\d+',
        minto_mode='MAG|refgenome'
    threads:
        config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" # BWA + samtools + msamtools + perl + parallel
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

            # Run multiple msamtools modes + samtools sort in parallel using GNU parallel.
            # We need 4 threads for msamtools and 'params.sort_threads' threads for samtools sort.
            # We cap it by 'threads' limit to parallel.
            common_args="--label {wildcards.omics}.{params.sample_alias} --total=$total_reads --mincount={params.mapped_reads_threshold} --pandas"
            parallel --jobs {threads} <<__EOM__
msamtools profile aligned.bam $common_args -o {output.raw_all_seq}     --multi=all  --unit=abund --nolen
msamtools profile aligned.bam $common_args -o {output.raw_prop_seq}    --multi=prop --unit=abund --nolen
msamtools profile aligned.bam $common_args -o {output.raw_prop_genome} --multi=prop --unit=abund --nolen --genome {input.genome_def}
msamtools profile aligned.bam $common_args -o {output.rel_prop_genome} --multi=prop --unit=rel           --genome {input.genome_def}
samtools sort aligned.bam -o sorted.bam -@ {params.sort_threads} -m {params.sort_memory}G --output-fmt=BAM
__EOM__

            # Index sorted.bam file, while also exporting it
            # Get mini-bed into shadowdir
            parallel --jobs {threads} <<__EOM__
samtools index sorted.bam sorted.bam.bai -@ {threads}
rsync -a sorted.bam {output.sorted}
rsync -a {input.bed_mini} in.bed
__EOM__

            # Use bedtools multicov to generate read-counts per gene
            (echo -e 'gene_length\\tID\\t{wildcards.omics}.{params.sample_alias}';
             bedtools multicov -bams sorted.bam -bed in.bed | cut -f4- | grep -v $'\\t0$') | gzip -c > out.bed.gz

            # Rsync output files
            parallel --jobs {threads} <<__EOM__
rsync -a sorted.bam.bai {output.index}
rsync -a out.bed.gz {output.absolute_counts}
__EOM__
        ) >& {log}
        """

rule old_merge_individual_msamtools_profiles:
    input:
        single=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.{type}.txt.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                minto_mode = wildcards.minto_mode,
                identity = wildcards.identity,
                type = wildcards.type,
                sample = ilmn_samples)
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{minto_mode}/Combined.p{identity}.profile.{type}.txt"
    run:
        import pandas as pd
        df = pd.read_csv(input.single[0], comment='#', header=0, index_col='ID', sep = "\t")
        for i in range(1, len(input.single)):
            df2 = pd.read_csv(input.single[i], comment='#', header=0, index_col='ID', sep = "\t")
            df  = df.join(df2, how='outer')
        df.to_csv(output.combined, sep = "\t", index = True)

def strip_commonpath(list_of_paths):
    topdir = os.path.commonpath(list_of_paths)
    stripped = [x.replace(topdir, '') for x in list_of_paths]
    stripped = [re.sub(r'^/+', '', x) for x in stripped]
    return(stripped)

# Get a list of msamtools profiles based on wildcards
def get_msamtools_profiles(wildcards):

    # Ensure right MINTO_MODE and filename combinations
    if     (wildcards.minto_mode in ['MAG', 'refgenome'] and wildcards.filename == 'gene_abundances') \
        or (wildcards.minto_mode == 'catalog' and wildcards.filename == 'genome_abundances'):
        raise Exception("MIntO mode {} cannot generate {} profiles".format(wildcards.minto_mode, wildcards.filename))

    # Ensure right MINTO_MODE and norm combination
    if (wildcards.minto_mode == 'catalog' and wildcards.type != 'TPM'):
        raise Exception("MIntO mode {} cannot generate {} profiles".format(wildcards.minto_mode, wildcards.type))

    profiles = expand("{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.{type}.txt.gz",
                            wd = wildcards.wd,
                            omics = wildcards.omics,
                            minto_mode = wildcards.minto_mode,
                            identity = wildcards.identity,
                            type = wildcards.type,
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
        files = lambda wildcards, input: ','.join(strip_commonpath(input.single))
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/merge_msamtools_profiles.{filename}.p{identity}.{type}.log"
    resources:
        mem=30
    threads: lambda wildcards,input: min(10, len(input.single))
    wildcard_constraints:
        identity='[0-9]+'
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
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

###############################################################################################
# Prepare genes for mapping to gene-database
## First, the index has to be generated
## Mapping, computation of read counts and TPM normalization is done in the same rule
## TPM normalization: sequence depth and genes’ length
###############################################################################################

# Memory is expected to be high, since gene catalogs are large.
# Assume the equivalent of 1000 genomes to begin with.
# It is roughly 50 MB per bacterial genome.
# So, start at 50 GB and increase by 50 GB each new attempt.
rule gene_catalog_bwaindex:
    input:
        genes="{gene_catalog_path}/{gene_catalog_name}"
    output:
        "{gene_catalog_path}/BWA_index/{gene_catalog_name}.0123",
        "{gene_catalog_path}/BWA_index/{gene_catalog_name}.amb",
        "{gene_catalog_path}/BWA_index/{gene_catalog_name}.ann",
        "{gene_catalog_path}/BWA_index/{gene_catalog_name}.bwt.2bit.64",
        "{gene_catalog_path}/BWA_index/{gene_catalog_name}.pac"
    shadow:
        "minimal"
    log:
        genes="{gene_catalog_path}/{gene_catalog_name}.bwaindex.log"
    threads: 2
    resources:
        mem = lambda wildcards, attempt: attempt*50
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        time (
            bwa-mem2 index {wildcards.gene_catalog_path}/{wildcards.gene_catalog_name} -p {wildcards.gene_catalog_name}
            ls {wildcards.gene_catalog_name}.*
            rsync -a {wildcards.gene_catalog_name}.* $(dirname {output[0]})/
        ) &> {log}
        """

rule gene_catalog_mapping_profiling:
    input:
        bwaindex=expand("{gene_catalog_path}/BWA_index/{gene_catalog_name}.{ext}",
                gene_catalog_path=gene_catalog_db,
                gene_catalog_name=gene_catalog_name,
                ext=['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']),
        fwd=get_fwd_files_only,
        rev=get_rev_files_only,
    output:
        filtered=   "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.bam",
        bwa_log=    "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.log",
        profile_tpm="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.TPM.txt.gz",
        map_profile="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz"
    shadow:
        "minimal"
    params:
        sample_alias=lambda wildcards: sample2alias[wildcards.sample],
        length=config["msamtools_filter_length"],
        mapped_reads_threshold=config["MIN_mapped_reads"],
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{sample}.p{identity}_bwa.log"
    wildcard_constraints:
        minto_mode='catalog'
    threads:
        config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #config["conda_env2_yml"] #BWA + samtools
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
rsync -a aligned.bam {output.filtered}
rsync -a profile.TPM.txt.gz {output.profile_tpm}
rsync -a profile.abund.all.txt.gz {output.map_profile}
__EOM__
        ) >& {log}
        """

###############################################################################################
# Computation of read counts to genes belonging to MAGs or publicly available genomes
###############################################################################################

# Ideally, bedtools multicov is run in the same rule as bwa-mapping to avoid reading BAM file over NFS.
# However, if the BAM file exists already, just running bedtools on the BAM files on NFS will be faster.
# I have currently disabled this competition between the rules for making 'bed.gz' by renaming output file.
# Could be turned back on in the future.

ruleorder: gene_abund_compute > genome_mapping_profiling

rule gene_abund_compute:
    input:
        bam=     "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.sorted.bam",
        index=   "{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.p{identity}.filtered.sorted.bam.bai",
        bed_mini="{wd}/DB/{minto_mode}/{minto_mode}-genes.bed.mini"
    output:
        absolute_counts="{wd}/{omics}/9-mapping-profiles/{minto_mode}/{sample}/{sample}.gene_abundances.p{identity}.bed.gz2"
    shadow:
        "minimal"
    params:
        sample_alias=lambda wildcards: sample2alias[wildcards.sample],
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/{sample}.gene_abundances.p{identity}.log"
    threads: 1
    resources:
        mem=5
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bedtools
    shell:
        """
        time (
            # Random sleep to avoid choking NFS
            sleep $((1 + $RANDOM % 120))

            # Rsync input files
            rsync -a {input.bam} in.bam
            rsync -a {input.index} in.bam.bai
            rsync -a {input.bed_mini} in.bed

            # Do the calculation
            (echo -e 'gene_length\\tID\\t{wildcards.omics}.{params.sample_alias}';
            bedtools multicov -bams in.bam -bed in.bed | cut -f4- | grep -v $'\\t0$') | gzip -c > out.bed.gz

            # Rsync output files
            rsync -a out.bed.gz {output.absolute_counts}
        ) >& {log}
        """

# Merge individual bedtools multicov BED files from genome mapping
# We don't set '--zeroes' because rule 'gene_abund_compute' removed all zero entries in individual BED files.
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
        files = lambda wildcards, input: ','.join(strip_commonpath(input.single))
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
        config["minto_dir"]+"/envs/r_pkgs.yml"
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
# Normalization of read counts to 10 marker genes (MG normalization)
## fetchMGs identifies 10 universal single-copy phylogenetic MGs
## (COG0012, COG0016, COG0018, COG0172, COG0215, COG0495, COG0525, COG0533, COG0541, and COG0552)
###############################################################################################

# Run fetchMGs
# fetchMG cannot handle '.' in gene names, so we replace '.' with '__MINTO_DOT__' in fasta headers; then change back in output.
rule fetchMG_genome_cds_faa:
    input:
        cds_faa      = "{wd}/DB/{minto_mode}/2-postprocessed/{genome}.faa",
        fetchMGs_dir = "{minto_dir}/data/fetchMGs-1.2".format(minto_dir = minto_dir)
    output: '{wd}/DB/{minto_mode}/fetchMGs/{genome}/{genome}.marker_genes.table'
    shadow:
        "minimal"
    log: '{wd}/logs/DB/{minto_mode}/fetchMGs/{genome}.log'
    threads: 8
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (
            sed 's/\\s.*//;s/\\./__MINTO_DOT__/g' {input.cds_faa} > HEADER_FIXED.faa
            {input.fetchMGs_dir}/fetchMGs.pl -outdir out -protein_only -threads {threads} -x {input.fetchMGs_dir}/bin -m extraction HEADER_FIXED.faa
            sed 's/__MINTO_DOT__/./g' out/HEADER_FIXED.all.marker_genes_scores.table > {output}
        ) >& {log}
        """

def get_genome_MG_tables(wildcards):
    #Collect the CDS faa files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    result = expand("{wd}/DB/{minto_mode}/fetchMGs/{genome}/{genome}.marker_genes.table",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(result)

rule merge_MG_tables:
    input: get_genome_MG_tables
    output: "{wd}/DB/{minto_mode}/genomes.marker_genes.table"
    log: "{wd}/logs/DB/{minto_mode}/merge_marker_genes_scores.table.log"
    shell:
        """
        time (
                head -n 1 {input[0]} > {output}
                for file in {input}; do
                    awk 'FNR>1' ${{file}} >> {output}
                done
        ) >& {log}
        """

###############################################################################################
# Normalization of read counts by sequence depth and genes’ length (TPM normalization)
# or marker genes (MG normalization)
# Normalize and add prefix to samples so that metaG and metaT do not clash when combined in future
###############################################################################################

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

rule gene_abund_normalization:
    input:
        unpack(get_gene_abund_normalization_input)
    output:
        norm_counts="{wd}/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.{norm}.tsv"
    params:
        mapped_reads_threshold=config["MIN_mapped_reads"],
        optional_arg_MG = lambda wildcards, input: "" if wildcards.norm == "TPM" else "--MG " + input.genomes_marker_genes
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{minto_mode}/gene_abundances.p{identity}.{norm}.log"
    wildcard_constraints:
        minto_mode='MAG|refgenome'
    threads: 4
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        time (
            Rscript {script_dir}/normalize_profiles.R --normalize {wildcards.norm} --threads {threads} --memory {resources.mem} --bed {input.absolute_counts} {params.optional_arg_MG} --out {output.norm_counts} --min-read-count {params.mapped_reads_threshold}
        ) &> {log}
        """

###############################################################################################
## Mappability rate
## Mappability stats
## Multimapping read count
###############################################################################################
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
        identity=r'\d+'
    threads: 1
    resources:
        mem=2
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
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
    input: lambda wildcards: "{wd}/{mag_omics}/mapping.yaml".format(wd=wildcards.wd, mag_omics=mag_omics)
    output:
        config_file="{wd}/data_integration.yaml"
    params:
        mapped_reads_threshold=config["MIN_mapped_reads"],
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

######################
# Program settings
######################

alignment_identity: {identity}
abundance_normalization: MG
MINTO_MODE: {MINTO_MODE}

MERGE_threads: 4
MERGE_memory: 5

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
