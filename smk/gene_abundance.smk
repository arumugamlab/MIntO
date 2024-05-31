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
        modify_genome_fasta_header, config_yml_integration, read_map_stats, merge_msamtools_genome_mapping_profiles

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'

# We prefer original non-error-corrected reads for profiling to preserve strain variation in population
read_dir=get_qc2_output_location(omics)

# If the original reads have been removed for whatever reasons (running out of space, may be?),
# then we look for error-corrected reads
if path.exists("{}/{}/{}/".format(working_dir, omics, read_dir)) is False:
    if path.exists("{}/{}/{}/".format(working_dir, omics, '6-corrected')) is True:
        read_dir='6-corrected'
    else:
        raise Error("ERROR in {}: One of {} or 6-corrected must exist. Please correct {}".format(config_path, get_qc2_output_location(omics), working_dir))

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
                    raise TypeError('ERROR in ', config_path, ': ILLUMINA sample', ilmn, 'does not exist. Please, complete ', config_path)
        n_samples=len(ilmn_samples)+3
else:
    print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)

if config['map_reference'] in ("MAG", "reference_genome","genes_db"):
    map_reference=config["map_reference"]
else:
    print('ERROR in ', config_path, ': map_reference variable is not correct. "map_reference" variable should be MAG, reference_genome or genes_db.')

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

if config['NAME_reference'] is None and map_reference == 'genes_db':
    print('ERROR in ', config_path, ': NAME_reference variable does not exit. Please, complete ', config_path)

mag_omics = 'metaG'
gene_catalog_db="None"
gene_catalog_name="None"
if map_reference == 'MAG':
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
        if map_reference == 'reference_genome':
            print('NOTE: MIntO is using "'+ config['PATH_reference']+'" as PATH_reference variable')
            reference_dir=config["PATH_reference"]
        elif map_reference == 'genes_db':
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

fetchMGs_dir=None
if 'MG' in normalization_modes and map_reference in ("MAG", "reference_genome"):
    if config['fetchMGs_dir'] is None:
        print('ERROR in ', config_path, ': fetchMGs_dir variable is empty. Please, complete ', config_path)
    elif path.exists(config['fetchMGs_dir']) is False:
        print('ERROR in ', config_path, ': fetchMGs_dir variable path does not exit. Please, complete ', config_path)
    elif path.exists(config['fetchMGs_dir']) is True:
        fetchMGs_dir=config["fetchMGs_dir"]
if 'MG' in normalization_modes and map_reference in ("genes_db"):
    print('WARNING in ', config_path, ': In "genes_db" mode, only TPM normalization is allowed.')


# Define all the outputs needed by target 'all'
if map_reference == 'MAG':
    post_analysis_out="MAG-genes"
    post_analysis_dir="9-MAG-genes-post-analysis"

elif map_reference == 'reference_genome':
    post_analysis_out="refgenome-genes"
    post_analysis_dir="9-refgenome-genes-post-analysis"

elif map_reference == 'genes_db':
    post_analysis_out="db-genes"
    post_analysis_dir="9-genes-db-post-analysis"

bwaindex_db="DB/{analysis_dir}/BWA_index/{analysis_name}".format(analysis_dir=post_analysis_dir, analysis_name=post_analysis_out)

def combined_genome_profiles():
    result = expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{label}.p{identity}.profile.{extension}",
                label = 'all',
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                identity = identity,
                extension = ['abund.prop.txt', 'abund.prop.genome.txt', 'relabund.prop.genome.txt']
                )
    return(result)

def combined_gene_abundance_profiles():
    result = expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.{norm}.csv",
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                norm = normalization_modes,
                identity = identity)
    return(result)

if map_reference == 'genes_db':
    reference_dir=config["PATH_reference"]
    def combined_genome_profiles():
        return()

    def gene_abundances_bwa_out(): # CHECK THIS PART - do not generate bwa index, normalization in a different way
        result = expand("{gene_catalog_path}/BWA_index/{gene_catalog_name}.pac",
                    gene_catalog_path=gene_catalog_db,
                    gene_catalog_name=gene_catalog_name),\
        expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.bam",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db-genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity),\
        expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.log",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db-genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity),\
        expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.TPM.txt.gz",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db-genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity)
        return(result)

def mapping_statistics():
    result = expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.{stats}.txt",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = post_analysis_out,
                    identity = identity,
                    stats = ['maprate', 'mapstats', 'multimap'])
    return(result)

def config_yaml():
    result = "{wd}/data_integration.yaml".format(
                wd = working_dir)
    return(result)

rule all:
    input:
        #gene_abundances_bwa_out(),
        combined_genome_profiles(),
        combined_gene_abundance_profiles(),
        mapping_statistics(),
        config_yaml()

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
    result = expand("{wd}/DB/9-{post_analysis_out}-post-analysis/fasta/{genome}.fna.hdr-mod",
                    wd=wildcards.wd,
                    post_analysis_out=wildcards.post_analysis_out,
                    genome=genomes)
    return(result)

rule modify_genome_fasta_header:
    input: "{wd}/DB/{post_analysis_dir}/genomes/{genome}/{genome}.fna",
    output: "{wd}/DB/{post_analysis_dir}/fasta/{genome}.fna.hdr-mod"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.reformat_fna.log"
    shell:
        """
        time (sed "s/^>gnl|X|/>{wildcards.genome}|/" {input} > {output}) >& {log}
        """

rule make_merged_genome_fna:
    input: get_genome_fna
    output:
        fasta_merge="{wd}/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}.fna"
    log:
        "{wd}/logs/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}.merge_genome.log"
    wildcard_constraints:
        post_analysis_out='MAG-genes|refgenome-genes'
    shell:
        """
        time (cat {input} > {output}) >& {log}
        """

rule make_genome_def:
    input: get_genome_fna
    output:
        genome_def="{wd}/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}.genome.def"
    shell:
        """
        grep '^>' {input} | sed 's^.*/^^' | sed 's/.fna.hdr-mod:>/\\t/' >> {output}
        """

# Memory is based on number of MAGs - 50 MB per genome; increase by 50 MB each new attempt.
rule genome_bwaindex:
    input:
        fna=get_genome_fna,
        fasta_merge=rules.make_merged_genome_fna.output.fasta_merge
    output:
        "{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.0123",
        "{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.amb",
        "{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.ann",
        "{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.bwt.2bit.64",
        "{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.pac",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}_bwaindex.log"
    threads: 2
    resources:
        mem = lambda wildcards, input, attempt: 0.05*len(input.fna)*attempt
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        time (\
                bwa-mem2 index {input.fasta_merge} -p {wildcards.post_analysis_out}
                ls {wildcards.post_analysis_out}.*
                rsync -a {wildcards.post_analysis_out}.* $(dirname {output[0]})/
            ) >& {log}
        """

rule genome_mapping_profiling:
    input:
        bwaindex=rules.genome_bwaindex.output,
        genome_def=rules.make_genome_def.output.genome_def,
        fwd=get_fwd_files_only,
        rev=get_rev_files_only,
    output:
        sorted="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam",
        index="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam.bai",
        bwa_log="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.log",
        raw_all_seq="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
        raw_prop_seq="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.txt.gz",
        raw_prop_genome="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.genome.txt.gz",
        rel_prop_genome="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.relabund.prop.genome.txt.gz",
    shadow:
        "minimal"
    params:
        length=config["msamtools_filter_length"],
        sort_threads=lambda wildcards, threads, resources: int(1+threads/4),
        sort_memory=lambda wildcards, threads, resources: int(resources.mem/int(1+threads/4)),
        mapped_reads_threshold=config["MIN_mapped_reads"],
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}.p{identity}_bwa.log"
    wildcard_constraints:
        identity=r'\d+'
    threads:
        config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" # BWA + samtools + msamtools + perl
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
        common_args="--label {wildcards.omics}.{wildcards.sample} --total=$total_reads --mincount={params.mapped_reads_threshold} --pandas"
        msamtools profile aligned.bam $common_args -o {output.raw_all_seq}     --multi=all  --unit=abund --nolen
        msamtools profile aligned.bam $common_args -o {output.raw_prop_seq}    --multi=prop --unit=abund --nolen
        msamtools profile aligned.bam $common_args -o {output.raw_prop_genome} --multi=prop --unit=abund --nolen --genome {input.genome_def}
        msamtools profile aligned.bam $common_args -o {output.rel_prop_genome} --multi=prop --unit=rel           --genome {input.genome_def}

        samtools sort aligned.bam -o sorted.bam -@ {params.sort_threads} -m {params.sort_memory}G --output-fmt=BAM
        samtools index sorted.bam sorted.bam.bai -@ {threads}

        rsync -a sorted.bam {output.sorted}
        rsync -a sorted.bam.bai {output.index}
        ) >& {log}
        """

rule old_merge_individual_msamtools_profiles:
    input:
        single=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.{type}.txt.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                post_analysis_out = wildcards.post_analysis_out,
                identity = wildcards.identity,
                type = wildcards.type,
                sample = ilmn_samples)
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/Combined.p{identity}.profile.{type}.txt"
    run:
        import pandas as pd
        df = pd.read_csv(input.single[0], comment='#', header=0, index_col='ID', sep = "\t")
        for i in range(1, len(input.single)):
            df2 = pd.read_csv(input.single[i], comment='#', header=0, index_col='ID', sep = "\t")
            df  = df.join(df2, how='outer')
        df.to_csv(output.combined, sep = "\t", index = True)

###############################################################################################
# The following two rules are identical. But first is for reference_genomes and MAGs. Second is for gene-catalogs
###############################################################################################

# Merge individual profiles from genome mapping

rule merge_msamtools_genome_mapping_profiles:
    input:
        single=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.{type}.txt.gz",
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            post_analysis_out = wildcards.post_analysis_out,
                                            identity = wildcards.identity,
                                            type = wildcards.type,
                                            sample = ilmn_samples)
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.profile.{type}.txt.raw"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/merge_msamtools_profiles.p{identity}.profile.{type}.log"
    resources:
        mem=30
    run:
        import shutil
        combine_profiles(input.single, 'combined.txt', log, key_columns=['ID'])
        shutil.copy2('combined.txt', output.combined)

# Merge normalized gene abundance or transcript profiles from gene catalog (TPM normalization)

rule merge_msamtools_gene_mapping_profiles:
    input:
        profile_tpm=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.{type}.txt.gz",
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            post_analysis_out = wildcards.post_analysis_out,
                                            identity = wildcards.identity,
                                            type = wildcards.type,
                                            sample=ilmn_samples)
    output:
        profile_tpm_all="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.{type}.csv.raw"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/merge_msamtools_profiles.p{identity}.profile.{type}.log"
    resources:
        mem=30
    wildcard_constraints:
        post_analysis_out='db-genes'
    run:
        import shutil
        combine_profiles(input.profile_tpm, 'combined.txt', log, key_columns=['ID'])
        shutil.copy2('combined.txt', output.profile_tpm_all)

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
        time (bwa-mem2 index {wildcards.gene_catalog_path}/{wildcards.gene_catalog_name} -p {wildcards.gene_catalog_name}
        ls {wildcards.gene_catalog_name}.*
        rsync -a {wildcards.gene_catalog_name}.* $(dirname {output[0]})/ ) &> {log}
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
        filtered="{wd}/{omics}/9-mapping-profiles/db-genes/{sample}/{sample}.p{identity}.filtered.bam",
        bwa_log="{wd}/{omics}/9-mapping-profiles/db-genes/{sample}/{sample}.p{identity}.filtered.log",
        profile_tpm="{wd}/{omics}/9-mapping-profiles/db-genes/{sample}/{sample}.p{identity}.filtered.profile.TPM.txt.gz",
        map_profile="{wd}/{omics}/9-mapping-profiles/db-genes/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz"
    shadow:
        "minimal"
    params:
        length=config["msamtools_filter_length"],
        mapped_reads_threshold=config["MIN_mapped_reads"],
        multiple_runs = lambda wildcards: "yes" if len(get_runs_for_sample(wildcards)) > 1 else "no",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/db-genes/{sample}.p{identity}_bwa.log"
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
        common_args="--label {wildcards.omics}.{wildcards.sample} --total=$total_reads --mincount={params.mapped_reads_threshold} --pandas"
        msamtools profile aligned.bam $common_args -o profile.TPM.txt.gz --multi prop --unit tpm
        msamtools profile aligned.bam $common_args -o profile.abund.all.txt.gz --multi all --unit abund --nolen

        rsync -a aligned.bam {output.filtered}
        rsync -a profile.TPM.txt.gz {output.profile_tpm}
        rsync -a profile.abund.all.txt.gz {output.map_profile}
        ) >& {log}
        """

###############################################################################################
# Computation of read counts to genes belonging to MAGs or publicly available genomes
###############################################################################################

rule gene_abund_compute:
    input:
        bam="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam",
        index="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam.bai",
        bed_mini="{wd}/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}.bed.mini"
    output:
        absolute_counts="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.genes_abundances.p{identity}.bed"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}.genes_abundances_counts.p{identity}.log"
    threads: 1
    resources:
        mem=5
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bedtools
    shell:
        """
        time (
            echo -e 'chr\\tstart\\tstop\\tstrand\\tfeature\\tID\\t{wildcards.sample}' > bed_file
            bedtools multicov -bams {input.bam} -bed {input.bed_mini} >> bed_file
            rsync bed_file {output.absolute_counts}
        ) >& {log}
        """

rule merge_gene_abund:
    input:
        single=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.genes_abundances.p{identity}.bed",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    post_analysis_out = wildcards.post_analysis_out,
                    identity = wildcards.identity,
                    sample=ilmn_samples),
    output:
        combined="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.bed.raw"
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances_counts.p{identity}.log"
    threads: 1
    resources:
        mem=30
    run:
        import shutil
        combine_profiles(input.single, 'combined.txt', log, key_columns=['chr','start','stop','strand','feature','ID'])
        shutil.copy2('combined.txt', output.combined)

# If this is a BED file, then sample names are just A, B, C.
# If this is a gene-catalog mode TPM.csv, then it is metaG.A, metaG.B, metaG.C.
# We use --prefix to prepend before relabeling.
rule relabel_merged_gene_abund:
    input:
        metadata=metadata,
        original="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{label}.p{identity}.{suffix}.raw"
    output:
        relabeled="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{label}.p{identity}.{suffix}"
    wildcard_constraints:
        label="all|genes_abundances",
        identity=r'\d+',
        suffix="|".join(['profile.abund.prop.txt', 'profile.abund.prop.genome.txt', 'profile.relabund.prop.genome.txt', 'bed'])
    shadow:
        "minimal"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/relabel_{label}.p{identity}.{suffix}.log"
    params:
        prefix = lambda wildcards: "" if (wildcards.suffix.endswith('bed')) else "--prefix {}.".format(omics)
    threads: 4
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    resources:
        mem=30
    shell:
        """
        time (Rscript {script_dir}/relabel_profiles_using_metadata.R --from sample --to sample_alias --threads {threads} --metadata {input.metadata} --input {input.original} --output relabeled.txt {params.prefix}) >& {log}
        rsync -a relabeled.txt {output.relabeled}
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
        cds_faa='{wd}/DB/{post_analysis_dir}/3-merged-features/{genome}.faa',
        fetchMGs_dir=str(fetchMGs_dir)
    output: '{wd}/DB/{post_analysis_dir}/fetchMGs/{genome}/{genome}.marker_genes.table'
    shadow:
        "minimal"
    log: '{wd}/logs/DB/{post_analysis_dir}/fetchMGs/{genome}.log'
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
    result = expand("{wd}/DB/{post_analysis_dir}/fetchMGs/{genome}/{genome}.marker_genes.table",
                    wd=wildcards.wd,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    genome=genomes)
    return(result)

rule merge_MG_tables:
    input: get_genome_MG_tables
    output: "{wd}/DB/{post_analysis_dir}/genomes.marker_genes.table"
    log: "{wd}/logs/DB/{post_analysis_dir}/merge_marker_genes_scores.table.log"
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
            'absolute_counts' : "{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.bed".format(
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    post_analysis_out = wildcards.post_analysis_out,
                    identity = wildcards.identity),
            'genomes_marker_genes' : "{wd}/DB/9-{post_analysis_out}-post-analysis/genomes.marker_genes.table".format(
                    wd = wildcards.wd,
                    post_analysis_out = wildcards.post_analysis_out)
            }
    else:
        return {
            'absolute_counts' : "{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.bed".format(
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    post_analysis_out = wildcards.post_analysis_out,
                    identity = wildcards.identity)
            }

rule gene_abund_normalization:
    input:
        unpack(get_gene_abund_normalization_input)
    output:
        norm_counts="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.{norm}.csv"
    params:
        mapped_reads_threshold=config["MIN_mapped_reads"],
        optional_arg_MG = lambda wildcards, input: "" if wildcards.norm == "TPM" else "--MG " + input.genomes_marker_genes
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/genes_abundances.p{identity}.{norm}.log"
    threads: 4
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        time (Rscript {script_dir}/normalize_profiles.R --normalize {wildcards.norm} --threads {threads} --memory {resources.mem} --bed {input.absolute_counts} {params.optional_arg_MG} --out {output.norm_counts} --sample-prefix {wildcards.omics}. --min-read-count {params.mapped_reads_threshold}) &> {log}
        """

###############################################################################################
## Mappability rate
## Mappability stats
## Multimapping read count
###############################################################################################
rule read_map_stats:
    input:
        map_profile=lambda wildcards: expand("{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            post_analysis_out = wildcards.post_analysis_out,
                                            identity = wildcards.identity,
                                            sample=ilmn_samples)
    output:
        maprate="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.maprate.txt",
        mapstats="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.mapstats.txt",
        multimap="{wd}/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.multimap.txt"
    params:
        map_profile_list=lambda wildcards, input: ",".join(input.map_profile),
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/{post_analysis_out}/all.p{identity}.read_map_stats.log"
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
map_reference: {map_reference}
MIN_mapped_reads: {params.mapped_reads_threshold}

MERGE_threads: 4
MERGE_memory: 5

MAG_omics: {mag_omics}

ANNOTATION_file:

# List annotation IDs to generate function profiles.
#
# If map_reference is 'MAG' or 'reference_genome', this list could contain elements from:
# 'eggNOG.OGs', 'eggNOG.KEGG_Pathway', 'eggNOG.KEGG_Module', 'eggNOG.KEGG_KO', 'eggNOG.PFAMs',
# 'kofam.KEGG_Pathway', 'kofam.KEGG_Module', 'kofam.KO',
# 'dbCAN.module', 'dbCAN.enzclass', 'dbCAN.subfamily', 'dbCAN.EC', 'dbCAN.eCAMI_subfamily', 'dbCAN.eCAMI_submodule'.
#
#   E.g.:
# - eggNOG.OGs
# - kofam.KEGG_Pathway
#
# If map_reference is 'gene_catalog', the names should match the ANNOTATION_file column names.

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
