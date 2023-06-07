#!/usr/bin/env python

'''
Alignment, normalization and integration step

Authors: Carmen Saenz
'''

# configuration yaml file
# import sys
from os import path
import pathlib

localrules: modify_cds_faa_header_for_fetchMG, make_merged_genome_fna, make_genome_def, merge_MG_tables,

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, local_dir, minto_dir, script_dir, metadata
include: 'config_parser.smk'

if omics == 'metaG':
    hq_dir="4-hostfree"
if omics == 'metaT':
    hq_dir="5-1-sortmerna"

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
    else:
        try:
            # Make list of illumina samples, if ILLUMINA in config
            ilmn_samples = list()
            if 'ILLUMINA' in config:
                #print("Samples:")
                for ilmn in config["ILLUMINA"]:
                    if path.exists(working_dir+'/'+omics+'/'+ hq_dir+'/'+ilmn+'/'+ilmn+'.1.fq.gz') is True:
                        #print(ilmn)
                        ilmn_samples.append(ilmn)
                    else:
                        raise TypeError('ERROR in ', config_path, ': ILLUMINA list of samples does not exist. Please, complete ', config_path)
            n_samples=len(ilmn_samples)+3
        except: 
            print('ERROR in ', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete ', config_path)
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

mag_omics = ''
if map_reference == 'MAG':
    mag_omics = 'metaG'
    if 'MAG_omics' in config and config['MAG_omics'] != None:
        mag_omics = config['MAG_omics']
    print('WARNING in ', config_path, ': MIntO is using "'+ working_dir+'/'+mag_omics+'/8-1-binning/mags_generation_pipeline/prokka" as PATH_reference variable')
else:
    if config['PATH_reference'] is None:
        print('ERROR in ', config_path, ': PATH_reference variable is empty. Please, complete ', config_path)
    elif path.exists(config['PATH_reference']) is False:
        print('ERROR in ', config_path, ': PATH_reference variable path does not exit. Please, complete ', config_path)
    else:
        if map_reference == 'reference_genome':
            print('WARNING in ', config_path, ': MIntO is using "'+ config['PATH_reference']+'" as PATH_reference variable')
            reference_dir=config["PATH_reference"]
        elif map_reference == 'genes_db':
            if path.exists(config['PATH_reference']+'/'+config['NAME_reference']) is True:
                print('WARNING in ', config_path, ': MIntO is using "'+ config['PATH_reference']+'/'+config['NAME_reference']+'" as PATH_reference and NAME_reference variables.')
                gene_catalog_db=config["PATH_reference"]
                gene_catalog_name=config["NAME_reference"]
            else:
                print('ERROR in ', config_path, ': NAME_reference variable does not exit. Please, complete ', config_path)

if config['BWAindex_threads'] is None:
    print('ERROR in ', config_path, ': BWAindex_threads variable is empty. Please, complete ', config_path)
elif type(config['BWAindex_threads']) != int:
    print('ERROR in ', config_path, ': BWAindex_threads variable is not an integer. Please, complete ', config_path)

if config['BWAindex_memory'] is None:
    print('ERROR in ', config_path, ': BWAindex_memory variable is empty. Please, complete ', config_path)
elif type(config['BWAindex_memory']) != int:
    print('ERROR in ', config_path, ': BWAindex_memory variable is not an integer. Please, complete ', config_path)

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
    reference_dir="{wd}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes".format(wd = working_dir, mag_omics = mag_omics)

if map_reference == 'reference_genome':
    post_analysis_out="refgenome-genes"
    post_analysis_dir="9-refgenome-genes-post-analysis"
    reference_dir=config["PATH_reference"]

gene_catalog_db="None"
gene_catalog_name="None"
bwaindex_db="DB/{analysis_dir}/BWA_index/{analysis_name}".format(analysis_dir=post_analysis_dir, analysis_name=post_analysis_out)

def gene_abundances_bwa_out():
    result = expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.{extension}",
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                identity = identity,
                extension = ['sorted.bam.bai', 'log', 'profile.abund.prop.txt.gz', 'profile.abund.prop.genome.txt.gz', 'profile.relabund.prop.genome.txt.gz'])
    return(result)

def gene_abundances_normalization_out():
    result = expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.{norm}.csv",
                wd = working_dir,
                omics = omics,
                post_analysis_out = post_analysis_out,
                norm = normalization_modes,
                identity = identity)
    return(result)

if map_reference == 'genes_db':
    post_analysis_dir="9-genes-db-post-analysis"
    post_analysis_out="db_genes"
    reference_dir="{wd}"
    def gene_abundances_bwa_out(): # CHECK THIS PART - do not generate bwa index, normalization in a different way
        result = expand("{gene_catalog_path}/BWA_index/{gene_catalog_name}.pac", 
                    gene_catalog_path=gene_catalog_db, 
                    gene_catalog_name=gene_catalog_name),\
        expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.bam",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db_genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity),\
        expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.log",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db_genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity),\
        expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile_TPM.txt.gz",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = "db_genes",
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
                    identity = identity)
        return(result)

def gene_abundances_map_prof():
    result = expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.txt",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = post_analysis_out,
                    identity = identity),\
    expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.maprate.txt",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = post_analysis_out,
                    identity = identity),\
    expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.mapstats.txt",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = post_analysis_out,
                    identity = identity),\
    expand("{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.multimap.txt",
                    wd = working_dir,
                    omics = omics,
                    post_analysis_out = post_analysis_out,
                    identity = identity)
    return(result)

def config_yaml():
    result = "{wd}/data_integration.yaml".format(
                wd = working_dir)
    return(result)

rule all:
    input: 
        #gene_abundances_bwa_out(),
        gene_abundances_normalization_out(),
        gene_abundances_map_prof(),
        config_yaml()


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

rule genome_bwaindex:
    input: 
        fasta_merge=rules.make_merged_genome_fna.output.fasta_merge
    output: 
        bwaindex="{wd}/DB/9-{post_analysis_out}-post-analysis/BWA_index/{post_analysis_out}.pac"
    params:
        tmp_bwaindex=lambda wildcards: "{local_dir}/{post_analysis_out}_bwaindex".format(local_dir=local_dir, post_analysis_out=post_analysis_out),
        indexprefix=post_analysis_out
    log:
        "{wd}/logs/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}_bwaindex.log"
    threads: config["BWAindex_threads"]
    resources:
        mem=config["BWAindex_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        rm -rf {params.tmp_bwaindex}
        mkdir -p {params.tmp_bwaindex}/BWA_index/
        time (\
                bwa-mem2 index {input} -p {params.tmp_bwaindex}/BWA_index/{params.indexprefix}
                rsync {params.tmp_bwaindex}/BWA_index/* $(dirname {output.bwaindex})/
                rm -rf {params.tmp_bwaindex} \
            ) &> {log}
        """

rule gene_abund_bwa_raw:
    input:
        bwaindex=rules.genome_bwaindex.output.bwaindex,
        genome_def=rules.make_genome_def.output.genome_def,
        hq_reads_fw=lambda wildcards: '{wd}/{omics}/{hq_dir}/{sample}/{sample}.1.fq.gz'.format(wd = wildcards.wd,omics=wildcards.omics, hq_dir=hq_dir, sample=wildcards.sample),
        hq_reads_rv=lambda wildcards: '{wd}/{omics}/{hq_dir}/{sample}/{sample}.2.fq.gz'.format(wd = wildcards.wd,omics=wildcards.omics, hq_dir=hq_dir, sample=wildcards.sample),
    output:
        sort="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam",
        index="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.sorted.bam.bai",
        bwa_log="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.log",
        profile_raw="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.txt.gz",
        map_profile="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz",
        map_mag_profile="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.abund.prop.genome.txt.gz",
        map_mag_profile_rel="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{identity}.filtered.profile.relabund.prop.genome.txt.gz",
    params:
        tmp_bwa=lambda wildcards: "{local_dir}/{omics}_{sample}_{post_analysis_out}_bwa_raw/".format(local_dir=local_dir, omics=wildcards.omics, sample=wildcards.sample, post_analysis_out=wildcards.post_analysis_out),
        length=config["msamtools_filter_length"],
        prefix="{sample}.p{identity}.filtered",
        memory=lambda wildcards, resources: resources.mem - 1,
        mapped_reads_threshold=config["MIN_mapped_reads"]
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}.p{identity}_bwa.log"
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" # BWA + samtools + msamtools + perl
    shell:
        """
        rm -rf {params.tmp_bwa}
        mkdir -p {params.tmp_bwa}
        bwaindex_dir=$(dirname {input.bwaindex})
        remote_dir=$(dirname {output.sort})
        (time (bwa-mem2 mem -a -t {threads} -v 3 ${{bwaindex_dir}}/{post_analysis_out} {input.hq_reads_fw} {input.hq_reads_rv}| \
                msamtools filter -S -b -l {params.length} -p {identity} -z 80 --besthit - > {params.tmp_bwa}{params.prefix}.bam) >& {params.tmp_bwa}{params.prefix}.log
        samtools sort {params.tmp_bwa}{params.prefix}.bam -o {params.tmp_bwa}{params.prefix}.sorted.bam -@ {threads} -m {params.memory}G --output-fmt=BAM 
        samtools index {params.tmp_bwa}{params.prefix}.sorted.bam {params.tmp_bwa}{params.prefix}.sorted.bam.bai -@ {threads}
        total_reads="$(grep Processed {params.tmp_bwa}{params.prefix}.log | perl -ne 'm/Processed (\\d+) reads/; $sum+=$1; END{{printf "%d\\n", $sum/2;}}')"
        echo $total_reads
        common_args="--label {wildcards.omics}.{wildcards.sample} --total=$total_reads --mincount={params.mapped_reads_threshold}"
        msamtools profile {params.tmp_bwa}{params.prefix}.bam $common_args -o {params.tmp_bwa}{params.prefix}.profile.abund.prop.txt.gz --multi=prop --unit=abund --nolen
        msamtools profile {params.tmp_bwa}{params.prefix}.bam $common_args -o {params.tmp_bwa}{params.prefix}.profile.abund.all.txt.gz  --multi=all  --unit=abund --nolen
        msamtools profile {params.tmp_bwa}{params.prefix}.bam $common_args -o {params.tmp_bwa}{params.prefix}.profile.abund.prop.genome.txt.gz --multi=prop --unit=abund --nolen --genome {input.genome_def}
        msamtools profile {params.tmp_bwa}{params.prefix}.bam $common_args -o {params.tmp_bwa}{params.prefix}.profile.relabund.prop.genome.txt.gz --multi=prop --unit=rel --genome {input.genome_def}
        rm {params.tmp_bwa}{params.prefix}.bam
        
        rsync {params.tmp_bwa}* ${{remote_dir}}
        rm -rf {params.tmp_bwa}) >& {log}
        """

###############################################################################################
# Prepare genes for mapping to gene-database
## First, the index has to be generated
## Mapping, computation of read counts and TPM normalization is done in the same rule
## TPM normalization: sequence depth and genes’ length 
###############################################################################################

rule gene_abund_bwaindex_gene_catalog:
    input:
        genes="{gene_catalog_path}/{gene_catalog_name}"
    output: 
        bwaindex="{gene_catalog_path}/BWA_index/{gene_catalog_name}.pac"
    params:
        tmp_bwaindex=lambda wildcards: "{local_dir}/{gene_catalog_name}_bwaindex/".format(local_dir=local_dir, gene_catalog_name=wildcards.gene_catalog_name),
    #log:
    #    "{wd}/logs/DB/{post_analysis_dir}/{gene_catalog_name}_bwaindex.log".format(wd = working_dir, post_analysis_dir=post_analysis_dir, gene_catalog_name=gene_catalog_name)
    threads: config["BWAindex_threads"]
    resources:
        mem=config["BWAindex_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        rm -rf {params.tmp_bwaindex}
        mkdir -p {params.tmp_bwaindex}
        time (bwa-mem2 index {wildcards.gene_catalog_path}/{wildcards.gene_catalog_name} -p {params.tmp_bwaindex}{wildcards.gene_catalog_name}
        rsync {params.tmp_bwaindex}* {wildcards.gene_catalog_path}/BWA_index/
        rm -rf {params.tmp_bwaindex}) &> {log}
        """
        
rule gene_abund_bwa_tpm:
    input:
        bwaindex="{gene_catalog_path}/BWA_index/{gene_catalog_name}.pac".format(gene_catalog_path=gene_catalog_db, gene_catalog_name=gene_catalog_name),
        hq_reads_fw=lambda wildcards: '{wd}/{omics}/{hq_dir}/{sample}/{sample}.1.fq.gz'.format(wd=wildcards.wd, omics=wildcards.omics,hq_dir=hq_dir, sample=wildcards.sample),
        hq_reads_rv=lambda wildcards: '{wd}/{omics}/{hq_dir}/{sample}/{sample}.2.fq.gz'.format(wd=wildcards.wd, omics=wildcards.omics,hq_dir=hq_dir, sample=wildcards.sample),
    output:
        filter="{wd}/{omics}/9-mapping-profiles/BWA_reads-db_genes/{sample}/{sample}.p{identity}.filtered.bam",
        bwa_log="{wd}/{omics}/9-mapping-profiles/BWA_reads-db_genes/{sample}/{sample}.p{identity}.filtered.log",
        profile_tpm="{wd}/{omics}/9-mapping-profiles/BWA_reads-db_genes/{sample}/{sample}.p{identity}.filtered.profile_TPM.txt.gz",
        map_profile="{wd}/{omics}/9-mapping-profiles/BWA_reads-db_genes/{sample}/{sample}.p{identity}.filtered.profile.abund.all.txt.gz"
    params:
        #extra_params=config["BWAalignm_parameters"],
        tmp_bwa=lambda wildcards: "{local_dir}/{omics}_{sample}_db_genes_bwa_tpm/".format(local_dir=local_dir, omics=omics, sample=wildcards.sample),
        length=config["msamtools_filter_length"],
        prefix="{sample}.p{identity}.filtered",
        memory=lambda wildcards, resources: resources.mem - 1,
        bwaindex="{gene_catalog_path}/BWA_index/{gene_catalog_name}".format(gene_catalog_path=gene_catalog_db, gene_catalog_name=gene_catalog_name),
        mapped_reads_threshold=config["MIN_mapped_reads"]
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-db_genes/{sample}.p{identity}_bwa.log"
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #config["conda_env2_yml"] #BWA + samtools
    shell:
        """
        rm -rf {params.tmp_bwa}
        mkdir -p {params.tmp_bwa}
        bwaindex_dir=$(dirname {input.bwaindex})
        remote_dir=$(dirname {output.filter})
        (time (bwa-mem2 mem -a -t {threads} -v 3 {params.bwaindex} {input.hq_reads_fw} {input.hq_reads_rv}| \
msamtools filter -S -b -l {params.length} -p {wildcards.identity} -z 80 --besthit - > {params.tmp_bwa}{params.prefix}.bam) >& {params.tmp_bwa}{params.prefix}.log
        total_reads="$(grep Processed {params.tmp_bwa}{params.prefix}.log | perl -ne 'm/Processed (\\d+) reads/; $sum+=$1; END{{printf "%d\\n", $sum/2;}}')"
        echo $total_reads
        msamtools profile {params.tmp_bwa}{params.prefix}.bam --label {wildcards.omics}.{wildcards.sample} -o {params.tmp_bwa}{params.prefix}.profile_TPM.txt.gz --total $total_reads --mincount {params.mapped_reads_threshold} --multi prop --unit tpm
        msamtools profile {params.tmp_bwa}{params.prefix}.bam --label {wildcards.omics}.{wildcards.sample} -o {params.tmp_bwa}{params.prefix}.profile.abund.all.txt.gz --total $total_reads --mincount {params.mapped_reads_threshold} --multi all --unit abund --nolen
        rsync {params.tmp_bwa}* ${{remote_dir}} ) >& {log}
        rm -rf {params.tmp_bwa}
        """

###############################################################################################
# Computation of read counts to genes belonging to MAGs or publicly available genomes
###############################################################################################

rule gene_abund_compute:
    input: 
        sort=expand("{{wd}}/{{omics}}/9-mapping-profiles/BWA_reads-{{post_analysis_out}}/{sample}/{sample}.p{{identity}}.filtered.sorted.bam", sample=ilmn_samples),
        index=expand("{{wd}}/{{omics}}/9-mapping-profiles/BWA_reads-{{post_analysis_out}}/{sample}/{sample}.p{{identity}}.filtered.sorted.bam.bai", sample=ilmn_samples),
        bed_subset="{wd}/DB/9-{post_analysis_out}-post-analysis/{post_analysis_out}_SUBSET.bed".format(wd = working_dir, post_analysis_out=post_analysis_out)
    output: 
        absolute_counts="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.bed"
    params:
        tmp_bwa=lambda wildcards: "{local_dir}/{omics}_{post_analysis_out}_abundances/".format(local_dir=local_dir, omics=omics, post_analysis_out=post_analysis_out),
        prefix="genes_abundances.p{identity}.bed"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances_counts.p{identity}.log"
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml" #bedtools
    shell:
        """
        rm -rf {params.tmp_bwa}
        mkdir -p {params.tmp_bwa}
        time (files='{ilmn_samples}'; echo ${{files}} | tr ' ' '\\t' > {params.tmp_bwa}filename_list
        echo -e 'chr\\tstart\\tstop\\tname\\tscore\\tstrand\\tsource\\tfeature\\tframe\\tinfo' > {params.tmp_bwa}column_names
        cat {params.tmp_bwa}filename_list >> {params.tmp_bwa}column_names; cat {params.tmp_bwa}column_names| tr '\\n' '\\t' > {params.tmp_bwa}column_names2
        sed 's/\t$//' {params.tmp_bwa}column_names2 >> {params.tmp_bwa}{params.prefix}; echo '' >> {params.tmp_bwa}{params.prefix}
        bedtools multicov -bams {input.sort} -bed {input.bed_subset} >> {params.tmp_bwa}{params.prefix}
        rsync {params.tmp_bwa}{params.prefix} {output.absolute_counts}) &> {log}
        rm -rf {params.tmp_bwa}
        """

###############################################################################################
# Normalization of read counts to 10 marker genes (MG normalization)
## fetchMGs identifies 10 universal single-copy phylogenetic MGs 
## (COG0012, COG0016, COG0018, COG0172, COG0215, COG0495, COG0525, COG0533, COG0541, and COG0552)
###############################################################################################

rule modify_cds_faa_header_for_fetchMG:
    input: '{something}.faa'
    output: '{something}_SUBSET.faa'
    shell:
        """
        sed 's/\\s.*//;s/\\./-/g' {input} > {output}
        """

# Run fetchMGs
rule fetchMG_genome_cds_faa:
    input: 
        cds_faa='{wd}/DB/{post_analysis_dir}/CD_transl/{genome}_SUBSET.faa',
        fetchMGs_dir=str(fetchMGs_dir)
    output: '{wd}/DB/{post_analysis_dir}/fetchMGs/{genome}/{genome}_SUBSET.all.marker_genes_scores.table'
    params:
        tmp_MG=lambda wildcards: "{local_dir}/{post_analysis_out}_marker_genes/{genome}".format(local_dir=local_dir, post_analysis_out=post_analysis_out, genome=wildcards.genome),
    log: '{wd}/logs/DB/{post_analysis_dir}/fetchMGs/{genome}.log'
    threads: 8
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell: 
        """
        rm -rf {params.tmp_MG}/{wildcards.genome}
        mkdir -p {params.tmp_MG}/
        {input.fetchMGs_dir}/fetchMGs.pl -outdir {params.tmp_MG}/{wildcards.genome} -protein_only -threads {threads} -x {input.fetchMGs_dir}/bin -m extraction {input.cds_faa} >& {log}
        rm -rf {params.tmp_MG}/{wildcards.genome}/temp
        rsync -a {params.tmp_MG}/{wildcards.genome}/* $(dirname {output})/
        rm -rf {params.tmp_MG}/{wildcards.genome}
        """

def get_genome_MG_tables(wildcards):
    #Collect the CDS faa files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    result = expand("{wd}/DB/{post_analysis_dir}/fetchMGs/{genome}/{genome}_SUBSET.all.marker_genes_scores.table",
                    wd=wildcards.wd,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    genome=genomes)
    return(result)

rule merge_MG_tables:
    input: get_genome_MG_tables
    output: "{wd}/DB/{post_analysis_dir}/all.marker_genes_scores.table"
    log: "{wd}/logs/DB/{post_analysis_dir}/merge_marker_genes_scores.table.log"
    shell:
        """
        time (\
                head -n 1 {input[0]} > {output}
                for file in {input}; do
                    awk 'FNR>1' ${{file}} >> {output}
                done\
            ) &> {log}
        """

rule gene_abund_normalization_MG:
    input:
        absolute_counts="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.bed",
        genomes_marker_genes="{wd}/DB/9-{post_analysis_out}-post-analysis/all.marker_genes_scores.table"
    output:
        norm_counts="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.MG.csv"
    params:
        mapped_reads_threshold=config["MIN_mapped_reads"]
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.MG.log"
    threads: 4
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell: 
        """
        time (Rscript {script_dir}/profile_MG.R --normalize MG --threads {threads} --memory {resources.mem} --bed {input.absolute_counts} --MG {input.genomes_marker_genes} --out {output.norm_counts} --omics {wildcards.omics} --min-read-count {params.mapped_reads_threshold}) &> {log}
        """

###############################################################################################
# Normalization of read counts by sequence depth and genes’ length (TPM normalization)
###############################################################################################

rule gene_abund_normalization_TPM:
    input:
        absolute_counts="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.bed",
    output:
        norm_counts="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.TPM.csv",
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.TPM.log"
    params:
        mapped_reads_threshold=config["MIN_mapped_reads"]
    #    tmp_norm=lambda wildcards: "{local_dir}/{omics}_{norm}_genes_abundances_normalization/".format(local_dir=local_dir, omics=omics, norm=wildcards.norm),
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml" #R
    shell:
        """
        time (Rscript {script_dir}/profile_MG.R --normalize TPM --threads {threads} --memory {resources.mem} --bed {input.absolute_counts} --out {output.norm_counts} --omics {wildcards.omics} --min-read-count {params.mapped_reads_threshold}) &> {log}
        """

###############################################################################################
# Merge normalized gene abundance or transcript profiles from gene catalog (TPM normalization)
###############################################################################################
rule gene_abund_tpm_merge:
    input:
        profile_tpm=expand("{{wd}}/{{omics}}/9-mapping-profiles/BWA_reads-{post_analysis_out}/{sample}/{sample}.p{{identity}}.filtered.profile_TPM.txt.gz", sample=ilmn_samples, post_analysis_out=post_analysis_out),
    output:
        profile_tpm_all="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.TPM.csv"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/gene_abund_merge.p{identity}.TPM_log"
    wildcard_constraints:
        post_analysis_out='db-genes'
    params:
        #profile_abund_all="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.abund.csv",
        tmp_tpm_merge=lambda wildcards: "{local_dir}/{omics}_{identity}.gene_abund_tpm_merge/".format(local_dir=local_dir, omics=omics, identity=identity),
        profile_tpm_list=lambda wildcards, input: ",".join(input.profile_tpm)
        #prefix_db="{gene_catalog_path}/{gene_catalog_name}",
        #mapped_reads_threshold=config["MIN_mapped_reads"]
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #config["conda_env_yml"] # base env
    shell:
        """ rm -rf {params.tmp_tpm_merge}
        mkdir -p {params.tmp_tpm_merge}
        prefix=$(basename {output.profile_tpm_all})
        time (sh {script_dir}/msamtools_merge_profiles.sh {input.profile_tpm[0]} '{params.profile_tpm_list}' db_genes {params.tmp_tpm_merge} ${{prefix}}
        rsync {params.tmp_tpm_merge}${{prefix}} {output.profile_tpm_all}) &> {log}
        rm -rf {params.tmp_tpm_merge}"""

###############################################################################################
# Merge map stats files
## Mappability rate
## Mappability stats
## Multimapping read count
###############################################################################################
rule gene_abund_profiling_merge:
    input: 
        map_profile=expand("{{wd}}/{{omics}}/9-mapping-profiles/BWA_reads-{{post_analysis_out}}/{sample}/{sample}.p{{identity}}.filtered.profile.abund.all.txt.gz", sample=ilmn_samples),
    output: 
        map_profile_all="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.txt",
        maprate="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.maprate.txt",
        mapstats="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.mapstats.txt",
        multimap="{wd}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.multimap.txt"
    params:
        tmp_prof_merge=lambda wildcards: "{local_dir}/{omics}_{identity}.gene_abund_profiling_merge/".format(local_dir=local_dir, omics=omics, identity=identity),
        map_profile_list=lambda wildcards, input: ",".join(input.map_profile),
        prefix="all.p{identity}.filtered.profile.abund.all"
    log:
        "{wd}/logs/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.gene_abund_profiling_merge_log"
    threads: config["BWA_threads"]
    resources:
        mem=config["BWA_memory"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """ rm -rf {params.tmp_prof_merge}
        mkdir -p {params.tmp_prof_merge}
        time (sh {script_dir}/msamtools_merge_profiles.sh {input.map_profile[0]} '{params.map_profile_list}' genome_abund {params.tmp_prof_merge} {params.prefix}.txt
        sh {script_dir}/msamtools_stats.sh '{params.map_profile_list}' {identity} {params.prefix} {params.tmp_prof_merge}
        rsync {params.tmp_prof_merge}{params.prefix}.maprate.txt {output.maprate}
        rsync {params.tmp_prof_merge}{params.prefix}.mapstats.txt {output.mapstats}
        rsync {params.tmp_prof_merge}{params.prefix}.multimap.txt {output.multimap}
        rsync {params.tmp_prof_merge}{params.prefix}.txt {output.map_profile_all}) &> {log}
        rm -rf {params.tmp_prof_merge}"""

###############################################################################################
# Generate configuration yml file for data integration - gene and function profiles
###############################################################################################
rule config_yml_integration:
    input: 
        map_profile_all=expand("{{wd}}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.txt", omics = omics, post_analysis_out=post_analysis_out, identity=identity),
        maprate=expand("{{wd}}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.maprate.txt", omics = omics, post_analysis_out=post_analysis_out, identity=identity),
        mapstats=expand("{{wd}}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.mapstats.txt", omics = omics, post_analysis_out=post_analysis_out, identity=identity),
        multimap=expand("{{wd}}/{omics}/9-mapping-profiles/BWA_reads-{post_analysis_out}/all.p{identity}.filtered.profile.abund.all.multimap.txt", omics = omics, post_analysis_out=post_analysis_out, identity=identity)
    output: 
        config_file="{wd}/data_integration.yaml"
    params: 
        tmp_integration_yaml=lambda wildcards: "{local_dir}{omics}_config_yml_integration/".format(local_dir=local_dir, omics = omics),
        mapped_reads_threshold=config["MIN_mapped_reads"],
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/config_yml_integration.log"
    shell: 
        """ rm -rf {params.tmp_integration_yaml}
        mkdir -p {params.tmp_integration_yaml}
        time (echo "######################
# General settings
######################
PROJECT: {project_id}
working_dir: {wildcards.wd}
omics: metaG_metaT
local_dir: {local_dir}
minto_dir: {minto_dir}
METADATA: {metadata}

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

# List annotation IDs matching to generate function profiles. 
# If map_reference= 'MAG' or 'reference_genome', this list correspond to:
# 'eggNOG_OGs','KEGG_Pathway','KEGG_Module','KEGG_KO','PFAMs','dbCAN.mod' and 'dbCAN.enzclass.
# The names should match the ANNOTATION_file column names.
#   E.g.:
# - eggNOG_OGs
# - KEGG_Pathway
ANNOTATION_ids:
 - eggNOG_OGs
 - KEGG_Pathway
 - KEGG_Module
 - KEGG_KO
 - PFAMs
 - dbCAN.mod
 - dbCAN.enzclass
" > {params.tmp_integration_yaml}data_integration.yaml

rsync {params.tmp_integration_yaml}data_integration.yaml {output.config_file}) >& {log}
rm -rf {params.tmp_integration_yaml}
        """
