#!/usr/bin/env python

'''
Binning preparation of contigs from assembly/coassembly.

Authors: Carmen Saenz, Mani Arumugam
'''

import snakemake

include: 'scripts/07-common-rules.smk'
include: 'scripts/08-common-rules.smk'

# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/binning-preparation_Mani.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --use-conda --conda-prefix /emc/cbmr/users/rzv923/ibdmdb_test/tmp_porus/ --configfile assemblies.smk.yaml --jobs 10
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/binning-preparation_Mani.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A"  --use-conda --conda-prefix /data/MIntO_snakemake_env/ --configfile assemblies.smk.yaml --jobs 10

localrules: filter_contigs_illumina_single, filter_contigs_illumina_coas, filter_contigs_nanopore, filter_contigs_illumina_single_nanopore, link_bam, combine_fasta, check_depths, combine_depth

# Scaffold type
SCAFFOLDS_type = list()
# Make list of illumina samples, if ILLUMINA in config
ilmn_samples = list()
if 'ILLUMINA' in config:
    SCAFFOLDS_type.append('illumina_single')
    #print("Samples:")
    for ilmn in config["ILLUMINA"]:
        #print(" "+ilmn)
        ilmn_samples.append(ilmn)
    ilmn_samples.sort()

# Make list of nanopore assemblies, if NANOPORE in config
nanopore_assemblies = list()
if 'NANOPORE' in config:
    #print("Nanopore assemblies:")
    SCAFFOLDS_type.append('nanopore')
    for nano in config["NANOPORE"]:
        #print(" "+nano)
        nanopore_assemblies.append(nano)

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config
hybrid_assemblies = list()
if 'HYBRID' in config:
    #print("Hybrid assemblies:")
    SCAFFOLDS_type.append('illumina_single_nanopore')
    for nano in config["HYBRID"]:
        for ilmn in config["HYBRID"][nano].split("+"):
            #print(" "+nano+"-"+ilmn)
            hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if COASSEMBLY in config
co_assemblies = list()
if 'COASSEMBLY' in config:
    #print("Coassemblies:")
    SCAFFOLDS_type.append('illumina_coas')
    for co in config["COASSEMBLY"]:
        #print(" "+co)
        co_assemblies.append(co)

#print(SCAFFOLDS_type)

# some variables

project_name = config['PROJECT']
local_dir = config['local_dir']
working_dir = config['working_dir']
omics = config['omics']
minto_dir=config["minto_dir"]
hybrid_max_k = config['METASPADES_hybrid_max_k']
illumina_max_k = config['METASPADES_illumina_max_k']
metadata=config["METADATA"]

# Define all the outputs needed by target 'all'

rule all:
    input: 
        fasta = "{wd}/{omics}/8-1-binning/{project}_scaffolds.2500.fasta".format(
                wd = working_dir, 
                omics = omics,
                project = project_name),
        depth = "{wd}/{omics}/8-1-binning/{project}_scaffolds.2500.depth.txt".format(
                wd = working_dir, 
                omics = omics,
                project = project_name),
        config_yaml = "{wd}/{omics}/mag_generation.yaml".format(
                wd = working_dir, 
                omics = omics)

###############################################################################################
# Filter contigs from 
#  1. MetaSPAdes-individual-assembled illumina samples,
#  2. MEGAHIT-co-assembled illumina samples, 
#  3. MetaSPAdes-hybrid-assembled illumina+nanopore data,
#  4. MetaFlye-assembled nanopore sequences
###############################################################################################

rule filter_contigs_nanopore:
    input: 
        #lambda wildcards: expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.polcirc.fasta", 
        lambda wildcards: expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.fasta", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                sample = nanopore_assemblies,
                assembly_preset = config['METAFLYE_presets'])
    output: 
        "{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta"
    wildcard_constraints:
        scaffolds_type = 'nanopore'
    resources:
        mem = 5
    params:
        min_length = lambda wildcards: wildcards.min_length
    run:
        filter_fasta_list_by_length(input, output[0], params.min_length)

rule filter_contigs_illumina_single_nanopore:
    input: 
        #lambda wildcards: expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.assembly.polcirc.fasta", 
        lambda wildcards: expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.scaffolds.fasta", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                assembly = hybrid_assemblies,
                kmer_dir = "k21-" + str(hybrid_max_k))
    output: 
        "{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta"
    wildcard_constraints:
        scaffolds_type = 'illumina_single_nanopore'
    resources:
        mem = 5
    params:
        min_length = lambda wildcards: wildcards.min_length
    run:
        filter_fasta_list_by_length(input, output[0], params.min_length)

rule filter_contigs_illumina_single:
    input: 
        lambda wildcards: expand("{wd}/{omics}/7-assembly/{illumina}/{kmer_dir}/{illumina}.scaffolds.fasta", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                illumina = ilmn_samples,
                kmer_dir = "k21-" + str(illumina_max_k))
    output: 
        "{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta"
    wildcard_constraints:
        scaffolds_type = 'illumina_single'
    resources:
        mem = 5
    params:
        min_length = lambda wildcards: wildcards.min_length
    run:
        filter_fasta_list_by_length(input, output[0], params.min_length)

rule filter_contigs_illumina_coas:
    input: 
        lambda wildcards: expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta",
                wd = wildcards.wd,
                omics = wildcards.omics,
                coassembly = co_assemblies,
                assembly_preset = config['MEGAHIT_presets'])
    output: 
        "{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta"
    wildcard_constraints:
        scaffolds_type = 'illumina_coas'
    resources:
        mem = 5
    params:
        min_length = lambda wildcards: wildcards.min_length
    run:
        filter_fasta_list_by_length(input, output[0], params.min_length)

################################################################################################
# Create BWA index for the filtered contigs
################################################################################################
rule BWA_index_contigs:
    input:
        scaffolds="{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta"
    output:
        bwaindex="{wd}/{omics}/6-mapping/BWA_index/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta.bwt.2bit.64"
    log:
        "{wd}/{omics}/logs/6-mapping/{project}_scaffolds_{scaffolds_type}.{min_length}.BWA_index.log"
    params:
        local_loc = lambda wildcards: "{local_dir}/{omics}/6-mapping/BWA_index".format(local_dir=local_dir, omics=wildcards.omics)
    resources:
        mem = 450
    threads: 4
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        outdir=$(dirname {output})
        mkdir -p $outdir {params.local_loc}
        time (bwa-mem2 index {input.scaffolds} -p {params.local_loc}/$(basename {input})) >& {log}
        ln -s {params.local_loc}/$(basename {input}).* $outdir/
        """

# Maps reads to contig-set using bwa2 
# Baseline memory usage is:
#   Main: 240GB
#   Sort: using config values: SAMTOOLS_sort_threads and SAMTOOLS_sort_memory_gb; 
#         e.g. 4 threads and 20GB = 80GB
# If an attempt fails, 40GB extra for every new attempt
rule map_contigs_BWA:
    input: 
        bwa_index = "{wd}/{omics}/6-mapping/BWA_index/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta.bwt.2bit.64",
        fwd='{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz',
        rev='{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz'
    output:
        bam='{wd}/{omics}/6-mapping/{illumina}/{illumina}.{project}_scaffolds_{scaffolds_type}.{min_length}.fasta.sorted.bam'
    params:
        map_threads = config['BWA_threads'],
        sort_threads = config['SAMTOOLS_sort_threads'],
        sort_mem = config['SAMTOOLS_sort_memory_gb'],
        local_loc = lambda wildcards: "{local_dir}/{omics}/6-mapping/{illumina}".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina)
    log:
        '{wd}/{omics}/6-mapping/{illumina}/{illumina}.{project}_scaffolds_{scaffolds_type}.{min_length}.bwa2.log'
    resources:
        mem = lambda wildcards, attempt: 240 + config['SAMTOOLS_sort_threads']*config['SAMTOOLS_sort_memory_gb'] + 40*attempt
    threads: 
        config['BWA_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
    shell:
        """
        mkdir -p $(dirname {output}) {params.local_loc}
        local_file={params.local_loc}/$(basename {output.bam})
        db_name=$(echo {input.bwa_index} | sed "s/.bwt.2bit.64//")
        time (bwa-mem2 mem -a -t {params.map_threads} $db_name {input.fwd} {input.rev} \
                | msamtools filter -buS -p 95 -l 45 - \
                | samtools sort -m {params.sort_mem}G --threads {params.sort_threads} - \
                > $local_file
        ) >& {log}
        ln --force -s $local_file $(dirname {output})/
        """

# symlink the bam files into a shorter name because the filename is used in the header of the output depth file.
# Since we map to different contig-sets, and the names of the bam files are thus different, the headers will not
# look alike for the different contig-sets. We need to ensure that the samples in columns are in the right order 
# in the depth files (see rule 'check_depths') and to make life easier we create symlinks with just sample names
# so that simple check for the first line in the header file will be enough. However, {project} and {min_length} 
# are added purely to propagate them as wildcards to earlier steps.
rule link_bam:
    input: 
        rules.map_contigs_BWA.output.bam
    output: 
        "{wd}/{omics}/8-1-binning/depth_{scaffolds_type}/{project}.{illumina}.{min_length}"
    shell:
        """
        mkdir -p $(dirname {output})
        ln --force -s {input} {output}
        """

# Run jgi contig depth script for each contig-set
rule contigs_depth:
    input: 
        bamlinks = lambda wildcards: expand("{wd}/{omics}/8-1-binning/depth_{scaffolds_type}/{project}.{illumina}.{min_length}", 
                                            wd = wildcards.wd,
                                            omics = wildcards.omics,
                                            scaffolds_type = wildcards.scaffolds_type,
                                            project = wildcards.project,
                                            min_length = wildcards.min_length,
                                            illumina=ilmn_samples)
    output: 
        depths="{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.depth.txt",
    params:
        #samples = lambda wildcards: " ".join(["{proj}.{i}.{len}".format(proj=wildcards.project, i=i, len=wildcards.min_length) for i in ilmn_samples])
        samples = lambda wildcards, input: [os.path.basename(i) for i in input.bamlinks]
    resources:
        mem = 450
    threads: 
        len(ilmn_samples)
    log:
        "{wd}/{omics}/8-1-binning/{project}_{scaffolds_type}.{min_length}_depth.log"
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml" #metabat2
    shell:
        """
        cd $(dirname {input[0]})
        time (jgi_summarize_bam_contig_depths --outputDepth {output} {params.samples}) >& {log}
        """

# Combine multiple fasta files into one
rule combine_fasta:
    input: 
        fasta=lambda wildcards: expand("{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.fasta", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                project = wildcards.project,
                min_length = wildcards.min_length,
                scaffolds_type = SCAFFOLDS_type)
    output: 
        fasta_combined="{wd}/{omics}/8-1-binning/{project}_scaffolds.{min_length}.fasta"
    resources:
       mem = 10
    threads: 
        1
    shell:
        """
        cat {input} > {output}
        """

# Sanity check to ensure that the order of sample-depths in the depth file is the same. Otherwise binning will be wrong!
rule check_depths:
    input: 
        depth=lambda wildcards: expand("{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaffolds_type}.{min_length}.depth.txt", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                project = wildcards.project,
                min_length = wildcards.min_length,
                scaffolds_type = SCAFFOLDS_type)
    output:
        ok="{wd}/{omics}/8-1-binning/{project}_scaffolds.{min_length}.ok"
    shell:
        """
        rm --force {output}
        uniq_headers=$(for file in {input.depth}; do head -1 $file; done | sort -u | wc -l)
        if [ "$uniq_headers" == "1" ]; then
            touch {output}
        else
            # Headers were not unique. So report and set exit code for this rule to non-zero
            >&2 echo "Headers in depth files were not unique - please check depth files"
            test "1" == "2"
        fi
        """

# Combine multiple depth files into one
rule combine_depth:
    input: 
        depth=lambda wildcards: expand("{wd}/{omics}/8-1-binning/{project}_scaffolds_{scaf_type}.{min_length}.depth.txt", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                project = wildcards.project,
                min_length = wildcards.min_length,
                scaf_type = SCAFFOLDS_type),
        depth_ok = "{wd}/{omics}/8-1-binning/{project}_scaffolds.{min_length}.ok"
    output: 
        depth_combined="{wd}/{omics}/8-1-binning/{project}_scaffolds.{min_length}.depth.txt"
    resources:
       mem = 10
    threads: 
        1
    shell:
        """
        head -1 {input.depth[0]} > {output}
        (for file in {input.depth}; do tail -n +2 $file; done) >> {output}
        """

###############################################################################################
# Generate configuration yml file for recovery of MAGs and taxonomic annotation step - binning 
###############################################################################################

rule config_yml_binning:
    input: 
        #depth_combined="{{wd}}/{{omics}}/8-1-binning/{project}_scaffolds.{min_length}.depth.txt".format(project=project_name, min_length = min_length),
        #depth_combined=lambda wildcards: expand("{{wd}}/{{omics}}/8-1-binning/{project}_scaffolds.{min_length}.depth.txt", 
        #        project = project_name,
        #        min_length = wildcards.min_length)
        depth=lambda wildcards: expand("{wd}/{omics}/8-1-binning/{project}_scaffolds.2500.depth.txt", 
                wd = wildcards.wd,
                omics = wildcards.omics,
                project = project_name),
    output: 
        config_file="{wd}/{omics}/mag_generation.yaml",
    params: 
        tmp_binning_yaml=lambda wildcards: "{local_dir}/{omics}_config_yml_binning/".format(local_dir=local_dir, omics = omics),
    resources:
        mem=2
    threads: 2
    log: 
        "{wd}/logs/{omics}/config_yml_mag_generation.log"
    shell: 
        """
        mkdir -p {params.tmp_binning_yaml}
        time (echo "######################
# General settings
######################
PROJECT: {project_name}
working_dir: {wildcards.wd}
omics: {wildcards.omics}
local_dir: {local_dir}
minto_dir: {minto_dir}
METADATA: {metadata}

######################
# Program settings
######################
# COMMON PARAMETERS
#
MIN_FASTA_LENGTH: "500000"

# VAMB settings
#
BINNERS:
- vamb_256
- vamb_384
- vamb_512
- vamb_768

VAMB_THREADS:
VAMB_memory:

# Use GPU in VAMB:
# could be "yes" or "not"
VAMB_GPU: yes


# CHECKM settings
#
CHECKM_THREADS: 
CHECKM_memory: 

# checkm threshold
CHECKM_COMPLETENESS: 90  #higher than this
CHECKM_CONTAMINATION: 5  # lower than this

# Clean up checkm files?
CLEAN_CHECKM: yes # can be yes or not

# COVERM settings
#
COVERM_THREADS: 
COVERM_memory: 

# SCORING THE BEST GENOMES settings
#
# this could be checkm or genome
SCORE_METHOD: "checkm" 


# PROKKA settings
#
RUN_PROKKA: yes
PROKKA_CPUS: 
PROKKA_memory: 

# PHYLOPHLAN METAGENOMICS settings
#
RUN_TAXONOMY: yes
TAXONOMY_DATABASE: SGB.Jan20
TAXONOMY_CPUS: 
TAXONOMY_memory: 
DATABASE_FOLDER:" > {params.tmp_binning_yaml}mag_generation.yaml

rsync {params.tmp_binning_yaml}mag_generation.yaml {output.config_file}) >& {log}
rm -rf {params.tmp_binning_yaml}
        """