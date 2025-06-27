#!/usr/bin/env python

'''
Assembling metagenomes from combinations of illumina/bgi-seq and nanopore sequencing data.

Authors: Carmen Saenz, Mani Arumugam
'''

import os.path

localrules: merge_runs, mark_circular_metaspades_contigs, mark_circular_flye_contigs, rename_megahit_contigs

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'
include: 'include/locations.smk'
include: 'include/fasta_bam_helpers.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule assembly_base from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# Variables from configuration yaml file

ilmn_samples            = list()
nanopore_samples        = list()
merged_illumina_samples = dict()

##############################################
# Register composite samples
##############################################

# Make list of illumina coassemblies, if MERGE_ILLUMINA_SAMPLES in config
if (x := validate_optional_key(config, 'MERGE_ILLUMINA_SAMPLES')):
    for m in x:
        merged_illumina_samples.append(m)

##############################################
# Get sample list
##############################################

# Make list of illumina samples, if ILLUMINA in config
if (x := validate_optional_key(config, 'ILLUMINA')):
    check_input_directory(x, locations = ['6-corrected', '5-corrected-runs', get_qc2_output_location(omics)])
    ilmn_samples = x

    # If it's composite sample, then don't need to see them until it gets merged later
    for i in x:
        if i in merged_illumina_samples.keys():
            ilmn_samples.remove(i)


# Make list of nanopore samples, if NANOPORE in config
if (x := validate_optional_key(config, 'NANOPORE')):
    check_input_directory(x, locations = ['6-corrected', get_qc2_output_location(omics)])
    nanopore_samples = x

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config
hybrid_assemblies = list()
if (x := validate_optional_key(config, 'HYBRID')):
    print(f"Found HYBRID in {config_path}. Enabling hybrid assembly.")
    for nano in x:
        for ilmn in x[nano].split("+"):
            hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if enable_COASSEMBLY is set to "yes" in config
co_assemblies = dict()
if validate_optional_key(config, 'enable_COASSEMBLY'):
    if (x := validate_optional_key(config, 'COASSEMBLY')):
        print(f"Found COASSEMBLY in {config_path}. Enabling co-assembly.")
        co_assemblies = x

###############################
# MetaSPAdes parameters
###############################

METASPADES_qoffset = validate_required_key(config, 'METASPADES_qoffset')
check_allowed_values('METASPADES_qoffset', METASPADES_qoffset, ('auto', '33', '64'))

METASPADES_threads = validate_required_key(config, 'METASPADES_threads')
METASPADES_memory  = validate_required_key(config, 'METASPADES_memory')

METASPADES_illumina_max_k = validate_required_key(config, 'METASPADES_illumina_max_k')
check_number_is_odd('METASPADES_illumina_max_k', METASPADES_illumina_max_k)

METASPADES_hybrid_max_k = validate_optional_key(config, 'METASPADES_hybrid_max_k')
if METASPADES_hybrid_max_k is not None:
    check_number_is_odd('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k)

# Figure out SPAdes build
# and verify kmer max

spades_script = 'spades.py' # from conda environment
if (x := validate_optional_key(config, 'METASPADES_custom_build')):
    spades_script = x # custom built, e.g. for higher K
    check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 300)
    if METASPADES_hybrid_max_k is not None:
        check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 300)
else:
    check_number_is_between('METASPADES_illumina_max_k', METASPADES_illumina_max_k, 19, 128)
    if METASPADES_hybrid_max_k is not None:
        check_number_is_between('METASPADES_hybrid_max_k', METASPADES_hybrid_max_k, 19, 128)

MEGAHIT_threads = validate_required_key(config, 'MEGAHIT_threads')
MEGAHIT_memory  = validate_required_key(config, 'MEGAHIT_memory')

###############################
# MEGAHIT parameter sets
###############################

MEGAHIT_presets = list()
mega_k_list = []

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

if co_assemblies and not MEGAHIT_presets:
    raise Exception(f"ERROR in {config_path}: MEGAHIT_presets list of MEGAHIT parameters to run per co-assembly is empty")

###############################
# MetaFlye parameter sets
###############################

METAFLYE_presets = dict()

# Check for METAFLYE_presets
if (x := validate_optional_key(config, 'METAFLYE_presets')):
    METAFLYE_presets = x

if nanopore_samples and not METAFLYE_presets:
    raise Exception(f"ERROR in {config_path}: METAFLYE_presets variable listing MetaFlye parameters is missing")

# Define all the outputs needed by target 'all'

def illumina_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{kmer_dir}/{sample}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = ilmn_samples,
                    kmer_dir = "k21-" + str(METASPADES_illumina_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

def illumina_co_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta",
                    wd = working_dir,
                    omics = omics,
                    coassembly = co_assemblies,
                    assembly_preset = MEGAHIT_presets)
    return(result)

def nanopore_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = nanopore_samples,
                    assembly_preset = METAFLYE_presets)
    return(result)

def hybrid_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    assembly = hybrid_assemblies,
                    kmer_dir = "k21-" + str(METASPADES_hybrid_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

# Figure out kmer details

def get_metaspades_kmer_option(kk):
    kmers = [21] # start at 21
    kmers.extend(list(range(33, kk, 22))) # add more kmers with 22 increment while < max_k
    kmers.extend([kk]) # add max_k
    kmer_option = ','.join([str(k) for k in kmers])
    return(kmer_option)

rule all:
    input:
        illumina_single_assembly_output(),
        illumina_co_assembly_output(),
        nanopore_single_assembly_output(),
        hybrid_assembly_output(),
        print_versions.get_version_output(snakefile_name)
    default_target: True

###############################################################################################
# Correct Illumina reads using SPAdes' spadeshammer
###############################################################################################

def get_hq_fastq_files(wildcards):
    return(expand("{wd}/{omics}/{location}/{illumina}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location=get_qc2_output_location(wildcards.omics),
                illumina = wildcards.illumina,
                run = wildcards.run,
                pair = [1, 2]))

rule correct_spadeshammer:
    input:
        reads=get_hq_fastq_files
    output:
        fwd=temp("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.1.fq.gz"),
        rev=temp("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.2.fq.gz"),
    shadow:
        "minimal"
    params:
        qoffset = METASPADES_qoffset
    resources:
        mem = lambda wildcards, attempt: attempt*METASPADES_memory
    log:
        "{wd}/logs/{omics}/5-corrected-runs/{illumina}/{run}_spadeshammer.log"
    threads:
        METASPADES_threads
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        mkdir -p $(dirname {output.fwd})
        time (
            {spades_script} --only-error-correction -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -m {resources.mem} -o {wildcards.run} --phred-offset {params.qoffset}
            rsync -a {wildcards.run}/corrected/{wildcards.run}.1.fq.00.0_0.cor.fastq.gz {output.fwd}; rsync -a {wildcards.run}/corrected/{wildcards.run}.2.fq.00.0_0.cor.fastq.gz {output.rev}
        ) >& {log}
        """

def get_corrected_runs_for_sample(wildcards):
    if wildcards.illumina in merged_illumina_samples:
        files = list()
        for r in merged_illumina_samples[wildcards.illumina].split('+'):
            rep = r.strip()
            files += expand("{wd}/{omics}/5-corrected-runs/{rep}/{run}.{pair}.fq.gz",
                                    wd = wildcards.wd,
                                    omics = wildcards.omics,
                                    rep = rep,
                                    run = get_runs_for_sample(wildcards.wd, wildcards.omics, rep),
                                    pair = wildcards.pair)
    else:
        files = expand("{wd}/{omics}/5-corrected-runs/{illumina}/{run}.{pair}.fq.gz",
                                    wd = wildcards.wd,
                                    omics = wildcards.omics,
                                    illumina = wildcards.illumina,
                                    run = get_runs_for_sample(wildcards.wd, wildcards.omics, wildcards.illumina),
                                    pair = wildcards.pair)
    return(files)

rule merge_runs:
    input:
        files=get_corrected_runs_for_sample
    output:
        combined="{wd}/{omics}/6-corrected/{illumina}/{illumina}.{pair}.fq.gz"
    wildcard_constraints:
        pair="1|2"
    shadow:
        "minimal"
    params:
        num_runs = lambda wildcards, input: len(input.files)
    shell:
        """
        if (( {params.num_runs} > 1 )); then
            cat {input} > combined.fq.gz
            rsync -a combined.fq.gz {output}
        else
            cd $(dirname {output})
            ln --force {input[0]} {output}
        fi
        """

ruleorder: hybrid_assembly_metaspades > illumina_assembly_metaspades

###############################################################################################
########  Individual assembly of illumina samples
###############################################################################################
rule illumina_assembly_metaspades:
    input:
        fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
        rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
    output:
        cont_fa    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/contigs.fasta",
        scaf_fa    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/scaffolds.fasta",
        scaf_gfa   = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/assembly_graph_with_scaffolds.gfa.gz",
        asm_log    = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/spades.log",
        asm_params = "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/params.txt",
    shadow:
        "minimal"
    params:
        qoffset  = METASPADES_qoffset,
        asm_mode = "--meta",
        kmer_option = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
    resources:
        mem = lambda wildcards, attempt: attempt*METASPADES_memory
    log:
        "{wd}/logs/{omics}/7-assembly/{illumina}/k21-{maxk}/{illumina}_metaspades.log"
    threads:
        METASPADES_threads
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        time (
            {spades_script} {params.asm_mode} --only-assembler -1 {input.fwd} -2 {input.rev} -t {threads} -m {resources.mem} -o outdir --tmp-dir tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a outdir/contigs.fasta {output.cont_fa}
            rsync -a outdir/scaffolds.fasta {output.scaf_fa}
            rsync -a outdir/spades.log {output.asm_log}
            rsync -a outdir/params.txt {output.asm_params}
            gzip -c outdir/assembly_graph_with_scaffolds.gfa > {output.scaf_gfa}
        ) >& {log}
        """

###############################################################################################
# Hybrid assembly of illumina and nanopore samples
#
# This is almost a replica of the rule above except that nanopore is added and assembly name is different.
# When I learn how to handle that in single rule, perhaps I can unify it.
# Until then, changes have to be made in both!
###############################################################################################
rule hybrid_assembly_metaspades:
    input:
        fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
        rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
        ont="{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz"
    output:
        cont_fa    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/contigs.fasta",
        scaf_fa    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/scaffolds.fasta",
        scaf_gfa   = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/assembly_graph_with_scaffolds.gfa.gz",
        asm_log    = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/spades.log",
        asm_params = "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/params.txt",
    shadow:
        "minimal"
    params:
        qoffset  = METASPADES_qoffset,
        asm_mode = "--meta",
        kmer_option = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
    resources:
        mem = lambda wildcards, attempt: attempt*METASPADES_memory
    log:
        "{wd}/logs/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/{nanopore}-{illumina}_metaspades.log"
    threads:
        METASPADES_threads
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        time (
            {spades_script} {params.asm_mode} --only-assembler -1 {input.fwd} -2 {input.rev} --nanopore {input.ont} -t {threads} -m {resources.mem} -o outdir --tmp-dir tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a outdir/contigs.fasta {output.cont_fa}
            rsync -a outdir/scaffolds.fasta {output.scaf_fa}
            rsync -a outdir/spades.log {output.asm_log}
            rsync -a outdir/params.txt {output.asm_params}
            gzip -c outdir/assembly_graph_with_scaffolds.gfa > {output.scaf_gfa}
        ) >& {log}
        """

###############################################################################################
########  Co-assembly
# This starts with 11G per sample, but if that fails, it increases by 5G per sample per repeated attempt
###############################################################################################

def get_megahit_parameters(wildcards, kk, k_list):
    if (wildcards.assembly_preset == 'rnaspades'):
        # We assume kk is about 1/2 max-length - see rnaSPAdes recommendation
        # We will then make kmers=[0.33max, ..., 0.5max] with step 22
        min_k = 2*int(kk/3)+1 # start at odd number near 1/3 max-length
        kmers = list(range(min_k, kk-1, 22)) # add more kmers with 22 increment while < max_k
        kmers.extend([kk])
        kmer_option = ','.join([str(k) for k in kmers])
        return "--k-list {}".format(kmer_option)
    elif (wildcards.assembly_preset.startswith('meta-custom')):
        k_list_i = int(wildcards.assembly_preset.rsplit("-", 1)[-1])
        return "--k-list {}".format(k_list[k_list_i - 1])
    elif (wildcards.assembly_preset in ['meta-large', 'meta-sensitive']):
        return "--presets {}".format(wildcards.assembly_preset)
    else:
        return ""

rule coassembly_megahit:
    input:
        fwd=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz', wd=wildcards.wd, omics=wildcards.omics, illumina=co_assemblies[wildcards.coassembly].split('+')),
        rev=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz', wd=wildcards.wd, omics=wildcards.omics, illumina=co_assemblies[wildcards.coassembly].split('+'))
    output:
        cont_fa    = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa",
        asm_log    = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/log.txt",
        asm_params = "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/options.json",
    shadow:
        "minimal"
    params:
        fwd_reads=lambda wildcards, input: ",".join(input.fwd),
        rev_reads=lambda wildcards, input: ",".join(input.rev),
        asm_params=lambda wildcards: get_megahit_parameters(wildcards, METASPADES_illumina_max_k, mega_k_list),
    resources:
        mem       = lambda wildcards, input, attempt: min(900, len(input.fwd)*(MEGAHIT_memory+6*(attempt-1))),
        mem_bytes = lambda wildcards, input, attempt: min(900, len(input.fwd)*(MEGAHIT_memory+6*(attempt-1)))*1024*1024*1024
    log:
        "{wd}/logs/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}_{assembly_preset}_coassembly_megahit.log"
    threads:
        MEGAHIT_threads
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        # Don't create the --out-dir directory as MEGAHIT wants it to not exist before
        time (
            megahit -1 {params.fwd_reads} -2 {params.rev_reads} -t {threads} -m {resources.mem_bytes} --out-dir outdir {params.asm_params}
        ) >& {log}
        rsync -a outdir/final.contigs.fa {output.cont_fa}
        rsync -a outdir/log {output.asm_log}
        rsync -a outdir/options.json {output.asm_params}
        """

###############################################################################################
# Assemble nanopore-raw sequences using flye 2.8.3
###############################################################################################
rule nanopore_assembly_metaflye:
    input:
        ont="{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz"
    output:
        "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly.fasta",
        "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_info.txt"
    log:
        "{wd}/logs/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}_{assembly_preset}_metaflye.log"
    resources:
        mem = lambda wildcards, attempt: 30*attempt
    threads: 16
    params:
        options = lambda wildcards: METAFLYE_presets[wildcards.assembly_preset]
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        """
        mkdir -p $(dirname {output[0]})
        time (
            flye --nano-raw {input} --out-dir $(dirname {output[0]}) --threads {threads} --meta {params.options}
        ) >& {log}
        """

###############################################################################################
# This checks whether the first and last k characters in the contig are the same.
# If yes, it removes the terminal k characters and appends '_circularA' to the fasta header.
# It also updates 'length' in the fasta header.
###############################################################################################
rule mark_circular_metaspades_contigs:
    input:
        "{wd}/{omics}/7-assembly/{sample}/k21-{maxk}/{which}.fasta"
    output:
        "{wd}/{omics}/7-assembly/{sample}/k21-{maxk}/{sample}.{which}.fasta"
    params:
        kmer = lambda wildcards: int(wildcards.maxk)
    run:
        import re
        regex = re.compile('length_[0-9]+_')
        # Open the output file
        with open(output[0], 'w') as out:
            # Go through the fasta file
            fiter = fasta_iter(input[0])
            for entry in fiter:
                header, seq = entry
                if len(seq) > params.kmer:
                    if seq[0:params.kmer] == seq[-params.kmer:]:
                        seq = seq[0:-params.kmer]
                        header = regex.sub("length_%s_" % len(seq), header) + '_circularA'
                out.write(">MetaSPAdes.k21-{maxk}.{sample}_{header}\n".format(maxk=wildcards.maxk, sample=wildcards.sample, header=header))
                out.write(seq+"\n")

###############################################################################################
# This rule identifies contigs marked circular in assembly_info.txt.
# It appends '_circularA' to the fasta header.
###############################################################################################
rule mark_circular_flye_contigs:
    input:
        "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly.fasta",
        "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/assembly_info.txt"
    output:
        "{wd}/{omics}/7-assembly/{nanopore}/{assembly_preset}/{nanopore}.assembly.fasta",
    run:

        # Go through the assembly_info file
        circular = {}
        for line in open(input[1]).readlines():
            words = line.split('\t')
            if (words[3] == "Y"):
                circular[words[0]] = 1

        # Open the output file
        with open(output[0], 'w') as out:
            # Go through the fasta file
            fiter = fasta_iter(input[0])
            for entry in fiter:
                header, seq = entry
                new_header = header.replace("contig", "NODE") + '_length_%s' % len(seq)
                if header in circular:
                    new_header = new_header+'_circularA'
                out.write(f'>MetaFlye.{wildcards.assembly_preset}.{wildcards.nanopore}_{new_header}\n')
                out.write(seq+"\n")

# Rename MEGAHIT contigs
rule rename_megahit_contigs:
    input:
        "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa"
    output:
        "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta"
    wildcard_constraints:
        assembly_preset = '|'.join(MEGAHIT_presets)
    conda:
        minto_dir + "/envs/MIntO_base.yml"
    shell:
        r"""
        perl -ne 's/^>k(\d+)_(\d+) (.*)len=(\d+)/>MEGAHIT.{wildcards.assembly_preset}.{wildcards.coassembly}_NODE_$2_length_$4_k_$1/ if m/^>/; print $_;' < {input} > {output}
        """
