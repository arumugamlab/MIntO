#!/usr/bin/env python

'''
Assembling metagenomes from combinations of illumina/bgi-seq and nanopore sequencing data.

Authors: Carmen Saenz, Mani Arumugam
'''

localrules: merge_runs, mark_circular_metaspades_contigs, mark_circular_flye_contigs, rename_megahit_contigs

#
# configuration yaml file
# import sys
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
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

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
    else:
        # Make list of illumina samples, if ILLUMINA in config
        ilmn_samples = list()
        #print("Samples:")
        for ilmn in config['ILLUMINA']:
            x = str(ilmn)
            if path.exists("{}/{}/{}/{}".format(working_dir, omics, '6-corrected', x)) is True:
                ilmn_samples.append(x)
            elif path.exists("{}/{}/{}/{}".format(working_dir, omics, '5-corrected-runs', x)) is True:
                ilmn_samples.append(x)
            elif path.exists("{}/{}/{}/{}".format(working_dir, omics, get_qc2_output_location(omics), x)) is True:
                ilmn_samples.append(x)
            else:
                raise Exception("ERROR in {}: ILLUMINA sequence does not exist for sample {}".format(config_path, x))
else:
    print('ERROR in', config_path, ': ILLUMINA list of samples is empty. Please, complete', config_path)

##############################################
# Register composite samples
##############################################

# Make list of illumina coassemblies, if MERGE_ILLUMINA_SAMPLES in config

merged_illumina_samples = dict()
if 'MERGE_ILLUMINA_SAMPLES' in config:
    if config['MERGE_ILLUMINA_SAMPLES'] is None:
        print('WARNING in ', config_path, ': MERGE_ILLUMINA_SAMPLES list of samples is empty. Skipping sample-merging.')
    else:
        merged_illumina_samples = config['MERGE_ILLUMINA_SAMPLES']
        #print(merged_illumina_samples)

# Figure out SPAdes version
spades_script = 'spades.py' # from conda environment
if 'METASPADES_custom_build' in config:
    spades_script = config['METASPADES_custom_build'] # custom built, e.g. for higher K

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config
hybrid_assemblies = list()
if 'HYBRID' in config:
    if config['HYBRID'] is None:
        print('WARNING in ', config_path, ': HYBRID list of samples is empty. Skipping hybrid assembly.')
    else:
        #print("Hybrid assemblies:")
        for nano in config["HYBRID"]:
            for ilmn in config["HYBRID"][nano].split("+"):
                #print(" "+nano+"-"+ilmn)
                hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if enable_COASSEMBLY is set to "yes" in config
co_assemblies = list()
if 'enable_COASSEMBLY' in config and config['enable_COASSEMBLY'] is not None and config['enable_COASSEMBLY'] is True:
    if 'COASSEMBLY' in config:
        if config['COASSEMBLY'] is None:
            print('WARNING in', config_path, ': COASSEMBLY list of samples is empty. Skipping co-assembly.')
        else:
            #print("Coassemblies:")
            for co in config["COASSEMBLY"]:
                #print(" "+co)
                co_assemblies.append(co)

if config['METASPADES_qoffset'] in ('auto', '33', '64'):
    pass
else:
    print('ERROR in ', config_path, ': METASPADES_qoffset variable is not correct. "METASPADES_qoffset" variable should be auto, 33 or 64.')

if config['METASPADES_threads'] is None:
    print('ERROR in ', config_path, ': METASPADES_threads variable is empty. Please, complete ', config_path)
elif type(config['METASPADES_threads']) != int:
    print('ERROR in ', config_path, ': METASPADES_threads variable is not an integer. Please, complete ', config_path)

if config['METASPADES_memory'] is None:
    print('ERROR in ', config_path, ': METASPADES_memory variable is empty. Please, complete ', config_path)
elif type(config['METASPADES_memory']) != int:
    print('ERROR in ', config_path, ': METASPADES_memory variable is not an integer. Please, complete ', config_path)

if type(config['METASPADES_hybrid_max_k']) != int or config['METASPADES_hybrid_max_k']%2==0:
    print('ERROR in ', config_path, ': METASPADES_hybrid_max_k variable must be an odd integer')
elif 'METASPADES_custom_build' in config:
    if config['METASPADES_hybrid_max_k'] < 300:
        hybrid_max_k = config['METASPADES_hybrid_max_k']
    else:
        print('ERROR in ', config_path, ': METASPADES_hybrid_max_k variable must be below 300.')
else:
    if config['METASPADES_hybrid_max_k'] < 128:
        hybrid_max_k = config['METASPADES_hybrid_max_k']
    else:
        print('ERROR in ', config_path, ': METASPADES_hybrid_max_k variable must be below 128.')

if type(config['METASPADES_illumina_max_k']) != int or config['METASPADES_illumina_max_k']%2==0:
    print('ERROR in ', config_path, ': METASPADES_illumina_max_k variable must be an odd integer')
elif 'METASPADES_custom_build' in config:
    if config['METASPADES_illumina_max_k'] < 300:
        illumina_max_k = config['METASPADES_illumina_max_k']
    else:
        print('ERROR in ', config_path, ': METASPADES_illumina_max_k variable must be below 300.')
else:
    if config['METASPADES_illumina_max_k'] < 128:
        illumina_max_k = config['METASPADES_illumina_max_k']
    else:
        print('ERROR in ', config_path, ': METASPADES_illumina_max_k variable must be below 128.')

if config['MEGAHIT_threads'] is None:
    print('ERROR in ', config_path, ': MEGAHIT_threads variable is empty. Please, complete ', config_path)
elif type(config['MEGAHIT_threads']) != int:
    print('ERROR in ', config_path, ': MEGAHIT_threads variable is not an integer. Please, complete ', config_path)

if config['MEGAHIT_presets'] is None and config['MEGAHIT_custom'] is None:
    print('ERROR in ', config_path, ': MEGAHIT_presets list of MEGAHIT parameters to run per co-assembly is empty. Please, complete ', config_path)

mega_k_list = []
if 'MEGAHIT_custom' not in config:
    config['MEGAHIT_custom'] = None
elif config['MEGAHIT_custom'] is not None:
    # if the custom k-s are set, that should be added to the MEGAHIT assembly types
    if isinstance(config['MEGAHIT_custom'], str):
        config['MEGAHIT_custom'] = [config['MEGAHIT_custom']]
    for i, k_list in enumerate(config['MEGAHIT_custom']):
        if k_list and not k_list.isspace():
            if config['MEGAHIT_presets'] is None:
                config['MEGAHIT_presets'] = [f'meta-custom-{i+1}']
            else:
                config['MEGAHIT_presets'].append(f'meta-custom-{i+1}')
            mega_k_list.append(k_list)

if config['METAFLYE_presets'] is None:
    print('ERROR in ', config_path, ': METAFLYE_presets list of METAFLYE parameters to run per long-read assembly is empty. Please, complete ', config_path)

if config['MEGAHIT_memory'] is None:
    print('ERROR in ', config_path, ': MEGAHIT_memory variable is empty. Please, complete ', config_path)
elif type(config['MEGAHIT_memory']) != int:
    print('ERROR in ', config_path, ': MEGAHIT_memory variable is not an integer. Please, complete ', config_path)
elif type(config['MEGAHIT_memory']) == int:
    memory_config=config['MEGAHIT_memory']

# Define all the outputs needed by target 'all'

def illumina_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{kmer_dir}/{sample}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = ilmn_samples,
                    kmer_dir = "k21-" + str(illumina_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

def illumina_co_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta",
                    wd = working_dir,
                    omics = omics,
                    coassembly = co_assemblies,
                    assembly_preset = config['MEGAHIT_presets'])
    return(result)

def nanopore_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.fasta",
                    wd = working_dir,
                    omics = omics,
                    sample = config["NANOPORE"] if "NANOPORE" in config else [],
                    assembly_preset = config['METAFLYE_presets'])
    return(result)

def hybrid_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{sequence}.fasta",
                    wd = working_dir,
                    omics = omics,
                    assembly = hybrid_assemblies,
                    kmer_dir = "k21-" + str(hybrid_max_k),
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
        fwd="{wd}/{omics}/5-corrected-runs/{illumina}/{run}.1.fq.gz",
        rev="{wd}/{omics}/5-corrected-runs/{illumina}/{run}.2.fq.gz",
    shadow:
        "minimal"
    params:
        qoffset=config["METASPADES_qoffset"]
    resources:
        mem = lambda wildcards, attempt: attempt*config["METASPADES_memory"]
    log:
        "{wd}/logs/{omics}/5-corrected-runs/{illumina}/{run}_spadeshammer.log"
    threads: config['METASPADES_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #METASPADES
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
        multiple_runs = lambda wildcards, input: "yes" if len(input.files) > 1 else "no"
    shell:
        """
        if [ "{params.multiple_runs}" == "yes" ]; then
            cat {input} > combined.fq.gz
            rsync -a combined.fq.gz {output}
        else
            cd $(dirname {output})
            ln --symbolic --relative --force {input[0]} {output}
        fi
        """

###############################################################################################
########  Individual assembly of illumina samples
###############################################################################################
rule illumina_assembly_metaspades:
    input:
        fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
        rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
    output:
        "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/contigs.fasta",
        "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/scaffolds.fasta",
    shadow:
        "minimal"
    params:
        qoffset=config["METASPADES_qoffset"],
        asm_mode = "--meta",
        kmer_option = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
        kmer_dir = lambda wildcards: "k21-" + wildcards.maxk
    resources:
        mem = lambda wildcards, attempt: attempt*config["METASPADES_memory"]
    log:
        "{wd}/logs/{omics}/7-assembly/{illumina}/k21-{maxk}/{illumina}_metaspades.log"
    threads:
        config['METASPADES_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        time (
            {spades_script} {params.asm_mode} --only-assembler -1 {input.fwd} -2 {input.rev} -t {threads} -m {resources.mem} -o {params.kmer_dir} --tmp-dir tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a {params.kmer_dir}/* $remote_dir/
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
        fwd=rules.correct_spadeshammer.output.fwd,
        rev=rules.correct_spadeshammer.output.rev,
        ont="{wd}/{omics}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz"
    output:
        "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/contigs.fasta",
        "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/scaffolds.fasta",
    shadow:
        "minimal"
    params:
        qoffset=config["METASPADES_qoffset"],
        asm_mode = "--meta",
        kmer_option = lambda wildcards: get_metaspades_kmer_option(int(wildcards.maxk)),
        kmer_dir = lambda wildcards: "k21-" + wildcards.maxk
    resources:
        mem = lambda wildcards, attempt: attempt*config["METASPADES_memory"]
    log:
        "{wd}/logs/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/{nanopore}-{illumina}_metaspades.log"
    threads:
        config['METASPADES_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        time (
            {spades_script} {params.asm_mode} --only-assembler -1 {input.fwd} -2 {input.rev} --nanopore {input.ont} -t {threads} -m {resources.mem} -o {params.kmer_dir} --tmp-dir tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a {params.kmer_dir}/* $remote_dir/
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
        fwd=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz', wd=working_dir, omics=omics, illumina=config["COASSEMBLY"][wildcards.coassembly].split('+')),
        rev=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz', wd=working_dir, omics=omics, illumina=config["COASSEMBLY"][wildcards.coassembly].split('+'))
    output:
        coassemblies= "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa"
    shadow:
        "minimal"
    params:
        fwd_reads=lambda wildcards, input: ",".join(input.fwd),
        rev_reads=lambda wildcards, input: ",".join(input.rev),
        asm_params=lambda wildcards: get_megahit_parameters(wildcards, illumina_max_k, mega_k_list),
        memory_config=config['MEGAHIT_memory']
    resources:
        mem = lambda wildcards, input, attempt: min(900, len(input.fwd)*(memory_config+6*(attempt-1))),
        mem_bytes=lambda wildcards, input, attempt: min(900, len(input.fwd)*(memory_config+6*(attempt-1)))*1024*1024*1024
    log:
        "{wd}/logs/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}_{assembly_preset}_coassembly_megahit.log"
    threads: config['MEGAHIT_threads']
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        # Don't create the --out-dir directory as MEGAHIT wants it to not exist before
        time (
            megahit -1 {params.fwd_reads} -2 {params.rev_reads} -t {threads} -m {resources.mem_bytes} --out-dir assembly {params.asm_params}
        ) >& {log}
        cd assembly
        tar cfz intermediate_contigs.tar.gz intermediate_contigs && rm -rf intermediate_contigs
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        rsync -a * $remote_dir/
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
        options = lambda wildcards: config['METAFLYE_presets'][wildcards.assembly_preset]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
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
        assembly_preset = '|'.join(config['MEGAHIT_presets'])
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        r"""
        perl -ne 's/^>k(\d+)_(\d+) (.*)len=(\d+)/>MEGAHIT.{wildcards.assembly_preset}.{wildcards.coassembly}_NODE_$2_length_$4_k_$1/ if m/^>/; print $_;' < {input} > {output}
        """
