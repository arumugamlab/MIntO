#!/usr/bin/env python

'''
Assembling metagenomes from combinations of illumina/bgi-seq and nanopore sequencing data.

Authors: Carmen Saenz, Mani Arumugam
'''

include: '../scripts/07-common-rules.smk'

localrules: mark_circular_metaspades_contigs, mark_circular_flye_contigs, rename_megahit_contigs

#
# configuration yaml file
# import sys
import os.path
from os import path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, local_dir, minto_dir, script_dir, metadata
include: 'config_parser.smk'

# Variables from configuration yaml file

# Make list of illumina samples, if ILLUMINA in config
if 'ILLUMINA' in config:
    if config['ILLUMINA'] is None:
        print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)
    else:
        try:
            # Make list of illumina samples, if ILLUMINA in config
            ilmn_samples = list()
            #print("Samples:")
            for ilmn in config['ILLUMINA']:
                x = str(ilmn)
                if path.exists(working_dir+'/'+omics+'/4-hostfree/'+ x +'/'+ x +'.1.fq.gz') is True:
                    ilmn_samples.append(x)
                else:
                    raise TypeError('ERROR in ', config_path, ': ILLUMINA list of samples does not exist. Please, complete ', config_path)
        except TypeError:
            print('ERROR in ', config_path, ': ILLUMINA list of samples does not exist or has an incorrect format. Please, complete ', config_path)
else:
    print('ERROR in ', config_path, ': ILLUMINA list of samples is empty. Please, complete ', config_path)

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

# Make list of illumina coassemblies, if COASSEMBLY in config
co_assemblies = list()
if 'COASSEMBLY' in config:
    if config['COASSEMBLY'] is None:
        print('WARNING in ', config_path, ': COASSEMBLY list of samples is empty. Skipping co-assembly.')
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

if config['MEGAHIT_presets'] is None:
    print('ERROR in ', config_path, ': MEGAHIT_presets list of MEGAHIT parameters to run per co-assembly is empty. Please, complete ', config_path)

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
    result = expand("{wd}/{omics}/7-assembly/{sample}/{kmer_dir}/{sample}.{sequence}.fasta.len", 
                    wd = working_dir,
                    omics = omics,
                    sample = ilmn_samples,
                    kmer_dir = "k21-" + str(illumina_max_k),
                    sequence = ["contigs", "scaffolds"])
    return(result)

def illumina_co_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}.contigs.fasta.len",
                    wd = working_dir,
                    omics = omics,
                    coassembly = co_assemblies,
                    assembly_preset = config['MEGAHIT_presets'])
    return(result)

def nanopore_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{assembly_preset}/{sample}.assembly.fasta.len", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["NANOPORE"] if "NANOPORE" in config else [],
                    assembly_preset = config['METAFLYE_presets'])
    return(result)

def hybrid_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{assembly}/{kmer_dir}/{assembly}.{sequence}.fasta.len", 
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
        hybrid_assembly_output()

###############################################################################################
# Correct Illumina reads using SPAdes' spadeshammer
###############################################################################################

def get_hq_fastq_files(wildcards):
    if wildcards.omics == "metaT":
        return(expand("{wd}/{omics}/5-1-sortmerna/{illumina}/{illumina}.{pair}.fq.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    illumina = wildcards.illumina,
                    pair = [1, 2]))
    else:
        return(expand("{wd}/{omics}/4-hostfree/{illumina}/{illumina}.{pair}.fq.gz",
                    wd = wildcards.wd,
                    omics = wildcards.omics,
                    illumina = wildcards.illumina,
                    pair = [1, 2]))

rule correct_spadeshammer:
    input: 
        reads=get_hq_fastq_files
    output: 
        fwd="{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz",
        rev="{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz",
    params:
        tmp_dir=lambda wildcards: "{local_dir}/{omics}_{illumina}_correction_metaspades/".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina),
        qoffset=config["METASPADES_qoffset"]
    resources:
        mem = lambda wildcards, attempt: attempt*config["METASPADES_memory"]
    log:
        "{wd}/logs/{omics}/6-corrected/{illumina}/{illumina}_spadeshammer.log"
    threads: config['METASPADES_threads']
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml" #METASPADES
    shell:
        """ 
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.fwd})
        time (\
            {spades_script} --only-error-correction -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} -m {resources.mem} -o {params.tmp_dir} --tmp-dir {params.tmp_dir}/tmp --phred-offset {params.qoffset}
            rsync -a {params.tmp_dir}/corrected/{wildcards.illumina}.1.fq.00.0_0.cor.fastq.gz {output.fwd}; rsync -a {params.tmp_dir}/corrected/{wildcards.illumina}.2.fq.00.0_0.cor.fastq.gz {output.rev}
        ) >& {log}
        rm -rf {params.tmp_dir}
        """

###############################################################################################
########  Individual assembly of illumina samples
###############################################################################################
rule illumina_assembly_metaspades:
    input: 
        fwd=rules.correct_spadeshammer.output.fwd,
        rev=rules.correct_spadeshammer.output.rev
    output: 
        "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/contigs.fasta",
        "{wd}/{omics}/7-assembly/{illumina}/k21-{maxk}/scaffolds.fasta",
    params:
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{illumina}_assembly_metaspades".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina),
        qoffset=config["METASPADES_qoffset"],
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
        mkdir -p {params.tmp_asm}
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        time (\
            {spades_script} --meta --only-assembler -1 {input.fwd} -2 {input.rev} -t {threads} -m {resources.mem} -o {params.tmp_asm}/{params.kmer_dir} --tmp-dir {params.tmp_asm}/tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a {params.tmp_asm}/{params.kmer_dir}/* $remote_dir/ 
        ) >& {log}
        rm -rf {params.tmp_asm}/{params.kmer_dir}
        rm -rf {params.tmp_asm}/tmp
        rmdir {params.tmp_asm}
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
    params:
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{nanopore}-{illumina}_assembly_metaspades".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina, nanopore=wildcards.nanopore),
        qoffset=config["METASPADES_qoffset"],
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
        mkdir -p {params.tmp_asm}
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        time (\
            {spades_script} --meta --only-assembler -1 {input.fwd} -2 {input.rev} --nanopore {input.ont} -t {threads} -m {resources.mem} -o {params.tmp_asm}/{params.kmer_dir} --tmp-dir {params.tmp_asm}/tmp --phred-offset {params.qoffset} -k {params.kmer_option}
            rsync -a {params.tmp_asm}/{params.kmer_dir}/* $remote_dir/ 
        ) >& {log}
        rm -rf {params.tmp_asm}/{params.kmer_dir}
        rm -rf {params.tmp_asm}/tmp
        rmdir {params.tmp_asm}
        """

###############################################################################################
########  Co-assembly 
# This starts with 11G per sample, but if that fails, it increases by 5G per sample per repeated attempt
###############################################################################################
rule coassembly_megahit:
    input: 
        fwd=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.1.fq.gz', wd=working_dir, omics=omics, illumina=config["COASSEMBLY"][wildcards.coassembly].split('+')),
        rev=lambda wildcards: expand('{wd}/{omics}/6-corrected/{illumina}/{illumina}.2.fq.gz', wd=working_dir, omics=omics, illumina=config["COASSEMBLY"][wildcards.coassembly].split('+'))
    output: 
        coassemblies= "{wd}/{omics}/7-assembly/{coassembly}/{assembly_preset}/final.contigs.fa" 
    params:
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{name}_coassembly_megahit_{assembly_preset}".format(local_dir=local_dir, omics=wildcards.omics, name=wildcards.coassembly, assembly_preset=wildcards.assembly_preset),
        fwd_reads=lambda wildcards, input: ",".join(input.fwd),
        rev_reads=lambda wildcards, input: ",".join(input.rev),
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
        local_dir="{params.tmp_asm}/{wildcards.assembly_preset}"
        tmp_dir="{params.tmp_asm}/tmp-{wildcards.assembly_preset}"
        mkdir -p $tmp_dir
        # Don't create the --out-dir directory as MEGAHIT wants it to not exist before
        time (\
            megahit -1 {params.fwd_reads} -2 {params.rev_reads} -t {threads} -m {resources.mem_bytes} --out-dir $local_dir --tmp-dir $tmp_dir --presets {wildcards.assembly_preset}
        ) >& {log}
        rm -rf $tmp_dir
        cd $local_dir
        tar cfz intermediate_contigs.tar.gz intermediate_contigs && rm -rf intermediate_contigs
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        rsync -a $local_dir/* $remote_dir/
        cd {params.tmp_asm}
        cd ..
        rm -rf $local_dir
        rmdir {params.tmp_asm}
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
    singularity: "docker://quay.io/biocontainers/flye:2.9--py38h69e0bdc_0"
    #"docker://quay.io/biocontainers/flye:2.8.3--py38h69e0bdc_1"
    shell:
        """
        mkdir -p $(dirname {output[0]})
        time (\
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
                out.write(f'>MetaSPAdes.k21-{wildcards.maxk}.{wildcards.sample}_{header}\n')
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
                new_header = header + '_length_%s' % len(seq)
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
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """ 
        perl -ne 's/^>k(\d+)_(\d+) (.*)len=(\d+)/>MEGAHIT.{wildcards.assembly_preset}.{wildcards.coassembly}_NODE_$2_length_$4_k_$1/ if m/^>/; print $_;' < {input} > {output}
        """ 
