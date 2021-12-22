#!/usr/bin/env python

'''
Assembling metagenomes from combinations of illumina/bgi-seq and nanopore sequencing data.

Authors: Carmen Saenz, Mani Arumugam
'''

include: 'scripts/07-common-rules.smk'

localrules: mark_circular_metaspades_contigs, mark_circular_flye_contigs, rename_megahit_contigs

# snakemake --use-singularity --snakefile /emc/cbmr/users/rzv923/MIntO/07-assembly.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --use-conda --conda-prefix /emc/cbmr/users/rzv923/ibdmdb_test/tmp_porus/ --configfile assemblies.smk.yaml --jobs 10
# snakemake --use-singularity --snakefile /emc/cbmr/users/rzv923/MIntO/07-assembly.smk --restart-times 1 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A"  --use-conda --conda-prefix /data/MIntO_snakemake_env/ --configfile assemblies.smk.yaml --jobs 10

# cobinning metaG - BWA
#configfile: "metaG_samples_HQ_reads.yaml"

# Figure out SPAdes version

spades_script = 'spades.py' # from conda environment
if 'METASPADES_custom_build' in config:
    spades_script = config['METASPADES_custom_build'] # custom built, e.g. for higher K

# Make list of nanopore-illumina hybrid assemblies, if HYBRID in config
hybrid_assemblies = list()
if 'HYBRID' in config:
    #print("Hybrid assemblies:")
    for nano in config["HYBRID"]:
        for ilmn in config["HYBRID"][nano].split("+"):
            #print(" "+nano+"-"+ilmn)
            hybrid_assemblies.append(nano+"-"+ilmn)

# Make list of illumina coassemblies, if COASSEMBLY in config
co_assemblies = list()
if 'COASSEMBLY' in config:
    #print("Coassemblies:")
    for co in config["COASSEMBLY"]:
            #print(" "+co)
            co_assemblies.append(co)

# some variables

local_dir = config['local_dir']
working_dir = config['working_dir']
omics = config['omics']
illumina_max_k = config['METASPADES_illumina_max_k']
hybrid_max_k = config['METASPADES_hybrid_max_k']


# Define all the outputs needed by target 'all'

def illumina_single_assembly_output():
    result = expand("{wd}/{omics}/7-assembly/{sample}/{kmer_dir}/{sample}.{sequence}.fasta.len", 
                    wd = working_dir,
                    omics = omics,
                    sample = config["ILLUMINA"] if "ILLUMINA" in config else [],
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
rule correct_spadeshammer:
    input: 
        fwd="{wd}/{omics}/4-hostfree/{illumina}/{illumina}.1.fq.gz",
        rev="{wd}/{omics}/4-hostfree/{illumina}/{illumina}.2.fq.gz"
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
        time ({spades_script} --only-error-correction -1 {input.fwd} -2 {input.rev} -t {threads} -m {resources.mem} -o {params.tmp_dir} --tmp-dir {params.tmp_dir}/tmp --phred-offset {params.qoffset}; \
        rsync -a {params.tmp_dir}/corrected/{wildcards.illumina}.1.fq.00.0_0.cor.fastq.gz {output.fwd}; rsync -a {params.tmp_dir}/corrected/{wildcards.illumina}.2.fq.00.0_0.cor.fastq.gz {output.rev};) >& {log}
        rm -rf {params.tmp_dir}
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
    params:
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{illumina}_assembly_metaspades/{illumina}".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina),
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
        time ({spades_script} --meta --only-assembler -1 {input.fwd} -2 {input.rev} -t {threads} -m {resources.mem} -o {params.tmp_asm}/{params.kmer_dir} --tmp-dir {params.tmp_asm}/tmp --phred-offset {params.qoffset} -k {params.kmer_option}
        rsync -a {params.tmp_asm}/{params.kmer_dir}/* $remote_dir/ ) >& {log}
        rm -rf {params.tmp_asm}/{params.kmer_dir}
        rm -rf {params.tmp_asm}/tmp
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
        "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/contigs.fasta",
        "{wd}/{omics}/7-assembly/{nanopore}-{illumina}/k21-{maxk}/scaffolds.fasta",
    params:
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{nanopore}-{illumina}_assembly_metaspades/{nanopore}-{illumina}".format(local_dir=local_dir, omics=wildcards.omics, illumina=wildcards.illumina, nanopore=wildcards.nanopore),
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
        time ({spades_script} --meta --only-assembler -1 {input.fwd} -2 {input.rev} --nanopore {input.ont} -t {threads} -m {resources.mem} -o {params.tmp_asm}/{params.kmer_dir} --tmp-dir {params.tmp_asm}/tmp --phred-offset {params.qoffset} -k {params.kmer_option}
        rsync -a {params.tmp_asm}/{params.kmer_dir}/* $remote_dir/ ) >& {log}
        rm -rf {params.tmp_asm}/{params.kmer_dir}
        rm -rf {params.tmp_asm}/tmp
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
        tmp_asm=lambda wildcards: "{local_dir}/{omics}_{name}_coassembly_megahit".format(local_dir=local_dir, omics=wildcards.omics, name=wildcards.coassembly),
        fwd_reads=lambda wildcards, input: ",".join(input.fwd),
        rev_reads=lambda wildcards, input: ",".join(input.rev),
    resources:
        mem = lambda wildcards, input, attempt: len(input.fwd)*(5+6*attempt),
        mem_bytes=lambda wildcards, input, attempt: len(input.fwd)*(5+6*attempt)*1024*1024*1024
    log:
        "{wd}/logs/{omics}/7-assembly/{coassembly}/{assembly_preset}/{coassembly}_{assembly_preset}_coassembly_megahit.log"
    threads: config['MEGAHIT_threads']
    conda: 
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        mkdir -p {params.tmp_asm}/tmp-{wildcards.assembly_preset}
        # Don't create the --out-dir directory as MEGAHIT wants it to not exist before
        time (megahit -1 {params.fwd_reads} -2 {params.rev_reads} -t {threads} -m {resources.mem_bytes} --out-dir {params.tmp_asm}/{wildcards.assembly_preset} --tmp-dir {params.tmp_asm}/tmp-{wildcards.assembly_preset} --presets {wildcards.assembly_preset}) >& {log}
        rm -rf {params.tmp_asm}/tmp-{wildcards.assembly_preset}
        remote_dir=$(dirname {output[0]})
        mkdir -p $remote_dir
        rsync -a {params.tmp_asm}/{wildcards.assembly_preset}/* $remote_dir/
        rm -rf {params.tmp_asm}/{wildcards.assembly_preset}
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
    singularity: "docker://quay.io/biocontainers/flye:2.8.3--py38h69e0bdc_1"
    shell:
        """
        mkdir -p $(dirname {output[0]})
        flye --nano-raw {input} --out-dir $(dirname {output[0]}) --threads {threads} --meta --genome-size 3.0m {params.options} >& {log}
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
    shell:
        """ 
        perl -ne 's/^>k(\d+)_(\d+) (.*)len=(\d+)/>MEGAHIT.{wildcards.assembly_preset}.{wildcards.coassembly}_NODE_$2_length_$4_k_$1/ if m/^>/; print $_;' < {input} > {output}
        """ 

# Generic rules
