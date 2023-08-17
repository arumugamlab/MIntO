#!/usr/bin/env python

localrules: bam_index, make_minimap_index, make_abi_blastdb, get_hq_mag_seqs, get_long_reads_for_mag, list_shortreads_for_mag

WORKDIR = config['working_dir']

if 'MAG_VERSION' in config: 
    MAG_VERSION = config['MAG_VERSION'] 
else:
    MAG_VERSION = 'undef'

if 'MAG_DB_DIR' in config: 
    MAG_DB_DIR = config['MAG_DB_DIR'] 
else:
    MAG_DB_DIR = 'undef'

# Make minimap index
rule make_minimap_index:
    input:
        f'{WORKDIR}/8-3-unique-MAGs/unique/{MAG_VERSION}.fasta'
    output:
        f'{MAG_DB_DIR}/{MAG_VERSION}.mmi'
    log:
        f'{MAG_DB_DIR}/{MAG_VERSION}.mmi.log'
    shell:
        """
        minimap2 -x map-ont -d {output} {input} >& {log}
        """

ruleorder: align_longreads_to_mag_db > bam_sort

# Make relevant long-read alignments if missing
rule align_longreads_to_mag_db:
    input: 
        db=f'{MAG_DB_DIR}/{MAG_VERSION}.mmi', 
        reads='{wd}/6-corrected/{nanopore}/{nanopore}.nanopore.fq.gz',
    output: 
        f'{{wd}}/6-mapping/{{nanopore}}/{{nanopore}}.{MAG_VERSION}.minimap.sorted.bam',
    log: 
        f'{{wd}}/6-mapping/{{nanopore}}/{{nanopore}}.{MAG_VERSION}.minimap.log',
    wildcard_constraints:
        nanopore = '[^-]+'
    params:
        sort_mem = lambda wildcards: 80 if wildcards.nanopore.endswith('X') else 10
    resources:
        mem = lambda wildcards: 4*80 if wildcards.nanopore.endswith('X') else 4*10
    threads: 16
    shell:
        """
        mkdir -p $(dirname {output})
        minimap2 -x map-ont -t {threads} -a {input.db} {input.reads} 2>{log} | samtools sort -@ 4 -m 80G > {output}
        """

rule align_shortreads_to_mag_db:
    input:
        fwd='{wd}/6-corrected/{illumina}/{illumina}.1.fq.gz',
        rev='{wd}/6-corrected/{illumina}/{illumina}.2.fq.gz',
        db=f'{MAG_DB_DIR}/{MAG_VERSION}.bwt.2bit.64'
    output:
        f'{{wd}}/6-mapping/{{illumina}}/{{illumina}}.{MAG_VERSION}.bam',
    log:
        f'{{wd}}/6-mapping/{{illumina}}/{{illumina}}.{MAG_VERSION}.bwa.log',
    params:
        db_name = '{location}/{name}'.format(location = MAG_DB_DIR, name = MAG_VERSION)
    resources:
        mem=20
    threads: 16
    shell:
        """
        bwa-mem2 mem -a -t {threads} {params.db_name} {input.fwd} {input.rev} 2>{log} \
                | msamtools filter -bS  -p 95 -l 45 - > {output}
        """

# Get sequences in this HQ MAG
rule get_hq_mag_seqs:
    input:
        '{wd}/8-3-unique-MAGs/unique/nr99/{mag}.fna'
    output:
        '{wd}/8-4-MAG-improvement/{mag}/{mag}.seqs.list'
    shell:
        """
        mkdir -p $(dirname {output})
        ln -s {input} $(dirname {output})/
        grep "^>" {input} | sed "s/>//" > {output}
        """

# Get the ONT reads mapping to this MAG
rule get_long_reads_for_mag:
    input:
        seqs='{wd}/8-4-MAG-improvement/{mag}/{mag}.seqs.list',
        bam=f'{{wd}}/6-mapping/{{nanopore}}/{{nanopore}}.{MAG_VERSION}.minimap.sorted.bam',
        bai=f'{{wd}}/6-mapping/{{nanopore}}/{{nanopore}}.{MAG_VERSION}.minimap.sorted.bam.bai',
    resources:
        mem = 4
    output:
        '{wd}/8-4-MAG-improvement/{mag}/{nanopore}.ONT.input.fq.gz'
    shell:
        """
        for i in $(cat {input.seqs}); do
            samtools view -h {input.bam} $i | samtools fastq - | paste - - - - | sort | uniq | tr '\\t' '\\n'
        done | bgzip > {output}
        """

# Get a list of the ILLUMINA reads mapping to this MAG
rule list_shortreads_for_mag:
    input:
        seqs='{wd}/8-4-MAG-improvement/{mag}/{mag}.seqs.list',
        bam=f'{{wd}}/6-mapping/{{illumina}}/{{illumina}}.{MAG_VERSION}.sorted.bam',
        bai=f'{{wd}}/6-mapping/{{illumina}}/{{illumina}}.{MAG_VERSION}.sorted.bam.bai',
    output:
        '{wd}/8-4-MAG-improvement/{mag}/{illumina}.ILN.input.list',
    resources:
        mem = 4
    shell:
        """
        for i in $(cat {input.seqs}); do
            samtools view {input.bam} $i | cut -f1
        done | sort -u > {output}
        """

# Get the ILLUMINA reads mapping to this MAG
rule retrieve_shortreads_for_mag:
    input: 
        names='{wd}/8-4-MAG-improvement/{mag}/{illumina}.ILN.input.list',
        fq='{wd}/6-corrected/{illumina}/{illumina}.{pair}.fq.gz',
    output:
        '{wd}/8-4-MAG-improvement/{mag}/{illumina}.ILN.input.{pair}.fq.gz',
    threads: 1
    resources:
        mem = 4
    shell:
        """
        filterSeq -t fastq --list {input.names} --gzip --input {input.fq} --output {output}
        """

# Utility functions

# Create blastdb

rule make_abi_blastdb:
    input:
        '{something}.f{asta}'
    output:
        '{something}.f{asta}.xnd'
    log:
        '{something}.f{asta}.xdformat.log'
    shell:
        """
        xdformat -n {input} >& {log}
        """

# Sort a bam file
# 40G per thread hardcoded below
rule bam_sort:
    input:
        '{something}.bam'
    output:
        '{something}.sorted.bam'
    threads: 4
    resources: 
        mem = lambda wildcards, threads: 40*threads
    shell:
        """
        samtools sort -@ {threads} -m 40G {input} > {output}
        """

# Index a bam file
rule bam_index:
    input:
        '{something}.bam'
    output:
        '{something}.bam.bai'
    threads: 8
    shell:
        """
        samtools index -@ {threads} {input}
        """

# Filter multiple files by length and output to one file
def filter_fasta_list_by_length(infile_list, outfile, min_length=2500):

    import os
    from pathlib import Path

    # Create the output folder
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))

    # Open the output file
    with open(outfile, 'w') as out:
        for infile in infile_list:
        # Go through the fasta file
            fiter = fasta_iter(infile)
            for entry in fiter:
                header, seq = entry
                if len(seq) >= int(min_length):
                    out.write(f'>{header}\n{seq}\n')
