#!/usr/bin/env python

localrules: bam_idx, faidx, get_fasta_length

#########################################
# SAM/BAM processing rules and functions
#########################################

# Index a bam file

rule bam_idx:
    #index a .bam file
    input:
        '{some}.bam'
    output:
        '{some}.bam.bai'
    threads: 8
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        if [ -s {input[0]} ] #if the input file is nonempty...
        then
            samtools index -@ {threads} {input}
        else
            touch {output}
        fi
        """

#########################################
# Fasta processing rules and functions
#########################################

rule faidx:
    #index a .fasta file
    input: 
        '{something}.f{asta}'
    output: 
        '{something}.f{asta}.fai'
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell: 
        """
        if [ -s {input[0]} ] #if the input is nonempty...
        then
            samtools faidx {input}
        else
            touch {output}
        fi
        """

rule get_fasta_length:
    input:
        '{something}.f{asta}.fai'
    output:
        '{something}.f{asta}.len'
    shell:
        """ 
        cut -f1,2 {input} | sort -k2,2nr > {output}
        """ 

def fasta_iter(infile):
    """
    Source: https://www.biostars.org/p/710/
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """

    from itertools import groupby

    #first open the file outside
    fh = open(infile)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

# Filter one fasta file by length and output to one file

def filter_fasta_by_length(infile, outfile, min_length=2500):

    import os

    # Create the output folder
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))

    # Open the output file
    with open(outfile, 'w') as out:
        # Go through the fasta file
        fiter = fasta_iter(infile)
        for entry in fiter:
            header, seq = entry
            if len(seq) >= min_length:
                out.write(f'>{header}\n{seq}\n')

# Filter multiple files by length and output to one file
def filter_fasta_list_by_length(infile_list, outfile, min_length=2500):

    import os

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
