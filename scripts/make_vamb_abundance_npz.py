#! /usr/bin/env python

import gzip
import argparse
import numpy as _np
import vamb

def validate_input_array(array):
    """
    Source: vamb codebase
    Returns array similar to input array but C-contiguous and with own data.
    """
    if not array.flags['C_CONTIGUOUS']:
        array = _np.ascontiguousarray(array)
    if not array.flags['OWNDATA']:
        array = array.copy()

    assert (array.flags['C_CONTIGUOUS'] and array.flags['OWNDATA'])
    return array

def load_jgi(filehandle):
    """
    Source: vamb codebase
    Load depths from the --outputDepth of jgi_summarize_bam_contig_depths.
    See https://bitbucket.org/berkeleylab/metabat for more info on that program.

    Usage:
        with open('/path/to/jgi_depths.tsv') as file:
            depths = load_jgi(file)
    Input:
        File handle of open output depth file
    Output:
        N_contigs x N_samples Numpy matrix of dtype float32
    """

    header = next(filehandle)
    fields = header.split('\t')
    if not fields[:3] == ["contigName", "contigLen", "totalAvgDepth"]:
        raise ValueError('Input file format error: First columns should be "contigName,"'
        '"contigLen" and "totalAvgDepth"')

    columns = tuple([i for i in range(3, len(fields)) if not fields[i].rstrip().endswith("-var")])
    array = _np.loadtxt(filehandle, dtype=_np.float32, usecols=columns)
    return validate_input_array(array)

def fasta_iter(infile):
    """
    Source: https://www.biostars.org/p/710/
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """

    from itertools import groupby

    #first open the file outside
    fh = gzip.open(infile, 'rt') if infile.endswith('.gz') else open(infile)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

# Main

parser = argparse.ArgumentParser(description="Make avamb-style npz file from jgi-style depth file")
parser.add_argument("--fasta", required=True, help="fasta file with all contigs (required)")
parser.add_argument("--jgi", required=True, help="jgi depth file")
parser.add_argument("--samples", required=True, nargs='+', help="list of samples")
parser.add_argument("--output", required=True, help="npz output file")
args = parser.parse_args()

fiter = fasta_iter(args.fasta)
headers = [header for header, seq in fiter]
refhash = vamb.vambtools.hash_refnames(headers)

jgipath = args.jgi
file = gzip.open(jgipath, 'rt') if jgipath.endswith('.gz') else open(jgipath)
rpkms = load_jgi(file)

abundance = vamb.parsebam.Abundance(rpkms, args.samples, 0.95, refhash)
abundance.save(args.output)
