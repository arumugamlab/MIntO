#! /usr/bin/env python3

# Borrows from vamb's code to make sure that we construct the abundance.npz file just like they would!

from pathlib import Path
import gzip
import argparse
import numpy as _np
import vamb
import vamb.vambtools as _vambtools
from vamb.parsecontigs import CompositionMetaData

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

# Parse command line arguments

parser = argparse.ArgumentParser(description="Make avamb-style npz file from jgi-style depth file")
parser.add_argument("--fasta",         required=True, type=Path, help="fasta file with all contigs (required)")
parser.add_argument("--abundance-tsv", required=True, type=Path, help="abundance TSV file (required)")
parser.add_argument("--output",        required=True, type=Path, help="npz output file (required)")
parser.add_argument("--minlength",     required=True, type=int,  help="minimum scaffold/contig length used for filtering input fasta (required)")
args = parser.parse_args()

# Start a fasta file iterator for input fasta

fiter = fasta_iter(str(args.fasta))

# Initiate arrays for contig names, lengths and mask status
# mask==True means sequences is included in analysis, so we mark them all True

contignames: list[str] = list()
lengths = _vambtools.PushArray(_np.int32)
mask    = bytearray()

# Iterate through all the fasta files in input

for header, seq in fiter:
    contignames.append(header)
    lengths.append(len(seq))
    mask.append(True)

# Call CompositionMetaData constructor
# It needs numpy arrays, so we convert them on the fly

comp_metadata = CompositionMetaData(
                                    _np.array(contignames, dtype=object),
                                    lengths.take(),
                                    _np.array(mask, dtype=bool),
                                    args.minlength,
                                    )

# Make abundance object and write it out

abundance = vamb.parsebam.Abundance.from_tsv(args.abundance_tsv, comp_metadata)
abundance.save(args.output)
