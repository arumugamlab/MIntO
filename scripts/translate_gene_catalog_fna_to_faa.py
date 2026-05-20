#!/usr/bin/env python3

"""
translate_gene_catalog_fna_to_faa.py

Robust, reproducible translation of bacterial gene nucleotide sequences (FNA)
into amino acid sequences (FAA), with explicit handling of:

  • Partial genes (missing start and/or stop codons)
  • Conservative vs. full start codon models
  • Automatic detection of alternate genetic code usage
    (TGA reassigned to tryptophan, e.g. Mycoplasma)
  • Explicit reporting of gene completeness via a compact partial code
  • Transparent warnings for rare start codons and alternate code switching

This script is designed for "infrastructure-grade" use in large-scale
genome annotation, metagenomics, and functional profiling pipelines where:

  • The original reference genome may be unavailable
  • Predicted genes may be truncated
  • Genetic code ambiguities must be handled safely and explicitly
  • Downstream proteomics, pangenomics, or functional annotation depends
    on deterministic, auditable translation behavior

Key design principles:
  - Each sequence is a potentially partial gene sequence
  - No silent assumptions
  - Deterministic translation
  - Explicit biological edge-case handling
  - Machine- and human-readable outputs
  - Safe defaults aligned with Prodigal and NCBI bacterial standards

Output:
  FASTA-formatted amino acid sequences written to stdout.
  Each record includes a 'partial=XY' flag where:
    X = 0 if start codon present, 1 if missing
    Y = 0 if stop codon present, 1 if missing

Genetic code logic:
  - Default: NCBI Table 11 (bacterial)
  - Automatically switches to TGA→W if TGA is observed internally
  - Final TGA is treated as stop unless internal TGA forces reassignment

Author: Mani Arumugam
"""

import sys
import re
from itertools import product
import argparse

# --- Codon Tables ---
STANDARD_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

ALT_CODON_TABLE = STANDARD_CODON_TABLE.copy()
ALT_CODON_TABLE['TGA'] = 'W'

STANDARD_STOP_CODONS = {'TAA', 'TAG', 'TGA'}
ALT_STOP_CODONS = {'TAA', 'TAG'}

# --- Start Codon Sets ---
PRODIGAL_STARTS = {'ATG', 'GTG', 'TTG'}
FULL_STARTS_TABLE11 = {'ATG', 'GTG', 'TTG', 'ATT', 'ATC', 'ATA', 'CTG'}
FULL_STARTS_TABLE4 = FULL_STARTS_TABLE11 | {'TTA'}

# --- Parse FASTA ---
def parse_fasta(handle):
    header = None
    seq_lines = []
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield header, ''.join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        yield header, ''.join(seq_lines)

# --- Format Codon Table Summary ---
def format_codon_table_summary(codon_table, start_codons, label):
    codons = [''.join(p) for p in product('TCAG', repeat=3)]
    AAs = []
    Starts = []
    Base1, Base2, Base3 = [], [], []

    for codon in codons:
        AAs.append(codon_table.get(codon, 'X'))
        Starts.append('M' if codon in start_codons else '-')
        Base1.append(codon[0])
        Base2.append(codon[1])
        Base3.append(codon[2])

    print(f"\n# {label}")
    print("    AAs  = " + ''.join(AAs))
    print("  Starts = " + ''.join(Starts))
    print("  Base1  = " + ''.join(Base1))
    print("  Base2  = " + ''.join(Base2))
    print("  Base3  = " + ''.join(Base3))

# --- Translate one sequence ---
def translate_sequence(nt_seq, header, start_codon_mode, prodigal_starts, full_starts_table11, full_starts_table4):
    seq = re.sub(r'\s+', '', nt_seq.upper())

    if len(seq) < 3:
        print(f"Warning: Sequence '{header}' is shorter than 3 bp. Skipping.", file=sys.stderr)
        return None

    if len(seq) % 3 != 0:
        print(f"Warning: Sequence '{header}' length not divisible by 3. Truncating.", file=sys.stderr)
        seq = seq[:len(seq) - len(seq) % 3]

    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if not codons:
        return None

    first_codon = codons[0]
    last_codon = codons[-1]
    internal_codons = codons[1:-1]

    uses_TGA_as_Trp = 'TGA' in codons[:-1]
    if uses_TGA_as_Trp:
        print(f"Warning: Using alternate genetic code (TGA → W) for sequence '{header}'", file=sys.stderr)
        codon_table = ALT_CODON_TABLE
        stop_codons = ALT_STOP_CODONS
        start_codons = full_starts_table4 if start_codon_mode == 'full' else prodigal_starts
    else:
        codon_table = STANDARD_CODON_TABLE
        stop_codons = STANDARD_STOP_CODONS
        start_codons = full_starts_table11 if start_codon_mode == 'full' else prodigal_starts

    if first_codon in start_codons:
        aa_seq = ['M']
        start_flag = '0'
        if start_codon_mode == 'full' and first_codon not in prodigal_starts:
            print(f"Warning: Rare start codon '{first_codon}' used in sequence '{header}'", file=sys.stderr)
    else:
        aa_seq = [codon_table.get(first_codon, 'X')]
        start_flag = '1'

    for codon in internal_codons:
        aa_seq.append(codon_table.get(codon, 'X'))

    if last_codon in stop_codons:
        stop_flag = '0'
    else:
        stop_flag = '1'
        aa_seq.append(codon_table.get(last_codon, 'X'))

    partial = f"{start_flag}{stop_flag}"
    return ''.join(aa_seq), partial

# --- Main ---
def main():
    parser = argparse.ArgumentParser(
        description="Translate bacterial gene nucleotide sequences (fna) to amino acid sequences (faa), detecting partial genes and alternate codon usage."
    )
    parser.add_argument('--fna', help='Input FASTA file of gene nucleotide sequences')
    parser.add_argument('--report-codon-table', action='store_true',
                        help='Print the codon table in NCBI 5-line summary format and exit')
    parser.add_argument('--start-codons', choices=['prodigal', 'full'], default='prodigal',
                        help='Set of start codons to use (default: prodigal)')

    args = parser.parse_args()

    if args.start_codons == 'prodigal':
        start_codons_std = PRODIGAL_STARTS
        start_codons_alt = PRODIGAL_STARTS
    else:
        start_codons_std = FULL_STARTS_TABLE11
        start_codons_alt = FULL_STARTS_TABLE4

    if args.report_codon_table:
        format_codon_table_summary(STANDARD_CODON_TABLE, start_codons_std, "Standard Codon Table (TGA = stop)")
        format_codon_table_summary(ALT_CODON_TABLE, start_codons_alt, "Alternate Codon Table (TGA = W)")
        sys.exit(0)

    if not args.fna:
        parser.error("the following arguments are required: --fna (unless --report-codon-table is used)")

    try:
        with open(args.fna, 'r') as handle:
            for header, seq in parse_fasta(handle):
                result = translate_sequence(
                    seq, header,
                    start_codon_mode=args.start_codons,
                    prodigal_starts=PRODIGAL_STARTS,
                    full_starts_table11=FULL_STARTS_TABLE11,
                    full_starts_table4=FULL_STARTS_TABLE4
                )
                if result is None:
                    continue
                aa_seq, partial = result
                print(f">{header} partial={partial}")
                print(aa_seq)
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: File not found: {args.fna}\n")
        sys.exit(1)

# --- Entry point ---
if __name__ == "__main__":
    if sys.version_info.major < 3:
        sys.stderr.write("ERROR: This script requires Python 3.\n")
        sys.exit(1)
    main()
