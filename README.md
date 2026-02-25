# motif-mark
Visualizing motifs from FASTA files using object oriented programming

This script reads a FASTA file and a list of motifs, finds all motif occurrences in the sequences, including overlapping matches, and produces a PNG image showing:
    gene backbone
    exon regions 
    motif locations 
    a legend of features

# Requirements
pycairo must be installed

-f: FASTA file
-m: file contiaining one motif per line

# Input files
FASTA file
    can contain multiple sequences
    exons must be uppercase letters 
    introns must be lowercase letters

    Example: 
        >chr1
        atcgGTCAGAatcg

# Motif file
one motif per line
motifs can incude IUPAC codes and can be DNA or RNA
    Example:
    YGCY
    catag
    GCAUG
    YYYYYYYYYY

# Running Script
python3 motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt

# Output
one PNG file using the same prefix as the input FASTA file
    Example:
    Figure_1.png