#!/usr/bin/env python3
"""
correct_annotation.py
This script calls gapped_gtf.py to produce the modified annotation file in .gtf format.
For more information about its use, please read the MANUAL.md
"""
__author__ = "Vincent Boivin, Gabrielle Deschamps-Francoeur, and Michelle Scott"
__email__ = "Michelle.Scott@Usherbrooke.ca"

import argparse
import tests
import gapped_gtf
import sys

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file in .gtf format.")

parser.add_argument("-o", "--output", help="Name of the correct_annotation output gtf. Default: add 'correct_annotation' suffix to original gtf name.", default='None')
parser.add_argument("-b", "--biotypes", help="List of gene biotypes to correct for. Must be given in a comma seperated list (no spaces). Default: snoRNA,scaRNA,tRNA,miRNA,snRNA ", default='snoRNA,scaRNA,tRNA,miRNA,snRNA')
parser.add_argument("-V", "--verbose", help="Print the progression of the gtf file reading",  action="store_true")
parser.add_argument("-f", "--fraction",
                    help="Minimal reciprocal overlap fraction to merge embedded genes to avoid conflict. "
                         "To disable, set to -1. Default: 0.85.",
                    default=0.85, type=float)
args = parser.parse_args()

gtf_file=args.annotation
output=args.output
biotypes_embedded=args.biotypes
biotypes_embedded=biotypes_embedded.split(',')
verbose = args.verbose
fraction = args.fraction

tests.check_gtf(gtf_file,check_gene_biotype=True)

if (fraction <= 0 or fraction > 1) and fraction != -1:
    print('-f/--fraction must be a value between 0 and 1, or -1 to disable. Exiting.', file=sys.stderr)
    sys.exit(1)

gapped_gtf.correct_annotation(gtf_file,output, verbose, fraction, biotypes_embedded=biotypes_embedded)
