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
args = parser.parse_args()

gtf_file=args.annotation
output=args.output
biotypes_embedded=args.biotypes
biotypes_embedded=biotypes_embedded.split(',')
verbose = args.verbose

tests.check_gtf(gtf_file,check_gene_biotype=True)

gapped_gtf.correct_annotation(gtf_file,output, verbose, biotypes_embedded=biotypes_embedded)
