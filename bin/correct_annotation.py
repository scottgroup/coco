#!/usr/bin/env python3
"""
CorrectAnnotation.py
This script calls gapped_gtf.py to produce the modified annotation file in .gtf format.
For more information about its use, please read the MANUAL.md
"""
__author__ = "Vincent Boivin, Gabrielle Deschamps-Francoeur, and Michelle Scott"
__email__ = "Michelle.Scott@Usherbrooke.ca"

import argparse
import sys, os
import tests
import gapped_gtf

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file in .gtf format.")

parser.add_argument("-o", "--output", help="Name of output gtf. Default: add 'CorrectAnnotation' suffix to original gtf name.", default='None')
args = parser.parse_args()

gtf_file=args.annotation
output=args.output
print('output', output)


tests.check_gtf(gtf_file,check_gene_biotype=True)
gapped_gtf.CorrectAnnotation(gtf_file,output)
