#!/usr/bin/env python3
"""
CorrectBedgraph.py
This script produces bedgraph from paired-end datasets.
For more information about its purpose and use, please read the MANUAL.md
"""
__author__ = "Vincent Boivin, Gabrielle Deschamps-Francoeur, and Michelle Scott"
__email__ = "Michelle.Scott@Usherbrooke.ca"

import argparse
import sys, os
import tests

parser = argparse.ArgumentParser()

parser.add_argument("bam_file", help="path to .bam file")
parser.add_argument("output", help="name of output bedgraph")
args = parser.parse_args()



output=args.output
bam_file=args.output

tests.check_bam(bam_file)
print('output', output)