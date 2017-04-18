#!/usr/bin/env python3
"""
CorrectCount.py
This script produces the abundance values (read counts, CPM and TPM) for each genes
by taking into account uniquely mapped and multi mapped reads (by default) through the
use of featureCounts from the Subread package.
For more information about its use, please read the MANUAL.md
"""
__author__ = "Vincent Boivin, Gabrielle Deschamps-Francoeur, and Michelle Scott"
__email__ = "Michelle.Scott@Usherbrooke.ca"

import argparse
import sys, os
import tests
import count_to_cpm

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file produced with CorrectAnnotation in .gtf format.")
parser.add_argument("bam_file", help="alignment file in .bam format.")
parser.add_argument("output", help="Name of the output file holding the counts per genes. ")

parser.add_argument("-c", "--countType", help="Decide whether to consider only uniquely mapped reads (uniqueOnly), only multi-mapped reads (multiOnly) or both (both) to produce read count values for genes. Default: both", choices=['uniqueOnly','multiOnly','both'], default='both')
parser.add_argument("-s", "--strand",help="Strandedness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default: 0", type=int, choices=[0, 1, 2],default=0)
parser.add_argument("-m", "--minOverlap", help="Minimum overlap a read must have with a feature to be assigned to its associated gene by featureCounts. Default: 10", type=int, default=10)
parser.add_argument("-p", "--paired", help="Use this option if you work on a paired-end dataset. Used by featureCounts.", action="store_true")
parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int, default=1)
parser.add_argument("-r", "--rawOnly", help="Use this option to output the featureCounts raw read counts only and not their calculated CPM and TPM values.", action="store_true")
args = parser.parse_args()

gtf_file=args.annotation
bam_file=args.bam_file
output=args.output
count_type=args.countType

strand=args.strand
minOverlap=args.minOverlap
thread=args.thread
paired=args.paired
rawonly=args.rawOnly

tests.check_gtf(gtf_file)
tests.check_bam(bam_file)
tests.check_output(output)
print('okay')
sys.exit()

if paired != True:
    paired=''
else:
    paired='-p '

if count_type =='uniqueOnly':
    command="featureCounts " \
            "--minOverlap %d " \
            "--largestOverlap " \
            "-s %d " \
            "-C " \
            "-T %s " \
            "%s" \
            "-a %s " \
            "-o %s " \
            "%s" %(minOverlap,strand,thread,paired,gtf_file,output_prefix,bam_file)
    os.system(command)
elif count_type=='multiOnly':
    print('Do the thing Gabe')
    print('...')
else:
    #For both, default.
    print('Do the thing Gabe')
    print('...')

if rawonly!=True:
    count_to_cpm.add_pm_counts(output_prefix,gtf_file,bam_file)

