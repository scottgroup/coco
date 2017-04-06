import argparse
import sys, os

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file produced with CorrectAnnotation in .gtf format.")
parser.add_argument("bam_file", help="alignment file in .bam format.")
parser.add_argument("output_prefix", help="prefix for the output.")

parser.add_argument("-c", "--count_type", help="Number of threads to be used by featureCounts. Default: 1", choices=['unique_only','mm_only','both'], default='both')
parser.add_argument("-s", "--strand",help="Strandness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default", choices=[0, 1, 2],default=0)
parser.add_argument("-p", "--paired", help="Use this option if you work on paired-end dataset. Used by featureCounts.", action="store_true")
parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int, default=1)
args = parser.parse_args()


gtf=args.gtf
bam_file=args.bam_file
output=args.output

print('you selected CorrectCount')
print(sys.argv)