import argparse
import sys, os

parser = argparse.ArgumentParser()

parser.add_argument("bam_file", help="path to .bam file")
parser.add_argument("output", help="name of output bedgraph")
args = parser.parse_args()

output=args.output
print('output', output)