#!/bin/env python2

"""
thanks to Kamil Slowikowski (https://github.com/slowkow)
for the original GTF.py lib: https://gist.github.com/slowkow/8101481
"""

from GTF import dataframe
import argparse
import os.path


def is_valid_file(argparser, arg):
    if not os.path.exists(arg):
        argparser.error("The file %s does not exist!" % arg)
    else:
        return arg


parser = argparse.ArgumentParser(description='Input GTF file and optional parameters')

parser.add_argument("-i", dest="infile", required=True,
                    help="input gtf file", metavar="INFILE",
                    type=lambda x: is_valid_file(parser, x))

parser.add_argument("-o", dest="outfile",
                    help="output csv file, defaults to `input file name`.csv", metavar="OUTFILE")

parser.add_argument("-columns", "-c", dest="columns", nargs="+",
                    help="list of columns to write to csv, e. g. `-columns column1 column2 column3`", metavar="COLUMN")

args = parser.parse_args()

args.outfile = args.outfile or os.path.splitext(args.infile)[0]+'.csv'

df = dataframe(args.infile)
print(" creating csv...")
df.to_csv(args.outfile, index=False, columns=args.columns)
print("SUCCESS: %s created." % args.outfile)
