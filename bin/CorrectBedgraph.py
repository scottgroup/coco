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
import correct_bg
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("bamfile", help="path to .bam file")
parser.add_argument("output", help="name of output bedgraph")
parser.add_argument("genomepath", help='complete path of the genome file for genomeCoverageBed for bedtools')
parser.add_argument("-m","--contains_multi", help='indicates that your file contains multimapped reads or not '
                                                  '(sorting will be longer, but accurate)', action="store_true")
parser.add_argument("-c", "--chunk_size", help="number of rows to treat at a time. Default: 2500000",
                    type=int, default=2500000 )
parser.add_argument("-u", "--ucsc_compatible", help="Add trackline and 'chr' prefix if needed.", action="store_true")
args = parser.parse_args()



output = args.output
bamfile = args.bamfile
chunk_size = args.chunk_size
ucsc = args.ucsc_compatible
genomepath = args.genomepath
multi = args.contains_multi

tests.check_bam(bamfile)
print('output:', output)

output_dir = os.path.dirname(output)
if output_dir == '':
    output_dir = os.getcwd()
print('output_dir:', output_dir)

x = correct_bg.prepare_bed12(bamfile, output_dir, multi)
if x !=0 :
    sys.exit('prepare_bed12 exit status: '+ str(x))

filelist = os.listdir(output_dir + '/chromo/')

for filename in filelist:

    if filename.endswith(".bed12") and 'corrected' not in filename:
        print('running chromosome: %s'%(filename.replace('.bed12','')))
        outfile = output_dir + '/chromo/corrected_' + os.path.basename(filename)

        for i, df_chunk in enumerate(pd.read_csv(output_dir + '/chromo/'+ filename, sep='\t',
                                    names=['chromo', 'start', 'end', 'name', 'score', 'strand', 'thick_start',
                                           'thick_end','RGB', 'nb_block', 'block_len', 'block_start'],
                                    dtype={'chromo': str}, chunksize=chunk_size)):
            if i == 0 :
                mode = 'w'
            else :
                mode = 'a'
            df_corrected = correct_bg.correct_bed12(df_chunk)
            df_corrected.to_csv(outfile, index=False, header=False, sep='\t', mode=mode)

os.remove('%s/%s'%(output_dir, bamfile.replace('.bam','.bed12')))
correct_bg.genome_cov(output_dir, output, genomepath, ucsc)