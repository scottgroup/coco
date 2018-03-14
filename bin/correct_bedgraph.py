#!/usr/bin/env python3
"""
correct_bedgraph.py
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
parser.add_argument("genomepath", help='complete path of the genome file for bedtools genomeCoverageBed.'
                                       'The file must contain two columns separated by a tabulation, the first column '
                                       'is the name of the chromosome and the second is the length of the chromosome.')
parser.add_argument("-m","--contains_multi", help='indicates that your file contains multimapped reads or not '
                                                  '(sorting will be longer, but accurate)', action="store_true")
parser.add_argument("-c", "--chunk_size", help="number of rows to treat at a time. Default: 2500000",
                    type=int, default=2500000 )
parser.add_argument("-u", "--ucsc_compatible", help="Add trackline and 'chr' prefix if needed.", action="store_true")
parser.add_argument("-t", "--thread", help="Number of threads to be used. Default: 1", type=int,
                    default=1)
args = parser.parse_args()



output = args.output
bamfile = args.bamfile
chunk_size = args.chunk_size
ucsc = args.ucsc_compatible
genomepath = args.genomepath
multi = args.contains_multi
thread = args.thread

tests.check_bam(bamfile)
print('output:', output)

output_dir = os.path.dirname(output)
if output_dir == '' or output_dir == '.':
    output_dir = os.getcwd()
print('output_dir:', output_dir)


x = correct_bg.prepare_bed12(bamfile, output_dir,''.join(os.path.basename(output).rsplit('.',1)[:-1]), multi)
temp_dir = ''.join(os.path.basename(output).rsplit('.',1)[:-1])
if x !=0 :
    sys.exit('prepare_bed12 exit status: '+ str(x))

filelist = os.listdir('%s/%s_chromo/'%(output_dir,temp_dir))

for filename in filelist:

    if filename.endswith(".bed12") and 'corrected' not in filename:
        print('running chromosome: %s'%(filename.replace('.bed12','')))
        outfile = output_dir + '/'+temp_dir+'_chromo/corrected_' + os.path.basename(filename)

        for i, df_chunk in enumerate(pd.read_csv(output_dir + '/'+temp_dir+'_chromo/'+ filename, sep='\t',
                                    names=['chromo', 'start', 'end', 'name', 'score', 'strand', 'thick_start',
                                           'thick_end','RGB', 'nb_block', 'block_len', 'block_start'],
                                    dtype={'chromo': str}, chunksize=chunk_size)):
            if i == 0 :
                mode = 'w'
            else :
                mode = 'a'
            df_corrected = correct_bg.correct_bed12(df_chunk, thread)
            df_corrected.to_csv(outfile, index=False, header=False, sep='\t', mode=mode)
os.remove('%s/%s.bed12'%(output_dir,temp_dir))
correct_bg.genome_cov(output_dir,temp_dir, output, genomepath, ucsc)
