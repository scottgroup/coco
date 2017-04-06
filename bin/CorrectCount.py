import argparse
import sys, os

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file produced with CorrectAnnotation in .gtf format.")
parser.add_argument("bam_file", help="alignment file in .bam format.")
parser.add_argument("output_prefix", help="prefix for the output.")

parser.add_argument("-c", "--countType", help="Decide whether to consider only uniquely mapped reads (uniqueOnly), only multi-mapped reads (multiOnly) or both (both) to produce read count values for genes. Default: both", choices=['uniqueOnly','multiOnly','both'], default='both')
parser.add_argument("-s", "--strand",help="Strandedness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default", choices=[0, 1, 2],default=0)
parser.add_argument("-p", "--paired", help="Use this option if you work on paired-end dataset. Used by featureCounts.", action="store_true")
parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int, default=1)
args = parser.parse_args()

gtf=args.gtf
bam_file=args.bam_file
output_prefix=args.output_prefix
count_type=args.count_type

strand=args.strand
thread=args.thread
paired=args.paired


if count_type =='uniqueOnly':
    command="python3 %sfeatureCounts.sh %s %s %s %s %s" %((os.path.realpath(__file__).replace('CorrectAnnotation.py','')),bam_file,gtf,output_prefix,str(strand),thread)
    if paired == True:
        command+=' -p'
    os.system(command)
elif count_type=='multiOnly':
    print('Do the thing Gabe')
    print('...')
else:
    #For both
    print('Do the thing Gabe')
    print('...')
