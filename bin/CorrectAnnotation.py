import argparse
import sys, os

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file in .gtf format.")

parser.add_argument("-o", "--output", help="Name of output gtf. Default: add 'CorrectAnnotation' suffix to original gtf name.", default='None')
parser.add_argument("-c", "--countType", help="Decide whether to consider only uniquely mapped reads (uniqueOnly), only multi-mapped reads (multiOnly) or both (both) to produce read count values for genes. Default: both", choices=['uniqueOnly','multiOnly','both'], default='both')
parser.add_argument("-s", "--strand",help="Strandness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default: 0", choices=[0, 1, 2],default=0)
parser.add_argument("-p", "--paired", help="Use this option if you work on paired-end dataset. Used by featureCounts.", action="store_true")
parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int, default=1)
args = parser.parse_args()

gtf=args.annotation
output=args.output
print('output', output)

command="python3 %sgapped_gtf.py %s %s" %((os.path.realpath(__file__).replace('CorrectAnnotation.py','')),gtf,output)
os.system(command)
