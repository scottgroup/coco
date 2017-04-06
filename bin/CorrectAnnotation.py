import argparse
import sys, os

parser = argparse.ArgumentParser()

parser.add_argument("annotation", help="annotation file in .gtf format.")

parser.add_argument("-o", "--output", help="Name of output gtf. Default: add 'CorrectAnnotation' suffix to original gtf name.", default='None')
args = parser.parse_args()

gtf=args.annotation
output=args.output
print('output', output)

command="python3 %sgapped_gtf.py %s %s" %((os.path.realpath(__file__).replace('CorrectAnnotation.py','')),gtf,output)
os.system(command)
