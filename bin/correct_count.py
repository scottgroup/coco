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
import distribute_counts
import distribute_multireads
import subprocess



def coco_unique(minOverlap, strand, thread, paired, gtf_file, output_prefix, bamfile, R_opt):
    command = "featureCounts " \
              "--minOverlap %d " \
              "--largestOverlap " \
              "-s %d " \
              "-C " \
              "-T %s " \
              "%s" \
              "-a %s " \
              "-o %s " \
              "%s " \
              "-B " \
              "%s" % (minOverlap, strand, thread, paired, gtf_file, output_prefix,R_opt, bamfile)
    x = os.system(command)
    return x


def coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile, unique_counts, output_dir, R_opt, v, chunksize):
    os.chdir(output_dir)
    if not unique_counts:
        sys.exit('multiOnly requires the option -u to be set (count matrix from uniquely mapped reads)')
    output_name = os.path.basename(output)
    tests.check_unique_count(unique_counts)
    print('Extracting multimapped reads')
    fetch_header = 'samtools view -H %s > %s/multi_%s.sam && '%(bamfile, output_dir, output_name)
    fetch_multi = "samtools view %s | grep 'NH:i:' | grep -vw 'NH:i:1' >> %s/multi_%s.sam && "%(bamfile, output_dir,
                                                                                                output_name)
    samtobam = 'samtools view -b %s/multi_%s.sam > %s/multi_%s.bam &&'%(output_dir,output_name, output_dir,output_name)
    rm_multisam = 'rm %s/multi_%s.sam'%(output_dir,output_name)
    x = os.system(fetch_header+fetch_multi+samtobam+ rm_multisam)
    if x !=0 :
        sys.exit('extracting multi exit status: %s'%(str(x)))
    command="featureCounts " \
            "--minOverlap %d " \
            "--largestOverlap " \
            "-s %d " \
            "-C " \
            "-T %s " \
            "%s" \
            "-a %s " \
            "-o %s "\
            "-M "\
            "-R %s " \
            "-B " \
            "%s" %(minOverlap, strand, thread, paired, gtf_file, 'multi_'+ output_name, R_opt,
                  'multi_'+output_name+'.bam')
    x = os.system(command)
    if x !=0:
        sys.exit(x)
    output_file = output_dir+'/'+output_name
    featurefile = '%s/multi_%s.bam.featureCounts'%(output_dir, output_name)
    distribute_multireads.distribute_multireads(featurefile, unique_counts, R_opt, gtf_file,
                                                output_file, v, chunksize, thread)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("annotation", help="annotation file produced with correct_annotation in .gtf format.")
    parser.add_argument("bamfile", help="alignment file in .bam format.")
    parser.add_argument("output", help="Name of the output file holding the counts per genes. ")

    parser.add_argument("-c", "--countType",
                        help="Decide whether to consider only uniquely mapped reads (uniqueOnly), only multi-mapped reads (multiOnly) or both (both) to produce read count values for genes. Default: both",
                        choices=['uniqueOnly', 'multiOnly', 'both'], default='both')
    parser.add_argument("-s", "--strand",
                        help="Strandedness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default: 0",
                        type=int, choices=[0, 1, 2], default=0)
    parser.add_argument("-m", "--minOverlap",
                        help="Minimum overlap a read must have with a feature to be assigned to its associated gene by featureCounts. Default: 10",
                        type=int, default=10)
    parser.add_argument("-p", "--paired",
                        help="Use this option if you work on a paired-end dataset. Used by featureCounts.",
                        action="store_true")
    parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int,
                        default=1)
    parser.add_argument("-r", "--rawOnly",
                        help="Use this option to output the featureCounts raw read counts only and not their calculated CPM and TPM values.",
                        action="store_true")
    parser.add_argument("-u", "--unique_counts", help="Path to uniquely mapped reads count matrix with following format: "
                                                      "gene_id  seqname start   end strand  length  accumulation (separated by a tabulation). "
                                                      "Required for multiOnly", type=str)
    parser.add_argument("-R", "--reportreads", help="featureCounts output format (SAM/BAM only available for >=v1.5.3)",
                                                      choices=['None','CORE','SAM', 'BAM'], type=str, default='None')
    parser.add_argument("-C", "--chunksize",
                        help="only with -R SAM/BAM, number of rows to read from sam/bam file at a time",
                        type=int, default=1000000)
    args = parser.parse_args()

    gtf_file = args.annotation
    bamfile = args.bamfile
    output = args.output
    count_type = args.countType
    unique_counts = args.unique_counts

    strand = args.strand
    minOverlap = args.minOverlap
    thread = args.thread
    paired = args.paired
    rawonly = args.rawOnly
    R_opt = args.reportreads
    chunksize = args.chunksize

    tests.check_gtf(gtf_file)
    tests.check_bam(bamfile)
    tests.check_output(output)

    v = subprocess.run(['featureCounts','-v'],
                       stderr=subprocess.PIPE).stderr.decode('utf-8').strip().split(' ')[-1]
    if v < 'v1.5.3':
        if R_opt == 'None':
            R_opt_unique = ''
        else:
            R_opt_unique = '-R'
        R_opt_multi = ''
    else:
        if R_opt == 'None':
            R_opt_unique = ''
            R_opt_multi = 'CORE'
        else:
            R_opt_unique = '-R %s'%R_opt
            R_opt_multi = R_opt
    print('Using, version %s, -R %s'%(str(v),R_opt))

    output_dir = os.path.dirname(output)
    if output_dir == '':
        output_dir = os.getcwd()

    if paired != True:
        paired = ''
    else:
        paired = '-p '

    if count_type =='uniqueOnly':
        output_file = output_dir + '/' + os.path.basename(output)
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, output_file, bamfile, R_opt_unique)
        if x!=0 :
            sys.exit(x)

    elif count_type=='multiOnly':
        coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile,
                   unique_counts, output_dir, R_opt_multi, v, chunksize)


    elif count_type == 'both':
        #For both, default.
        unique_output = output_dir+'/unique_'+os.path.basename(output)
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, unique_output, bamfile, R_opt_unique)
        if x != 0:
            sys.exit(x)
        coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile,
                   unique_output, output_dir, R_opt_multi, v, chunksize)

    else:
        sys.exit(1)

    if count_type !='uniqueOnly' and R_opt == 'None':
        featurefile = '%s/multi_%s.bam.featureCounts'%(output_dir, os.path.basename(output))
        os.remove(featurefile)

    if rawonly!=True:
        count_to_cpm.add_pm_counts(output, gtf_file, bamfile, count_type)
    print('coco cc finished successfully')

if __name__ == '__main__':
    main()
