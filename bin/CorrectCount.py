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



def coco_unique(minOverlap, strand, thread, paired, gtf_file, output_prefix, bamfile):
    command = "featureCounts " \
              "--minOverlap %d " \
              "--largestOverlap " \
              "-s %d " \
              "-C " \
              "-T %s " \
              "%s" \
              "-a %s " \
              "-o %s " \
              "-B " \
              "%s" % (minOverlap, strand, thread, paired, gtf_file, output_prefix, bamfile)
    x = os.system(command)
    return x


def coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile, unique_counts, output_dir):
    if not unique_counts:
        sys.exit('multiOnly requires the option -u to be set (count matrix from uniquely mapped reads)')
    output_name = os.path.basename(output)
    tests.check_unique_count(unique_counts)
    print('Extracting multimapped reads')
    fetch_header = 'samtools view -H %s > %s/multi_%s.sam && '%(bamfile, output_dir, output_name)
    fetch_multi = "samtools view %s | grep 'NH:i:' | grep -vw 'NH:i:1' >> %s/multi_%s.sam && "%(bamfile, output_dir, output_name)
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
            "-R " \
            "-B " \
            "%s" %(minOverlap, strand, thread, paired, gtf_file, 'multi_'+ output_name,
                   output_dir+'/multi_'+output_name+'.bam')
    x = os.system(command)
    if x !=0:
        sys.exit(x)
    group_reads = 'perl %s/count_from_bed.pl %s/%s.multi_%s.bam.featureCounts > %s/%s_grouped.txt'%(os.path.dirname(__file__),
                                                                                        output_dir,
                                                                                        output_dir.strip('/').replace('/','.'),
                                                                                        output_name,
                                                                                        output_dir,
                                                                                        output_name)
    os.system(group_reads)
    output_file = output_dir+'/'+output_name
    distribute_counts.ratio_mmg(output_name +'_grouped.txt', gtf_file, unique_counts, output_file)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("annotation", help="annotation file produced with CorrectAnnotation in .gtf format.")
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

    tests.check_gtf(gtf_file)
    tests.check_bam(bamfile)
    tests.check_output(output)
    print('okay')
    output_dir = os.path.dirname(output)
    if output_dir == '':
        output_dir = os.getcwd()

    if paired != True:
        paired = ''
    else:
        paired = '-p '

    if count_type =='uniqueOnly':
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, output, bamfile)
        if x!=0 :
            sys.exit(x)

    elif count_type=='multiOnly':
        coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile, unique_counts, output_dir)

    else:
        #For both, default.
        unique_output = output_dir+'/unique_'+os.path.basename(output)
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, unique_output, bamfile)
        if x != 0:
            sys.exit(x)
        coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile, unique_output, output_dir)


    if rawonly!=True:
        count_to_cpm.add_pm_counts(output, gtf_file, bamfile)


if __name__ == '__main__':
    main()