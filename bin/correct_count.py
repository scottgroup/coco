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
import distribute_multireads
import subprocess
import distribute_embedded_counts as dist_emb
import pandas as pd
import gc

def read_gtf(gtf_file):
    df_gtf = pd.read_csv(gtf_file, sep='\t',
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score',
                                'strand', 'frame', 'attributes'],
                         dtype={'seqname': str, 'start': int, 'end': int},
                         usecols=['seqname', 'source', 'feature', 'start', 'end', 'attributes'])
    df_gtf['gene_id'] = df_gtf['attributes'].str.extract('gene_id\s"([^;]+);?"', expand=True)
    df_gtf['transcript_id'] = df_gtf['attributes'].str.extract('transcript_id\s"([^;]+);?"', expand=True)
    df_gtf['gene_name'] = df_gtf['attributes'].str.extract('gene_name\s"([^;]+);?"', expand=True)
    df_gtf['transcript_name'] = df_gtf['attributes'].str.extract('transcript_name\s"([^;]+);?"', expand=True)
    df_gtf['gene_biotype'] = df_gtf['attributes'].str.extract('gene_biotype\s"([^;]+);?"', expand=True)
    df_gtf['transcript_biotype'] = df_gtf['attributes'].str.extract('transcript_biotype\s"([^;]+);?"', expand=True)
    df_gtf['transcript_support_level'] = df_gtf['attributes'].str.extract('transcript_support_level\s"([^;]+);?"',
                                                                          expand=True)
    df_gtf['transcript_support_level'] = pd.to_numeric(df_gtf['transcript_support_level'], errors='coerce')
    df_gtf = df_gtf.drop(['attributes'], axis=1)
    return df_gtf


def coco_unique(minOverlap, strand, thread, paired, gtf_file, output_prefix, bamfile, R_opt, feature):
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
              "%s " \
              "-B " \
              "%s" % (minOverlap, strand, thread, paired, gtf_file, output_prefix,R_opt, feature, bamfile)
    x = os.system(command)
    return x

def extract_multi(output_dir, output, bamfile, thread):
    output_name = os.path.basename(output)
    print('Extracting multimapped reads')
    fetch_header = 'samtools view -H %s > %s/multi_%s.sam'%(bamfile, output_dir, output_name)
    x = os.system(fetch_header)
    if x!=0:
        sys.exit('fetch_header exit status %d'%x)
    fetch_multi = "samtools view %s | grep 'NH:i:' | grep -vwE 'NH:i:[10]' >> %s/multi_%s.sam"%(bamfile, output_dir,
                                                                                                output_name)
    x = os.system(fetch_multi)
    if x!=0:
        if x == 256:
            print('Input file does not have multimapped reads (based on NH tag), skipping coco multi')
            return x
        else:
            sys.exit('fetch_multi exit status %d' % x)
    samtobam = 'samtools view -bS -@ %d %s/multi_%s.sam > %s/multi_%s.bam'%(thread, output_dir, output_name, output_dir,
                                                                            output_name)
    x = os.system(samtobam)
    if x!=0:
        sys.exit('samtobam exit status %d'%x)
    os.remove('%s/multi_%s.sam' % (output_dir, output_name))
    return x


def coco_multi(minOverlap, strand, thread, paired, gtf_file, output, infile,
               unique_counts, output_dir, R_opt, v, chunksize, feature,ftype, df_gtf, df_gtf_intron):
    os.chdir(output_dir)
    output_name = os.path.basename(output)
    tests.check_unique_count(unique_counts)
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
            "%s " \
            "-B " \
            "%s" %(minOverlap, strand, thread, paired, gtf_file, 'multi_'+ output_name, R_opt, feature,
                  infile)
    x = os.system(command)
    if x !=0:
        sys.exit(1)
    featurefile = '%s/%s.featureCounts'%(output_dir, os.path.basename(infile))
    distribute_multireads.distribute_multireads(featurefile, unique_counts, R_opt, df_gtf,
                                                output, v, chunksize, thread, ftype,df_gtf_intron)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("annotation", help="annotation file produced with correct_annotation in .gtf format.")
    parser.add_argument("bamfile", help="alignment file in .bam format.")
    parser.add_argument("output", help="Name of the output file holding the counts per genes. ")

    parser.add_argument("-c", "--countType",
                        help="Decide whether to consider only uniquely mapped reads (uniqueOnly) or both uniquely and multimapped reads (both) to produce read count values for genes. Default: both",
                        choices=['uniqueOnly', 'both'], default='both')
    parser.add_argument("-s", "--strand",
                        help="Strandedness: strand-specific read counting. Acceptable values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default: 0",
                        type=int, choices=[0, 1, 2], default=0)
    parser.add_argument("-m", "--minOverlap",
                        help="Minimum overlap a read must have with a feature to be assigned to its associated gene by featureCounts. Default: 10",
                        type=int, default=10)
    parser.add_argument("-i", "--meanInsertSize",
                        help="Mean insert size of the alignment file. Used to compute TPM if the insert size is greater than read length. Default: 0",
                        type=int, default=0)
    parser.add_argument("-p", "--paired",
                        help="Use this option if you work on a paired-end dataset. Used by featureCounts.",
                        action="store_true")
    parser.add_argument("-t", "--thread", help="Number of threads to be used by featureCounts. Default: 1", type=int,
                        default=1)
    parser.add_argument("-r", "--rawOnly",
                        help="Use this option to output the featureCounts raw read counts only and not their calculated CPM and TPM values.",
                        action="store_true")
    parser.add_argument("-R", "--reportreads", help="featureCounts output format (SAM/BAM only available for >=v1.5.3). Default: None",
                                                      choices=['None','CORE','SAM', 'BAM'], type=str, default='None')
    parser.add_argument("-C", "--chunksize",
                        help="only with -R SAM/BAM, number of rows to read from sam/bam file at a time. Default: 1000000",
                        type=int, default=1000000)
    args = parser.parse_args()

    gtf_file = args.annotation
    gtf_file = os.path.abspath(gtf_file)
    bamfile = args.bamfile
    output = args.output
    count_type = args.countType

    strand = args.strand
    minOverlap = args.minOverlap
    meanInsertSize = args.meanInsertSize
    thread = args.thread
    paired = args.paired
    rawonly = args.rawOnly
    R_opt = args.reportreads
    chunksize = args.chunksize

    tests.check_gtf(gtf_file)
    tests.check_bam(bamfile)
    tests.check_output(output)

    vfc = subprocess.Popen(['featureCounts','-v'], stderr=subprocess.PIPE)
    std_out, err = vfc.communicate()
    v = err.decode('utf-8').strip().split(' ')[-1]
    if v < 'v1.5.3':
        if R_opt == 'None':
            R_opt_unique = ''
        else:
            R_opt_unique = '-R'
        R_opt_multi = ''
        R_opt_emb =''
    else:
        if R_opt == 'None':
            R_opt_unique = ''
            R_opt_multi = 'CORE'
        else:
            R_opt_unique = '-R %s'%R_opt
            R_opt_multi = R_opt
        R_opt_emb = 'CORE'
    print('Using featureCounts, version %s'%(str(v)))

    output = os.path.abspath(output)
    output_dir = os.path.dirname(output)

    bamfile = os.path.abspath(bamfile)

    if paired != True:
        paired = ''
    else:
        paired = '-p '

    gtf_file_intron = gtf_file.replace('.gtf','.introns.gtf')
    df_gtf_full = read_gtf(gtf_file)
    df_gtf_intron = read_gtf(gtf_file_intron)

    if count_type =='uniqueOnly':
        output_file = os.path.join(output_dir, os.path.basename(output))
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, output_file, bamfile, R_opt_unique,'')
        if x!=0 :
            sys.exit(1)
        # embedded correction
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file_intron,
                        output_file+'.intron', bamfile, '', '-g transcript_id')
        if x!=0 :
            sys.exit(1)
        dist_emb.correct_embedded(df_gtf_intron, output_file, output_file + '.intron',
                                  output_file + '_final',count_type)


    elif count_type == 'both':
        #For both, default.
        unique_output = output_dir+'/unique_'+os.path.basename(output)
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, unique_output, bamfile, R_opt_unique,'')
        if x != 0:
            sys.exit(1)
        # embedded correction
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file.replace('.gtf','.introns.gtf'),
                        unique_output+'.intron', bamfile, '', '-g transcript_id')
        if x != 0:
            sys.exit(1)

        x = extract_multi(output_dir, output, bamfile, thread)
        multibam = '%s/multi_%s.bam' % (output_dir, os.path.basename(output))
        if x == 0:
            coco_multi(minOverlap, strand, thread, paired, gtf_file_intron,
                       output + '.intron', multibam, unique_output, output_dir, R_opt_emb, v, chunksize,
                       '-g transcript_id', 'intron', df_gtf_full, df_gtf_intron)
            gc.collect()
            intronfile = '%s/multi_%s.bam.featureCounts' % (output_dir, os.path.basename(output))
            newintronfile = '%s/multi_%s.intron.bam.featureCounts' % (output_dir, os.path.basename(output))
            os.rename(intronfile, newintronfile)
            coco_multi(minOverlap, strand, thread, paired, gtf_file, output, multibam,
                       unique_output, output_dir, R_opt_multi, v, chunksize,'','gene', df_gtf_full, None)
            gc.collect()

            os.remove('%s/multi_%s.bam' % (output_dir, os.path.basename(output)))
        else:
            os.rename(unique_output, output)
            os.rename(unique_output+'.summary', output+'.summary')
            os.rename(unique_output+'.intron', output+'.intron')
            os.rename(unique_output+'.intron.summary', output+'.intron.summary')
            count_type = 'uniqueOnly'
        dist_emb.correct_embedded(df_gtf_intron, output, output + '.intron',
                                  output+'_final', count_type)

    else:
        sys.exit(1)

    # doing some cleaning
    os.remove(output + '.intron')
    os.remove(output)
    os.rename(output+'_final', output)

    if count_type == 'uniqueOnly':
        os.remove(output+'.intron.summary')
        os.remove(output + '.summary')

    elif  count_type == 'both':
        os.remove(os.path.join(output_dir,'unique_'+os.path.basename(output))+'.intron.summary')
        os.remove(os.path.join(output_dir, 'unique_'+os.path.basename(output)) + '.summary')
        os.remove(os.path.join(output_dir,'multi_'+os.path.basename(output))+'.intron.summary')
        os.remove(os.path.join(output_dir,'multi_'+ os.path.basename(output)) + '.summary')
        os.remove(os.path.join(output_dir, 'multi_' + os.path.basename(output)) + '.intron')
        os.remove(os.path.join(output_dir, 'unique_' + os.path.basename(output)) + '.intron')
        os.remove(os.path.join(output_dir, 'unique_' + os.path.basename(output)))
        os.remove(os.path.join(output_dir, 'multi_' + os.path.basename(output)))
        if R_opt != 'None':
            os.remove(os.path.join(output_dir, 'multi_' + os.path.basename(output))+'.intron.bam.featureCounts')


    if count_type !='uniqueOnly' and R_opt == 'None':
        genefile = '%s/multi_%s.bam.featureCounts'%(output_dir, os.path.basename(output))
        intronfile = '%s/multi_%s.intron.bam.featureCounts' % (output_dir, os.path.basename(output))
        os.remove(genefile)
        os.remove(intronfile)

    if rawonly!=True:
        count_to_cpm.add_pm_counts(output, df_gtf_full, bamfile,
                                   mean_insert_size=meanInsertSize)
    print('coco cc finished successfully')

if __name__ == '__main__':
    main()
