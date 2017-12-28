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
import distribute_embedded_counts as dist_emb
import gtf


def read_gtf(gtf_file):
    df_gtf = gtf.dataframe(gtf_file)
    df_gtf = df_gtf[
        ['seqname', 'source', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'exon_number',
         'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype', 'transcript_support_level']]
    df_gtf['seqname'] = df_gtf['seqname'].map(str)
    df_gtf['start'] = df_gtf['start'].map(int)
    df_gtf['end'] = df_gtf['end'].map(int)
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

def extract_multi(output_dir, output, bamfile):
    output_name = os.path.basename(output)
    print('Extracting multimapped reads')
    fetch_header = 'samtools view -H %s > %s/multi_%s.sam'%(bamfile, output_dir, output_name)
    x = os.system(fetch_header)
    if x!=0:
        sys.exit('fetch_header exit status %d'%x)
    fetch_multi = "samtools view %s | grep 'NH:i:' | grep -vw 'NH:i:1' >> %s/multi_%s.sam"%(bamfile, output_dir,
                                                                                                output_name)
    x = os.system(fetch_multi)
    if x!=0:
        if x == 256:
            print('Input file does not have multimapping read, skipping coco multi')
            return x
        else:
            sys.exit('fetch_multi exit status %d' % x)
    samtobam = 'samtools view -b %s/multi_%s.sam > %s/multi_%s.bam'%(output_dir,output_name, output_dir,output_name)
    print(samtobam)
    x = os.system(samtobam)
    if x!=0:
        sys.exit('samtobam exit status %d'%x)
    os.remove('%s/multi_%s.sam' % (output_dir, output_name))
    return x


def coco_multi(minOverlap, strand, thread, paired, gtf_file, output, infile,
               unique_counts, output_dir, R_opt, v, chunksize, feature,ftype, df_gtf, df_gtf_intron):
    os.chdir(output_dir)
    output_name = os.path.basename(output)
    if not unique_counts:
        sys.exit('multiOnly requires the option -u to be set (count matrix from uniquely mapped reads)')
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
        sys.exit(x)
    output_file = output_dir+'/'+output_name
    featurefile = '%s/%s.featureCounts'%(output_dir, os.path.basename(infile))
    distribute_multireads.distribute_multireads(featurefile, unique_counts, R_opt, df_gtf,
                                                output_file, v, chunksize, thread, ftype,df_gtf_intron)


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
        R_opt_emb =''
    else:
        if R_opt == 'None':
            R_opt_unique = ''
            R_opt_multi = 'CORE'
        else:
            R_opt_unique = '-R %s'%R_opt
            R_opt_multi = R_opt
        R_opt_emb = 'CORE'
    print('Using, version %s, -R %s'%(str(v),R_opt))

    output_dir = os.path.dirname(output)
    if output_dir == '' or output_dir == '.':
        output_dir = os.getcwd()

    bam_dir = os.path.dirname(bamfile)
    if bam_dir == '' or bam_dir == '.':
        bam_dir = os.getcwd()
        bamfile = os.path.join(bam_dir, bamfile)

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
            sys.exit(x)
        # embedded correction
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file_intron,
                        output_file+'.intron', bamfile, '', '-g transcript_id')
        if x!=0 :
            sys.exit(x)
        #dist_emb....
        dist_emb.correct_embedded(df_gtf_intron, output_file, output_file + '.intron',
                                  output_file + '_final',count_type)


    elif count_type=='multiOnly':
        x = extract_multi(output_dir, output, bamfile)
        if x == 0:
            coco_multi(minOverlap, strand, thread, paired, gtf_file, output, bamfile,
                       unique_counts, output_dir, R_opt_multi, v, chunksize,'','gene', df_gtf_full, None)
            coco_multi(minOverlap, strand, thread, paired, gtf_file.replace('.gtf','.introns.gtf'),
                       output+'.intron', bamfile, unique_counts, output_dir, R_opt_emb, v, chunksize,'-g transcript_id',
                       'intron',df_gtf_full, df_gtf_intron)
            os.remove('%s/multi_%s.bam' % (output_dir, os.path.basename(output)))
            dist_emb.correct_embedded(df_gtf_intron, output, output + '.intron',
                                      output_dir + os.path.basename(output) + '_final', count_type)


    elif count_type == 'both':
        #For both, default.
        unique_output = output_dir+'/unique_'+os.path.basename(output)
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file, unique_output, bamfile, R_opt_unique,'')
        if x != 0:
            sys.exit(x)
        # embedded correction
        x = coco_unique(minOverlap, strand, thread, paired, gtf_file.replace('.gtf','.introns.gtf'),
                        unique_output+'.intron', bamfile, '', '-g transcript_id')
        if x != 0:
            sys.exit(x)

        x = extract_multi(output_dir, output, bamfile)
        multibam = '%s/multi_%s.bam' % (output_dir, os.path.basename(output))
        if x == 0:
            coco_multi(minOverlap, strand, thread, paired, gtf_file_intron,
                       output + '.intron', multibam, unique_output, output_dir, R_opt_emb, v, chunksize,
                       '-g transcript_id', 'intron', df_gtf_full, df_gtf_intron)
            intronfile = '%s/multi_%s.bam.featureCounts' % (output_dir, os.path.basename(output))
            newintronfile = '%s/multi_%s.intron.bam.featureCounts' % (output_dir, os.path.basename(output))
            os.rename(intronfile, newintronfile)
            coco_multi(minOverlap, strand, thread, paired, gtf_file, output, multibam,
                       unique_output, output_dir, R_opt_multi, v, chunksize,'','gene', df_gtf_full, None)

            os.remove('%s/multi_%s.bam' % (output_dir, os.path.basename(output)))
        else:
            os.rename(unique_output, os.path.join(output_dir,os.path.basename(output)))
            os.rename(unique_output+'.intron', os.path.join(output_dir, os.path.basename(output)+'.intron'))
            count_type = 'uniqueOnly'
        dist_emb.correct_embedded(df_gtf_intron, output, output + '.intron',
                                  os.path.join(output_dir,os.path.basename(output)+'_final'), count_type)

    else:
        sys.exit(1)

    if count_type !='uniqueOnly' and R_opt == 'None':
        genefile = '%s/multi_%s.bam.featureCounts'%(output_dir, os.path.basename(output))
        intronfile = '%s/multi_%s.intron.bam.featureCounts' % (output_dir, os.path.basename(output))
        os.remove(genefile)
        os.remove(intronfile)

    if rawonly!=True:
        count_to_cpm.add_pm_counts(output, df_gtf_full, bamfile, count_type)
    print('coco cc finished successfully')

if __name__ == '__main__':
    main()
