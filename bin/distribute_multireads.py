import pandas as pd
import sys
import math
import os
import gc


def distribute_counts(df, df_unique):
    '''
    Distributes the multimapped counts based on the uniquely mapped reads ratio,
    if none of genes have uniquely mappped reads, the counts are equally distributed within the group
    :param df: contains the multimapped reads read assignments having at least columns 'gene_id' and 'QNAME'
    :param df_unique: the count matrix of uniquely mapped reads, having at least columns 'gene_id' and 'accumulation'
    :return: a merged matrix having the distributed counts
    '''

    df_merged = pd.merge(df, df_unique, on='gene_id')
    df_merged['group_tot'] = df_merged.groupby('QNAME')['accumulation'].transform('sum')
    df_merged['group_len'] = df_merged.groupby('QNAME')['accumulation'].transform('size')
    df_merged.loc[df_merged.group_tot == 0, 'dist_counts'] = 1.0 / df_merged.group_len
    df_merged.loc[df_merged.group_tot != 0, 'dist_counts'] = df_merged.accumulation / df_merged.group_tot
    df_merged = df_merged.drop(['group_tot','group_len'],axis=1)
    return df_merged


def read_count_matrix(count_file, ftype):
    '''
    Reads the uniquely mapped reads matrix obtained by featureCounts
    :param count_file: path to the featureCounts matrix
    :return: df having the uniquely mapped counts
    '''
    if ftype=='intron':
        use_columns = ['gene_id', 'length', 'accumulation']
    else:
        use_columns = ['gene_id', 'accumulation']
    df= pd.read_csv(count_file, sep='\t', comment='#',
                    names=['gene_id', 'seqname', 'start', 'end', 'strand', 'length', 'accumulation'],
                    usecols = use_columns)
    df['accumulation'] = pd.to_numeric(df['accumulation'],
                                       errors='coerce')  # if header in file, column name will be changed to NaN,
    df = df.dropna(axis=0)  # and row will be removed
    return df


def read_corefile(corefile, v_subread):
    if v_subread >= 'v1.5.3':
        colnames = ['QNAME', 'status', 'nb_targets', 'gene_id']
    else:
        colnames = ['QNAME', 'status', 'gene_id', 'nb_targets']
    df = pd.read_csv(corefile, names=colnames,dtype={'status':str, 'gene_id':str}, sep='\t',
                     usecols=['QNAME', 'status', 'gene_id'])
    df = df[df['status']=='Assigned']
    df = df.drop(['status'], axis=1)
    return df


def unique_counts(df_gtf, unique_file, ftype):
    df_gtf = df_gtf[df_gtf.feature == 'gene'][['gene_id']]
    df_unique = read_count_matrix(unique_file, ftype)
    df_unique = df_unique.merge(df_gtf, on='gene_id', how='outer')
    df_unique = df_unique.fillna(0.0)
    df_unique = df_unique.set_index('gene_id')
    df_unique = df_unique.reset_index()
    return df_unique


def distribute_samfile(samfile, chunksize, df_unique, nb_threads, R_opt):
    '''
    Distribute the multimapped counts based on the uniquely mapped reads ratio and report the ratio in the
    alignment file (SAM or BAM)
    :param samfile: name of the sam/bam file
    :param chunksize: number of rows to treat at a time
    :param df_unique: matrix containing the uniquely mapped counts
    :param nb_threads: number of threads to use for sorting
    :param R_opt: either SAM or BAM
    :return: a matrix having the distributed counts. Writes the distributed ratio in the SAM/BAM file
    '''
    ext = R_opt.lower()
    fetch_header = 'samtools view -H %s.%s > %s.out.sam && '%(samfile,ext,samfile)
    sort_sam = 'samtools sort -n -o %s.sorted.%s -@ %s %s.%s && '%(samfile, ext, nb_threads, samfile, ext)
    view_sam = 'samtools view %s.sorted.%s > %s.noheader.sam'%(samfile, ext, samfile)
    x = os.system(fetch_header + sort_sam + view_sam)
    if x !=0:
        sys.exit(x)
    sortedsam = '%s.sorted.%s'%(samfile,ext)
    os.remove(sortedsam)
    sam_noheader = samfile+'.noheader.sam'
    with open(sam_noheader, 'rb') as f:
        file_len = sum((1 for line in f))
    print(file_len)

    total_chunk = math.floor(float(file_len) / float(chunksize))
    print('total_chunk %d' % total_chunk)
    df_last = pd.DataFrame()
    df_count = pd.DataFrame()
    for n, df_chunk in enumerate(pd.read_csv(sam_noheader, sep='\t',
                                             names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
                                                    'TLEN', 'SEQ', 'QUAL', 'TAG1', 'TAG2', 'TAG3', 'TAG4', 'TAG5',
                                                    'TAG6','TAG7','TAG8','TAG9','TAG10','TAG11'],
                                             dtype={'RNAME':str}, chunksize=chunksize)):

        # samfile must be sorted by name, keep the last name for the next chunk, if the read has multiple alignments
        # that are in different chunks. Ensures all reads are only counted once.
        last_name = df_chunk.tail(1).QNAME.values[0]
        df_chunk = pd.concat([df_chunk, df_last])
        if n != total_chunk:
            df_last = df_chunk[df_chunk.QNAME == last_name].copy(deep=True)
            df_chunk = df_chunk[df_chunk.QNAME != last_name]

        print('chunk %d' % n)
        df_chunk = df_chunk.fillna('')
        df_chunk['TAGS'] = df_chunk['TAG4'] + '|' + df_chunk['TAG5'] + '|' + df_chunk['TAG6'] + '|' + df_chunk['TAG7'] \
                           + '|' + df_chunk['TAG8'] + '|' + df_chunk['TAG9'] + '|' + df_chunk['TAG10'] + '|' \
                           + df_chunk['TAG11']
        df_chunk['TAGS'] = df_chunk['TAGS'].str.rstrip('|')
        df_chunk = df_chunk.drop(['TAG4', 'TAG5', 'TAG6', 'TAG7', 'TAG8', 'TAG9', 'TAG10', 'TAG11'], axis=1)
        df_chunk['gene_id'] = df_chunk['TAGS'].str.split('|').str[-1].str.replace('XT:Z:','')
        df_chunk.loc[df_chunk.gene_id.str.contains('XS:Z:'),'gene_id'] = 'NaN'
        df_merged = distribute_counts(df_chunk[['QNAME','FLAG','TAG2','gene_id']], df_unique)
        df_chunk = df_chunk.drop(['gene_id'], axis=1)
        df_merged = df_merged[['QNAME','FLAG','TAG2','gene_id', 'dist_counts']]
        df_count = pd.concat([df_count, df_merged[['gene_id','dist_counts']]])
        df_merged['dist_counts'] = 'XC:f:'+df_merged['dist_counts'].round(decimals=3).map(str)
        df_sam = pd.merge(df_chunk, df_merged[['QNAME','FLAG','TAG2','dist_counts']], on=['QNAME','FLAG','TAG2'],
                          how='outer')
        del df_chunk, df_merged
        df_sam.loc[df_sam.dist_counts.isnull(), 'dist_counts']='XC:f:0.000'
        df_sam.to_csv('%s.out.sam'%samfile, sep='\t',
                         index=False, header=False, mode='a')
        del df_sam
    sed = "sed -i 's/|/\t/g' %s.out.sam"%samfile
    x = os.system(sed)
    if x!=0:
        sys.exit(x)
    os.remove(samfile+'.'+ext)
    os.remove(sam_noheader)
    if R_opt == 'BAM':
        samtobam = 'samtools view -b -@ %d %s.out.sam > %s.out.bam && rm %s.out.sam'%(nb_threads, samfile, samfile,
                                                                                     samfile)
        x = os.system(samtobam)
        if x!=0:
            sys.exit('sam to bam exited with status : %d'%x)
    return df_count

def distribute_multireads(featurefile, uniquefile, R_opt, df_gtf, output_file, v_subread, chunksize, nb_threads, ftype,
                          gtf_intron):
    df_unique = unique_counts(df_gtf, uniquefile, ftype)
    output_col = ['gene_id','tot']
    if R_opt=='SAM' or R_opt =='BAM':
        df_dist = distribute_samfile(featurefile, chunksize, df_unique, nb_threads, R_opt)
    else:
        df_multi = read_corefile(featurefile, v_subread)
        if ftype == 'intron':
            output_col = ['transcript_id','tot','length']
            df_multi = df_multi.rename(columns={'gene_id':'transcript_id'})
            gtf_intron = gtf_intron[gtf_intron.feature=='transcript']
            df_multi = df_multi.merge(gtf_intron[['gene_id','transcript_id']], on='transcript_id')
        df_dist = distribute_counts(df_multi, df_unique)
        del df_multi
        gc.collect()
    df_dist = df_dist.fillna(0)
    if ftype=='intron':
        del df_unique
        df_group = df_dist.groupby('transcript_id').dist_counts.sum().reset_index()
        df_unique_intron = read_count_matrix(uniquefile+'.intron',)
        df_unique_intron = df_unique_intron.rename(columns={'gene_id':'transcript_id'})
        df_tot = pd.merge(df_group, df_unique_intron, on='transcript_id', how='outer')
    else :
        df_group = df_dist.groupby('gene_id').dist_counts.sum().reset_index()
        df_tot = pd.merge(df_group,df_unique, on='gene_id', how='outer')
    df_tot = df_tot.fillna(0)
    df_tot['tot'] = df_tot.accumulation + df_tot.dist_counts
    df_tot[output_col].to_csv(output_file, header=False, index=False, sep='\t')
