import pandas as pd
import gtf
import sys
import math
import os


def distribute_counts(df, df_unique):
    # df contains the multimapped reads read assignments having at least columns 'gene_id' and 'QNAME'
    # df_unique is the count matrix of uniquely mapped reads, having at least columns 'gene_id' and 'accumulation'
    df_merged = pd.merge(df, df_unique, on='gene_id')
    df_merged['dist_counts'] = -1
    df_merged['group_tot'] = df_merged.groupby('QNAME')['accumulation'].transform('sum')
    df_merged['group_len'] = df_merged.groupby('QNAME')['accumulation'].transform('size')
    df_merged.loc[df_merged.group_tot == 0, 'dist_counts'] = 1.0 / df_merged.group_len
    df_merged.loc[df_merged.group_tot != 0, 'dist_counts'] = df_merged.accumulation / df_merged.group_tot
    return df_merged


def read_count_matrix(count_file):
    df= pd.read_csv(count_file, sep='\t', comment='#',
                    names=['gene_id', 'seqname', 'start', 'end', 'strand', 'length', 'accumulation'])
    df = df.drop(['seqname', 'start', 'end', 'strand', 'length'], axis=1)
    df['accumulation'] = pd.to_numeric(df['accumulation'],
                                       errors='coerce')  # if header in file, column name will be changed to NaN,
    df = df.dropna(axis=0)  # and row will be removed
    return df


def read_corefile(corefile, colnames):
    # colnames must contain 'QNAME','gene_id','status' and'nb_targets'
    df = pd.read_csv(corefile, names=colnames, sep='\t')
    df = df[df['status']=='Assigned']
    df = df.drop(['status','nb_targets'], axis=1)
    return df


def unique_counts(gtf_file, unique_file):
    if gtf_file.endswith('.gtf'):
        print('Reading gtf')
    else:
        print('Annotation file:',gtf_file)
        print('error: Wrong annotation format. Only .gtf files are accepted. Filename must contain ".gtf" extension.')
        sys.exit(1)
    try:
        df_gtf=gtf.dataframe(gtf_file)
        df_gtf = df_gtf[
            ['seqname', 'source', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'exon_number',
             'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype', 'transcript_support_level']]
        df_gtf['seqname'] = df_gtf['seqname'].map(str)
        df_gtf['start'] = df_gtf['start'].map(int)
        df_gtf['end'] = df_gtf['end'].map(int)
    except:
        print("Error: gtf file cannot be converted to dataframe. Make sure the annotation file provided is in gene"
              "transfer format (.gtf) and respects Ensembl's format.")
        sys.exit(1)
    df_gtf = df_gtf[df_gtf.feature == 'gene']
    df_gtf = df_gtf.drop(['seqname', 'source', 'feature', 'start', 'end', 'strand', 'transcript_id', 'exon_number',
                          'transcript_name', 'transcript_biotype', 'transcript_support_level','gene_biotype',
                          'gene_name'], axis=1)
    df_unique = read_count_matrix(unique_file)
    df_unique = df_unique.merge(df_gtf, on='gene_id', how='outer')
    df_unique = df_unique.fillna(0.0)
    df_unique = df_unique.set_index('gene_id')
    df_unique = df_unique.reset_index()
    return df_unique


def distribute_samfile(samfile, chunksize, df_unique, nb_threads, R_opt):
    ext = R_opt.lower()
    fetch_header = 'samtools view -H %s.%s > %s.out.sam && '%(samfile,ext,samfile)
    sort_sam = 'samtools sort -n -o %s.sorted.%s -@ %s %s.%s && '%(samfile, ext, nb_threads, samfile, ext)
    view_sam = 'samtools view %s.sorted.%s > %s.noheader.sam'%(samfile, ext, samfile)
    x = os.system(fetch_header + sort_sam + view_sam)
    if x !=0:
        sys.exit(x)
    sortedsam = '%s.sorted.%s'%(samfile,ext)
    os.remove(sortedsam)
    sam = samfile+'.noheader.sam'
    with open(sam, 'rb') as f:
        file_len = sum((1 for line in f))
    print(file_len)

    total_chunk = math.floor(float(file_len) / float(chunksize))
    print('total_chunk %d' % total_chunk)
    df_last = pd.DataFrame()
    df_count = pd.DataFrame()
    for n, df_chunk in enumerate(pd.read_csv(sam, sep='\t',
                                             names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
                                                    'TLEN', 'SEQ', 'QUAL', 'TAG1', 'TAG2', 'TAG3', 'TAG4', 'TAG5',
                                                    'TAG6','TAG7','TAG8','TAG9','TAG10','TAG11'],
                                             dtype={'RNAME':str}, chunksize=chunksize)):

        print('chunk %d' % n)
        # merge optional tags in one column
        df_chunk['TAGS'] = df_chunk[df_chunk.columns[14:]].apply(lambda x: '|'.join(x.dropna().astype(str)), axis=1)
        df_chunk = df_chunk.drop(['TAG4', 'TAG5', 'TAG6', 'TAG7', 'TAG8', 'TAG9', 'TAG10', 'TAG11'], axis=1)
        # extract XT:Z: tag if present
        df_chunk['gene_id'] = df_chunk['TAGS'].apply(lambda x: x.split('|')[-1].strip().replace('XT:Z:', '') if 'XT:Z:' in x else 'NaN')

        # samfile must be sorted by name, keep the last name for the next chunk, if the read has multiple alignments
        # that are in different chunks. Ensures all reads are only counted once.
        last_name = df_chunk.tail(1).QNAME.values[0]
        df_chunk = pd.concat([df_chunk, df_last])
        if n != total_chunk:
            df_last = df_chunk[df_chunk.QNAME == last_name].copy(deep=True)
            df_chunk = df_chunk[df_chunk.QNAME != last_name]

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
    if R_opt == 'BAM':
        samtobam = 'samtools view -b -@ %d %s.out.sam > %s.out.bam && rm %s.out.sam'%(nb_threads, samfile, samfile,
                                                                                     samfile)
        x = os.system(samtobam)
        if x!=0:
            sys.exit('sam to bam exited with status : %d'%x)
    return df_count

def distribute_multireads(featurefile, uniquefile, R_opt, gtf_file, output_file, v_subread, chunksize, nb_threads):
    df_unique = unique_counts(gtf_file, uniquefile)
    if R_opt=='SAM' or R_opt =='BAM':
        df_dist = distribute_samfile(featurefile, chunksize, df_unique, nb_threads, R_opt)
    else:
        if v_subread >= 'v1.5.3':
            colnames = ['QNAME','status','nb_targets','gene_id']
        else:
            colnames = ['QNAME', 'status', 'gene_id', 'nb_targets']
        print(colnames)
        df_multi = read_corefile(featurefile, colnames)
        df_dist = distribute_counts(df_multi, df_unique)
    df_dist = df_dist.fillna(0)

    df_group = df_dist.groupby('gene_id').dist_counts.sum().reset_index()
    print('merge df multi unique')
    df_tot = pd.merge(df_group,df_unique, on='gene_id')
    df_tot['tot'] = df_tot.accumulation + df_tot.dist_counts
    df_tot[['gene_id','tot']].to_csv(output_file, header=False, index=False, sep='\t')
