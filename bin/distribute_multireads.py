import pandas as pd
import gtf
import sys


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
                                              errors='coerce')  # if header is present in file, column name will be changed to NaN,
    df = df.dropna(axis=0)  # and row will be removed
    return df


def read_corefile(corefile, colnames):
    # colnames must contain 'QNAME','gene_id','status' and'nb_targets'
    df = pd.read_csv(corefile, names=colnames, sep='\t')
    df = df[df['status']=='Assigned']
    df = df.drop(['status','nb_targets'], axis=1)
    return df


def read_samfile(samfile):
    df = pd.read_csv(samfile, sep='\t',
                     names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL',
                            'TAG1', 'TAG2', 'TAG3', 'TAG4', 'TAG5', 'TAG6','TAG7','TAG8','TAG9','TAG10','TAG11'],
                     dtype={'RNAME':str})
    # merge optional tags in one column
    df['TAGS'] = df[df.columns[14:]].apply(lambda x: '\t'.join(x.dropna().astype(str)), axis=1)
    # extract XT:Z: tag if present
    df['gene_id'] = df['TAGS'].apply(lambda x: x.split('\t')[-1].strip().replace('XT:Z:', '') if 'XT:Z:' in x else 'NaN')
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
        print("Error: gtf file cannot be converted to dataframe. Make sure the annotation file provided is in gene \
        transfer format (.gtf) and respects Ensembl's format.")
        sys.exit(1)
    df_gtf = df_gtf[df_gtf.feature == 'gene']
    df_gtf = df_gtf.drop(['seqname', 'source', 'feature', 'start', 'end', 'strand', 'transcript_id', 'exon_number',
                          'transcript_name', 'transcript_biotype', 'transcript_support_level'], axis=1)
    df_unique = read_count_matrix(unique_file)
    df_unique = df_unique.merge(df_gtf, on='gene_id', how='outer')
    df_unique = df_unique.fillna(0.0)
    df_unique = df_unique.set_index('gene_id')
    df_unique = df_unique.reset_index()
    return df_unique


def distribute_multireads(featurefile, uniquefile, R_opt, gtf_file, output_file):
    print('unique_counts')
    df_unique = unique_counts(gtf_file, uniquefile)
    print('read multi featurefile')
    if R_opt=='SAM':
        df_multi = read_samfile(featurefile)
    elif R_opt == 'CORE':
        colnames = ['QNAME','status','nb_targets','gene_id']
        df_multi = read_corefile(featurefile, colnames)
    elif R_opt == '':
        colnames = ['QNAME', 'status', 'gene_id', 'nb_targets']
        df_multi = read_corefile(featurefile, colnames)
    else:
        sys.exit(1)
    print('distribute counts')
    df_dist = distribute_counts(df_multi, df_unique)
    if R_opt == 'SAM':
        print('write sam output')
        outfile = output_file + '_distributed_multi.sam'
        df_dist.to_csv(outfile, header=False, index=False, sep='\t')
    print('groupby gen_id')
    df_group = df_dist.groupby('gene_id').dist_counts.sum().reset_index()
    print('merge df multi unique')
    df_tot = pd.merge(df_group,df_unique, on='gene_id')
    df_tot['tot'] = df_tot.accumulation + df_tot.dist_counts
    df_tot[['gene_id','tot']].to_csv(output_file, header=False, index=False, sep='\t')
