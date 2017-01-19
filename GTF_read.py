import pandas as pd
import sys
import re

def parse_gtf(row):
    attribute=row['attribute']
    try:
        if 'transcript_version' in attribute:
            if 'exon_number' in attribute:
                att_line=re.match('gene_id "(.*?)"; .*?transcript_id "(.*?)"; .*?exon_number "(.*?)"; .*?gene_name "(.*?)"; .*?gene_biotype "(.*?)";.*?transcript_name "(.*?)";.*?transcript_biotype "(.*?)";.*?', attribute)
                row['gene_id']=att_line.group(1)
                row['transcript_id']=att_line.group(2)
                row['exon_number']=att_line.group(3)
                row['gene_name']=att_line.group(4)
                row['gene_biotype']=att_line.group(5)
                row['transcript_name']=att_line.group(6)
                row['transcript_biotype']=att_line.group(7)
            else:
                att_line=re.match('gene_id "(.*?)"; .*?transcript_id "(.*?)"; .*?gene_name "(.*?)"; .*?gene_biotype "(.*?)";.*?', attribute)
                row['gene_id']=att_line.group(1)
                row['transcript_id']=att_line.group(2)
                row['gene_name']=att_line.group(3)
                row['gene_biotype']=att_line.group(4)
                # row['transcript_name']=att_line.group(8)
                # row['transcript_source']=att_line.group(9)
                # row['transcript_biotype']=att_line.group(10)
        else:
            att_line=re.match('gene_id "(.*?)"; .*?gene_name "(.*?)"; .*?gene_biotype "(.*?)";.*?', attribute)
            row['gene_id']=att_line.group(1)
            row['gene_name']=att_line.group(2)
            print(row['gene_name'])
            row['gene_biotype']=att_line.group(3)

    except(AttributeError, TypeError):
        print(row)
        print(attribute)
        sys.exit("Error: unexpected gtf format made parsing impossible")
    #print(row['gene_name'])
    return row


def fetch_genes(gtf_file,feature_to_keep='GTE'):
    df_gtf = pd.read_csv(gtf_file, sep='\t',names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame','attribute'])
    df_gtf = df_gtf.dropna(axis='index', how='any', subset=['start'])
    #df_gtf=df_gtf[-500:]
    if feature_to_keep == 'GTE':
        keep=['gene','transcript','exon']
        print(len(df_gtf))
        df_gtf = df_gtf[df_gtf.feature.isin(keep) == True]
        print(len(df_gtf))
    elif feature_to_keep != 'all' and feature_to_keep != 'GTE':
        df_gtf = df_gtf[df_gtf.feature == feature_to_keep]

    #add method to return normal gtf
    #df_complete_gtf = df_gtf.copy(deep=True)
    #df_complete_gtf['attribute'] = [x.split(';')[0].rstrip('" ').lstrip('gene_id "') for x in df_complete_gtf['attribute']]

    df_gtf=df_gtf.apply(lambda row: parse_gtf(row), axis=1)
    # df_gtf['start']=df_gtf['start'].map(int)
    # df_gtf['end']=df_gtf['end'].map(int)
    # df_gtf['length']=df_gtf['end']-df_gtf['start']
    del df_gtf['attribute']
    return df_gtf


if __name__=='__main__':
    gtf_file='/home/gabrielle/genome/hg38_87/human_ensembl_87.gtf'
    #feature_to_keep='gene'
    df_gtf=fetch_genes(gtf_file)
    #print(df_gtf[:5])
    # df_gtf=df_gtf.rename(columns={'attribute':'ens_id'})
    # #df_gtf['length']=df_gtf['length'].map(int)
    # #df_gtf=df_gtf[['ens_id','length']]
    # df_gtf.to_csv(path_or_buf='/home/vincent/Desktop/Results/Seq/hg38/FabDupSan_human_ensembl_83_length.csv',
    #                            index=False, sep='\t', header=False)
    df_gtf=df_gtf[['chr', 'source', 'feature', 'start', 'end', 'strand', 'gene_id','transcript_id','exon_number','gene_name','gene_biotype','transcript_name','transcript_biotype']]
    dexon=df_gtf[df_gtf.feature == 'exon']
    dexon['feature']='transcript'
    dexon=dexon[['transcript_id','transcript_name','transcript_biotype']]
    dexon=dexon.drop_duplicates()
    del df_gtf['transcript_name'],df_gtf['transcript_biotype']
    df_gtf=pd.merge(df_gtf,dexon,how='left',on=['transcript_id'])
    df_gtf.to_csv(path_or_buf='/home/gabrielle/genome/hg38_87/human_ensembl_87.csv',
                        index=False, sep='\t', header=True)
    # print(df_gtf[:5])

