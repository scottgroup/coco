import pandas as pd
pd.options.mode.chained_assignment = None
import sys
import subprocess
from gtf import dataframe
from distribute_counts import read_count_matrix


def get_main_transcript(dtranscript):
    """
    Keeps only the main transcript for each genes. The main transcript is determined by all these parameters,ranging from most important to least:
    transcript_support_level, source (ensembl_havana first, havana second, all the others after) and transcript_name.
    :param dtranscript: dataframe containing all transcripts
    :return:dataframe containing one main transcript per gene_id
    """
    gene_list = set(dtranscript.gene_id.unique())
    dtranscript2=dtranscript[dtranscript['gene_biotype']==dtranscript['transcript_biotype']]
    short_list = set(dtranscript2.gene_id.unique())
    gene_lost = gene_list-short_list
    dtranscript2 = dtranscript2.append(dtranscript[dtranscript.gene_id.isin(gene_lost)])
    dtranscript2.loc[dtranscript2['source']=='ensembl_havana','source']='aa' #source with highest level of confidence.
    dtranscript2.loc[dtranscript2['source']=='havana','source']='ab' #source with second level of confidence.
    if 'transcript_support_level' in dtranscript2.columns:
        dtranscript2=dtranscript2.sort_values(by=['gene_id','transcript_support_level','source','transcript_name'],
                                            ascending=[True,True,True,True])
    else:
        dtranscript2=dtranscript2.sort_values(by=['gene_id','source','transcript_name'],ascending=[True,True,True])
    dtranscript2=dtranscript2.drop_duplicates(subset=['gene_id'])
    return dtranscript2


def get_true_length_from_gtf(df_gtf):
    #Gives true length of transcripts without introns
    dtranscript=df_gtf[df_gtf['feature']=='transcript']
    dtranscript=get_main_transcript(dtranscript)
    dtranscript=dtranscript[['gene_id','transcript_id']]
    dexon=df_gtf[df_gtf['feature']=='exon']
    dexon=dexon[dexon['transcript_id'].isin(dtranscript['transcript_id'])==True]
    dexon['length']=dexon['end']-dexon['start']
    df_gtf=df_gtf[df_gtf['feature'] == 'gene']
    dexon_groupby=dexon.groupby('gene_id', as_index=False)['length'].sum()
    df_gtf=pd.merge(df_gtf,dexon_groupby,how='left',on='gene_id')
    df_gtf=df_gtf.sort_values(by='gene_id')
    return df_gtf


def add_pm_counts(count_file,df_gtf,bam_file, count_type):
    """
    Takes the input featureCounts output count file and modifies it to add CPM and TPM values. (adds gene_name as well).
    Uses the gtf.py script to read the gtf in a dataframe format. Takes about a minute to do so. You may skip the use of
    that script by specifying the --rawOnly option in CorrectCount.

    :param count_file: featureCounts count file.
    :param df_gtf: dataframe containing the gtf information.
    :param bam_file: bam file, used to calculate max read size.
    """
    read_bam = ('samtools','view', bam_file)
    head_file = ('head','-n','100000')
    select_col = ('cut','-f','10')
    wc = ('wc','-L')
    c1 = subprocess.Popen(read_bam, stdout=subprocess.PIPE)
    c2 = subprocess.Popen(head_file, stdin=c1.stdout, stdout=subprocess.PIPE)
    c3 = subprocess.Popen(select_col, stdin=c2.stdout, stdout=subprocess.PIPE)
    max_read_size=subprocess.check_output(wc, stdin=c3.stdout)

    try:
        max_read_size = int(str(max_read_size, 'utf-8'))
    except:
        print("error: There was a problem while reading the bam file to get the max read size.\n"
              "command performed: 'samtools view %s | head -n 10| cut -f 10 | wc -L'\n"
              "output obtained: %s" %(bam_file,max_read_size))
        sys.exit(1)
    df_gtf = get_true_length_from_gtf(df_gtf)
    df_gtf = df_gtf[['gene_id', 'gene_name', 'gene_biotype', 'length']]
    if count_type == 'uniqueOnly':
        dcount = read_count_matrix(count_file)
        dcount = dcount.rename(columns={'accumulation': 'count'})
    else:
        dcount = pd.read_csv(count_file, sep='\t', header=None, names=['gene_id', 'count'])
    dcount[['count']] = dcount[['count']].astype(float)
    Assigned = dcount['count'].sum()
    dcount['cpm'] = (dcount['count'].map(float) / Assigned) * 1E6
    dcount = pd.merge(dcount, df_gtf, how='right', on='gene_id')
    dcount['temp'] = dcount['count'] / dcount['length']
    dcount.loc[dcount['length'] < max_read_size, 'temp'] = dcount['count'] / max_read_size
    sum_temp = dcount['temp'].sum()
    dcount['tpm'] = (dcount['temp'] / sum_temp) * 1E6
    # dcount.loc[dcount['length'] < read_length_dict[exp], 'tpm'] = dcount['cpm']
    del dcount['temp']
    # dcount=dcount.sort_values(by=['tpm'],ascending=False)
    dcount = dcount[['gene_id', 'gene_name', 'count', 'cpm', 'tpm']]
    dcount.to_csv(path_or_buf=count_file, index=False, sep='\t', header=True)

