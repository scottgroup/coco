import pandas as pd
pd.options.mode.chained_assignment = None
import os,sys
import re
from gtf import dataframe
from distribute_counts import read_count_matrix


def get_main_transcript(dtranscript):
    """
    Keeps only the main transcript for each genes. The main transcript is determined by all these parameters,ranging from most important to least:
    transcript_support_level, source (ensembl_havana first, havana second, all the others after) and transcript_name.
    :param dtranscript: dataframe containing all transcripts
    :return:dataframe containing one main transcript per gene_id
    """
    dtranscript=dtranscript[dtranscript['gene_biotype']==dtranscript['transcript_biotype']]
    dtranscript.loc[dtranscript['source']=='ensembl_havana','source']='aa' #source with highest level of confidence.
    dtranscript.loc[dtranscript['source']=='havana','source']='ab' #source with second level of confidence.
    if 'transcript_support_level' in dtranscript.columns:
        dtranscript=dtranscript.sort_values(by=['gene_id','transcript_support_level','source','transcript_name'],ascending=[True,True,True,True])
    else:
        dtranscript=dtranscript.sort_values(by=['gene_id','source','transcript_name'],ascending=[True,True,True])
    dtranscript=dtranscript.drop_duplicates(subset=['gene_id'])
    return dtranscript


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


def add_pm_counts(count_file,gtf_file,bam_file, count_type):
    """
    Takes the input featureCounts output count file and modifies it to add CPM and TPM values. (adds gene_name as well).
    Uses the gtf.py script to read the gtf in a dataframe format. Takes about a minute to do so. You may skip the use of that script
    by specifying the --rawOnly option in CorrectCount.

    :param count_file: featureCounts count file.
    :param gtf_file: annotation file in gtf format.
    :param bam_file: bam file, used to calculate max read size.
    """
    command='samtools view %s | head -n 100000 | cut -f 10 | wc -L' %(bam_file)
    max_read_size=os.system(command)
    try:
        max_read_size=int(max_read_size)
    except:
        print("error: There was a problem while reading the bam file to get the max read size.\n"
              "command performed: 'samtools view %s | head -n 100000 | cut -f 10 | wc -L'\n"
              "output obtained: %s" %(bam_file,max_read_size))
        sys.exit(1)
    try:
        df_gtf=dataframe(gtf_file)
        df_gtf = df_gtf[
            ['seqname', 'source', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'exon_number',
             'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype', 'transcript_support_level']]
        df_gtf['seqname'] = df_gtf['seqname'].map(str)
        df_gtf['start'] = df_gtf['start'].map(int)
        df_gtf['end'] = df_gtf['end'].map(int)
    except:
        raise
        #print('error: gtf file %s could not be read. Please specify a proper annotation file in gene transfert format (.gtf) obtained from Ensembl.' %(gtf_file))
        #sys.exit(1)
    df_gtf=get_true_length_from_gtf(df_gtf)
    df_gtf=df_gtf[['gene_id','gene_name','gene_biotype','length']]
    if count_type  == 'uniqueOnly':
        dcount = read_count_matrix(count_file)
        dcount = dcount.rename(columns={'accumulation':'count'})
    else :
        dcount = pd.read_csv(count_file, sep=' ', header=None, names=['gene_id','count'])
    Assigned=dcount['count'].sum()
    dcount['cpm']=(dcount['count'].map(float)/Assigned)*1E6
    dcount=pd.merge(dcount,df_gtf,how='right',on='gene_id')
    dcount['temp']=dcount['count']/dcount['length']
    dcount.loc[dcount['length'] < max_read_size,'temp']=dcount['count']/max_read_size
    sum_temp=dcount['temp'].sum()
    dcount['tpm']=(dcount['temp']/sum_temp)*1E6
    #dcount.loc[dcount['length'] < read_length_dict[exp], 'tpm'] = dcount['cpm']
    del dcount['temp']
    #dcount=dcount.sort_values(by=['tpm'],ascending=False)
    dcount=dcount[['gene_id','gene_name','count','cpm','tpm']]
    dcount.to_csv(path_or_buf=count_file,
                      index=False, sep='\t', header=True)



#if __name__ == '__main__':
    #main(file='/home/vincent/Desktop/Sequencing/Methods_Compare/Total/raw_count/Rsubread_CoCo/%s.CorrectCount.Rsubread.count')
    #main(file='/home/vincent/Desktop/Sequencing/Methods_Compare/Total/raw_count/before_correction/Notcoco/%s.CorrectCount_NOTCOCO.Rsubread.count')
    #main(file='/home/vincent/Desktop/Sequencing/Methods_Compare/Total/raw_count/before_correction/Initial/%s.CorrectCount.Rsubread.count')
    #main(file='/home/vincent/Desktop/Sequencing/Methods_Compare/Total/raw_count/old_CoCo/%s.CorrectCount_wo_dup.Rsubread.count')
