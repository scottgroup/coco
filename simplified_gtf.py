import pandas as pd
pd.options.mode.chained_assignment = None
import GTF_read
import sys, time, os
import csv
pd.set_option('display.width', 400)


def Intersect(dataf1,dataf2,output='./Intersect',name='name',score=0,keep='all',r=0,join='loj',option=''):
    intersect_output=output+'.temp.intersect.bed'
    dataf_1_output=output+'.temp.1.bed'
    dataf_2_output=output+'.temp.2.bed'
    if score == 0:
        dataf1['score']=0
        dataf2['score']=0
        score='score'
    dataf1=dataf1[['chr','start','end',name,score,'strand']]
    dataf2=dataf2[['chr','start','end',name,score,'strand']]

    dataf1['start']=dataf1['start'].map(int)
    dataf1['end']=dataf1['end'].map(int)
    dataf2['start']=dataf2['start'].map(int)
    dataf2['end']=dataf2['end'].map(int)
    dataf1.to_csv(path_or_buf=dataf_1_output,
                          index=False, sep='\t', header=False)
    dataf2.to_csv(path_or_buf=dataf_2_output,
                          index=False, sep='\t', header=False)
    command= "bedtools intersect \
    -a %s \
    -b %s \
    -s " \
    %(dataf_1_output, dataf_2_output)
    command+="-"+join+" "
    if r > 0 :
        command+= "-f %s -r " \
        %(str(r))
    if option != '':
        command+= option
    command+="| sed 's/\t\.\t/\tx\t/g' \
    >  %s" \
    %(intersect_output)
    os.system(command)
    Intersect_dataf=pd.read_csv(filepath_or_buffer=intersect_output,
                       index_col=False, sep='\t', header=None,
                         names=['chr','start','end',name,score,'strand',
                                'dataf2_file','start_2','end_2',name+'_2','overlap','strand_2'])
    del Intersect_dataf['dataf2_file'], Intersect_dataf['strand_2']
    if keep=='name_only':
        Intersect_dataf=Intersect_dataf[[name,name+'_2']]
    command= "rm %s &&\
    rm %s && \
    rm %s" \
    %(dataf_1_output, dataf_2_output, intersect_output)
    os.system(command)
    return Intersect_dataf


def get_simplified_exons(dIntersect,dgene):
    """
    Takes the dataframe holding the gaps corresponding to snoRNAs and
    :param dIntersect:
    :param dgene:
    :return:
    """
    dIntersect_one=dIntersect[dIntersect['count'] == 1]
    dIntersect_one=dIntersect_one[['chr','gene_id','gap_start','gap_end','strand']]
    dIntersect_more_than_one=dIntersect[dIntersect['count'] >= 1]
    dIntersect_more_than_one=dIntersect_more_than_one[['chr','gene_id','gap_start','gap_end','strand']]
    dIntersect_more_than_one=dIntersect_more_than_one.sort_values(by=['gene_id','gap_start'])

    ###For the genes that have only one gap, no iterative part, makes it much faster.
    print('one gaps')
    start = time.time()
    dexon_one=dgene[dgene['gene_id'].isin(list(dIntersect_one['gene_id'])) == True]
    print(len(dexon_one))
    dexon_one=pd.merge(dexon_one,dIntersect_one,on=['chr','strand','gene_id'],how='left')
    dexon_one_1=dexon_one[['chr','gene_id','start','gap_start','strand']]
    dexon_one_1=dexon_one_1.rename(columns={'gap_start':'end'})
    dexon_one_2=dexon_one[['chr','gene_id','gap_end','end','strand']]
    dexon_one_2=dexon_one_2.rename(columns={'gap_end':'start'})
    dexon_one=pd.concat([dexon_one_1,dexon_one_2])
    del dexon_one_1, dexon_one_2
    end = time.time()
    print(end - start)

    ###For the genes that have more than one gap, done in an iterative manner.
    print('more gaps')
    start = time.time()
    dgene_host_more_than_one=dgene[dgene['gene_id'].isin(list(dIntersect_more_than_one['gene_id'])) == True]
    print(len(dgene_host_more_than_one))
    dexon_two=pd.DataFrame(data={'chr':[],'gene_id':[],'gap_end':[],'end':[],'strand':[]})
    #print(dgene_host_more_than_one[dgene_host_more_than_one['gene_id']=='ENSG00000178971'][:1])
    for index,row in dgene_host_more_than_one.iterrows():
        dgene_temp=dIntersect_more_than_one[dIntersect_more_than_one['gene_id'] == row['gene_id']]
        gap_list=list(zip(dgene_temp['gap_start'].map(int), dgene_temp['gap_end'].map(int)))
        gap_list_gap_ends=[ value[0]-1 for value in gap_list ]+[row['end']]
        gap_list_gap_starts=[row['start']]+[ value[1]+1 for value in gap_list ]
        dexon_temp=pd.DataFrame(data={'chr':row['chr'],'gene_id':row['gene_id'],'start':gap_list_gap_starts,'end':gap_list_gap_ends,'strand':row['strand']})
        dexon_two=pd.concat([dexon_two,dexon_temp])
    end = time.time()
    print(end - start)

    #Merging the exons from single and many gap(s) host genes
    dexon=pd.concat([dexon_one,dexon_two])
    dexon=dexon[dexon['start'] < dexon['end']]
    dexon=dexon[['chr','start','end','strand','gene_id']]
    dexon=dexon.drop_duplicates()
    dexon=dexon.sort_values(by=['gene_id','start','chr'])
    dexon=dexon.reset_index(drop=True)

    #Adding exon number...
    dupplicate_serie=dexon.groupby(dexon['gene_id'],as_index=False).size()
    dupplicate_quants=list(dupplicate_serie.values)
    unique_keys=list(dupplicate_serie.keys())
    dataf_dup=pd.DataFrame(data={'gene_id':unique_keys,
                         'count':dupplicate_quants})
    dexon=pd.merge(dexon,dataf_dup, on='gene_id', how='left')
    dexon_temp=dexon[:]
    dexon_temp['old_index']=dexon_temp.index
    dexon_temp=dexon_temp.drop_duplicates(subset='gene_id')
    dexon=pd.merge(dexon,dexon_temp[['gene_id','old_index']],how='left',on='gene_id')
    dexon['exon_number']=dexon.index-dexon['old_index']+1
    del dexon['count'], dexon['old_index']
    return dexon


def fix_positions(dgene_within):
    dgene_within_starts=dgene_within.copy(deep=True)
    dgene_within_starts=dgene_within_starts[dgene_within_starts['overlap'] != -1]
    dgene_within_starts=dgene_within_starts[dgene_within_starts['end_2'] <= dgene_within_starts['old_end']]
    dgene_within_starts=dgene_within_starts.sort_values(by=['gene_id','end_2'], ascending=[True,False])    #We want the overlaping exon with the biggest end at the top of the list
    dgene_within_starts=dgene_within_starts.drop_duplicates(subset=['gene_id']) #Drops all the overlapping exons but the one with biggest end for each gene
    dgene_within_starts['start']=dgene_within_starts['end_2']
    dgene_within_starts.loc[dgene_within_starts['start']>dgene_within_starts['old_start'],'start']=dgene_within_starts['old_start']
    dgene_within_starts['fix']='starts'
    dgene_within_starts=dgene_within_starts[['gene_id','start','fix']]

    dgene_within_ends=dgene_within.copy(deep=True)
    dgene_within_ends=dgene_within_ends[dgene_within_ends['overlap'] != -1]
    dgene_within_ends=dgene_within_ends[dgene_within_ends['start_2'] >= dgene_within_ends['old_start']]
    dgene_within_ends=dgene_within_ends.sort_values(by=['gene_id','start_2'], ascending=[True,True])    #We want the overlaping exon with the smallest start at the top of the list
    dgene_within_ends=dgene_within_ends.drop_duplicates(subset=['gene_id']) #Drops all the overlapping exons but the one with smallest start for each gene
    dgene_within_ends['end']=dgene_within_ends['start_2']
    dgene_within_ends.loc[dgene_within_ends['end']<dgene_within_ends['old_end'],'end']=dgene_within_ends['old_end']
    dgene_within_ends['fix']='ends'
    dgene_within_ends=dgene_within_ends[['gene_id','end','fix']]

    dgene_within_bigger_smaller=dgene_within.copy(deep=True)
    dgene_within_bigger_smaller=dgene_within_bigger_smaller[dgene_within_bigger_smaller['overlap'] != -1]
    dgene_within_bigger_smaller.loc[(dgene_within_bigger_smaller['start_2'] <= dgene_within_bigger_smaller['old_start'])&(dgene_within_bigger_smaller['end_2'] >= dgene_within_bigger_smaller['old_end']),'keep']='yes'
    dgene_within_bigger_smaller.loc[(dgene_within_bigger_smaller['start_2'] >= dgene_within_bigger_smaller['old_start'])&(dgene_within_bigger_smaller['end_2'] <= dgene_within_bigger_smaller['old_end']),'keep']='yes'
    dgene_within_bigger_smaller=dgene_within_bigger_smaller[dgene_within_bigger_smaller['keep']=='yes']
    del dgene_within_bigger_smaller['keep']
    dgene_within_bigger_smaller=dgene_within_bigger_smaller.drop_duplicates(subset=['gene_id'])
    dgene_within_bigger_smaller['end']=dgene_within_bigger_smaller['old_end']
    dgene_within_bigger_smaller['start']=dgene_within_bigger_smaller['old_start']
    #Reduce all genes that are either entirely covered by an exon or that englobes an exon entirely to their original size.
    dgene_within_bigger_smaller['fix']='bigsmall'
    dgene_within_bigger_smaller=dgene_within_bigger_smaller[['gene_id','start','end','fix']]

    dgene_within_fix=pd.concat([dgene_within_starts,dgene_within_ends,dgene_within_bigger_smaller])
    dgene_within_fix=dgene_within_fix.rename(columns={'start':'fixed_start','end':'fixed_end'})
    dgene_within_fix=dgene_within_fix.drop_duplicates()

    dgene_within_fix_start=dgene_within_fix[['gene_id','fixed_start']][:]
    dgene_within_fix_start=dgene_within_fix_start.sort_values(by=['gene_id','fixed_start'], ascending=[True,False]) #We want to keep the biggest fixed start
    dgene_within_fix_start=dgene_within_fix_start.drop_duplicates(subset='gene_id')

    dgene_within_fix_end=dgene_within_fix[['gene_id','fixed_end']][:]
    dgene_within_fix_end=dgene_within_fix_end.sort_values(by=['gene_id','fixed_end'], ascending=[True,True]) #We want to keep the smallest fixed end
    dgene_within_fix_end=dgene_within_fix_end.drop_duplicates(subset='gene_id')

    dgene_within_fix=pd.merge(dgene_within_fix_start,dgene_within_fix_end,how='outer',on='gene_id')

    del dgene_within['start_2'],dgene_within['end_2'],dgene_within['score'],dgene_within['overlap'],dgene_within['gene_id_2']
    dgene_within=dgene_within.drop_duplicates()
    dgene_within=pd.merge(dgene_within,dgene_within_fix,how='left',on='gene_id')

    dgene_within.loc[(dgene_within['fixed_start'].isnull() == False), 'start']=dgene_within['fixed_start']
    dgene_within.loc[(dgene_within['fixed_end'].isnull() == False), 'end']=dgene_within['fixed_end']
    del dgene_within['fixed_start'],dgene_within['fixed_end'],dgene_within['old_start'],dgene_within['old_end']
    return dgene_within



def enlarge_gene_within(dgene_within,df_gtf,window):
    dgene_within_old=dgene_within.copy(deep=True)
    dgene_within_old['old_start']=dgene_within_old['start']
    dgene_within_old['old_end']=dgene_within_old['end']
    dgene_within['length']=dgene_within['end']-dgene_within['start']
    dgene_within['start']=dgene_within['start']-(dgene_within['length']*(window-1)/2).round(0)
    dgene_within['end']=dgene_within['end']+(dgene_within['length']*(window-1)/2).round(0)

    dexon_canonical=df_gtf[df_gtf.feature == 'exon'][['chr','start','end','strand','gene_id']]
    dexon_canonical=dexon_canonical.drop_duplicates()
    dexon_canonical=dexon_canonical[dexon_canonical['gene_id'].isin(dgene_within['gene_id'])==False]    #To keep only exons from other genes
    dgene_within=Intersect(dgene_within,dexon_canonical,name='gene_id')
    dgene_within=pd.merge(dgene_within,dgene_within_old[['gene_id','old_start','old_end']],how='left',on='gene_id')
    del dexon_canonical     #Just to free a bit of RAM...

    dgene_within=fix_positions(dgene_within)
    dgene_within=pd.merge(dgene_within,dgene_within_old[['gene_id','source', 'feature', 'transcript_id', 'exon_number', 'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype']],how='left',on='gene_id')
    return dgene_within


def build_simplified_gtf(dgene,dexon_overlap,output,dgene_within=None):
    if dgene_within is not None:
        dgene_within=dgene_within.rename(columns={'start':'new_start','end':'new_end'})
        dgene_within=dgene_within[['gene_id','new_start','new_end']]
        dgene=pd.merge(dgene,dgene_within,how='left',on='gene_id')
        dgene.loc[dgene['new_start'].isnull()==False,'start']=dgene['new_start']
        dgene.loc[dgene['new_end'].isnull()==False,'end']=dgene['new_end']
        del dgene['new_start'],dgene['new_end']
    del dgene['exon_number']
    dtranscript=dgene.copy(deep=True)
    dtranscript['feature']='transcript'
    dtranscript['transcript_id']=dtranscript['gene_id']+'T'
    dtranscript['transcript_name']=dtranscript['gene_name']
    dtranscript['transcript_biotype']=dtranscript['gene_biotype']
    dtranscript['transcript_source']=dtranscript['source']
    #For all the genes that did not have an overlap with any snoRNA, miRNA, etc., make a single exon that has the whole length of the gene
    dexon=dgene[dgene['gene_id'].isin(dexon_overlap['gene_id']) == False]
    dexon['exon_number']=1

    dexon_overlap=pd.merge(dexon_overlap,dgene[['gene_id', 'source', 'feature', 'transcript_id', 'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype']], how='left', on='gene_id')
    dexon=pd.concat([dexon,dexon_overlap])
    dexon['exon_number']=dexon['exon_number'].map(str)
    dexon['feature']='zexon'
    dexon['transcript_id']=dexon['gene_id']+'T'
    dexon['exon_id']=dexon['gene_id']+'E'+'.'+dexon['exon_number']
    dexon['transcript_name']=dexon['gene_name']
    dexon['transcript_biotype']=dexon['gene_biotype']
    dexon['transcript_source']=dexon['source']
    dgene=pd.concat([dgene,dtranscript,dexon])
    dgene=dgene.sort_values(by=['gene_id','feature'])
    dgene.loc[dgene['feature'] == 'zexon','feature']='exon'
    dgene=dgene.reset_index(drop=True)
    dgene['exon_number']=dgene['exon_number'].map(str)
    dgene.loc[dgene['feature'] == 'gene','attribute']='gene_id "'+dgene['gene_id']+'"; gene_name "'+dgene['gene_name']+'"; gene_source "'+dgene['source']+'"; gene_biotype "'+dgene['gene_biotype']+'";'
    dgene.loc[dgene['feature'] == 'transcript','attribute']='gene_id "'+dgene['gene_id']+'"; transcript_id "'+dgene['transcript_id']+'"; gene_name "'+dgene['gene_name']+'"; gene_source "'+dgene['source']+'"; gene_biotype "'+dgene['gene_biotype']+'"; transcript_name "'+dgene['transcript_name']+'"; transcript_source "'+dgene['transcript_source']+'"; transcript_biotype "'+dgene['transcript_biotype']+'";'
    dgene.loc[dgene['feature'] == 'exon','attribute']='gene_id "'+dgene['gene_id']+'"; transcript_id "'+dgene['transcript_id']+'"; exon_number "'+dgene['exon_number']+'"; gene_name "'+dgene['gene_name']+'"; gene_source "'+dgene['source']+'"; gene_biotype "'+dgene['gene_biotype']+'"; transcript_name "'+dgene['transcript_name']+'"; transcript_source "'+dgene['transcript_source']+'"; transcript_biotype "'+dgene['transcript_biotype']+'"; exon_id "'+dgene['exon_id']+'";'
    dgene['score']='.'
    dgene['frame']='.'
    dgene=dgene[['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame','attribute']]
    dgene['start']=dgene['start'].map(int)
    dgene['end']=dgene['end'].map(int)
    dgene.to_csv(path_or_buf=output,
                          index=False, sep='\t', header=False, quoting=csv.QUOTE_NONE)



def main(window=3, csvgtf_provided=False, gtf_file='/home/gabrielle/genome/hg38_87/human_ensembl_87.gtf'):

    if csvgtf_provided == True:
        df_gtf = pd.read_csv(gtf_file, sep='\t')
    else:
        df_gtf=GTF_read.fetch_genes(gtf_file, feature_to_keep='all')
    dgene=df_gtf[df_gtf.feature == 'gene']
    dgene['start']=dgene['start'].map(int)
    dgene['end']=dgene['end'].map(int)
    dgene['chr']=dgene['chr'].map(str)
    biotypes_within=['snoRNA', 'scaRNA', 'tRNA', 'miRNA', 'snRNA']
    dgene_within=dgene[dgene['gene_biotype'].isin(biotypes_within) == True]
    dgene_big=dgene[dgene['gene_biotype'].isin(biotypes_within) == False]
    if window != 1:
        dgene_within=enlarge_gene_within(dgene_within,df_gtf,window)

    dIntersect=Intersect(dgene_big,dgene_within,name='gene_id')
    dIntersect=dIntersect[dIntersect['overlap'] != -1]

    dIntersect=dIntersect.rename(columns={'start_2':'gap_start','end_2':'gap_end'})
    del dIntersect['gene_id_2'],dIntersect['overlap']
    dIntersect=dIntersect.drop_duplicates()

    dupplicate_serie=dIntersect.groupby(dIntersect['gene_id'],as_index=False).size()
    dupplicate_quants=list(dupplicate_serie.values)
    unique_keys=list(dupplicate_serie.keys())
    dataf_dup=pd.DataFrame(data={'gene_id':unique_keys,
                         'count':dupplicate_quants})
    dIntersect=pd.merge(dIntersect,dataf_dup, on='gene_id', how='left')

    dexon=get_simplified_exons(dIntersect,dgene)

    gtf_file=gtf_file.replace('.csv','.Simplified_gap.gtf')
    if window > 1:
        build_simplified_gtf(dgene,dexon,gtf_file,dgene_within)
    else:
        build_simplified_gtf(dgene,dexon,gtf_file)

if __name__=='__main__':
    main()