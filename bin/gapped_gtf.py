import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.width', 400)
import sys, os
import csv
from gtf import dataframe

def getoverlap(a, b):
    if max(0, min(a[1], b[1]) - max(a[0], b[0])) > 0:
        return (min(a[0], b[0]), max(a[1], b[1]))
    else:
        return a

def unique_list(list):
    seen = set()
    seen_add = seen.add
    return [x for x in list if not (x in seen or seen_add(x))]

def intersect(dataf1,dataf2,output='./Intersect',name=('name','name'),score=0,keep='all',r=0,join='loj',option=''):
    intersect_output=output+'.temp.intersect.bed'
    dataf_1_output=output+'.temp.1.bed'
    dataf_2_output=output+'.temp.2.bed'
    if score == 0:
        dataf1['score']=0
        dataf2['score']=0
        score='score'
    dataf1=dataf1[['seqname','start','end',name[0],score,'strand']]
    dataf2=dataf2[['seqname','start','end',name[1],score,'strand']]

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
                         names=['seqname','start','end',name[0],score,'strand',
                                'dataf2_file','start_2','end_2',name[1],'overlap','strand_2'])
    del Intersect_dataf['dataf2_file'], Intersect_dataf['strand_2']
    if keep=='name_only':
        Intersect_dataf=Intersect_dataf[[name,name+'_2']]
    command= "rm %s &&\
    rm %s && \
    rm %s" \
    %(dataf_1_output, dataf_2_output, intersect_output)
    os.system(command)
    return Intersect_dataf


def make_group_biotype(biotype_dataf):
    biotype_group_dict={'protein_coding':['IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IGLV_gene', 'IGM_gene', 'IG_V_gene', 'IGZ_gene', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'polymorphic_pseudogene', 'protein_coding', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene'],
                        'pseudogene':['disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene', 'processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene'],
                        'Long_noncoding':['3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense', 'antisense_RNA', 'lincRNA', 'ncrna_host', 'processed_transcript', 'sense_intronic', 'sense_overlapping'],
                        'snoRNA':['snoRNA'],'snRNA':['snRNA'],'scaRNA':['scaRNA'],'miRNA':['miRNA'],'tRNA':['tRNA','Mt_tRNA'],'7SL':['7SL'],'7SK':['7SK']}
    biotype_dataf.loc[biotype_dataf['gene_name'].str.contains('7SL') == True,'gene_biotype']='7SL'
    biotype_dataf.loc[biotype_dataf['gene_name'].str.contains('Metazoa_SRP') == True,'gene_biotype']='7SL'
    biotype_dataf.loc[biotype_dataf['gene_name'].str.contains('7SK') == True,'gene_biotype']='7SK'
    biotype_dataf.loc[(biotype_dataf['gene_name'].str.contains('SCA') == True) & (biotype_dataf['gene_biotype'] == 'snoRNA'),'gene_biotype']='scaRNA'
    for biotype_group in biotype_group_dict:
        biotype_dataf.loc[biotype_dataf['gene_biotype'].isin(biotype_group_dict[biotype_group]) == True,'gene_biotype_group']=biotype_group
    biotype_dataf.loc[biotype_dataf['gene_biotype_group'].isnull() == True,'gene_biotype_group']='Other'
    return biotype_dataf


def drill_an_exon(dIntersect,dexon_fix):
    print('Splitting exons overlapping embedded genes...')
    dIntersect=dIntersect.drop_duplicates()
    dupplicate_serie=dIntersect.groupby(dIntersect['exon_id'],as_index=False).size()
    dupplicate_quants=list(dupplicate_serie.values)
    unique_keys=list(dupplicate_serie.keys())
    dataf_dup=pd.DataFrame(data={'exon_id':unique_keys,
                         'count':dupplicate_quants})
    dIntersect=pd.merge(dIntersect,dataf_dup, on='exon_id', how='left')

    dIntersect_first_slice=dIntersect[dIntersect['count']==1].copy(deep=True)
    dIntersect_second_slice=dIntersect[dIntersect['count']==1].copy(deep=True)

    dIntersect_first_slice['end']=dIntersect_first_slice['start_2']
    dIntersect_first_slice=dIntersect_first_slice[dIntersect_first_slice['end'] > dIntersect_first_slice['start']]
    dIntersect_second_slice['start']=dIntersect_second_slice['end_2']
    dIntersect_second_slice=dIntersect_second_slice[dIntersect_second_slice['start'] < dIntersect_second_slice['end']]
    dexon_slice=pd.concat([dIntersect_first_slice,dIntersect_second_slice])
    dexon_slice=dexon_slice.sort_values(by=['exon_id'])

    ### Fixing for exons that have more than one embeded gene.
    dIntersect_many_slice=dIntersect[dIntersect['count']>1].copy(deep=True)
    dIntersect_many_slice.loc[dIntersect_many_slice['start_2']< dIntersect_many_slice['start'],'start_2']=dIntersect_many_slice['start']
    dIntersect_many_slice.loc[dIntersect_many_slice['end_2']> dIntersect_many_slice['end'],'end_2']=dIntersect_many_slice['end']
    dIntersect_many_slice=dIntersect_many_slice.sort_values(by=['exon_id','start_2'])
    dexon_more_than_one=dexon_fix[dexon_fix['exon_id'].isin(list(dIntersect_many_slice['exon_id'])) == True]

    dexon_more_than_one_fixed_exons=dexon_more_than_one[dexon_more_than_one['gene_name']=='empy_dataf'].copy(deep=True)
    for index,row in dexon_more_than_one.iterrows():
        dexon_temp=dIntersect_many_slice[dIntersect_many_slice['exon_id'] == row['exon_id']]
        old_gap_list=list(zip(dexon_temp['start_2'].map(int), dexon_temp['end_2'].map(int)))
        gap_list=old_gap_list.copy()
        for i, gap in enumerate(old_gap_list):
            for j, gap_2 in enumerate(old_gap_list[:i] + old_gap_list[i+1:]):
                gap_list[i]=getoverlap(gap, gap_2)
        gap_list=sorted(unique_list(gap_list), key=lambda tup: tup[0])
        gap_list_gap_ends=[ value[0]-1 for value in gap_list ]+[row['end']]
        gap_list_gap_starts=[row['start']]+[ value[1]+1 for value in gap_list ]
        dexon_temp=pd.DataFrame(data={'seqname':row['seqname'],'exon_id':row['exon_id'],'start':gap_list_gap_starts,'end':gap_list_gap_ends,'strand':row['strand'],'source':row['source'],'feature':row['feature'],'gene_id':row['gene_id'],'transcript_id':row['transcript_id'],'exon_number':row['exon_number'],'gene_name':row['gene_name'],'gene_biotype':row['gene_biotype'],'transcript_name':row['transcript_name'],'transcript_biotype':row['transcript_biotype'],'transcript_support_level':row['transcript_support_level'],'score':row['score']})
        dexon_more_than_one_fixed_exons=pd.concat([dexon_more_than_one_fixed_exons,dexon_temp])
        dexon_more_than_one_fixed_exons=dexon_more_than_one_fixed_exons[dexon_more_than_one_fixed_exons['end']-dexon_more_than_one_fixed_exons['start']>1]
    dexon_more_than_one_fixed_exons=dexon_more_than_one_fixed_exons[dexon_more_than_one_fixed_exons['end']-dexon_more_than_one_fixed_exons['start']>1]
    dexon_slice=dexon_slice[['seqname','start','end','exon_id','strand']]
    dexon_fix=dexon_fix[['gene_id', 'transcript_id', 'exon_number', 'feature', 'gene_name', 'gene_biotype', 'source', 'transcript_name', 'transcript_biotype', 'transcript_support_level', 'exon_id', 'score']]
    dexon_slice=pd.merge(dexon_slice,dexon_fix,how='left',on='exon_id')
    dexon_slice=pd.concat([dexon_slice,dexon_more_than_one_fixed_exons])
    return dexon_slice


def fix_exon_number(dexon):
    dexon=dexon.sort_values(by=['transcript_id','exon_id','start']).reset_index()
    dupplicate_serie=dexon.groupby(dexon['transcript_id'],as_index=False).size()
    dupplicate_quants=list(dupplicate_serie.values)
    unique_keys=list(dupplicate_serie.keys())
    dataf_dup=pd.DataFrame(data={'transcript_id':unique_keys,
                         'count':dupplicate_quants})
    dexon=pd.merge(dexon,dataf_dup, on='transcript_id', how='left')
    dexon_temp=dexon[:]
    dexon_temp=dexon.copy(deep=True)
    dexon_temp['old_index']=dexon_temp.index
    dexon_temp=dexon_temp.drop_duplicates(subset='transcript_id')
    dexon=pd.merge(dexon,dexon_temp[['transcript_id','old_index']],how='left',on='transcript_id')
    dexon['exon_number']=dexon.index-dexon['old_index']+1
    del dexon['index'],dexon['old_index'],dexon['count']
    dexon['exon_id']=dexon['transcript_id']+'.'+dexon['exon_number'].map(str)
    return dexon


def fix_gene_and_transcript_size(df_gtf,dexon):
    print('Fixing gene and transcript size...')
    df_gtf=df_gtf[df_gtf['feature']!='exon']
    dexon_groupby_transcript_min=dexon.groupby(dexon['transcript_id'])['start'].min().to_frame().rename(columns={'start':'min_start'})
    dexon_groupby_transcript_max=dexon.groupby(dexon['transcript_id'])['end'].max().to_frame().rename(columns={'end':'max_end'})
    dexon_groupby_gene_min=dexon.groupby(dexon['gene_id'])['start'].min().to_frame().rename(columns={'start':'min_start'})
    dexon_groupby_gene_max=dexon.groupby(dexon['gene_id'])['end'].max().to_frame().rename(columns={'end':'max_end'})
    ##Fixing transcript positions
    df_gtf=pd.merge(df_gtf,dexon_groupby_transcript_min,left_on='transcript_id',right_index=True, how='left')
    df_gtf=pd.merge(df_gtf,dexon_groupby_transcript_max,left_on='transcript_id',right_index=True, how='left')
    df_gtf.loc[(df_gtf['feature']=='transcript')&(df_gtf['start'] < df_gtf['min_start']),'start']=df_gtf['min_start']
    df_gtf.loc[(df_gtf['feature']=='transcript')&(df_gtf['end'] > df_gtf['max_end']),'end']=df_gtf['max_end']
    del df_gtf['min_start'],df_gtf['max_end']
    #Fixing gene positions
    df_gtf=pd.merge(df_gtf,dexon_groupby_gene_min,left_on='gene_id',right_index=True, how='left')
    df_gtf=pd.merge(df_gtf,dexon_groupby_gene_max,left_on='gene_id',right_index=True, how='left')
    df_gtf.loc[(df_gtf['feature']=='gene')&(df_gtf['start'] < df_gtf['min_start']),'start']=df_gtf['min_start']
    df_gtf.loc[(df_gtf['feature']=='gene')& (df_gtf['end'] > df_gtf['max_end']),'end']=df_gtf['max_end']
    df_gtf['start']=df_gtf['start'].map(int)
    df_gtf['end']=df_gtf['end'].map(int)
    del df_gtf['min_start'],df_gtf['max_end']
    return df_gtf

def build_gapped_gtf(df_gtf,dexon,output):
    df_gtf=fix_gene_and_transcript_size(df_gtf,dexon)
    print('Building correct_annotation gtf...')
    df_gtf=pd.concat([df_gtf,dexon])
    df_gtf['exon_number']=df_gtf['exon_number'].fillna(value=0,axis=0).map(int).map(str)
    df_gtf.loc[df_gtf['feature'] == 'exon','feature']='zexon'   #Makes proper gtf sorting simpler
    df_gtf.loc[df_gtf['feature'] == 'gene','transcript_id']='A' #Makes proper gtf sorting simpler
    df_gtf=df_gtf.sort_values(by=['gene_id','transcript_id','feature'])
    df_gtf.loc[df_gtf['feature'] == 'zexon','feature']='exon'
    df_gtf=df_gtf.reset_index(drop=True)
    df_gtf['transcript_support_level']=df_gtf['transcript_support_level'].fillna(value='NA')
    df_gtf.loc[df_gtf['feature'] == 'gene','attribute']='gene_id "'+df_gtf['gene_id']+'"; gene_name "'+df_gtf['gene_name']+'"; gene_source "'+df_gtf['source']+'"; gene_biotype "'+df_gtf['gene_biotype']+'";'
    df_gtf.loc[df_gtf['feature'] == 'transcript','attribute']='gene_id "'+df_gtf['gene_id']+'"; transcript_id "'+df_gtf['transcript_id']+'"; gene_name "'+df_gtf['gene_name']+'"; gene_source "'+df_gtf['source']+'"; gene_biotype "'+df_gtf['gene_biotype']+'"; transcript_name "'+df_gtf['transcript_name']+'"; transcript_biotype "'+df_gtf['transcript_biotype']+'"; transcript_support_level "'+df_gtf['transcript_support_level']+'";'
    df_gtf.loc[df_gtf['feature'] == 'exon','attribute']='gene_id "'+df_gtf['gene_id']+'"; transcript_id "'+df_gtf['transcript_id']+'"; exon_number "'+df_gtf['exon_number']+'"; gene_name "'+df_gtf['gene_name']+'"; gene_source "'+df_gtf['source']+'"; gene_biotype "'+df_gtf['gene_biotype']+'"; transcript_name "'+df_gtf['transcript_name']+'"; transcript_biotype "'+df_gtf['transcript_biotype']+'"; transcript_support_level "'+df_gtf['transcript_support_level']+'"; exon_id "'+df_gtf['exon_id']+'";'
    df_gtf['score']='.'
    df_gtf['frame']='.'
    df_gtf=df_gtf[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame','attribute']]
    df_gtf.to_csv(path_or_buf=output,
                          index=False, sep='\t', header=False, quoting=csv.QUOTE_NONE)
    print('All done!')



def correct_annotation(gtf_file, output, biotypes_embedded=('snoRNA', 'scaRNA', 'tRNA', 'miRNA', 'snRNA')):
    """
    correct_annotation builds a new gtf from the one provided, with holes on exons from genes that overlap the specified embedded biotypes.
    Read the MANUAL.md for extensive description.
    :param gtf_file: orifinal gtf file
    :param biotypes_embedded: list of the embedded biotypes. Default: 'snoRNA','scaRNA','tRNA','miRNA' and 'snRNA'
    :output: modified gtf file with the .correct_annotation.gtf prefix.
    """
    if gtf_file.endswith('.gtf')==True:
        print('Reading gtf')
    else:
        print('Annotation file:',gtf_file)
        print('error: Wrong annotation format. Only .gtf files are accepted')
        sys.exit(1)
    try:
        df_gtf=dataframe(gtf_file)
        df_gtf=df_gtf[['seqname', 'source', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id',
                       'exon_number', 'gene_name', 'gene_biotype', 'transcript_name', 'transcript_biotype',
                       'transcript_support_level']]
        df_gtf['seqname']=df_gtf['seqname'].map(str)
        df_gtf['start']=df_gtf['start'].map(int)
        df_gtf['end']=df_gtf['end'].map(int)
    except:
        print("error: gtf file cannot be converted to dataframe. Make sure the annotation file provided is in gene transfer format (.gtf) and comes from Ensembl.")
        sys.exit(1)
    df_gtf.loc[df_gtf['transcript_name'].isnull()==True,'transcript_name']=df_gtf['gene_name']
    df_gtf.loc[df_gtf['transcript_biotype'].isnull()==True,'transcript_biotype']=df_gtf['gene_biotype']
    df_gtf=df_gtf[df_gtf['feature'].isin(['gene','transcript','exon'])==True]
    dgene=df_gtf[df_gtf.feature == 'gene']
    dgene['start']=dgene['start'].map(int)
    dgene['end']=dgene['end'].map(int)
    dgene['seqname']=dgene['seqname'].map(str)
    #dgene=make_group_biotype(dgene)
    dgene_embedded=dgene[dgene['gene_biotype'].isin(biotypes_embedded) == True]
    dgene_host=dgene[dgene['gene_biotype'].isin(biotypes_embedded) == False]

    dexon_host=df_gtf[(df_gtf.feature == 'exon') & (df_gtf['gene_id'].isin(dgene_host['gene_id'])==True)]
    dexon_host['exon_id']=dexon_host['transcript_id']+'.'+dexon_host['exon_number'].map(int).map(str)
    dexon_not_host=df_gtf[(df_gtf.feature == 'exon') & (df_gtf['gene_id'].isin(dgene_host['gene_id'])==False)]
    dexon_not_host['exon_id']=dexon_not_host['transcript_id']+'.'+dexon_not_host['exon_number'].map(int).map(str)
    dIntersect=intersect(dexon_host,dgene_embedded,name=('exon_id','gene_id'))
    dIntersect=dIntersect[dIntersect['overlap'] != -1]
    del dIntersect['overlap']
    dexon_slice=drill_an_exon(dIntersect,dexon_host)
    dexon_host=dexon_host[dexon_host['exon_id'].isin(dIntersect['exon_id']) == False]
    dexon_host=pd.concat([dexon_host,dexon_slice])
    dexon=pd.concat([dexon_host,dexon_not_host])
    dexon=fix_exon_number(dexon)
    if output =='None':
        output=gtf_file.replace('.gtf','.correct_annotation.gtf')
    build_gapped_gtf(df_gtf,dexon,output)

if __name__=='__main__':
    if len(sys.argv)>1:
        correct_annotation(gtf_file=sys.argv[1])
    else:
        print('error: Not enough arguments!')
        sys.exit(1)