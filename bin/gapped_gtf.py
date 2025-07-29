import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.width', 400)
import sys, os
import csv
import re
import gzip

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
                                'dataf2_file','start_2','end_2',name[1],'overlap','strand_2'], dtype={'seqname':str})
    del Intersect_dataf['dataf2_file'], Intersect_dataf['strand_2']
    if keep=='name_only':
        Intersect_dataf=Intersect_dataf[[name,name+'_2']]
    os.remove(dataf_1_output)
    os.remove(dataf_2_output)
    os.remove(intersect_output)
    return Intersect_dataf


def drill_an_exon(dIntersect,dexon_fix):
    print('Splitting exons overlapping embedded genes...')
    dIntersect=dIntersect.drop_duplicates()
    dataf_dup = dIntersect.groupby(dIntersect['exon_id'],as_index=False).size().reset_index()
    dataf_dup = dataf_dup.rename(columns={0: 'size'})
    dIntersect=pd.merge(dIntersect,dataf_dup, on='exon_id', how='left')

    dIntersect_first_slice=dIntersect[dIntersect['size']==1].copy(deep=True)
    dIntersect_second_slice=dIntersect[dIntersect['size']==1].copy(deep=True)

    dIntersect_first_slice['end']=dIntersect_first_slice['start_2']
    dIntersect_first_slice=dIntersect_first_slice[dIntersect_first_slice['end'] > dIntersect_first_slice['start']]
    dIntersect_second_slice['start']=dIntersect_second_slice['end_2']
    dIntersect_second_slice=dIntersect_second_slice[dIntersect_second_slice['start'] < dIntersect_second_slice['end']]
    if pd.__version__ >= '0.23.0':
        dexon_slice=pd.concat([dIntersect_first_slice,dIntersect_second_slice], sort=False)
    else:
        dexon_slice = pd.concat([dIntersect_first_slice, dIntersect_second_slice])
    dexon_slice=dexon_slice.sort_values(by=['exon_id'])

    ### Fixing for exons that have more than one embeded gene.
    dIntersect_many_slice=dIntersect[dIntersect['size']>1].copy(deep=True)
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
        dexon_temp = pd.DataFrame(data={'seqname':row['seqname'],'exon_id':row['exon_id'],'start':gap_list_gap_starts,
                                        'end':gap_list_gap_ends,'strand':row['strand'],'source':row['source'],
                                        'feature':row['feature'],'gene_id':row['gene_id'],'transcript_id':row['transcript_id'],
                                        'exon_number':row['exon_number'],'gene_name':row['gene_name'],
                                        'gene_biotype':row['gene_biotype'],'transcript_name':row['transcript_name'],
                                        'transcript_biotype':row['transcript_biotype'],
                                        'transcript_support_level':row['transcript_support_level'],'score':row['score']})
        if pd.__version__ >= '0.23.0':
            dexon_more_than_one_fixed_exons=pd.concat([dexon_more_than_one_fixed_exons,dexon_temp], sort=False)
        else:
            dexon_more_than_one_fixed_exons = pd.concat([dexon_more_than_one_fixed_exons, dexon_temp])
        dexon_more_than_one_fixed_exons=dexon_more_than_one_fixed_exons[dexon_more_than_one_fixed_exons['end']-dexon_more_than_one_fixed_exons['start']>1]
    dexon_more_than_one_fixed_exons=dexon_more_than_one_fixed_exons[dexon_more_than_one_fixed_exons['end']-dexon_more_than_one_fixed_exons['start']>1]
    dexon_slice=dexon_slice[['seqname','start','end','exon_id','strand']]
    dexon_fix = dexon_fix[['gene_id', 'transcript_id', 'exon_number', 'feature', 'gene_name', 'gene_biotype', 'source',
                         'transcript_name', 'transcript_biotype', 'transcript_support_level', 'exon_id', 'score']]
    dexon_slice=pd.merge(dexon_slice,dexon_fix,how='left',on='exon_id')
    if pd.__version__ >= '0.23.0':
        dexon_slice = pd.concat([dexon_slice, dexon_more_than_one_fixed_exons], sort=False)
    else:
        dexon_slice = pd.concat([dexon_slice, dexon_more_than_one_fixed_exons])
    return dexon_slice


def fix_exon_number(dexon):
    if 'exon_number' in dexon.columns:
        dexon[['exon_number']] = dexon[['exon_number']].astype(int)
        exon_col = 'exon_number'
    else:
        exon_col = 'exon_id'
        dexon = dexon.sort_values(by=['transcript_id', 'exon_id']).reset_index()
    dexon_minus = dexon[dexon.strand =='-'].copy(deep=True)
    dexon_plus = dexon[dexon.strand == '+'].copy(deep=True)
    dexon_minus = dexon_minus.sort_values(by=['transcript_id', exon_col, 'start'], ascending=[True, True, False])
    dexon_plus = dexon_plus.sort_values(by=['transcript_id', exon_col, 'start'], ascending=[True, True, True])
    if pd.__version__ >= '0.23.0':
        dexon = pd.concat([dexon_minus, dexon_plus], sort=False)
    else:
        dexon = pd.concat([dexon_minus, dexon_plus])
    if 'exon_number' in dexon.columns:
        dexon[['exon_number']] = dexon[['exon_number']].astype(str)
    del dexon_minus, dexon_plus
    dataf_dup = dexon.groupby(dexon['transcript_id'],as_index=False).size().reset_index()
    dataf_dup = dataf_dup.rename(columns={0: 'size'})
    dexon=pd.merge(dexon,dataf_dup, on='transcript_id', how='left')
    dexon_temp=dexon.copy(deep=True)
    dexon_temp['old_index']=dexon_temp.index
    dexon_temp=dexon_temp.drop_duplicates(subset='transcript_id')
    dexon=pd.merge(dexon,dexon_temp[['transcript_id','old_index']],how='left',on='transcript_id')
    dexon['exon_number']=dexon.index-dexon['old_index']+1
    del dexon['old_index'],dexon['size']
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

def build_gapped_gtf(df_gtf,output):
    df_gtf['exon_number'] = df_gtf['exon_number'].astype(float)
    df_gtf['exon_number'] = df_gtf['exon_number'].fillna(value=0, axis=0)
    df_gtf['exon_number'] = df_gtf['exon_number'].astype(int).astype(str)
    df_gtf.loc[df_gtf['feature'] == 'exon','feature']='zexon'   #Makes proper gtf sorting simpler
    df_gtf.loc[df_gtf['feature'] == 'gene','transcript_id']='A' #Makes proper gtf sorting simpler
    df_gtf=df_gtf.sort_values(by=['gene_id','transcript_id','feature'])
    df_gtf.loc[df_gtf['feature'] == 'zexon','feature']='exon'
    df_gtf=df_gtf.reset_index(drop=True)
    df_gtf['transcript_support_level']=df_gtf['transcript_support_level'].fillna(value='NA')
    df_gtf['attribute'] = ''
    for att in ['gene_id', 'gene_name', 'gene_biotype']:
        df_gtf.loc[~df_gtf[att].isnull(), 'attribute'] = df_gtf.loc[~df_gtf[att].isnull(), 'attribute'] \
                                                         + att + ' "' + df_gtf.loc[~df_gtf[att].isnull(), att] + '"; '
    for att in ['transcript_id', 'transcript_name', 'transcript_biotype', 'transcript_support_level']:
        df_gtf.loc[(df_gtf['feature'] != 'gene') & (~df_gtf[att].isnull()), 'attribute'] = \
            df_gtf.loc[(df_gtf['feature'] != 'gene') & (~df_gtf[att].isnull()), 'attribute'] + att + ' "' \
            + df_gtf.loc[(df_gtf['feature'] != 'gene') & (~df_gtf[att].isnull()), att] + '"; '
    for att in ['exon_number', 'exon_id']:
        df_gtf.loc[(df_gtf['feature'] == 'exon') & (~df_gtf[att].isnull()), 'attribute'] = \
            df_gtf.loc[(df_gtf['feature'] == 'exon') & (~df_gtf[att].isnull()), 'attribute'] + att + ' "' \
            + df_gtf.loc[(df_gtf['feature'] == 'exon') & (~df_gtf[att].isnull()), att] + '"; '
    df_gtf['attribute'] = df_gtf['attribute'].str.strip()
    df_gtf['score']='.'
    df_gtf['frame']='.'
    df_gtf=df_gtf[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame','attribute']]
    df_gtf.to_csv(path_or_buf=output,
                          index=False, sep='\t', header=False, quoting=csv.QUOTE_NONE)


def get_closest_feature(df_gtf):
    df = df_gtf[['seqname','gene_id','start','end']].copy(deep=True)
    df = df.sort_values(['seqname','end'])
    df['prev_end'] = df.end.shift(1)
    df = df.sort_values(['seqname','start'])
    df['next_start'] = df.start.shift(-1)
    return df

def fetch_overlapping_intron(df_intersect, df_gtf):
    list_emb = df_intersect.gene_id.unique()
    df_closest_plus = get_closest_feature(df_gtf[(df_gtf.feature=='exon') & (df_gtf.strand=='+')])
    df_closest_minus =  get_closest_feature(df_gtf[(df_gtf.feature=='exon') & (df_gtf.strand=='-')])
    if pd.__version__ >= '0.23.0':
        df_closest = pd.concat([df_closest_plus,df_closest_minus], ignore_index=True, sort=False)
    else:
        df_closest = pd.concat([df_closest_plus, df_closest_minus], ignore_index=True)
    df_closest = df_closest[df_closest.gene_id.isin(list_emb)]
    df_closest[['prev_end','next_start']] = df_closest[['prev_end','next_start']].astype(int)
    df_closest = df_closest.rename(columns={'gene_id':'gene_id_emb'})
    df_closest.drop(['seqname','start','end'],axis=1, inplace=True)
    df_intersect = df_intersect.rename(columns={'start_2': 'start_emb', 'end_2': 'end_emb'})
    # remove overlapping features having the same coordinates
    df_intersect = df_intersect[(df_intersect.start != df_intersect.start_emb) |
                                (df_intersect.end != df_intersect.end_emb)]
    df_merged = pd.merge(df_intersect[['gene_id', 'exon_id', 'start_emb', 'end_emb']], df_gtf, on='exon_id',
                         suffixes=['_emb', '_host'])

    df_merged['over5p'] = 0
    df_merged['over3p'] = 0
    df_merged.loc[(df_merged.start_emb > df_merged.start) & (df_merged.start_emb <= df_merged.end), 'over5p'] = 1
    df_merged.loc[(df_merged.end_emb >= df_merged.start) & (df_merged.end_emb < df_merged.end), 'over3p'] = 1
    df_contained = df_merged[(df_merged.over5p == 0) & (df_merged.over3p == 0)].copy(deep=True)
    df_contained=df_contained.rename(columns={'gene_id':'gene_id_host'})
    df_contained.drop(['over3p', 'over5p'],axis=1, inplace=True)
    df_temp = df_merged.groupby(['gene_id_emb', 'gene_id_host'], as_index=False)[[
        'start_emb', 'end_emb', 'over5p', 'over3p']].max()

    df_temp = df_temp.merge(df_temp[['gene_id_emb', 'gene_id_host']], on='gene_id_host', suffixes=['', '_other'])
    df_overlap = df_temp[['gene_id_emb', 'gene_id_host', 'over5p', 'over3p']].copy(deep=True)
    df_temp.drop(['over5p', 'over3p'], axis=1, inplace=True)
    df_temp = df_temp[df_temp['gene_id_emb'] != df_temp['gene_id_emb_other']]
    df_temp = df_temp.merge(df_gtf[df_gtf.feature == 'exon'], left_on='gene_id_emb_other', right_on='gene_id')
    df_temp.drop(['gene_id_emb_other'],axis=1, inplace=True)
    df_merged = df_merged[(df_merged.over5p != 0) | (df_merged.over3p != 0)]
    df_paired = df_merged.groupby(['gene_id_emb', 'gene_id_host'], as_index=False)[['start_emb', 'end_emb']].last()
    df_paired = df_paired.merge(df_gtf[df_gtf.feature == 'exon'], left_on='gene_id_host', right_on='gene_id')
    if pd.__version__ >= '0.23.0':
        df_paired = pd.concat([df_paired, df_temp, df_contained], sort=False)
    else:
        df_paired = pd.concat([df_paired, df_temp, df_contained])
    del df_temp, df_contained
    df_paired = df_paired.merge(df_overlap, on=['gene_id_emb', 'gene_id_host'])
    del df_overlap
    df_paired = df_paired.merge(df_closest, on='gene_id_emb')
    del df_closest
    df_paired['diff5p_intron_all'] = df_paired.start_emb - df_paired.prev_end
    df_paired['diff3p_intron_all'] = df_paired.next_start - df_paired.end_emb
    df_paired.drop(['prev_end','next_start'],axis=1,inplace=True)
    df_paired['diff5p_intron_host'] = df_paired.start_emb - df_paired.end
    df_paired['diff3p_intron_host'] = df_paired.start - df_paired.end_emb
    df_paired['diff5p_exon'] = df_paired.start_emb - df_paired.start
    df_paired['diff3p_exon'] = df_paired.end - df_paired.end_emb
    df_paired.loc[(df_paired['diff5p_intron_host'] <= 0) & (df_paired['over5p'] != 0), 'diff5p_intron_host'] = 10E7
    df_paired.loc[(df_paired['diff3p_intron_host'] <= 0) & (df_paired['over3p'] != 0), 'diff3p_intron_host'] = 10E7
    df_paired.loc[(df_paired['diff5p_intron_all'] <= 0) & (df_paired['over5p'] != 0), 'diff5p_intron_all'] = 10E7
    df_paired.loc[(df_paired['diff3p_intron_all'] <= 0) & (df_paired['over3p'] != 0), 'diff3p_intron_all'] = 10E7
    df_paired.loc[(df_paired['diff5p_exon'] <= 0) & (df_paired['over5p'] != 0), 'diff5p_exon'] = 10E7
    df_paired.loc[(df_paired['diff3p_exon'] <= 0) & (df_paired['over3p'] != 0), 'diff3p_exon'] = 10E7

    df_paired['min5p'] = df_paired[['diff5p_intron_host', 'diff5p_exon','diff5p_intron_all']].min(axis=1)
    df_paired['min3p'] = df_paired[['diff3p_intron_host', 'diff3p_exon','diff3p_intron_all']].min(axis=1)
    df_grouped = df_paired.groupby(['gene_id_emb', 'start_emb', 'end_emb', 'over5p', 'over3p'])[[
         'gene_id_host', 'min5p', 'min3p']].min().reset_index()
    del df_paired
    df_grouped['left_start'] = df_grouped.start_emb - df_grouped.min5p.astype(int)
    df_grouped['left_end'] = df_grouped.start_emb - 1
    df_grouped['right_start'] = df_grouped.end_emb + 1
    df_grouped['right_end'] = df_grouped.end_emb + df_grouped.min3p.astype(int)
    df_grouped['intron_type'] = '.retained_intron'
    df_grouped.loc[(df_grouped['over5p'] != 0) & (df_grouped['over3p'] == 0), ['right_start', 'right_end']] = -1
    df_grouped.loc[(df_grouped['over5p'] == 0) & (df_grouped['over3p'] != 0), ['left_start', 'left_end']] = -1
    df_grouped.loc[(df_grouped['over5p'] == 0) & (df_grouped['over3p'] == 0), 'left_start'] = df_grouped.end_emb \
                                                                                            + df_grouped.min3p.astype(int)
    df_grouped.loc[(df_grouped['over5p'] == 0) & (df_grouped['over3p'] == 0), 'left_end'] = df_grouped.start_emb \
                                                                                              - df_grouped.min5p.astype(int)
    df_grouped.loc[(df_grouped['over5p'] == 0) & (df_grouped['over3p'] == 0), ['right_end','right_start']] = -1
    df_grouped.loc[(df_grouped['over5p'] == 0) & (df_grouped['over3p'] == 0), 'intron_type'] = '.contained_exon'
    df_fakehost = df_grouped[df_grouped.intron_type == '.contained_exon'].copy(deep=True)
    if not df_fakehost.empty:
        df_fakehost['left_start'] = df_fakehost.start_emb + 1
        df_fakehost['left_end'] = df_fakehost.end_emb + df_fakehost.min3p.astype(int) - 1
        df_fakehost['right_end'] = df_fakehost.end_emb - 1
        df_fakehost['right_start'] = df_grouped.start_emb - df_grouped.min5p.astype(int) + 1
        df_fakehost['intron_type'] = '.fakehost'
        if pd.__version__ >= '0.23.0':
            df_grouped = pd.concat([df_grouped,df_fakehost],ignore_index=True, sort=False)
        else:
            df_grouped = pd.concat([df_grouped, df_fakehost], ignore_index=True)
    df_grouped.drop(['over5p', 'over3p', 'min5p', 'min3p', 'start_emb', 'end_emb'], axis=1, inplace=True)
    df_grouped.loc[df_grouped.left_start >= df_grouped.left_end,['left_start','left_end']] = -1
    df_grouped.loc[df_grouped.right_start >= df_grouped.right_end, ['right_start', 'right_end']] = -1
    df_grouped['trx_start'] = df_grouped.left_start
    df_grouped['trx_end'] = df_grouped.right_end
    df_grouped.loc[df_grouped.left_start ==-1, 'trx_start'] = df_grouped.right_start
    df_grouped.loc[df_grouped.right_end ==-1, 'trx_end'] = df_grouped.left_end
    df_left = df_grouped[['gene_id_emb', 'gene_id_host', 'left_start', 'left_end', 'intron_type','trx_start','trx_end']].copy(deep=True)
    df_right = df_grouped[['gene_id_emb', 'gene_id_host', 'right_start', 'right_end', 'intron_type','trx_start','trx_end']].copy(
        deep=True)
    del df_grouped
    df_left['exon_type'] = '.left'
    df_left = df_left.rename(columns={'left_start': 'exon_start', 'left_end': 'exon_end'})
    df_right['exon_type'] = '.right'
    df_right = df_right.rename(columns={'right_start': 'exon_start', 'right_end': 'exon_end'})
    if pd.__version__ >= '0.23.0':
        df_final = pd.concat([df_right, df_left], sort=False)
    else:
        df_final = pd.concat([df_right, df_left])
    df_final = df_final[(df_final.exon_start!=-1)]
    df_weird = df_final.groupby(['gene_id_emb'])['gene_id_host'].nunique().reset_index()
    weird_list = df_weird[df_weird.gene_id_host >1].gene_id_emb.tolist()
    if len(weird_list)>0:
        print('Warning! The intron correction will not be done on the following genes because of their particular situation:\n%s'
          %('\n'.join(weird_list)), file=sys.stderr)
    df_final = df_final[~df_final.gene_id_emb.isin(weird_list)]
    return df_final, weird_list


def create_gtf(df_intron, df_gtf, output):
    df_intron['temp_col'] = df_intron.gene_id_host
    df_intron.loc[(df_intron.intron_type == '.fakehost'),'gene_id_host'] = df_intron.gene_id_emb
    df_intron.loc[(df_intron.intron_type == '.contained_exon') |
                  (df_intron.intron_type == '.fakehost'), 'gene_id_emb'] = df_intron.temp_col
    df_intron.drop(['temp_col'], axis=1, inplace=True)

    gene_cols = ['gene_id','gene_name','gene_biotype','seqname','strand','source']
    gene_gtf = df_gtf[(df_gtf.feature=='gene') &
                        (df_gtf.gene_id.isin(df_intron.gene_id_host.unique()))][gene_cols]
    df_new_gtf = df_gtf[(df_gtf.feature=='gene') &
                        (df_gtf.gene_id.isin(df_intron.gene_id_host.unique()))].copy(deep=True)
    df_intron['gene_start'] = df_intron.groupby(['gene_id_host'])['trx_start'].transform('min')
    df_intron['gene_end'] = df_intron.groupby(['gene_id_host'])['trx_end'].transform('max')
    df_new_gtf = df_new_gtf.merge(df_intron[['gene_id_host','gene_start','gene_end']],
                                  left_on='gene_id', right_on='gene_id_host')
    df_new_gtf['start'] = df_new_gtf['gene_start']
    df_new_gtf['end'] = df_new_gtf['gene_end']
    df_new_gtf = df_new_gtf.drop(['gene_id_host','gene_start','gene_end'],axis=1)
    emb_gtf = df_gtf[df_gtf.gene_id.isin(df_intron[df_intron.intron_type.str.contains('intron')].gene_id_emb.unique())].copy(deep=True)
    emb_trx  = emb_gtf[emb_gtf.feature=='transcript'].copy(deep=True)
    emb_gene = emb_gtf[emb_gtf.feature=='gene']
    emb_trx['length'] = emb_trx.end - emb_trx.start
    emb_trx = emb_trx.sort_values('length',ascending=False)
    emb_trx.drop_duplicates(subset='gene_id', keep='first', inplace=True) #keep longest transcript only
    emb_exon = emb_gtf[(emb_gtf.feature=='exon') & (emb_gtf.transcript_id.isin(emb_trx.transcript_id.unique()))]
    if pd.__version__ >= '0.23.0':
        df_new_gtf = pd.concat([df_new_gtf, emb_gene, emb_trx, emb_exon], sort=False)
    else:
        df_new_gtf = pd.concat([df_new_gtf, emb_gene, emb_trx, emb_exon])
    del emb_trx, emb_exon, emb_gene, emb_gtf
    df_intron = df_intron.merge(df_gtf[df_gtf.feature == 'gene'][['gene_id', 'gene_name']], left_on='gene_id_emb',
                                right_on='gene_id')
    df_intron = df_intron.rename(columns={'gene_name': 'gene_name_emb'})

    df_intron['transcript_id'] = df_intron.gene_id_emb + df_intron.intron_type
    df_intron['transcript_name'] = df_intron.gene_name_emb + df_intron.intron_type
    df_intron['exon_id'] = df_intron.gene_id_emb + df_intron.exon_type
    df_intron['transcript_support_level'] = '1'
    df_intron['transcript_biotype'] = 'processed_transcript'
    df_trx = df_intron[['transcript_id','transcript_name','trx_start','trx_end','gene_id_host','transcript_biotype',
                        'transcript_support_level']]
    df_trx.drop_duplicates(inplace=True)
    df_trx = df_trx.rename(columns={'gene_id_host':'gene_id','trx_start':'start','trx_end':'end'})
    df_trx['feature'] = 'transcript'
    df_trx = df_trx.merge(gene_gtf, on='gene_id')
    df_exon = df_intron[['transcript_id','transcript_name','exon_start','exon_end','gene_id_host','transcript_biotype',
                        'transcript_support_level','exon_id']]
    del df_intron
    df_exon = df_exon.rename(columns={'gene_id_host': 'gene_id', 'exon_start': 'start', 'exon_end': 'end'})
    df_exon = df_exon.merge(gene_gtf, on='gene_id')
    df_exon['feature'] = 'exon'
    df_exon = fix_exon_number(df_exon)
    if pd.__version__ >= '0.23.0':
        df_new_gtf = pd.concat([df_new_gtf,df_trx,df_exon], sort=False)
    else:
        df_new_gtf = pd.concat([df_new_gtf, df_trx, df_exon])
    df_new_gtf.drop_duplicates(inplace=True)
    build_gapped_gtf(df_new_gtf, output)


def fetch_closest_exons(df_intersect, df_gtf):
    # remove overlapping features having the same coordinates
    df_intersect = df_intersect[(df_intersect.start != df_intersect.start_2) |
                                (df_intersect.end != df_intersect.end_2)]
    df_intersect['gene_length'] = df_intersect.end - df_intersect.start
    df_selected = df_intersect[df_intersect.groupby(['gene_id_emb']).gene_length.transform('min') == df_intersect.gene_length][['gene_id','gene_id_emb']]
    df_selected = df_selected.rename(columns={'gene_id':'gene_id_host'})
    list_emb = df_intersect.gene_id_emb.unique()
    df = get_closest_feature(df_gtf)
    df = df[df.gene_id.isin(list_emb)]
    df['left_start'] = df.prev_end + 1
    df['left_end'] = df.start - 1
    df['right_start'] = df.end + 1
    df['right_end'] = df.next_start - 1
    df = df.drop(['start','end','seqname','prev_end','next_start'],axis=1)
    df = df.merge(df_selected, left_on='gene_id', right_on='gene_id_emb')
    df[['left_start', 'left_end','right_start','right_end']] = df[['left_start', 'left_end','right_start','right_end']].astype(int)
    df.loc[df.right_start >= df.right_end, ['right_start','right_end']] = -1
    df.loc[df.left_start >= df.left_end, ['left_start','left_end']] = -1
    df['trx_start'] = df.left_start
    df['trx_end'] = df.right_end
    df.loc[df.trx_start == -1, 'trx_start'] = df.right_start
    df.loc[df.trx_end == -1, 'trx_end'] = df.left_end
    df['intron_type'] = '.spliced_intron'
    df_left = df[['gene_id_emb', 'gene_id_host', 'left_start', 'left_end', 'intron_type','trx_start','trx_end']].copy(deep=True)
    df_right = df[['gene_id_emb', 'gene_id_host', 'right_start', 'right_end', 'intron_type','trx_start','trx_end']].copy(
        deep=True)
    del df
    df_left['exon_type'] = '.left'
    df_left = df_left.rename(columns={'left_start': 'exon_start', 'left_end': 'exon_end'})
    df_right['exon_type'] = '.right'
    df_right = df_right.rename(columns={'right_start': 'exon_start', 'right_end': 'exon_end'})
    if pd.__version__ >= '0.23.0':
        df_final = pd.concat([df_right, df_left], sort=False)
    else:
        df_final = pd.concat([df_right, df_left])
    df_final = df_final[(df_final.exon_start != -1)]
    return df_final

def check_biotypes(df_gtf,biotypes_embedded):
    biotypes_in_gtf=df_gtf['gene_biotype'].unique()
    for biotype in biotypes_embedded:
        if biotype not in biotypes_in_gtf:
            print('Warning! gene_biotype %s is not present in the gtf.' %(biotype), file=sys.stderr)
    if len(set(biotypes_in_gtf) & set(biotypes_embedded)) == 0:
        print('No embedded genes are present in the gtf file. No need to perform any correction.')
        return False
    return True


def read_gtf(gtf_file, verbose=False):
    if verbose:
        print('Reading gtf')
    gtf_list = []
    # Open an optionally gzipped file
    fn_open = gzip.open if gtf_file.endswith('.gz') else open
    with fn_open(gtf_file) as f:
        for line in f:
            try:
                line = line.decode('utf-8')
            except AttributeError:
                pass
            if line[0] == '#':
                continue
            line = line.strip()
            line = line.split('\t')
            if line[2] not in ['gene', 'transcript', 'exon']:
                continue
            cols_to_keep = [0, 1, 2, 3, 4, 5, 6]
            attribute = line[8]
            gtf_entry = [item for idx, item in enumerate(line) if idx in cols_to_keep]
            att_list = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_biotype', 'transcript_name',
                        'exon_id', 'transcript_biotype', 'transcript_support_level', 'gene']
            for att in att_list:
                try:
                    res = re.search(att + '\s"?([^;"]+)"?;?' , attribute).group(1)
                except AttributeError:
                    if att == 'gene_id':
                        print("Some entries are missing a gene_id. Exiting.", file=sys.stderr)
                        sys.exit(1)
                    res = None
                gtf_entry.append(res)
            gtf_list.append(gtf_entry)
    df_gtf = pd.DataFrame(gtf_list,
                          columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                                   'gene_id', 'transcript_id', 'exon_number', 'gene_name',
                                   'gene_biotype', 'transcript_name', 'exon_id',
                                   'transcript_biotype', 'transcript_support_level', 'gene'])
    if len(df_gtf[df_gtf.gene_biotype.notnull()]) == 0:
        print('There is no "gene_biotype" in your gtf file, can\'t select genes to correct. Exiting.', file=sys.stderr)
        sys.exit(1)

    df_gtf.loc[(df_gtf.gene_name.isnull()) &
               (df_gtf.gene.notnull()), 'gene_name'] = df_gtf.loc[(df_gtf.gene_name.isnull()) &
                                                                  (df_gtf.gene.notnull()), 'gene']
    df_gtf[['start', 'end']] = df_gtf[['start', 'end']].astype(int)
    return df_gtf


def add_entry(merged_entries, d_merged, prev_seqname, prev_start, prev_end, prev_strand,
              gene_ids, gene_names, gene_biotype, sources):
    merged_entries.append([prev_seqname, prev_start, prev_end, prev_strand,
                           '-'.join(gene_ids), '-'.join(gene_names), gene_biotype, '-'.join(sources)])
    for i in gene_ids:
        d_merged[i] = {
            'ids': '-'.join(gene_ids),
            'names': '-'.join(gene_names)
        }
    print(','.join(gene_ids), file=sys.stderr)
    return merged_entries, d_merged

def merge_overlapping_emb_genes(df_gtf, emb_biotypes, output, min_overlap):
    # identify overlapping genes
    df_emb = df_gtf[df_gtf.gene_biotype.isin(emb_biotypes)]
    emb_bed = os.path.join(os.path.dirname(output), 'emb_genes.' + os.path.basename(output) + '.bed')
    overlapping_bed = os.path.join(os.path.dirname(output), 'overlapping_emb_genes.' + os.path.basename(output) + '.bed')
    df_emb[df_emb.feature == 'gene'][['seqname', 'start', 'end', 'gene_id', 'source', 'strand']].to_csv(emb_bed,
                                                                             sep='\t', index=False, header=False)
    intersect_command = "bedtools intersect -a %s -b %s -wa -wb -s -f %0.3f -r | awk ' $4 != $10' > %s" % (
        emb_bed, emb_bed, min_overlap, overlapping_bed)
    x = os.system(intersect_command)
    if x != 0:
        print("Merge overlapping embedded genes exited with : %d" % x, file=sys.stderr)
        sys.exit(1)
    os.remove(emb_bed)
    d_biotypes = df_emb[['gene_id', 'gene_biotype']].set_index('gene_id').gene_biotype.to_dict()
    d_names = df_emb[['gene_id', 'gene_name']].set_index('gene_id').gene_name.to_dict()
    df_overlap = pd.read_csv(overlapping_bed, sep='\t',
                             names=['seqname1', 'start1', 'end1', 'gene_id1', 'source1', 'strand1',
                                    'seqname2', 'start2', 'end2', 'gene_id2', 'source2', 'strand2'])
    os.remove(overlapping_bed)
    if df_overlap.empty:
        return df_gtf
    df_overlap['gene_biotype1'] = df_overlap.gene_id1.map(d_biotypes)
    df_overlap['gene_biotype2'] = df_overlap.gene_id2.map(d_biotypes)
    df_overlap['gene_name1'] = df_overlap.gene_id1.map(d_names)
    # only merge overlapping genes having the same gene_biotype
    df_overlap = df_overlap[df_overlap.gene_biotype1 == df_overlap.gene_biotype2]
    overlapping_genes = df_overlap.gene_id1.unique().tolist()

    bed = df_overlap[[
        'seqname1', 'start1', 'end1', 'gene_id1', 'strand1', 'gene_name1', 'gene_biotype1', 'source1'
    ]].drop_duplicates().sort_values(['strand1', 'gene_biotype1', 'seqname1', 'start1', 'end1']).values
    if len(bed) == 0:
        return df_gtf
    prev_seqname = ''
    prev_start = -1
    prev_end = -1
    prev_strand = ''
    merged_entries = []
    d_merged = {}
    prev_biotype = ''
    print('Warning! The following genes were merged because they share more than %0.2f of their sequence '
          'with respect to the -f/--fraction argument:' % min_overlap, file=sys.stderr)
    for entry in bed:
        seqname, start, end, gene_id, strand, gene_name, gene_biotype, source = entry
        start = int(start)
        end = int(end)
        if seqname != prev_seqname or strand != prev_strand or start > prev_end or prev_biotype != gene_biotype:
            if prev_seqname != '':
                merged_entries, d_merged = add_entry(merged_entries, d_merged, prev_seqname, prev_start, prev_end,
                                                     prev_strand, gene_ids, gene_names, prev_biotype, sources)
            prev_seqname = seqname
            prev_start = start
            prev_end = end
            prev_strand = strand
            prev_biotype = gene_biotype
            gene_ids = [gene_id]
            gene_names = [gene_name]
            sources = [source]
        else:
            prev_end = max([prev_end, end])
            gene_ids.append(gene_id)
            gene_names.append(gene_name)
            if source not in sources:
                sources.append(source)
    merged_entries, d_merged = add_entry(merged_entries, d_merged, prev_seqname, prev_start, prev_end,
                                         prev_strand, gene_ids, gene_names, prev_biotype, sources)
    df_merged = pd.DataFrame(merged_entries,
                             columns=['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_biotype', 'source'])
    df_merged['feature'] = 'gene'
    df_merged['score'] = '.'
    # update gtf
    df_gtf = df_gtf[~((df_gtf.gene_id.isin(overlapping_genes)) & (df_gtf.feature == 'gene'))]
    df_gtf = pd.concat([df_gtf, df_merged], ignore_index=True)
    for k in d_merged.keys():
        # k = original gene_id
        df_gtf.loc[df_gtf.gene_id == k, 'gene_name'] = d_merged[k]['names']
        df_gtf.loc[df_gtf.gene_id == k, 'gene_id'] = d_merged[k]['ids']
    return df_gtf


def correct_annotation(gtf_file, output, verbose, fraction, biotypes_embedded=('snoRNA', 'scaRNA', 'tRNA', 'miRNA', 'snRNA')):
    """
    correct_annotation builds a new gtf from the one provided, with holes on exons from genes that overlap the specified embedded biotypes.
    Read the MANUAL.md for extensive description.
    :param gtf_file: original gtf file
    :param biotypes_embedded: list of the embedded biotypes. Default: 'snoRNA','scaRNA','tRNA','miRNA' and 'snRNA'
    :output: modified gtf file with the .correct_annotation.gtf prefix.
    """
    if gtf_file.endswith('.gtf') or gtf_file.endswith('.gtf.gz'):
        print('Reading gtf')
    else:
        print('Annotation file:', gtf_file, file=sys.stderr)
        print('error: Wrong annotation format. Only .gtf and .gtf.gz files are accepted', file=sys.stderr)
        sys.exit(1)
    df_gtf = read_gtf(gtf_file, verbose)
    df_gtf.loc[(df_gtf.feature != 'gene') & (df_gtf.transcript_support_level.isnull()), 'transcript_support_level'] = '5'
    df_gtf = df_gtf[[
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'strand',
        'gene_id',
        'transcript_id',
        'exon_number',
        'gene_name',
        'gene_biotype',
        'transcript_name',
        'transcript_biotype',
        'transcript_support_level'
    ]]
    df_gtf['seqname'] = df_gtf['seqname'].map(str)
    df_gtf['start'] = df_gtf['start'].map(int)
    df_gtf['end'] = df_gtf['end'].map(int)
    gene_biotype_dict = df_gtf[df_gtf.feature == 'gene'][[
        'gene_id', 'gene_biotype']].copy(deep=True).set_index('gene_id').gene_biotype.to_dict()
    df_gtf.loc[df_gtf.gene_biotype.isnull(), 'gene_biotype'] = df_gtf.loc[df_gtf.gene_biotype.isnull(), 'gene_id'].map(
        gene_biotype_dict)
    if output == 'None':
        output = gtf_file.rsplit('.gtf', 1)[0] + '.correct_annotation.gtf'
    do_corr = check_biotypes(df_gtf,biotypes_embedded)
    df_gtf.loc[df_gtf['gene_name'].isnull(), 'gene_name'] = df_gtf['gene_id']
    df_gtf.loc[df_gtf['transcript_name'].isnull(),'transcript_name']=df_gtf['gene_name']
    df_gtf.loc[df_gtf['transcript_biotype'].isnull(),'transcript_biotype']=df_gtf['gene_biotype']
    df_gtf=df_gtf[df_gtf['feature'].isin(['gene','transcript','exon'])==True]
    df_gtf['exon_id'] = 'NaN'
    df_gtf.loc[df_gtf['feature']=='exon','exon_id'] = df_gtf['transcript_id'] + '.' + df_gtf['exon_number'].map(str).str.split('.').str.get(0)
    if do_corr is True:
        if fraction != -1:
            # merge similar overlapping annotations
            df_gtf = merge_overlapping_emb_genes(df_gtf, biotypes_embedded, output, fraction)
        dgene=df_gtf[df_gtf.feature == 'gene']
        dgene['start']=dgene['start'].map(int)
        dgene['end']=dgene['end'].map(int)
        dgene['seqname']=dgene['seqname'].map(str)
        dgene_embedded=dgene[dgene['gene_biotype'].isin(biotypes_embedded) == True]
        dgene_host=dgene[dgene['gene_biotype'].isin(biotypes_embedded) == False]

        dexon_host=df_gtf[(df_gtf.feature == 'exon') & (df_gtf['gene_id'].isin(dgene_host['gene_id'])==True)]
        dexon_not_host=df_gtf[(df_gtf.feature == 'exon') & (df_gtf['gene_id'].isin(dgene_host['gene_id'])==False)]
        dIntersect=intersect(dexon_host,dgene_embedded,output=output  + '.' + 'Intersect', name=('exon_id','gene_id'))
        dIntersect=dIntersect[dIntersect['overlap'] != -1]
        del dIntersect['overlap']

        if dIntersect.empty is False:
            df_overlapping_intron, emb_genes = fetch_overlapping_intron(dIntersect, df_gtf)
            emb_genes.extend(df_overlapping_intron.gene_id_emb.unique().tolist())
            other_emb = dgene_embedded[~dgene_embedded['gene_id'].isin(emb_genes)]
        else:
            other_emb = dgene_embedded
        other_emb = other_emb.rename(columns={'gene_id':'gene_id_emb'})
        dIntersect_gene = intersect(dgene_host,other_emb,output=output + '.' + 'Intersect',
                                    name=('gene_id','gene_id_emb'))
        dIntersect_gene = dIntersect_gene[dIntersect_gene['overlap'] != -1]
        del dIntersect_gene['overlap']
        if pd.__version__ >= '0.23.0' :
            df_minus = pd.concat([dexon_host[(dexon_host.strand == '-')], dgene_embedded[(dgene_embedded.strand == '-')]],
                                 sort=False)
            df_plus = pd.concat([dexon_host[(dexon_host.strand == '+')], dgene_embedded[(dgene_embedded.strand == '+')]],
                                sort=False)
        else:
            df_minus = pd.concat([dexon_host[(dexon_host.strand == '-')], dgene_embedded[(dgene_embedded.strand == '-')]])
            df_plus = pd.concat([dexon_host[(dexon_host.strand == '+')], dgene_embedded[(dgene_embedded.strand == '+')]])
        df_intron_minus = fetch_closest_exons(dIntersect_gene, df_minus)
        df_intron_plus = fetch_closest_exons(dIntersect_gene, df_plus)
        if dIntersect.empty is False:
            if pd.__version__ >= '0.23.0':
                df_intron = pd.concat([df_intron_minus,df_intron_plus,df_overlapping_intron],ignore_index=True, sort=False)
            else:
                df_intron = pd.concat([df_intron_minus, df_intron_plus, df_overlapping_intron], ignore_index=True)
        else:
            if pd.__version__ >= '0.23.0':
                df_intron = pd.concat([df_intron_minus, df_intron_plus], ignore_index=True, sort=False)
            else:
                df_intron = pd.concat([df_intron_minus, df_intron_plus], ignore_index=True)
        create_gtf(df_intron,df_gtf, output.replace('.gtf','.introns.gtf'))
        if dIntersect.empty is False:
            dexon_slice=drill_an_exon(dIntersect,dexon_host)
            dexon_host=dexon_host[dexon_host['exon_id'].isin(dIntersect['exon_id']) == False]
            if pd.__version__ >= '0.23.0':
                dexon_host=pd.concat([dexon_host,dexon_slice], sort=False)
                dexon=pd.concat([dexon_host,dexon_not_host], sort=False)
            else:
                dexon_host = pd.concat([dexon_host, dexon_slice])
                dexon = pd.concat([dexon_host, dexon_not_host])
            dexon=fix_exon_number(dexon)

            df_gtf=fix_gene_and_transcript_size(df_gtf,dexon)
            if pd.__version__ >= '0.23.0':
                df_gtf=pd.concat([df_gtf,dexon], sort=False)
            else:
                df_gtf = pd.concat([df_gtf, dexon])
        df_gtf = df_gtf.reset_index()
        df_gtf.loc[df_gtf.gene_name.isnull(),'gene_name'] = df_gtf.gene_id
    build_gapped_gtf(df_gtf,output)
    if do_corr is False:
        open(output.replace('.gtf','.introns.gtf'), 'w').close()
    print('All done!')


