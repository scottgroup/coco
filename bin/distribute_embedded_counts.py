import pandas as pd
import distribute_multireads

pd.set_option('expand_frame_repr',False)


def correct_embedded(df_gtf, gene_file, intron_file, outfile,count_type):
    df_gtf = df_gtf[df_gtf.feature=='transcript']
    if count_type == 'uniqueOnly':
        df_gene = distribute_multireads.read_count_matrix(gene_file, 'gene')
        df_counts = distribute_multireads.read_count_matrix(intron_file, 'intron')
        df_counts = df_counts.rename(columns={'gene_id':'transcript_id'})
    else:
        df_gene = pd.read_csv(gene_file,sep='\t',names=['gene_id','accumulation'])
        df_counts = pd.read_csv(intron_file, sep='\t',
                                names=['transcript_id', 'accumulation','length'])
    df_counts[['length']] = df_counts[['length']].astype(float)
    df_counts = df_counts.merge(df_gtf[['gene_id','transcript_id','gene_name']], on='transcript_id')
    df_intron = df_counts[(df_counts.transcript_id.str.contains('intron')) |
                          (df_counts.transcript_id.str.contains('.fakehost'))].copy(deep=True)
    df_intron['reads_per_nt'] = df_intron['accumulation'] / df_intron['length']
    df_embedded = df_counts[~(df_counts.transcript_id.str.contains('intron')) &
                            ~(df_counts.transcript_id.str.contains('.fakehost')) &
                            ~(df_counts.transcript_id.str.contains('.contained_exon'))].copy(deep=True)
    df_contained = df_counts[(df_counts.transcript_id.str.contains('.contained_exon'))].copy(deep=True)
    df_embedded = df_embedded.drop(['accumulation'],axis=1)
    df_embedded = df_embedded.merge(df_gene[['gene_id','accumulation']], on='gene_id')
    if pd.__version__ >= '0.23.0':
        df_embedded = pd.concat([df_embedded, df_contained], sort=False)
    else:
        df_embedded = pd.concat([df_embedded, df_contained])
    df_embedded['reads_per_nt'] = df_embedded['accumulation'] / df_embedded['length']
    del df_counts
    df_intron['embedded_id'] = df_intron.transcript_id.str.split('.').str.get(0)
    df_merged = df_intron.merge(df_embedded, left_on='embedded_id', right_on='gene_id',
                                suffixes=['_intron','_emb'])
    df_merged['host_diff'] = 0.0
    df_merged['emb_diff'] = 0.0
    del df_intron, df_embedded
    # if the read distribution in the intron is higher than or equal to the embedded gene, the counts go to the host
    # if the embedded gene is in a spliced intron, only correct the embedded, because intron counts are not counted
    df_merged.loc[(df_merged.transcript_id_intron.str.contains('.retained_intron')) &
                  (df_merged.reads_per_nt_intron >= df_merged.reads_per_nt_emb), 'host_diff'] = df_merged.accumulation_emb
    df_merged.loc[(df_merged.transcript_id_intron.str.contains('intron'))
                  & (df_merged.reads_per_nt_intron >= df_merged.reads_per_nt_emb), 'emb_diff'] = -df_merged.accumulation_emb

    # if the host exon is totally included in the embedded gene(ex: ENSG00000142002), the embedded gene becomes a
    # fakehost and its counts are redistributed to the contained exon of the "real host" gene
    df_merged.loc[(df_merged.transcript_id_intron.str.contains('.fakehost')) &
                  (df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb), 'host_diff'] = -df_merged.accumulation_emb
    df_merged.loc[(df_merged.transcript_id_intron.str.contains('.fakehost'))
                  & (df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb), 'emb_diff'] = df_merged.accumulation_emb

    # if the read distribution is higher in the embedded gene, remove the baseline equal to the intron distribution
    df_merged.loc[(df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb) &
                  (df_merged.transcript_id_intron.str.contains('intron')), 'emb_diff'] = -df_merged.reads_per_nt_intron * df_merged.length_emb
    # if the overlapping feature is a retained_intron, the counts go to the host gene
    df_merged.loc[(df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb) &
                  (df_merged.transcript_id_intron.str.contains('.retained_intron')), 'host_diff'] = df_merged.reads_per_nt_intron*df_merged.length_emb
    # if the overlapping feature is a contained exon,
    df_merged.loc[(df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb) &
                  (df_merged.transcript_id_intron.str.contains('.fakehost')), 'emb_diff'] = df_merged.reads_per_nt_intron * df_merged.length_emb
    df_merged.loc[(df_merged.reads_per_nt_intron < df_merged.reads_per_nt_emb) &
                  (df_merged.transcript_id_intron.str.contains('.fakehost')), 'host_diff'] = -df_merged.reads_per_nt_intron * df_merged.length_emb
    intron_diff = df_merged.groupby('gene_id_intron')['host_diff'].sum().reset_index()

    df_gene = df_gene.merge(intron_diff[['gene_id_intron','host_diff']], left_on='gene_id',right_on='gene_id_intron',
                            how='left')
    df_gene = df_gene.merge(df_merged[['gene_id_emb','emb_diff']], left_on='gene_id',right_on='gene_id_emb',
                            how='left')
    del df_merged, intron_diff
    df_gene.loc[df_gene.host_diff.notnull(),'accumulation'] = df_gene.accumulation + df_gene.host_diff
    df_gene.loc[df_gene.emb_diff.notnull(),'accumulation'] = df_gene.accumulation + df_gene.emb_diff
    # Since the intron part is an approximation, the accumulation resulting may be under 0 (very rarely),
    # so we change it to zero to avoid problems afterwards.
    df_gene.loc[df_gene.accumulation < 0, 'accumulation'] = 0
    df_gene = df_gene.sort_values('gene_id')

    df_gene[['gene_id','accumulation']].to_csv(outfile, index=False, sep='\t', header=False)
    del df_gene
