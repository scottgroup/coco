import pandas as pd
import GTF
import sys


def read_count_matrix(count_file):
    df= pd.read_csv(count_file, sep='\t', comment='#',
                            names=['gene_id', 'seqname', 'start', 'end', 'strand', 'length', 'accumulation'])
    df = df.drop(['seqname', 'start', 'end', 'strand', 'length'], axis=1)
    df['accumulation'] = pd.to_numeric(df['accumulation'],
                                              errors='coerce')  # if header is present in file, column name will be changed to NaN,
    df = df.dropna(axis=0)  # and row will be removed
    return df

def ratio_mmg(input_file, gtf_file, unique_counts, output_file):
    if gtf_file.endswith('.gtf'):
        print('Reading gtf')
    else:
        print('Annotation file:',gtf_file)
        print('error: Wrong annotation format. Only .gtf files are accepted')
        sys.exit(1)
    try:
        df_gtf=GTF.dataframe(gtf_file)
    except:
        print("Error: gtf file cannot be converted to dataframe. Make sure the annotation file provided is in gene \
        transfer format (.gtf) and comes from Ensembl.")
        sys.exit()
    df_gtf = df_gtf[df_gtf.feature=='gene']
    df_unique = read_count_matrix(unique_counts)
    df_unique = df_unique.merge(df_gtf, on='gene_id', how='outer')
    df_unique = df_unique.fillna(0.0)
    df_unique = df_unique.set_index('gene_id')
    dict_unique = df_unique['accumulation'].to_dict()
    dict_multi = dict.fromkeys(dict_unique, 0.0)
    with open(input_file, 'r') as f:
        for line in f:
            line = line.split('\t')
            group_genes = line[0].split(',')
            group_count = float(line[1])
            if len(group_genes) > 1:
                unique_tot = 0.0
                for g in group_genes:
                    unique_tot += dict_unique[g]
                if unique_tot == 0:
                    for g in group_genes:
                        dict_multi[g] += group_count * 1 / len(group_genes)
                else:
                    for g in group_genes:
                        dict_multi[g] += group_count * dict_unique[g]/unique_tot
            else:
                dict_multi[group_genes[0]] += group_count
    df_multi = pd.DataFrame.from_dict(dict_multi, orient='index')
    df_multi.columns = ['count_multi']
    df_merged = pd.merge(df_unique, df_multi, left_index=True, right_index=True)
    df_merged['tot'] = df_merged.accumulation + df_merged.count_multi
    outfile = output_file
    df_merged['tot'].to_csv(outfile, sep='\t', header=False)


