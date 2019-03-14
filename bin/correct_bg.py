import pandas as pd
import os
import numpy as np
import sys
import multiprocessing as mp
import subprocess
import shutil

def prepare_bed12(bamfile, output_dir, output_name, thread):
    command = 'bash %s/prepare_bed12.sh %s %s %s %s'%(os.path.dirname(__file__),bamfile,
                                               output_dir, output_name, thread)
    x = os.system(command)
    return x


def calc_reads(row):
    row_len = [int(i.strip()) for i in row.block_len.split(',')]
    row_start = [int(i.strip()) for i in row.block_start.split(',')]
    for i in range(0, len(row_start)-1):
        if row_start[i] + row_len[i] >= row_start[i+1]:
            row['r1_start'] = ','.join([str(j) for j in row_start[:i+1]])
            row['r1_len'] = ','.join([str(j) for j in row_len[:i+1]])
            row['r2_start'] = ','.join([str(j) for j in row_start[i+1:]])
            row['r2_len'] = ','.join([str(j) for j in row_len[i+1:]])
    return row


def apply_func(args):
    df, func = args
    df = df.apply(lambda row : func(row),axis=1)
    return df


def select_longest(row):
    read1_len = [int(i) for i in row.r1_len.split(',')]
    read2_len = [int(i) for i in row.r2_len.split(',')]
    row.block_len = ','.join([str(i) for i in np.max([read1_len, read2_len], axis=0)])
    return row


def select_blocks(row):
    read1_len = [int(i) for i in row.r1_len.split(',')]
    read2_len = [int(i) for i in row.r2_len.split(',')]
    read1_start = [int(i) for i in row.r1_start.split(',')]
    read2_start = [int(i) for i in row.r2_start.split(',')]
    shared_exons = list(set(read1_start) & set(read2_start))
    if len(shared_exons)!=0:
        shared_exons_ind1 = [read1_start.index(i1) for i1 in shared_exons]
        shared_exons_ind2 = [read2_start.index(i2) for i2 in shared_exons]
        kept_block_len = read1_len[:min(shared_exons_ind1)]
        kept_block_len.extend(
            np.max([[read1_len[i] for i in shared_exons_ind1],[read2_len[j] for j in shared_exons_ind2]], axis=0))
        kept_block_len.extend(read2_len[max(shared_exons_ind2)+1:])
        kept_block_start = read1_start[:max(shared_exons_ind1)]
        kept_block_start.extend(read2_start[max(shared_exons_ind2):])
    else:
        kept_block_start = read1_start + read2_start[1:]
        kept_block_len = read1_len[:-1] + [read2_start[0] + read2_len[0]-read1_start[-1]] + read2_len[1:]
    row['block_len'] = ','.join([str(i) for i in kept_block_len])
    row['block_start'] = ','.join([str(i) for i in kept_block_start])
    return row


def correct_bed12(df, nb_threads):
    df['r1_start'] = -1
    df['r1_len'] = -1
    df['r2_start'] = -1
    df['r2_len'] = -1

    pool = mp.Pool(processes=nb_threads)
    results = pool.map(apply_func, [[df_pool, calc_reads] for df_pool in np.array_split(df, nb_threads)])
    pool.close()
    if pd.__version__ >= '0.23.0':
        df = pd.concat(list(results), sort=True)
    else:
        df = pd.concat(list(results))
    del results

    # Keep reads with no overlaps unchanged and remove from main df
    df_wo_overlap = df[df.r1_start == -1]
    df = df[df.r1_start != -1]
    df_wo_overlap = df_wo_overlap.drop(['r1_start','r2_start','r1_len','r2_len'], axis=1)

    # Keep only R1 for identical pair and remove from main df
    df_identical = df[(df.r1_start == df.r2_start) & (df.r1_len == df.r2_len)]
    df = df[(df.r1_start != df.r2_start) | (df.r1_len != df.r2_len)]
    df_identical['block_len'] = df_identical['r1_len']
    df_identical['block_start'] = df_identical['r1_start']
    df_identical = df_identical.drop(['r1_start','r2_start','r1_len','r2_len'], axis=1)

    # if the block starts are the same, keep the longest block and read1 starts, and remove from main df
    df_same_start = df[df.r1_start == df.r2_start]
    df = df[df.r1_start != df.r2_start]
    pool = mp.Pool(processes=nb_threads)
    results = pool.map(apply_func, [[df_pool, select_longest] for df_pool in np.array_split(df_same_start, nb_threads)])
    pool.close()
    if pd.__version__ >= '0.23.0':
        df_same_start = pd.concat(list(results), sort=True)
    else:
        df_same_start = pd.concat(list(results))
    del results
    df_same_start['block_start'] = df.r1_start
    df_same_start = df_same_start.drop(['r1_start','r2_start','r1_len','r2_len'], axis=1)

    # the remaining reads have partial exon overlap, so select and regroup the overlapping exon while keeping the ones
    # only represented by one read
    pool = mp.Pool(processes=nb_threads)
    results = pool.map(apply_func, [[df_pool, select_blocks] for df_pool in np.array_split(df, nb_threads)])
    pool.close()
    if pd.__version__ >= '0.23.0':
        df = pd.concat(list(results), sort=True)
    else:
        df = pd.concat(list(results))
    del results
    df = df.drop(['r1_start','r2_start','r1_len','r2_len'], axis=1)
    if pd.__version__ >= '0.23.0':
        df_bed12 = pd.concat([df_wo_overlap, df_identical, df_same_start, df], sort=True)
    else:
        df_bed12 = pd.concat([df_wo_overlap, df_identical, df_same_start, df])
    del df_wo_overlap, df_identical, df_same_start, df
    df_bed12['nb_block'] = df_bed12.block_len.str.count(',')+1
    df_bed12['nb_start'] = df_bed12.block_start.str.count(',')+1
    df_bed12 = df_bed12[df_bed12.nb_block == df_bed12.nb_start]
    df_bed12 = df_bed12.drop(['nb_start'], axis=1)
    return df_bed12


def genome_cov(output_dir, temp_dir, output, genomepath, ucsc):
    cat_file = 'cat %s/%s_chromo/corrected*.bed12 > %s/corrected_%s.bed12'%(output_dir, temp_dir, output_dir, os.path.basename(output))
    x = os.system(cat_file)
    if x !=0 :
        sys.exit('cat exit status: ' + str(x))
    shutil.rmtree('%s/%s_chromo/' % (output_dir, temp_dir))
    if ucsc is True:
        track_opt = "-trackline -trackopts 'name=\"%s\"'"%(os.path.basename(output))
    else :
        track_opt = ''
    genomecov = 'bedtools genomecov %s -split -bg -i %s/corrected_%s.bed12 -g %s > %s/%s'%(track_opt, output_dir,
                                                                                     os.path.basename(output),
                                                                                     genomepath, output_dir,
                                                                                     os.path.basename(output))
    x = os.system(genomecov)
    if x !=0 :
        sys.exit('genomecov exit status: ' + str(x))

    if ucsc is True:
        with open(genomepath) as f:
            # skips the trackline
            f.readline()
            # reads the first line of the bedgraph
            line = f.readline()
        if line[0:3] != 'chr':
            v = subprocess.check_output(['awk', '--version']).decode('utf-8').strip().split()[2].strip(',')
            if v >= '4.1.0':
                add_chr = "awk -i inplace 'NR == 1 { print; OFS=\"\t\" } NR > 1 {$1 = \"chr\" $1; print }' %s "%(output)
            else:
                add_chr = "awk 'NR == 1 { print; OFS=\"\t\" } NR > 1 {$1 = \"chr\" $1; print }' %s > %s.tmp && mv %s.tmp %s"%(output, output, output, output)                         
            os.system(add_chr)
        sedM = "sed -i 's/chrM\t/chrMT\t/' %s"%output
        os.system(sedM)
