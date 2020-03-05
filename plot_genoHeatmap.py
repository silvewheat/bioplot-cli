# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 16:04:22 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import re
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter



def loadgroup(groupfile):
    id2groups = {}
    with open(groupfile) as f:
        for line in f:
            tline = line.strip().split()
            id2groups[tline[0]] = tline[1]
    return id2groups


def findmajor(genos):
    "genos = ['A|A', 'A|T', 'T|T']"
    alleles = [allele for geno in genos for allele in re.split(r'\||/', geno)]
    majorallele = Counter(alleles).most_common()[0][0]
    return majorallele


def geno2count(genos, allele, mask='N'):
    """
    转换genos列表为目标allele计数的列表
    genos = np.array(['N/N', 'T|T', 'A|T', 'A|A'])
    allele = 'A'
    mask = 'N'
    return [-2, 0, 1, 2]
    """
    counts = np.array([x.count(allele) for x in genos])
    # 分别对 N/N N|N N 三种可能的缺失进行屏蔽
    counts[genos.values == f'{mask}/{mask}'] = -2
    counts[genos.values == f'{mask}|{mask}'] = -2
    counts[genos.values == f'{mask}'] = -2
    return counts

def replaceGeno2Counts(genos, samples, mask='N'):
    """
    将原来dataframe中的geno信息转换为samples中的major allele的计数
    samples一般为外群
    """
    allele = findmajor(genos[samples])
    return geno2count(genos, allele, mask)





@click.command()
@click.option('--genofile', help='input genotype file')
@click.option('--outgroup', help='samples list for outgroup (用于确定derived allele)')
@click.option('--indlist', help='最终画出这些个体(文本文件, 一行一个个体)')
@click.option('--chrom', help='画图区域所在染色体', type=str)
@click.option('--start', help='画图区域的染色体起始位置', type=int)
@click.option('--end', help='画图区域的染色体终止位置', type=int)
@click.option('--maf', help='只画最小等位基因频率高于这个值的, default=0.2', type=float, default=0.2)
@click.option('--groupfile', help='样本分群文件, 两列, [样本ID 组名]')
@click.option('--mask', help='miss使用的字母(对miss进行屏蔽), default=N', default='N')
@click.option('--outfile', help='输出图片')
@click.option('--figsize', nargs=2, type=float, help='图像长宽, 默认15 5', default=(15, 5))
def main(genofile, outgroup, indlist, chrom, start, end, maf, groupfile, mask, outfile):
    """
    从genofile中读取基因分型画单倍型热图
    genofile使用genomics_general/VCF_processing/parseVCF.py生成
    把outgroup中频率最高的定为祖先allele
    """
    samples = [x.strip() for x in open(indlist)]
    outgroup = [x.strip() for x in open(outgroup)]
    cols = ['#CHROM', 'POS'] + samples
    df = pd.read_csv(genofile, sep='\t', dtype={'#CHROM': str}, usecols=cols)
    print(f'{genofile}:\n{df.shape}')
    id2groups = loadgroup(groupfile)
    df = df.loc[df['#CHROM']==chrom, :]
    df = df.drop(labels='#CHROM', axis=1).set_index('POS')[samples]
    #df = df.drop(columns='#CHROM').set_index('POS')[samples]  New in version 0.21.0.
    df = df.loc[(df.index>=start)&(df.index<=end), :]
    print(f'filter {chrom}:{start}-{end}:\n{df.shape}')
    if int(pd.__version__.split('.')[1]) >= 23:
        df[samples] = df[samples].apply(replaceGeno2Counts, args=(outgroup, mask), axis=1, result_type='broadcast')
    else:
        df[samples] = df[samples].apply(replaceGeno2Counts, args=(outgroup, mask), axis=1, broadcast=True)
    df = df.astype(int)
    freqs = df.apply(lambda x: np.sum(x[x!=-2]), axis=1) / (df.shape[1]*2)
    df = df.replace(-2, np.nan)
    df.columns = [f'{id2groups[x]}_{x}' for x in df.columns]
    fig, ax = plt.subplots(1, 1, figsize=(100, 18))
    ax.set_facecolor("grey")
    df = df.loc[((1-maf)>=freqs)&(freqs>=maf), :]
    print(f'filter maf({maf}):\n{df.shape}')
    sns.heatmap(df.T, yticklabels=2, cmap='OrRd', ax=ax)
    plt.savefig(outfile, dpi=300)
    plt.close()


if __name__ == '__main__':
    main()



