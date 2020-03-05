# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 17:28:50 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns




def loadgroup(groupfile):
    id2groups = {}
    with open(groupfile) as f:
        for line in f:
            tline = line.strip().split()
            id2groups[tline[0]] = tline[1]
    return id2groups



@click.command()
@click.option('--vcffile', help='input genotype file')
@click.option('--outlist', help='samples list for outgroup (用于确定derived allele)', default=None)
@click.option('--querylist', help='最终画出这些个体(文本文件, 一行一个个体)')
@click.option('--region', help='要画的区域，如12:1000-2000', default=None)
@click.option('--maf', help='只画最小等位基因频率高于这个值的, default=0.2', type=float, default=0.2)
@click.option('--groupfile', help='样本分群文件, 两列, [样本ID 组名]')
@click.option('--outfile', help='输出图片')
@click.option('--figsize', nargs=2, type=float, help='图像长宽, 默认20 5', default=(20, 5))
@click.option('--ticklabelsize', help='刻度文字大小, 默认1', default=1, type=int)
@click.option('--dpi', help='图片分辨率, 默认500', default=500)
def main(vcffile, outlist, querylist, region, maf, groupfile, outfile, figsize, ticklabelsize, dpi):
    """
    从vcf文件中读取基因分型画单倍型热图
    把outgroup中频率最高的定为祖先allele
    不给outgroup就按vcf的ref和alt去画了
    """
    querysamples = [x.strip() for x in open(querylist)]
    if outlist:
        outsamples = [x.strip() for x in open(outlist)]
        vcf_outgroup = VCF(vcffile, gts012=True, samples=outsamples)
    vcf_query = VCF(vcffile, gts012=True, samples=querysamples)
    if len(querysamples) > len(vcf_query.samples):
        miss = set(querysamples) - set(vcf_query.samples)
        print(f'query sample miss: {miss}')
        for ind in miss:
            querysamples.remove(ind)
    df = []
    index = []
    if outlist:
        # 定outgroup的allele
        for variant_outgroup, variant_query in zip(vcf_outgroup(region),
                                                   vcf_query(region)):
            counts = np.bincount(variant_outgroup.gt_types) # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
            try:
                major_gt = np.argmax([counts[0], counts[2]]) # 比较0和2哪个多
            except IndexError: # 没有HOM_ALT
                major_gt = 0
            arr = variant_query.gt_types
            if major_gt == 0:
                # 如果major allele是alt，则对换ref和alr
                arr[arr==2] = -9
                arr[arr==0] = 2
                arr[arr==-9] = 0
            df.append(arr.tolist())
            index.append(variant_query.POS)
    else:
        for variant_query in vcf_query(region):
            arr = variant_query.gt_types
            df.append(arr.tolist())
            index.append(variant_query.POS)

    df = pd.DataFrame(df, columns=vcf_query.samples, index=index)
    df = df[querysamples] # 排序
    df = df.replace(3, np.nan) # 自动 int to float
    print(f'{os.path.basename(vcffile)} {region}:\n{df.shape}')

    # 改名
    id2groups = loadgroup(groupfile)
    df.columns = [f'{id2groups[x]}_{x}' for x in df.columns]

    # MAF筛选
    freqs = df.sum(axis=1).values / (df.count(axis=1).values * 2)
    df = df.loc[((1-maf)>=freqs)&(freqs>=maf), :]
    print(f'filter maf({maf}):\n{df.shape}')

    # 画图
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_facecolor("grey")
    sns.heatmap(df.T, yticklabels=1, cmap='OrRd', ax=ax)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(ticklabelsize)
    plt.savefig(outfile, dpi=dpi)
    plt.close()


if __name__ == '__main__':
    main()

