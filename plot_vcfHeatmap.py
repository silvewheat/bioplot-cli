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




def loadrename(rename):
    old2new = {}
    with open(rename) as f:
        for line in f:
            tline = line.strip().split()
            old2new[tline[0]] = tline[1]
    return old2new



@click.command()
@click.option('--vcffile', help='input genotype file')
@click.option('--outlist', help='samples list for outgroup (用于确定derived allele)', default=None)
@click.option('--querylist', help='最终画出这些个体(文本文件, 一行一个个体)')
@click.option('--region', help='要画的区域，如12:1000-2000', default=None)
@click.option('--regionfile', help='要画的多个区域，每行会输出一个文件，三列, chr\\tstart\\tend', default=None)
@click.option('--maf', help='只画最小等位基因频率高于这个值的, default=0.2', type=float, default=0.2)
@click.option('--rename', help='修改输出的样本ID, 两列, [旧ID 新ID]', default=None)
@click.option('--outfile', help='输出图片')
@click.option('--figsize', nargs=2, type=float, help='图像长宽, 默认20 5', default=(20, 5))
@click.option('--ticklabelsize', help='刻度文字大小, 默认1', default=1, type=int)
@click.option('--dpi', help='图片分辨率, 默认500', default=500)
@click.option('--outprefix', help='输出文件前缀，与regionfile配合使用', default=None)
@click.option('--outsuffix', help='输出文件后缀（格式），与regionfile配合使用', default=None)
def main(vcffile, outlist, querylist, region, regionfile, maf, rename, outfile, figsize, ticklabelsize, dpi, outprefix, outsuffix):
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
    regions = []
    if regionfile:
        with open(regionfile) as f:
            for line in f:
                tline = line.strip().split()
                regions.append(f'{tline[0]}:{tline[1]}-{tline[2]}')
    else:
        regions.append(region)
    nregion = 0
    for region in regions:
        df = []
        index = []
        nregion += 1
        if outlist:
            # 定outgroup的allele
            for variant_outgroup, variant_query in zip(vcf_outgroup(region),
                                                       vcf_query(region)):
                if (len(variant_query.ALT) == 1) and (len(variant_query.REF) == 1) and (len(variant_query.ALT[0]) == 1): # 只保留双等位SNP
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
                if (len(variant_query.ALT) == 1) and (len(variant_query.REF) == 1) and (len(variant_query.ALT[0]) == 1): # 只保留双等位SNP
                    arr = variant_query.gt_types
                    df.append(arr.tolist())
                    index.append(variant_query.POS)

        df = pd.DataFrame(df, columns=vcf_query.samples, index=index)
        df = df[querysamples] # 排序
        df = df.replace(3, np.nan) # 自动 int to float
        print(f'{os.path.basename(vcffile)} {region}:\n{df.shape}')

        # 改名
        if rename:
            old2new = loadrename(rename)
            df.columns = [old2new[x] for x in df.columns]

        # MAF筛选
        freqs = df.sum(axis=1).values / (df.count(axis=1).values * 2)
        df = df.loc[((1-maf)>=freqs)&(freqs>=maf), :]
        print(f'filter maf({maf}):\n{df.shape}')

        # 画图
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_facecolor("grey")
        sns.heatmap(df.T, yticklabels=1, cmap='OrRd', ax=ax)
        ax.set_title(region)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(ticklabelsize)
        if not outprefix:
            plt.savefig(outfile, dpi=dpi)
        else:
            plt.savefig(f'{outprefix}_{nregion}.{outsuffix}')
        plt.close()


if __name__ == '__main__':
    main()

