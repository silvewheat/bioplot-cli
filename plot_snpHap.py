# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:13:00 2017

@author: Caiyd
"""


import os
import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def load_bcf2df(bcffile, region, regions_file):
    snp_list = []
    index_list = []
    sw_dict = {r'1|1': 2,
               r'1|0': 1,
               r'0|1': 1,
               r'0|0': 0,
               r'1/1': 2,
               r'0/1': 1,
               r'0/0': 0,
               r'./.': np.nan}
    if regions_file:
        cmd = f'''bcftools view -i 'TYPE=="snp" & N_ALT==1' -R {regions_file} {bcffile}'''
    elif region:
        cmd = f'''bcftools view -i 'TYPE=="snp" & N_ALT==1' {bcffile} {region}'''
    else:
        cmd = f'''bcftools view -i 'TYPE=="snp" & N_ALT==1' {bcffile}'''
    for line in os.popen(cmd):
        if line[0] != '#':
            line = line.strip().split()
            snp_list.append([sw_dict[x.split(':')[0]] for x in line[9:]])
            index_list.append('%s:%s' % (line[0], line[1]))
        elif line[:6] == '#CHROM':
            name_list = line.strip().split('\t')[9:]
    return pd.DataFrame(snp_list, columns=name_list, index=index_list)


def load_sample_order(order_file):
    sample_order_list = []
    with open(order_file) as f:
        for line in f:
            if line[0] != '#':
                line = line.strip()
                sample_order_list.append(line)
    return sample_order_list

def filter_af(df, afcutoff):
    """
    allele freq
    过滤掉ALT频率低于cutoff的位点
    """
    print(f'filter alt allele frequence, cutoff is {afcutoff}')
    print(f'before: {df.shape}')
    afarray = np.nansum(df.values, axis=1) / (np.sum(~np.isnan(df.values), axis=1) * 2)
    newdf = df.loc[afarray>=afcutoff, :]
    print(f'after: {newdf.shape}')
    return newdf



def filter_maf(df, afcutoff):
    """
    allele freq
    过滤掉minor allele frequency频率低于cutoff的位点
    """
    print(f'filter minor allele frequence, cutoff is {afcutoff}')
    print(f'before: {df.shape}')
    afarray = np.nansum(df.values, axis=1) / (np.sum(~np.isnan(df.values), axis=1) * 2)
    afarray[afarray>0.5] = 1 - afarray[afarray>0.5]
    print(f'max MAF: {max(afarray)}')
    print(f'min MAF: {min(afarray)}')
    newdf = df.loc[afarray>=afcutoff, :]
    print(f'after: {newdf.shape}')
    return newdf



def plot(df, font_scale, outfile):
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(30, 10))
    ax.set_facecolor('grey')
#    sns.set(font_scale=font_scale)
    sns.heatmap(df.T,
                cmap='Oranges',
                square=False,
                cbar=False,
                linewidths=0,
                linecolor='k',
                ax=ax)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(font_scale)
    plt.savefig(outfile, dpi=500)


@click.command()
@click.option('--bcffile', help='输入的bcf文件')
@click.option('--region', help='要画的区域，如12:1000-2000', default=None)
@click.option('--regions-file', help='区域文件，默认None，覆盖--region', default=None)
@click.option('--orderfile', help='个体ID顺序，画出的图会按这个排序，一行一个个体ID')
@click.option('--afcutoff', type=float, help='alt的频率，低于这个频率的snp不画', default=None)
@click.option('--mafcutoff', type=float, help='minor allele的频率，低于这个频率的snp不画', default=None)
@click.option('--font-scale', help='坐标轴文字大小尺度, 默认是1', default=1, type=float)
@click.option('--outfile', help='输出文件')
def main(bcffile, region, regions_file, orderfile, afcutoff, mafcutoff, font_scale, outfile):
    """
    从bcf文件中画snp的单倍型热图
    """
    print('---plot_snpHap---')
    df = load_bcf2df(bcffile, region, regions_file)
    print('load snp done.')
    sample_order_list = load_sample_order(orderfile)
    df = df[sample_order_list]
    print('sort order done.')
    if afcutoff:
        df = filter_af(df, afcutoff)
    if mafcutoff:
        df = filter_maf(df, mafcutoff)
    plot(df, font_scale, outfile)


if __name__ == '__main__':
    main()
