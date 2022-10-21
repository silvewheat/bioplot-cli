# -*- coding: utf-8 -*-
"""
Created on Wed May 12 16:28:00 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import os
import typer
import allel
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from itertools import cycle
from matplotlib.colors import LinearSegmentedColormap


def select_samples(samples_all, samples_query):
    """
    从samples_all里面筛选samples_query，产生布尔掩码
    """
    return [True if sample in samples_query else False for sample in samples_all]


def load_groupfile(groupfile):
    group_ordered = []
    samples_queried = []
    haps_queried = []
    group2haps = defaultdict(list)
    hap2group = {}
    group2samples = defaultdict(list)
    with open(groupfile) as f:
        for line in f:
            sample, group = line.strip().split()
            if group not in group_ordered:
                group_ordered.append(group)  
            samples_queried.append(sample)
            group2samples[group].append(sample)
            haps_queried.append(f'{sample}_1')
            haps_queried.append(f'{sample}_2')
            group2haps[group].append(f'{sample}_1')
            group2haps[group].append(f'{sample}_2')
            hap2group[f'{sample}_1'] = group
            hap2group[f'{sample}_2'] = group
    return samples_queried, haps_queried, group_ordered, group2haps, hap2group, group2samples


def load_vcf2array(vcffile, region=None, samples=None, outsamples=None):
    callset = allel.read_vcf(vcffile, region=region, samples=samples,
                             fields=['samples', 'calldata/GT', 'variants/POS'])
    gt_array = callset['calldata/GT'] # 三维array
    n_sites, n_samples, n_hap = gt_array.shape
    print(f'{n_sites} sites were loaded.')
    samples_all = callset['samples']
    pos_array = callset['variants/POS']
    
    # 只保留双等位，因为对于多等位杂合后续没法区分，如0/2和 1/1加了之后都会转化为2
    selection_biallelic = np.max(np.max(gt_array, axis=2), axis=1) < 2
    gt_array = gt_array[selection_biallelic, :, :]
    pos_array = pos_array[selection_biallelic]
    n_sites, n_samples, n_hap = gt_array.shape
    print(f'{n_sites} biallelic sites were remained.')

    # 把outsamples中最高频率的allele设置为alt
    if outsamples:
        selection_outsamples = select_samples(samples_all, outsamples)
        print(f'{np.sum(selection_outsamples)} outgroup samples in vcf file.')
        gt_array_out = gt_array[:, selection_outsamples, :].reshape(n_sites, np.sum(selection_outsamples)*n_hap)
        selection_swtich = np.sum(gt_array_out==0, axis=1) > np.sum(gt_array_out>0, axis=1) # ref(0)的数量比非ref(!=0)但不是miss(-1)的数量多
        print(f'Swtich REF and ALT in {np.sum(selection_swtich)} sites.')
        assert gt_array.min() >= -1
        assert gt_array.max() <= 1
        gt_swtich = gt_array[selection_swtich, :, :]
        gt_swtich[gt_swtich == 1] = 9
        gt_swtich[gt_swtich == 0] = 1
        gt_swtich[gt_swtich == 9] = 0
        gt_array[selection_swtich, :, :] = gt_swtich
    
    # 按单倍体展开为2维矩阵, mat_haplo第一个维度是单倍体，第二个维度是位点。
    haplotypes_1 = gt_array[:,:,0]
    haplotypes_2 = gt_array[:,:,1]
    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T
    
    hapIDs = []
    for sample in callset['samples']:
        hapIDs.append(f'{sample}_1')
        hapIDs.append(f'{sample}_2')
    
    return mat_haplo, hapIDs


def main(vcffile: str = typer.Option(..., help="vcf文件"),
         regionfile: str = typer.Option(..., help="要画区域，两列，第一列为chrom:start-end,第二列为该区域输出的名字"),
         groupfile: str = typer.Option(..., help="两列，第一列个体ID，第二列对应的分组。样本会根据这个文件筛选，分组的排列顺序与这个文件一致"),
         outprefix: str = typer.Option(..., help="输出文件前缀"),
         mafcutoff: float = typer.Option(0.05, help="全局maf阈值"),
         groupcolor: str = typer.Option(..., help="分组的颜色，两列，一列组名，一列#开头的HEX"),
         outgroup: str = typer.Option(None, help="将ALT指定为outgroup中频率最高的allele, 参数为groupfile中的分组名称，多个分组的话用英文逗号(,)隔开, 如pop1,pop2")):
    """
    按单倍体画热图，分组内进行聚类

    """
    samples_queried, haps_queried, group_ordered, group2haps, hap2group, group2samples = load_groupfile(groupfile)
    
    regions = []
    names = []
    with open(regionfile) as f:
        for line in f:
            region, name = line.strip().split()
            regions.append(region)
            names.append(name)
    
    group2color = {x.split()[0]: x.strip().split()[1] for x in open(groupcolor)}     
    
    outsamples = []
    if outgroup:
        for group in outgroup.split(','):
            outsamples.extend(group2samples[group])

    for region, name in zip(regions, names):
        print(name, region)
        hap_array, hapIDs = load_vcf2array(vcffile, region, samples_queried, outsamples)
        df = pd.DataFrame(hap_array.T, columns=hapIDs)
        # 排序
        haps_sorted = []
        for group in group_ordered:
            haps = df[group2haps[group]].sum().sort_values(ascending=False).index
            haps_sorted.extend(haps)
        df = df[haps_sorted]
        # 一定要用series，不能用dict
        hap2color = pd.Series({x: group2color[hap2group[x]] for x in df.columns})
        # maf
        maf = df.sum(axis=1) / len(haps_sorted)
        # 画图
        figsize=(20, 5)
        myColors = ((179/255, 179/255, 179/255, 1.0), (255/255, 247/255, 236/255, 1.0),  (127/255, 0, 0, 1.0))
        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
        
        g = sns.clustermap(df.loc[maf>=mafcutoff, haps_sorted].T, 
                           row_cluster=False, col_cluster=False,
                           yticklabels=False, xticklabels=False,
                           cmap=cmap, figsize=figsize, row_colors=hap2color,
                           vmin=-1, vmax=1, cbar_pos=(0.5, 0.8, 0.03, 0.1),
                           cbar_kws={"ticks":[-0.7, 0, 0.7]})
        
        # legend
        for label, color in group2color.items(): # 3.6以后默认都是有序字典了
            g.ax_row_dendrogram.bar(0, 0, color=color,
                                    label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc='lower left', ncol=3, bbox_to_anchor=(1, 1))
        g.cax.set_yticklabels(['miss', 'REF', 'ALT'])
        g.fig.suptitle(f'{region}\n{name}', x=0.8, y=0.9)
        plt.savefig(f'{outprefix}_{name}.jpg', dpi=300)
        plt.close()

if __name__ == '__main__':
    typer.run(main)
