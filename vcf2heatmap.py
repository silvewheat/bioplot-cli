# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 22:16:40 2020

Last edit on Fri Nov 20 17:14:05 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

from matplotlib.colors import LinearSegmentedColormap
from itertools import cycle
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import allel
import os
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
matplotlib.use('Agg')


def select_samples(samples_all, samples_query):
    """
    从samples_all里面筛选samples_query，产生布尔掩码
    """
    return [True if sample in samples_query else False for sample in samples_all]


def load_vcf2array(vcffile, region: None, chrom2sites: None, samples, outsamples: None):
    callset = allel.read_vcf(vcffile, region=region, samples=samples,
                             fields=['samples', 'calldata/GT', 'variants/POS'])
    gt_array = callset['calldata/GT']  # 三维array
    samples_all = callset['samples']
    pos_array = callset['variants/POS']
    chrom = region.split(':')[0]

    # 只保留sitefile中的位点
    if chrom2sites:
        selection_requiredsites = [
            True if x in chrom2sites[chrom] else False for x in pos_array]
        gt_array = gt_array[selection_requiredsites, :, :]
        pos_array = pos_array[selection_requiredsites]
        print(f'{np.sum(selection_requiredsites)} remained according to the sitefile.')

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
        gt_array_out = gt_array[:, selection_outsamples, :].reshape(
            n_sites, np.sum(selection_outsamples)*n_hap)
        selection_swtich = np.sum(gt_array_out == 0, axis=1) > np.sum(
            gt_array_out > 0, axis=1)  # ref(0)的数量比非ref(!=0)但不是miss(-1)的数量多
        print(f'Swtich REF and ALT in {np.sum(selection_swtich)} sites.')
        assert gt_array.min() >= -1
        assert gt_array.max() <= 1
        gt_swtich = gt_array[selection_swtich, :, :]
        gt_swtich[gt_swtich == 1] = 9
        gt_swtich[gt_swtich == 0] = 1
        gt_swtich[gt_swtich == 9] = 0
        gt_array[selection_swtich, :, :] = gt_swtich
    return gt_array, callset['samples'], pos_array


def cal_ref_freq(gt_array):
    """
    gt_array为load_vcf2array产生的3维ndarray
    注意，返回的是ref的frequency
    """
    return allel.GenotypeArray(gt_array).count_alleles().to_frequencies()[:, 0]


def load_groupfile(groupfile):
    group2samples = defaultdict(list)
    sample2group = {}
    samples_queried = []
    with open(groupfile) as f:
        for line in f:
            sample, group = line.strip().split()
            group2samples[group].append(sample)
            samples_queried.append(sample)
            sample2group[sample] = group
    return group2samples, sample2group, samples_queried


@click.command()
@click.option('--vcffile', required=True)
@click.option('--region', help='要输出的区域, 如12:1000-2000', default=None)
@click.option('--regionfile', help='要输出的多个区域，三列, chr\\tstart\\tend，可增加第四列作为每个区间的名字', default=None)
@click.option('--sitefile', help='前两列为chrom\\tpos。只输出包含在sitefile中的位点，可以搭配region，regionfile以及别的频率选项进一步筛选', default=None)
@click.option('--groupfile', required=True, help='两列，第一列个体ID，第二列对应的分组。结果中只保留存在于这个文件中的个体，输出图片中分组顺序会按照这个来排，如果没有使用聚类选项的话个体顺序也会按这个排列, #开头的忽略')
@click.option('--groupcolor', help='指定每个group的分组颜色，两列，一列组名，一列16进制颜色ID或者常见颜色名,不指定的话固定一组颜色循环', default=None)
@click.option('--outgroup', help='将ALT指定为outgroup中频率最高的allele, 参数为groupfile中的分组名称，多个分组的话用英文逗号(,)隔开, 如pop1,pop2', default=None)
@click.option('--min-maf-global', 'minMafGlobal', type=float, default=None, help='计算所有groupfile中的个体的总体maf, 只保留达到这一阈值的位点, 默认不过滤')
@click.option('--max-altfreq-pop', 'maxAltFreqs', nargs=2, type=(str, float), multiple=True, default=None, help='过滤指定群体中ALT的最大频率，参数格式为：pop1 0.01，该选项可以多次使用')
@click.option('--min-freqdif', 'minFreqDifs', nargs=3, type=(str, str, float), multiple=True, default=None, help='过滤指定两个群体中ref allele的最低频率差, 参数格式为：pop1 pop2 0.8')
@click.option('--outfile', help='输出文件名，格式随后缀变化，如.jpg')
@click.option('--dpi', default=300, show_default=True, help='输出图片的dpi')
@click.option('--figsize', nargs=2, default=(10, 5), show_default=True, type=float, help='图像长宽, 默认15 5')
@click.option('--outprefix', help='输出文件前缀，与regionfile配合使用', default=None)
@click.option('--outsuffix', help='输出文件后缀（格式），与regionfile配合使用', default=None)
def main(vcffile, region, regionfile, sitefile, groupfile, groupcolor, outgroup, minMafGlobal, maxAltFreqs, minFreqDifs, outfile, dpi, figsize, outprefix, outsuffix):
    """
    目前只写了除了双等位的情况，多等位被过滤掉了。多等位的画图以后再写吧。
    """
    # 读取个体群体对应信息及outgroup
    group2samples, sample2group, samples_queried = load_groupfile(groupfile)
    if outgroup:
        samples_outgroup = []
        for group in outgroup.split(','):
            samples = group2samples[group]
            if len(samples) == 0:
                print(
                    f'Warning: The outgroup "{group}" is not exist in the groupfile')
            samples_outgroup.extend(group2samples[group])
    else:
        samples_outgroup = []
    print(f'{len(samples_outgroup)} samples in outgroup.')
    # 读取需要保留的位点
    if sitefile:
        chrom2sites = defaultdict(set)
        with open(sitefile) as f:
            for line in f:
                chrom, pos = line.strip().split()[:2]
                try:
                    chrom2sites[chrom].add(int(pos))
                except ValueError:
                    pass
    else:
        chrom2sites = None
    # 读取待处理区域
    regions = []
    region_names = []
    if regionfile:
        with open(regionfile) as f:
            for line in f:
                tline = line.strip().split()
                regions.append(f'{tline[0]}:{tline[1]}-{tline[2]}')
                if len(tline) == 4:
                    region_names.append(tline[3])
                else:
                    region_names.append(f'{tline[0]}:{tline[1]}-{tline[2]}')
    else:
        regions.append(region)
        region_names.append(region)
    # 分区域处理
    nregion = 0
    for region, region_name in zip(regions, region_names):
        nregion += 1
        print(f'region{nregion}')
        # 读取vcf
        gt_array, samples_all, sites = load_vcf2array(
            vcffile, region, chrom2sites, samples_queried, samples_outgroup)
        nsites, nsamples, nhap = gt_array.shape
        print(f'loaded {gt_array.shape}')
        if len(samples_all) != len(samples_queried):
            samples_notfound = [
                x for x in samples_queried if x not in samples_all]
            print(
                f'Warning: {len(samples_notfound)} were not found in vcffile.')
            print(', '.join(samples_notfound))

        selection_sites = np.ones(nsites, dtype=bool)
        # 根据全局最小等位基因频率过滤
        if minMafGlobal:
            af_global = cal_ref_freq(gt_array)
            selection_afGlobal_1 = af_global >= minMafGlobal
            selection_afGlobal_2 = af_global < (1-minMafGlobal)
            selection_sites = np.logical_and.reduce(
                [selection_sites, selection_afGlobal_1, selection_afGlobal_2])

        # 根据部分群体的ALT频率过滤
        if maxAltFreqs:
            selection_altfreq = []
            for group, cutoff in maxAltFreqs:
                print(group, cutoff)
                samples_groupX = group2samples[group]
                selection_samples = select_samples(samples_all, samples_groupX)
                af_groupX = cal_ref_freq(
                    gt_array[:, selection_samples, :])  # 这个返回的是ref的frequency
                selection_altfreq.append(af_groupX >= (1-cutoff))
            selection_altfreq = np.logical_and.reduce(selection_altfreq)
            selection_sites = np.logical_and(
                selection_sites, selection_altfreq)

        # 根据部分群体间ref频率差过滤
        if minFreqDifs:
            selection_freqdiff = []
            for group1, group2, cutoff in minFreqDifs:
                selection_group1 = select_samples(
                    samples_all, group2samples[group1])
                selection_group2 = select_samples(
                    samples_all, group2samples[group2])
                af_group1 = cal_ref_freq(gt_array[:, selection_group1, :])
                af_group2 = cal_ref_freq(gt_array[:, selection_group2, :])
                af_diff = np.abs(af_group1 - af_group2)
                selection_freqdiff.append(af_diff >= cutoff)
            selection_freqdiff = np.logical_and.reduce(selection_freqdiff)
            selection_sites = np.logical_and(
                selection_sites, selection_freqdiff)

        # 过滤得到最终的GT矩阵
        # 第三维度加和变二维度矩阵, 便于后面画图。./. -> -2, 0/1 -> 1, 1/1 -> 2 !按单倍型画的话，直接reshape展开就行了
        gt_array = np.sum(gt_array[selection_sites, :, :], axis=2)
        sites = sites[selection_sites]
        print(f'{len(sites)} sites are remains.')

        df = pd.DataFrame(gt_array, columns=samples_all)
#        df = df.replace(-2, -1) # 之前./.两个-1相加变-2了，缺失最终还是用-1表示吧，其实不用换，因为画图的时候设定了vmax和vmin，那就不换了吧
        sample_ordered = [x for x in samples_queried if x in samples_all]
        df = df[sample_ordered]

        # 画图
        if groupcolor:
            group2color = {x.split()[0]: x.strip().split()[1]
                           for x in open(groupcolor)}
        else:
            group2color = {x: y for x, y in zip(group2samples.keys(),
                                                cycle(['#4c72b0', '#dd8452', '#55a868', '#c44e52',
                                                      '#8172b3', '#da8bc3', '#ccb974', '#64b5cd']))}
        sample2color = pd.Series(
            {x: group2color[sample2group[x]] for x in df.columns})
        # 自定义离散camp
        myColors = ((179/255, 179/255, 179/255, 1.0), (255/255, 247/255, 236 /
                    255, 1.0), (252/255, 140/255, 89/255, 1.0), (127/255, 0, 0, 1.0))
        cmap = LinearSegmentedColormap.from_list(
            'Custom', myColors, len(myColors))

        g = sns.clustermap(df.T, row_cluster=False, col_cluster=False, yticklabels=False, xticklabels=False,
                           cmap=cmap, figsize=figsize, row_colors=sample2color,
                           vmin=-1, vmax=2,
                           cbar_kws={"ticks": [-0.65, 0.1, 0.9, 1.6]})
        g.cax.set_yticklabels(['miss', 'hom-REF', 'het', 'hom-ALT'])
        # 分组row_colors对应的legend
        for label, color in group2color.items():  # 3.6以后默认都是有序字典了
            g.ax_row_dendrogram.bar(0, 0, color=color,
                                    label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc="best", ncol=2)
#        g.fig.suptitle(f'region{nregion}: {region}', x=0.05, y=0.05)
        g.fig.suptitle(f'{region_name}', x=0.05, y=0.05)
        if not outprefix:
            g.savefig(outfile, dpi=dpi)
        else:
            #g.savefig(f'{outprefix}_{nregion}.{outsuffix}', dpi=dpi)
            g.savefig(f'{outprefix}_{region_name}.{outsuffix}', dpi=dpi)
        plt.close()


if __name__ == '__main__':
    main()
