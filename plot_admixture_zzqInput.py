# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 19:53:53 2017

@author: Caiyd
"""


import json
import click
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns




@click.command()
@click.option('--infile', help='输入文件，前三列为sampleID, group1, group2, 之后为admixture的结果，一列一个成分。文件为tab分割')
@click.option('--colorfile', help='颜色配置，一行一个#开头的RGB值(16进制)')
@click.option('--outfile', help='输出文件名')
def main(infile, colorfile, outfile):
    colors = [x.strip() for x in open(colorfile)]
    df = pd.read_csv(infile, sep='\t', header=None)
    indlabels = df.iloc[:, 0].values.tolist()
    df = df.iloc[:, 3:].T
    fig, ax = plt.subplots(1, 1, figsize=(40, 10))
    sns.set_style('white')
    ind = np.arange(df.shape[1])
    wd = 1
    bt = [0] * len(ind)
    k = df.shape[0]
    for i in range(0, k):
        ax.bar(x=ind, height=df.iloc[i, :].values, width=wd, bottom=bt, color=colors[i], linewidth=0, align='center')
        bt = df.iloc[i, :].values + bt
    ax.set_xlim(-0.5, len(ind)-0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks(ind)
    ax.set_xticklabels(indlabels, rotation='vertical')
    ax.yaxis.set_visible(False)
    plt.savefig(outfile)

if __name__ == '__main__':
    main()
