# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:39:34 2017

@author: Caiyd
"""


import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def load_mismatch(mismatchfile):
    df = pd.read_csv(mismatchfile, sep='\t')
    return df


def c2t_ratio(df):
    return np.sum(df['T'].values) / np.sum(df[['A', 'C', 'G', 'T']].values)


def g2a_ratio(df):
    return np.sum(df['A'].values) / np.sum(df[['A', 'C', 'G', 'T']].values)


def allmu_ratio(df):
    n_right = 0 # 正确的碱基数
    for idx, base in zip(range(4), ('A', 'C', 'G', 'T')):
        n_right += np.sum(df.loc[df['Ref'] == idx, base].values)
    return 1 - n_right / np.sum(df[['A', 'C', 'G', 'T']].values)


def cal_freq(df):
    p5_c2t = []
    p5_g2a = []
    p5_other = []
    p3_c2t = []
    p3_g2a = []
    p3_other = []
    for i in range(25):
        p5df = df.loc[((df['posi'] == i) & (df['strand']==0))]
        p5_c2t.append(c2t_ratio(p5df.loc[p5df['Ref'] == 1, :]))
        p5_g2a.append(g2a_ratio(p5df.loc[p5df['Ref'] == 2, :]))
        p5_other.append(allmu_ratio(p5df.loc[(p5df['Ref'] == 0) | (p5df['Ref'] == 3), :]))

        p3df = df.loc[((df['posi'] == i) & (df['strand']==1))]
        p3_c2t.append(c2t_ratio(p3df.loc[p3df['Ref'] == 1, :]))
        p3_g2a.append(g2a_ratio(p3df.loc[p3df['Ref'] == 2, :]))
        p3_other.append(allmu_ratio(p3df.loc[(p3df['Ref'] == 0) | (p3df['Ref'] == 3), :]))
    p3_c2t.reverse()
    p3_g2a.reverse()
    p3_other.reverse()
    return p5_c2t, p5_g2a, p5_other, p3_c2t, p3_g2a, p3_other


def plot(p5_c2t, p5_g2a, p5_other, p3_c2t, p3_g2a, p3_other, outprefix):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3))
    plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None,
                    wspace=0.06, hspace=None)
    for ay, label, color in zip((p5_c2t, p5_g2a, p5_other),
                                 ('C>T', 'G>A', 'others'),
                                 ('red', 'blue', 'grey')):

        ax[0].plot(np.arange(1, 26, 1), ay, linestyle='solid', label=label, lw=2, color=color)
    ax[0].set_ylim([0, 0.3])
    ax[0].set_xlabel("Distance from the 5' end")
    for ay, label, color in zip((p3_c2t, p3_g2a, p3_other),
                                 ('C>T', 'G>A', 'others'),
                                 ('red', 'blue', 'grey')):
        ax[1].plot(np.arange(1, 26, 1), ay, linestyle='solid', label=label, lw=2, color=color)
    ax[1].set_ylim([0, 0.3])
    ax[1].set_xticks(np.arange(0, 26, 5))
    ax[1].set_xticklabels([str(x) for x in np.arange(25, -1, -5)])
    ax[1].get_yaxis().tick_right()
    ax[1].set_xlabel("Distance from the 3' end")
    plt.savefig(f'{outprefix}.pdf', dpi=500)


@click.command()
@click.option('--mismatchfile', help='angsd -doMisMatch 1 的输出文件')
@click.option('--outprefix', help='输出图片文件的前缀')
def main(mismatchfile, outprefix):
    df = load_mismatch(mismatchfile)
    p5_c2t, p5_g2a, p5_other, p3_c2t, p3_g2a, p3_other = cal_freq(df)
    plot(p5_c2t, p5_g2a, p5_other, p3_c2t, p3_g2a, p3_other, outprefix)


if __name__ == '__main__':
    main()
