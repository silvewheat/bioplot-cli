# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 22:10:10 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import json
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns




def load_data(infile, group_col, val_col, color_col):
    df = pd.read_csv(infile,
                     sep='\t',
                     usecols = [group_col, val_col, color_col],
                     dtype={group_col: str,
                            val_col: float,
                            color_col: str})
    return df





def plot(df, group_col, val_col, color_col, outfile):
    fig, ax = plt.subplots(1, 1, figsize=(5,4))
    lefts = np.arange(df[group_col].unique().shape[0])
    grouplabels = df[group_col].unique()
    bt = [0] * len(lefts)
    wd = 0.8


    for group in
    for hap, color in zip(('norm', 'del_het', 'del_hom'),
                          ('#b3b3b3', '#abc9ea', '#597dbf')):
        print(hap, color)
        ax.bar(left=lefts, height=df.loc[df['type']==hap, 'perc'].values, width=wd, bottom=bt, color=color, label=hap)
        bt = bt + df.loc[df['type']==hap, 'perc'].values


    for ngroup, (group, num) in enumerate(df.loc[df['type'] == 'norm', ['group', 'sum']].values):
        left = lefts[ngroup]
        # Create annotation
        ax.text(
            left,
            1.02,
            f'n={num}',
            ha='center', va='bottom')


    ax.set_xlim(-0.5, len(lefts)-0.5)
    ax.set_xticks(lefts)
    ax.set_xticklabels(indlabels, rotation='vertical')
    ax.set_ylim(-0, 1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #plt.legend(bbox_to_anchor=(1.0,0.5))
    plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
    plt.xticks(rotation=40)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    plt.legend(bbox_to_anchor=(1.0,0.5))
    plt.subplots_adjust(left=0.1, bottom=0.3, right=0.8, top=0.9)
    plt.savefig(outfile, dpi=300, transparent=True)




@click.command()
@click.option('--infile')
@click.option('--group-col')
@click.option('--val-col')
@click.option('--color-col')
@click.option('--outfile')
def main(infile, group_col, val_col, color_col, outfile):
    df = load_data(infile)
