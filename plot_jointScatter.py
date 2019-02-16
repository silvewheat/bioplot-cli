# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 15:28:26 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt



def load_infile(infile, x_col, y_col):
    """
    infile contain several chromosomes
    """
    df = pd.read_csv(infile,
                     sep='\t',
                     usecols = [x_col, y_col],
                     dtype={
                            x_col: float,
                            y_col: float})
    df.dropna(inplace=True)
    return df



def plot(df, x_col, y_col, cutoff_x, cutoff_y, markersize, outfile):
    fig, axScatter = plt.subplots(figsize=(5.5, 5.5))

    # scatter plot
    axScatter.scatter(df[x_col], df[y_col],
                      s=markersize,
                      color='k')


    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(axScatter)
    axHistx = divider.append_axes("top", 0.8, pad=0.1, sharex=axScatter)
    axHisty = divider.append_axes("right", 0.8, pad=0.1, sharey=axScatter)

    # make some labels invisible
    _ = plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible=False)

    # make some spines invisible
    axHistx.spines['right'].set_visible(False)
    axHistx.spines['top'].set_visible(False)
    axHisty.spines['right'].set_visible(False)
    axHisty.spines['top'].set_visible(False)
    axScatter.spines['right'].set_visible(False)
    axScatter.spines['top'].set_visible(False)

    # hist plot
    _ = axHistx.hist(df[x_col], bins=100, color='#CCC0DA')
    _ = axHisty.hist(df[y_col], bins=100, orientation='horizontal', color='#D8E4BC')


    # fix xlim ylim
    scatter_xlim = axScatter.get_xlim()
    scatter_ylim = axScatter.get_ylim()
    axScatter.set_ylim(scatter_ylim)
    axScatter.set_xlim(scatter_xlim)
    histx_xlim = axHistx.get_xlim()
    histx_ylim = axHistx.get_ylim()
    axHistx.set_xlim(histx_xlim)
    axHistx.set_ylim(histx_ylim)
    histy_xlim = axHisty.get_xlim()
    histy_ylim = axHisty.get_ylim()
    axHisty.set_xlim(histy_xlim)
    axHisty.set_ylim(histy_ylim)

    # cutoff
    cutcolor = '#808080'
    cutstyle = 'dashed'
    axScatter.vlines(cutoff_x, scatter_ylim[0], scatter_ylim[1], linestyle=cutstyle, color=cutcolor)
    axHistx.vlines(cutoff_x, histx_ylim[0], histx_ylim[1], linestyle=cutstyle, color=cutcolor)
    axScatter.hlines(cutoff_y, scatter_xlim[0], scatter_xlim[1], linestyle=cutstyle, color=cutcolor)
    axHisty.hlines(cutoff_y, histy_xlim[0], histy_xlim[1], linestyle=cutstyle, color=cutcolor)


    # 高亮过cutoff的点
    tdf = df.loc[(df[x_col]>=cutoff_x) & (df[y_col]>=cutoff_y), :]
    axScatter.scatter(tdf[x_col], tdf[y_col],
                      s=markersize,
                      color='#E74C3C')

    plt.savefig(outfile, dpi=300, transparent=True)





@click.command()
@click.option('--infile', help='tsv文件,包含header')
@click.option('--x-col', help='x轴值列名')
@click.option('--y-col', help='y轴值列名')
@click.option('--cutoff-x', help='highlight cutoff in x axis', type=float)
@click.option('--cutoff-y', help='hightlight cutoff in y axis', type=float)
@click.option('--markersize', default=3, help='散点大小, default is 3', type=float)
@click.option('--outfile', help='输出文件,根据拓展名判断输出格式')
def main(infile, x_col, y_col, cutoff_x, cutoff_y, markersize, outfile):
    df = load_infile(infile, x_col, y_col)
    plot(df, x_col, y_col, cutoff_x, cutoff_y, markersize, outfile)


if __name__ == '__main__':
    main()
