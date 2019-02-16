# -*- coding: utf-8 -*-
"""
Created on Thu May 10 16:12:12 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns



def cal_parameter(data, leftgaplen, rightgaplen, colorfile):
    """
    根据读入文件计算图片大小等
    """
    n_left = data['leftindex'].unique().shape[0]
    leftgap = (data['leftlabel'].unique().shape[0] - 1) * leftgaplen
    leftlen = n_left + leftgap # 每个record长度都是1
    n_right = data['rightindex'].unique().shape[0]
    rightgap = (data['rightlabel'].unique().shape[0] - 1) * rightgaplen
    rightlen = n_right + rightgap
    len_max = max(leftlen, rightlen)
    ax_ymax = len_max * 1.2 # 两遍空10%
    ax_xmax = ax_ymax * 0.8
    barwidth = ax_xmax / 10
    left_top = (ax_ymax + leftlen) / 2
    leftbottom = left_top - leftlen
    right_top = (ax_ymax + rightlen) / 2
    rightbottom = right_top - rightlen
    all_labels = np.r_[data['leftlabel'].unique(), data['rightlabel'].unique()]
    if colorfile is None:
        colordict = {}
        pal = "hls"
        cls = sns.color_palette(pal, len(all_labels))
        for i, l in enumerate(all_labels):
            colordict[l] = cls[i]
    return ax_ymax, ax_xmax, barwidth, leftbottom, rightbottom, colordict




def cal_bar_pos(data, leftgaplen, rightgaplen, ax_ymax, ax_xmax, barwidth, leftbottom, rightbottom):
    """
    输入读入文件的dataframe
    计算左右两个bar的绘图位置
    """
    # 左侧
    leftbars = []
    left_x1 = 0
    left_x2 = barwidth
    offset = leftbottom
    for label in data['leftlabel'].unique():
        high = data.loc[data['leftlabel']==label, 'leftindex'].unique().shape[0]
        left_y1 = offset
        left_y2 = left_y1 + high
        leftbars.append([label, left_x1, left_x2, left_y1, left_y2])
        offset = left_y2 + leftgaplen
    leftbars = pd.DataFrame(leftbars, columns=['label', 'x1', 'x2', 'y1', 'y2'])

    # 右侧
    rightbars = []
    right_x2 = ax_xmax
    right_x1 = right_x2 - barwidth
    offset = rightbottom
    for label in data['rightlabel'].unique():
        high = data.loc[data['rightlabel']==label, 'rightindex'].unique().shape[0]
        right_y1 = offset
        right_y2 = right_y1 + high
        rightbars.append([label, right_x1, right_x2, right_y1, right_y2])
        offset = right_y2 + rightgaplen
    rightbars = pd.DataFrame(rightbars, columns=['label', 'x1', 'x2', 'y1', 'y2'])
    return leftbars, rightbars


def convolve(left_y, right_y):
    ay = np.array(50 * [left_y] + 50 * [right_y])
    ay = np.convolve(ay, 0.05 * np.ones(20), mode='valid')
    ay = np.convolve(ay, 0.05 * np.ones(20), mode='valid')
    return ay


def cal_strips(data, leftbars, rightbars):
    left_label2y1 = {x:y for x,y in leftbars[['label', 'y1']].values}
    right_label2y1 = {x:y for x,y in rightbars[['label', 'y1']].values}
    left_indexbase = {}
    for label in data['leftlabel'].values:
        left_indexbase[label] = data.loc[data['leftlabel']==label, 'leftindex'].min()
    right_indexbase = {}
    for label in data['rightlabel'].values:
        right_indexbase[label] = data.loc[data['rightlabel']==label, 'rightindex'].min()
    data['left_y1'] =  data['leftlabel'].map(left_label2y1)
    data['left_y1'] = data['left_y1'] + (data['leftindex'] - data['leftlabel'].map(left_indexbase)) * 1
    data['right_y1'] = data['rightlabel'].map(right_label2y1)
    data['right_y1'] = data['right_y1'] + (data['rightindex'] - data['rightlabel'].map(right_indexbase)) * 1

    strips = []
    left_x = leftbars['x2'][0]
    right_x = rightbars['x1'][0]
    for left_y1, right_y1 in data[['left_y1', 'right_y1']].values:
        left_y2 = left_y1 + 1
        right_y2 = right_y1 + 1
        strip_y1 = convolve(left_y1, right_y1)
        strip_y2 = convolve(left_y2, right_y2)
        strips.append([strip_y1, strip_y2])
    return strips, left_x, right_x



def plot(ax_xmax, ax_ymax, leftbars, rightbars, strips, colordict, outfile):
    fig, ax = plt.subplots(1, 1, figsize=(8, 10))
    for label, x1, x2, y1, y2 in leftbars.values:
        ax.fill_between([x1, x2], y1,  y2, color=colordict[label])
        ax.text(x1-0.2, y1+0.5, label, {'ha': 'right', 'va': 'center'})
    for label, x1, x2, y1, y2 in rightbars.values:
        ax.fill_between([x1, x2], y1,  y2, color=colordict[label])
        ax.text(x2+0.5, y1+0.5, label, {'ha': 'right', 'va': 'center'})
    for y1, y2 in strips:
        ax.fill_between(
                        np.linspace(leftbars['x2'][0], rightbars['x1'][0], len(y1)),
                        y1, y2, alpha=0.6,
                        color='k')
    ax.set_xlim([-0.1, ax_xmax+0.1])
    ax.set_ylim([-0.1, ax_ymax+0.1])
    ax.axis('off')
    plt.savefig(outfile)




@click.command()
@click.option('--infile')
@click.option('--rightgaplen', default=0.5, help='block之间的距离, default is 0.5')
@click.option('--leftgaplen', default=0.5, help='block之间的距离, default is 0.5')
@click.option('--colorfile', default=None, help='左右两边label的颜色, default is None, tsv文件, 两列')
@click.option('--outfile', help='输出图片文件,如out.pdf')
def main(infile, rightgaplen, leftgaplen, colorfile, outfile):
    """
    \b
    infile:
    # index从0开始，连续整数, 索引值相同的会处于同一位置, 下面这个示例为source1的一个位置同时连接到target1和target2
    leftindex   leftlabel    rightindex    rightlabel
    0    source1    0    target1
    0    source1    1    target2
    1    source2    2    target2
    """
    data = pd.read_csv(infile, sep='\t', low_memory=False)
    ax_ymax, ax_xmax, barwidth, leftbottom, rightbottom, colordict = cal_parameter(data, rightgaplen, leftgaplen, colorfile)
    leftbars, rightbars = cal_bar_pos(data, rightgaplen, leftgaplen, ax_ymax, ax_xmax, barwidth, leftbottom, rightbottom)
    strips, left_x, right_x = cal_strips(data, leftbars, rightbars)
    plot(ax_xmax, ax_ymax, leftbars, rightbars, strips, colordict, outfile)

if __name__ == '__main__':
    main()


