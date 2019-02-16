# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 01:32:47 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def load_blastout(blastoutfile):
    """
    '7 std qcovs qcovhsp qcovus'
    """
    result = []
    with open(blastoutfile) as f:
        for line in f:
            if line[0] != '#':
                tline = line.strip().split()
                result.append([tline[0], tline[1], float(tline[2]), int(tline[3]), int(tline[6]), int(tline[7]), int(tline[8]), int(tline[9]), int(tline[-3])])
    return pd.DataFrame(result, columns=['qid', 'sid', 'ident', 'alnlen', 'qstart', 'qend', 'sstart', 'send', 'qcov'])


def load_querylen(querylen):
    lendict = {x.split()[0]: int(x.strip().split()[1]) for x in open(querylen).readlines()}
    return lendict


def plot(tmpdf, sid, lendict, outfile):
    fig, ax = plt.subplots(1, 1, figsize=(25,3))
    # 画subject
    sstart = tmpdf[['sstart', 'send']].min().min()
    send = tmpdf[['sstart', 'send']].max().max()
    slen = send - sstart + 1
    ax.fill_between([sstart, send], 7, 8)

    # 计算query的scale以及之间的gap长度
    query_totallen = 0
    n_query = 0
    for query in tmpdf.sort_values('sstart')['qid'].unique():
        query_totallen += lendict[query]
        n_query += 1
    qscale = (slen * 0.95) / query_totallen # query和subject对齐，并且总gap为subject长度的5%
    gap = (slen - query_totallen * qscale) / (n_query-1) if n_query > 1 else 0

    # 画query
    offset = sstart
    qoffset = {}
    for query in tmpdf.sort_values('sstart')['qid'].unique():
        qoffset[query] = offset
        xrange = [offset, offset+(lendict[query]*qscale)]
        ax.fill_between(xrange, 4, 5)
        plt.text(offset+(lendict[query]*qscale/2), 4.5, query, horizontalalignment='center')
        qcov = np.unique(tmpdf.loc[tmpdf['qid']==query, 'qcov'])[0]
        plt.text(offset+(lendict[query]*qscale/2), 4, f'{qcov}%', horizontalalignment='center')
        offset = offset + (lendict[query]*qscale) + gap


    # 画比对连线
    patches = []
    idents = []
    for qid, qstart, qend, sstart, send, ident, alnlen in tmpdf[['qid', 'qstart', 'qend', 'sstart', 'send', 'ident', 'alnlen']].values:
        offset = qoffset[qid]
        qalnlen = (qend-qstart+1)*qscale
        qstart_scaled = (qstart-1)*qscale
        polygon = Polygon([[offset+qstart_scaled, 5], [offset+qstart_scaled+qalnlen, 5], [send, 7], [sstart, 7]], closed=True)
        patches.append(polygon)
        idents.append(ident)
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.9)
    colors = idents
    p.set_array(np.array(colors))
    ax.add_collection(p)
    fig.colorbar(p, ax=ax)
    ax.set_xlabel(sid)
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(outfile, dpi=300)
    plt.close()


@click.command()
@click.option('--blastoutfile', help='blast结果')
@click.option('--querylen', help='两列，queryid\tquerylen')
@click.option('--minqcov', help='某个subject至少有一个query比对上去的总%coverage达到这个值才被画出来, 默认30', default=30, type=int)
@click.option('--minident', help='%相似度达到这个值的比对结果才会被画出来, 默认0', default=0, type=int)
@click.option('--minaln', help='长度达到这个值的比对结果才被画出来, 默认100', default=100, type=int)
@click.option('--outprefix', help='输出图片前缀')
def main(blastoutfile, querylen, minqcov, minident, minaln, outprefix):
    """
    画blast输出结果
    '7 std qcovs qcovhsp qcovus'
    """
    print(__doc__)
    print(blastoutfile)
    df = load_blastout(blastoutfile)
    lendict = load_querylen(querylen)
    for sid in df['sid'].unique():
        outfile = f'{outprefix}_{sid}.pdf'
        tmpdf = df.loc[(df['sid']==sid) &
                       (df['ident']>=minident) &
                       (df['alnlen']>=minaln), :]
        if (tmpdf.shape[0] > 0) and (tmpdf['qcov'].max()) >= minqcov:
            plot(tmpdf, sid, lendict, outfile)


if __name__ == '__main__':
    main()

