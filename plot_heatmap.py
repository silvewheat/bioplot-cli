import typer
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns



def main(matrixfile: str = typer.Argument(..., help="输入的矩阵文件，包含header，第一列为index，tab分割"),
         outfile: str = typer.Argument(..., help="输出图片的文件名"),
         width: float = typer.Option(10, help="输出图片宽度."),
         height: float = typer.Option(10, help="输出图片高度."),
         transpose: bool = typer.Option(False, help="是否对原始矩阵进行转置."),
         row_cluster: bool = typer.Option(False, help="是否对行进行聚类."),
         col_cluster: bool = typer.Option(False, help="是否对列进行聚类."),
         metric: str = typer.Option('euclidean', help="聚类时使用的距离度量,详见scipy.spatial.distance.pdist."),
         yticklabels: bool = typer.Option(True, help="是否显示Y轴刻度名."),
         xticklabels: bool = typer.Option(True, help="是否显示X轴刻度名."),
         annot: bool = typer.Option(False, help="是否在格子里显示数值.")
         ):
    # 要显示就显示全部的
    yticklabels = 1 if yticklabels else yticklabels
    xticklabels = 1 if xticklabels else xticklabels
    annot = None if not annot else annot # False不管用，得改成None
    
    df = pd.read_csv(matrixfile, sep='\t', index_col=0)
    
    if transpose:
        df = df.T
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#00008B", "#F0FFFF", 'red'])
    g = sns.clustermap(df, row_cluster=row_cluster, col_cluster=col_cluster, yticklabels=yticklabels, xticklabels=xticklabels, metric=metric,
                       cmap=cmap, figsize=(width, height), annot=annot)
    g.ax_heatmap.set_facecolor('#BFBFBF')
    plt.savefig(outfile)
    

if __name__ == '__main__':
    typer.run(main)