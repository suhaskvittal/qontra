"""
    author: Suhas Vittal

    PERSONAL USE ONLY.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as lines

def light_theme(context='paper'):
    sns.set_theme(context=context, style='whitegrid', palette='dark', font_scale=1.0,\
            rc={'lines.linewidth': 1.1, 'axes.linewidth': 1.3, 'axes.edgecolor': 'black',\
            'axes.facecolor': '#fafafa', 'axes.spines.left': True, 'axes.spines.right': True,\
            'axes.spines.top': True, 'axes.spines.bottom': True})

def dark_theme(context='paper'):
    sns.set_theme(context=context, style='darkgrid', palette='pastel', font_scale=1.0,\
            rc={'lines.linewidth': 1.1, 'axes.linewidth': 1.3, 'axes.edgecolor': 'white',\
            'axes.facecolor': '#121212', 'axes.spines.left': True, 'axes.spines.right': True,\
            'axes.spines.top': True, 'axes.spines.bottom': True, 'grid.color': '#fafafa',\
            'text.color': '.15'})

light_theme()

def mkhist(data, x_name, y_name, category_key=None, palette=None, use_log=False, show_shape=False, binwidth=1, multiple='stack'):
    return sns.histplot(data=data, x=x_name, stat=y_name, hue=category_key,\
                    palette=palette, log_scale=use_log, kde=show_shape,
                    binwidth=binwidth, multiple=multiple)

def mkline(data, x_name, y_name, category_key=None, style_key=None, palette=None, alpha=1.0,\
            ci='sd', estimator=np.mean, marker='o', ticks=None):
    fig = sns.lineplot(data=data, x=x_name, y=y_name, hue=category_key, style=style_key,\
                    palette=palette, alpha=alpha, ci=ci,\
                    sort=True, estimator=estimator, markers=True, marker=marker)
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.3)
    if ticks is not None:
        fig.set_xlim(min(ticks), max(ticks))
        fig.set_xticks(ticks)
    return fig

def mkline2ax(data, x_name, y1_name, y2_name, category_key=None, style_key=None, palette1='pastel',\
        palette2='dark', alpha=1.0, ci1='sd', ci2='sd', est1=np.mean, est2=np.mean):
    fig = mkline(data, x_name, y1_name, category_key=category_key, style_key=style_key,\
                palette=palette1, alpha=alpha, ci=ci1, estimator=est1)
    ax = sns.lineplot(data=data, x=x_name, y=y2_name, hue=category_key, style=style_key,\
                    palette=palette2, alpha=alpha, ci=ci2,\
                    sort=True, estimator=est2, markers=True,\
                    marker='x', ax=fig.axes.twinx(), legend=False)
    ax.grid(False)
    return fig

def mkscatter(data, x_name, y_name, category_key=None, style_key=None, palette=None,\
                alpha=1.0, ci='sd', estimator=np.mean):
    fig = sns.scatterplot(data=data, x=x_name, y=y_name, hue=category_key, style=style_key,\
                        palette=palette, alpha=alpha, ci=ci, estimator=estimator)
    fig.patch.set_edgecolor('black')
    fig.patch.set_linewidth(1.3)
    return fig

def mkheatmap(data, cbar=True, linewidths=1.0, linecolor='black', cmap=None):
    fig = sns.heatmap(data=data, linewidths=linewidths, linecolor=linecolor, cmap=cmap)
    return fig

def add_baseline(fig, xmin, xmax, y):
    fig.plot([xmin, xmax], [y, y], linewidth=1.5, color='#545454', alpha=0.6, linestyle='dashed')

def add_title(fig, title):
    fig.set_title(title)

