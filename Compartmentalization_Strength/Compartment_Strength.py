import cooler
import cooltools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from cooltools.api.saddle import saddle_strength
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


# This script calculates the strength of compartmentalization in each tissue using E1 values for 10 chromosomes.
# We then present the computational results using saddle plots and scatter plots.


def saddleplot(track, saddledata, ax, n_bins, vrange=None, qrange=(0.0, 1.0),
               vmin=0.5, vmax=2, xlabel=None, title=None):

    digitized_track, binedges = cooltools.digitize(track, n_bins, vrange=vrange, qrange=qrange)
    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    lo, hi = qrange
    binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Heatmap
    norm = LogNorm(vmin=vmin, vmax=vmax)
    img = ax.pcolormesh(X, Y, C, norm=norm, cmap="coolwarm")
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_xticks([lo, 0.2, 0.4, 0.6, 0.8, hi])
    ax.set_xticklabels(['B','0.2','0.4','0.6', '0.8', 'A'], fontsize=10)
    ax.yaxis.set_visible(False)

    divider = make_axes_locatable(ax)

    # Margins
    # left margin hist
    ax_left = divider.append_axes("left", size="20%", pad=0.1, sharey=ax)
    ax_left.barh(binedges, height=1/len(binedges), width=groupmean, align="edge", edgecolor="k", linewidth=1)
    ax_left.set_xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    ax_left.set_ylim(hi, lo)
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["bottom"].set_visible(False)
    ax_left.spines["right"].set_visible(False)
    ax_left.spines["left"].set_visible(False)
    ax_left.xaxis.set_visible(False)
    ax_left.yaxis.set_visible(False)
    # top margin hist
    ax_top = divider.append_axes("top", size="20%", pad=0.1, sharex=ax)
    ax_top.bar(binedges, width=1/len(binedges), height=groupmean, align="edge", edgecolor="k", linewidth=1)
    ax_top.set_xlim(lo, hi)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.set_ylabel('E1       ', fontsize=12, rotation=0)
    ax_top.set_yticks([])
    ax_top.xaxis.set_visible(False)
    ax_top.set_title(title, fontsize=14)

    # Colorbar
    ax_bar = divider.append_axes("right", size="5%", pad=0.1)
    ax_bar.set_ylabel('average observed/expected contact frequency')
    plt.colorbar(img, cax=ax_bar)


def saddle_plot_strength(clr, cvd, e1_data, ax, strength_data, label, title=None, 
                         n_groups=None, q_lo=None, q_hi=None, view_df=None):
    # Saddle
    interaction_sum, interaction_count =  cooltools.saddle(
        clr, cvd, e1_data, 'cis', n_bins=n_groups, qrange=(q_lo,q_hi), view_df=view_df)
    # plot
    saddleplot(e1_data, interaction_sum/interaction_count, ax, n_groups, qrange=(q_lo,q_hi),
               xlabel='Label {0}'.format(label), title=title)
    # Strength
    strength_data[label] = saddle_strength(interaction_sum, interaction_count)
    return strength_data


def strengthplot_scatter(data_path, ax, corner_extent, label_data, fig, title=None, vmax=None, vmin=None):
    
    data = pd.read_csv(data_path, header=0, index_col=0, sep='\t')
    strength_e = data.iloc[:, corner_extent+1]
    strength_e = dict(strength_e)

    label_set = set(label_data['label'].values)
    for label_i in label_set:
        if label_i not in strength_e.keys():
            strength_e[label_i] = None
    
    label_data['Compartment Strength'] = label_data['label'].apply(lambda x: strength_e[x])
    label_ndata = label_data[~(label_data['Compartment Strength']>0)].copy()

    im = ax.scatter(label_data['x'], label_data['y'], c=label_data['Compartment Strength'], 
                    s=4, marker='s', cmap='inferno', vmax=vmax, vmin=vmin)
    ax.scatter(label_ndata['x'], label_ndata['y'], c='#d5d0cf', s=3, marker='s')
    fig.colorbar(im, ax=ax, orientation='vertical', ticks=[1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
    ax.set_title(title)



sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
chrom_id_list = [2, 4, 5, 6, 8, 9, 11, 15, 17, 19]
chrom_idx = [int(x) for x in np.array(chrom_id_list) - 1]
chrom_name = ['chr{0}'.format(x) for x in chrom_id_list]
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
embryo_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Heatmap/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartmentalization_Strength/'

for sample_id in sampleid_list:

    embryo_cool = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.cool'.format(embryo_dir, sample_id)
    embryo_clr = cooler.Cooler(embryo_cool)
    cooler.balance_cooler(embryo_clr, cis_only=True, store=True)

    view_df = pd.DataFrame({'chrom': chrom_name, 
                            'start': 0, 
                            'end': embryo_clr.chromsizes.values[chrom_idx],
                            'name': chrom_name})

    embryo_cvd = cooltools.expected_cis(clr=embryo_clr, view_df=view_df)
    embryo_e1_path = '{0}E1_Value/{1}.e1.value'.format(save_dir, sample_id)
    embryo_e1_data = pd.read_csv(embryo_e1_path, header=0, index_col=None, sep='\t')
    embryo_e1_data = embryo_e1_data.loc[:, ['chrom', 'start', 'end', 'embryo_e1']]
    embryo_e1_data.columns = ['chrom', 'start', 'end', 'E1']

    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    # Saddle
    col_num = 6
    label_num = len(label_set)
    if label_num % col_num == 0:
        row_num = label_num // col_num
    else:
        row_num = label_num // col_num + 1
    Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
    Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
    N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each

    f, axs = plt.subplots(figsize=(4.5 * col_num, 4 * row_num), 
                          nrows=row_num, ncols=col_num,
                          gridspec_kw={"hspace": 0.3, "wspace": 0.3})
    
    strength_data = pd.DataFrame()
    for i, label_i in enumerate(label_set):
        col_idx = i % col_num
        row_idx = i // col_num
        if row_num == 1:
            ax_embryo = axs[col_idx]
        else:
            ax_embryo = axs[row_idx, col_idx]

        tissue_cool = '{0}Data_Cool/{1}_C{2}.cool'.format(save_dir, sample_id, label_i)
        tissue_clr = cooler.Cooler(tissue_cool)
        cooler.balance_cooler(tissue_clr, cis_only=True, store=True)
        strength_data = saddle_plot_strength(tissue_clr, embryo_cvd, embryo_e1_data, 
                                             ax_embryo, strength_data, label_i, n_groups=N_GROUPS, 
                                             q_lo=Q_LO, q_hi=Q_HI, view_df=view_df)
        
    strength_path = '{0}Strength_{1}.tsv'.format(save_dir, sample_id)
    strength_data = strength_data.T
    strength_data.to_csv(strength_path, header=True, index=True, sep='\t')
    plot_path = '{0}Saddle.Strength_{1}.pdf'.format(save_dir, sample_id)
    plt.savefig(plot_path)
    plt.close()

    # Scatter
    corner_extent = 5
    if sample_id == 'MS0612-5':
        vmax = 2.0
        vmin = 1.0
        figsize = (4, 6)
    elif sample_id == 'MS0612-3.odd70':
        vmax = 1.8
        vmin = 1.0
        figsize = (4, 4.5)
    fig = plt.figure(figsize=figsize)
    ax_embryo = fig.add_subplot(1, 1, 1)
    strengthplot_scatter(strength_path, ax_embryo, corner_extent, label_data, 
                         fig, title='Compartmentalization Strength', vmax=vmax, vmin=vmin)

    plot_path = '{0}Scatter.Strength_{1}.pdf'.format(save_dir, sample_id)
    plt.savefig(plot_path)
    plt.close()
