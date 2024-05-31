import cooler
import cooltools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
mpl.rcParams['pdf.fonttype'] = 42


# The script uses cooltools to calculate the E1 value of chromosomes and infer AB compartments.


def delete_matrix(matrix, delete_num):
    delete_idx = [range(delete_num)]
    delete_idx_ = [-x for x in range(1, delete_num + 1)]
    matrix = np.delete(matrix, delete_idx, 0)
    matrix = np.delete(matrix, delete_idx, 1)
    matrix = np.delete(matrix, delete_idx_, 0)
    matrix = np.delete(matrix, delete_idx_, 1)
    return matrix


def plot_ab(matrix, eigenvectors: np.ndarray, gc_cov: np.ndarray, 
            chrom_id: int, xmin: int, xmax: int, norm, ax):
    
    ax.matshow(matrix, norm=norm, cmap=fruitpunch)

    ymax = xmin
    ymin = xmax
    ax.set_ylabel('Contact Frequency', fontsize=16)
    ax.set_yticks(range(ymax + 10, ymin, 10))
    ax.xaxis.set_visible(False)

    divider = make_axes_locatable(ax)

    # E1
    ax1 = divider.append_axes("top", size="15%", pad=0.20, sharex=ax)
    ax1.plot([xmin, xmax], [0, 0], 'k', lw=0.25)
    ax1.fill_between(range(cgi.shape[0]), eigenvectors, 0, where=(eigenvectors > 0), facecolor='red', alpha=0.5)
    ax1.fill_between(range(cgi.shape[0]), eigenvectors, 0, where=(eigenvectors <= 0), facecolor='blue', alpha=0.5)
    ax1.plot(eigenvectors, label='E1', color='black', alpha=1, lw=0.2)
    ax1.set_ylabel('E1', fontsize=16)
    ax1.set_xticks([])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

    # GC
    ax2 = divider.append_axes("top", size="15%", pad=0.20, sharex=ax)
    ax2.plot(gc_cov, label='GC')
    ax2.set_ylabel('GC', fontsize=16)
    ax2.set_xticks([])
    ax2.set_title('Chrom {0}'.format(chrom_id.split('chr')[-1]), fontsize=16)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)


delete_bin = 3
resolution = 1_000_000
norm = LogNorm(vmax=0.05, vmin=0.003)
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Heatmap/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/AB_Compartment/'
chrom_id_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']

gc_path = '{0}mm10_gc_cov_1MB.tsv'.format(save_dir)
gc_cov = pd.read_csv(gc_path, header=0, index_col=None, sep='\t')

for sample_id in sampleid_list:

    plot_path = '{0}Lineplot_AB_Compartment_{1}.pdf'.format(save_dir, sample_id)
    cool_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.cool'.format(data_dir, sample_id)
    clr = cooler.Cooler(cool_path)
    cooler.balance_cooler(clr, store=True)

    row_num = 4
    col_num = 5
    f, axs = plt.subplots(figsize=(7 * col_num, 8 * row_num), 
                          nrows=row_num, ncols=col_num, 
                          gridspec_kw={"hspace": 0.2, "wspace": 0})

    for i, chr_id in enumerate(chrom_id_list):

        # Smooth and Imputation
        cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(chr_id),
                                  clr.matrix(balance=False).fetch(chr_id),
                                  cutoff=3, max_levels=8)
        cgi = interp_nan(cg, pad_zeros=True)
        cgi = delete_matrix(cgi, delete_bin)

        # GC
        gc_cov_chr = gc_cov.loc[gc_cov['chrom'] == chr_id, 'GC'].values
        gc_cov_chr = gc_cov_chr[delete_bin: -delete_bin]

        # E1
        eigenvalues, eigenvectors  = cooltools.api.eigdecomp.cis_eig(
            cgi, phasing_track=gc_cov_chr, ignore_diags=0, clip_percentile=99, n_eigs=30, sort_metric='pearsonr')
        impu_E1 = eigenvectors[0, :]
        
        # Plot
        col_idx = i % col_num 
        row_idx = i // col_num
        ax_heatmap = axs[row_idx, col_idx]
        view_df = pd.DataFrame({'chrom': [chr_id], 
                                'start': 0, 
                                'end': clr.chromsizes[chr_id], 
                                'name': [chr_id]})
        xmin = int(view_df['start'].values[0]/resolution)
        xmax = int(view_df['end'].values[0]/resolution) - 2 * delete_bin
        plot_ab(cgi, impu_E1, gc_cov_chr, chr_id, xmin, xmax, norm, ax_heatmap)
        
    plt.savefig(plot_path)