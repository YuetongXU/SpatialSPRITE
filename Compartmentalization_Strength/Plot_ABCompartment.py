import cooler
import cooltools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
mpl.rcParams['pdf.fonttype'] = 42


# The script uses cooltools to calculate the E1 value of chromosomes and infer AB compartments.
# The script then calculated the pearson correlation coefficient between tissue E1 and embryo E1


def delete_matrix(matrix, delete_num):
    delete_idx = [range(delete_num)]
    delete_idx_ = [-x for x in range(1, delete_num + 1)]
    matrix = np.delete(matrix, delete_idx, 0)
    matrix = np.delete(matrix, delete_idx, 1)
    matrix = np.delete(matrix, delete_idx_, 0)
    matrix = np.delete(matrix, delete_idx_, 1)
    return matrix


def imputation(clr, chr_id: str, delete_bin: int = 3) -> np.ndarray:
    # Smooth and Imputation
    cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(chr_id), 
                              clr.matrix(balance=False).fetch(chr_id), 
                              cutoff=3, max_levels=8)
    cgi = interp_nan(cg, pad_zeros=True)
    cgi = delete_matrix(cgi, delete_bin)
    return cgi


def calculate_e1(matrix: np.ndarray, gc_cov_chr: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors  = cooltools.api.eigdecomp.cis_eig(
        matrix, phasing_track=gc_cov_chr, ignore_diags=0, clip_percentile=99, n_eigs=30, sort_metric='pearsonr')
    e1_array = eigenvectors[0, :]
    return e1_array


def plot_ab(matrix, eigenvectors: np.ndarray, chrom_id: int, xmin: int, xmax: int, norm, ax):
    
    ax.matshow(matrix, norm=norm, cmap=fruitpunch)

    ymax = xmin
    ymin = xmax
    ax.set_ylabel('Contact Frequency', fontsize=16)
    ax.set_yticks(range(ymax + 10, ymin, 10))
    ax.xaxis.set_visible(False)

    divider = make_axes_locatable(ax)
    ax1 = divider.append_axes("top", size="15%", pad=0.20, sharex=ax)
    ax1.plot([xmin, xmax], [0, 0], 'k', lw=0.25)
    ax1.fill_between(range(matrix.shape[0]), eigenvectors, 0, where=(eigenvectors > 0), facecolor='red', alpha=0.5)
    ax1.fill_between(range(matrix.shape[0]), eigenvectors, 0, where=(eigenvectors <= 0), facecolor='blue', alpha=0.5)
    ax1.plot(eigenvectors, label='E1', color='black', alpha=1, lw=0.2)
    ax1.set_ylabel('E1', fontsize=16)
    ax1.set_xticks([])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_title('Chrom {0}'.format(chrom_id.split('chr')[-1]), fontsize=16)


delete_bin = 3
resolution = 1_000_000
norm = LogNorm(vmax=0.05, vmin=0.003)
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
embryo_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Heatmap/'
gc_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/AB_Compartment/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartmentalization_Strength/'
chrom_id_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']

gc_path = '{0}mm10_gc_cov_1MB.tsv'.format(gc_dir)
gc_cov = pd.read_csv(gc_path, header=0, index_col=None, sep='\t')

for sample_id in sampleid_list:

    # embryo E1
    embryo_cool = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.cool'.format(embryo_dir, sample_id)
    embryo_clr = cooler.Cooler(embryo_cool)
    cooler.balance_cooler(embryo_clr, cis_only=True, store=True)

    gc_e1_data = pd.DataFrame()
    for i, chr_id in enumerate(chrom_id_list):
        gc_cov_chr = gc_cov[gc_cov['chrom'] == chr_id].copy()
        gc_cov_chr = gc_cov_chr.iloc[delete_bin: -delete_bin, :]
        gc_value_chr = gc_cov_chr['GC'].values
        embryo_cgi = imputation(embryo_clr, chr_id)
        embryo_e1 = calculate_e1(embryo_cgi, gc_value_chr)
        gc_cov_chr['embryo_e1'] = embryo_e1
        gc_e1_data = pd.concat([gc_e1_data, gc_cov_chr], ignore_index=True)

    # Tissue
    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    for label_i in label_set:

        plot_path = '{0}AB_Compartment/{1}.C{2}.pdf'.format(save_dir, sample_id, label_i)

        tissue_cool = '{0}Data_Cool/{1}_C{2}.cool'.format(save_dir, sample_id, label_i)
        tissue_clr = cooler.Cooler(tissue_cool)
        cooler.balance_cooler(tissue_clr, cis_only=True, store=True)

        row_num = 4
        col_num = 5
        f, axs = plt.subplots(figsize=(7 * col_num, 8 * row_num), 
                              nrows=row_num, ncols=col_num, 
                              gridspec_kw={"hspace": 0.2, "wspace": 0})
        
        tissue_e1_data = pd.DataFrame()
        for i, chr_id in enumerate(chrom_id_list):
            gc_cov_chr = gc_cov[gc_cov['chrom'] == chr_id]
            gc_cov_chr = gc_cov_chr.iloc[delete_bin: -delete_bin, :]
            gc_value_chr = gc_cov_chr['GC'].values
            tissue_cgi = imputation(tissue_clr, chr_id, delete_bin)
            tissue_e1 = calculate_e1(tissue_cgi, gc_value_chr)
            gc_cov_chr['C{0}_e1'.format(label_i)] = tissue_e1
            tissue_e1_data = pd.concat([tissue_e1_data, gc_cov_chr], ignore_index=True)

            # Plot
            col_idx = i % col_num 
            row_idx = i // col_num
            ax_heatmap = axs[row_idx, col_idx]
            view_df = pd.DataFrame({'chrom': [chr_id], 
                                    'start': 0, 
                                    'end': tissue_clr.chromsizes[chr_id], 
                                    'name': [chr_id]})
            xmin = int(view_df['start'].values[0]/resolution)
            xmax = int(view_df['end'].values[0]/resolution) - 2 * delete_bin
            plot_ab(tissue_cgi, tissue_e1, chr_id, xmin, xmax, norm, ax_heatmap)
            
        plt.savefig(plot_path)

        gc_e1_data = gc_e1_data.merge(tissue_e1_data, how='inner', on=['chrom', 'start', 'end', 'GC'])

    save_e1_path = '{0}E1_Value/{1}.e1.value'.format(save_dir, sample_id, label_i)
    gc_e1_data.to_csv(save_e1_path, header=True, index=False, sep='\t')
    
    # pearson
    pearson_data = pd.DataFrame()
    for i, chr_id in enumerate(chrom_id_list):
        e1_chr = gc_e1_data[gc_e1_data['chrom'] == chr_id]
        embryo_e1_chr = e1_chr['embryo_e1']
        pearson_chrom_list = []
        for label_i in label_set:
            tissue_e1_chr = e1_chr['C{0}_e1'.format(label_i)]
            pearson_ij = round(stats.pearsonr(tissue_e1_chr, embryo_e1_chr)[0], 4)
            pearson_chrom_list.append(pearson_ij)
        pearson_chrom_data = pd.DataFrame(pearson_chrom_list, columns=[chr_id], 
                                          index=['C'+str(x) for x in label_set])
        if i == 0:
            pearson_data = pearson_chrom_data
        else:
            pearson_data = pearson_data.merge(pearson_chrom_data, left_index=True, right_index=True)
            
    pearson_data.loc['mean', :] = pearson_data.mean(axis=0)
    pearson_data = pearson_data.apply(lambda x: round(x, 4))
    pearson_path = '{0}E1_Value/{1}.e1.chrom.pearson'.format(save_dir, sample_id)
    pearson_data.to_csv(pearson_path, header=True, index=True, sep='\t')

    # Pearson Boxplot
    pearson_data = pd.read_csv(pearson_path, header=0, index_col=0, sep='\t')
    com_cmap = mpl.colormaps.get_cmap('seismic')
    com_cmap_list = list(com_cmap(pearson_data.loc['mean', :].values))
    pearson_data.drop(index=['C999', 'mean'], inplace=True)

    plot_data = pd.DataFrame()
    for i, col_i in enumerate(pearson_data.columns.values):
        data_i = pearson_data[col_i].to_frame()
        data_i.columns = ['pearson']
        data_i['chr_id'] = [col_i] * data_i.shape[0]
        data_i['cate'] = data_i.index.values
        if i == 0:
            plot_data = data_i
        else:
            plot_data = pd.concat([plot_data, data_i])

    fig = plt.figure(figsize=(13, 4))
    gs = GridSpec(1, 40, figure=fig)

    ax_box = fig.add_subplot(gs[0, 0:39])
    sns.boxplot(x='chr_id', y='pearson', data=plot_data, palette=com_cmap_list, 
                hue='chr_id', width=0.5, ax=ax_box)
    ax_box.set_xlabel('Chrom ID')
    ax_box.set_ylabel('Pearson')

    ax_cbar = fig.add_subplot(gs[0, 39])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    im = mpl.cm.ScalarMappable(norm=norm, cmap='seismic')
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')

    plot_path = '{0}E1_Value/Boxplot_{1}.e1.chrom.pearson.pdf'.format(save_dir, sample_id)
    plt.savefig(plot_path)

