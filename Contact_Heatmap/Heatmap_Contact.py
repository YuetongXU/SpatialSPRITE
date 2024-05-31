import cooler
import cooltools
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
mpl.rcParams['pdf.fonttype'] = 42


# The script uses cooltools software to plot the contact heatmap of chromosomes。


def chrom_info(resolution):
    chrom_data = pd.read_csv('/home/xuyuetong/Reference/mm10_reference/mm10_chrom_sizes.txt', header=None, index_col=None)
    chrom_data.columns = ['chrID', 'chrLength']
    chrom_data['chrID'] = chrom_data['chrID'].apply(lambda x: x.strip('chr'))
    chrom_data['binNum'] = chrom_data['chrLength'].apply(lambda x: int(x / resolution) + 1)

    binNum_acc = 0
    binNum_acc_list = []
    for i in chrom_data['binNum']:
        binNum_acc += i
        binNum_acc_list.append(binNum_acc)
    chrom_data['binNum_acc'] = binNum_acc_list
    return chrom_data


fig_nrow = 4
fig_ncol = 5
resolution = 1_000_000
norm = LogNorm(vmin=0.0035, vmax=0.1)
chrom_size = chrom_info(resolution)
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
chr_id_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
               'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
               'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Heatmap/'

for sample_id in sampleid_list:

    plot_path = '{0}Heatmap_Contact_{1}.pdf'.format(data_dir, sample_id)
    cool_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.cool'.format(data_dir, sample_id)
    clr = cooler.Cooler(cool_path)
    cooler.balance_cooler(clr, store=True)

    f, ax = plt.subplots(figsize=(6*fig_ncol, 5*fig_nrow), ncols=fig_ncol, nrows=fig_nrow)

    for i, chr_id in enumerate(chr_id_list):

        # Smooth and Imputation
        cg = adaptive_coarsegrain(clr.matrix(balance=True).fetch(chr_id),
                                  clr.matrix(balance=False).fetch(chr_id),
                                  cutoff=3, max_levels=8)
        cgi = interp_nan(cg, pad_zeros=True)

        # Calculate Coverage
        clr_cover, _ = cooltools.coverage(clr)
        chrom_df = chrom_size[chrom_size['chrID'] == chr_id.split('chr')[-1]]
        chr_accbinstart = chrom_df['binNum_acc'].values[0] - chrom_df['binNum'].values[0]
        chr_accbinend = chrom_df['binNum_acc'].values[0]
        chr_cover = clr_cover[chr_accbinstart: chr_accbinend]

        # Plot
        col_idx, row_idx = i % fig_ncol, i//fig_ncol
        ax_chr = ax[row_idx, col_idx]

        im = ax_chr.matshow(cgi, cmap=fruitpunch, norm=norm)
        ax_chr.set_title('Chrom {0}'.format(chr_id.split('chr')[-1]), fontsize=16)
        ax_chr.xaxis.set_visible(False)

        divider = make_axes_locatable(ax_chr)

        ax_colorbar = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(im, cax=ax_colorbar)

        ax_cover = divider.append_axes("bottom", size="20%", pad=0.2, sharex=ax_chr)
        ax_cover.plot(chr_cover, label='coverage')

    plt.subplots_adjust(hspace=0.2, wspace=0)
    plt.savefig(plot_path)

