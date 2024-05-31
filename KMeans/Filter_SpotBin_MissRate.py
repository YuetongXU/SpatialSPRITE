import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42


# This script first calculates the distribution of missing rate for bin and spot.
# The results show that some bins have no interaction with any spot, which may be caused by systematic errors. 
# After excluding the bin without interaction, we also removed the bin and spot with high missing rate (top 5%).


def chrom_info(resolution, threshold):
    chrom_size = pd.read_csv('/home/xuyuetong/Reference/mm10_reference/mm10_chrom_sizes.txt', header=None, index_col=None)
    chrom_size.columns = ['chrID', 'chrLength']
    chrom_size['chrID'] = chrom_size['chrID'].apply(lambda x: x.strip('chr'))
    chrom_size['binNum'] = chrom_size['chrLength'].apply(lambda x: int(x/resolution) + 1)
    chrom_size['binThreshold'] = chrom_size['binNum'].apply(lambda x: int(x*threshold))

    binNum_acc = 0
    binNum_acc_list = []
    for i in chrom_size['binNum']:
        binNum_acc += i
        binNum_acc_list.append(binNum_acc)
    chrom_size['binNum_acc'] = binNum_acc_list
    return chrom_size


def sprite_matrix(sprite_path: str, matrix_path: str, chrom_size: pd.DataFrame, resolution: int) -> None:

    spot_dict = {}
    chrom_length_list = chrom_size['binNum_acc'].values

    with open(sprite_path, 'r') as f:

        for line in f.readlines():
            lines = line.strip().split('\t')
            spot_id = '.'.join(lines[0].split('.')[4:-1])
            sites = lines[1:]

            eigenvector  = np.array([0]*chrom_length_list[-1])
            for site_a in sites:

                site_a_chr, site_a_loc = site_a.split('_')[0:2]

                if site_a_chr == 'X':
                    site_a_chr = 20
                elif site_a_chr == 'Y':
                    site_a_chr = 21
                else:
                    site_a_chr = int(site_a_chr)

                site_a_bin = int(int(site_a_loc)/resolution)

                if site_a_chr == 1:
                    site_a_index = site_a_bin
                else:
                    site_a_index = chrom_length_list[site_a_chr - 2] + site_a_bin

                eigenvector[site_a_index] = eigenvector[site_a_index] + 1

            if spot_id in spot_dict.keys():
                spot_dict[spot_id] = spot_dict[spot_id] + eigenvector
            else:
                spot_dict[spot_id] = eigenvector

    matrix_data = pd.DataFrame(spot_dict).T
    matrix_data.to_csv(matrix_path, header=True, index=True, sep='\t')


def plot_missrate(rowmiss_data, colmiss_data, save_path):

    f, ax = plt.subplots(figsize=(8, 3), ncols=2, nrows=1)

    ax_rowmiss = ax[0]
    sns.histplot(rowmiss_data, ax=ax_rowmiss, binwidth=0.02)
    ax_rowmiss.set_xlim(0, 1)
    ax_rowmiss.set_title('Spot MissRate')

    ax_colmiss = ax[1]
    sns.histplot(colmiss_data, ax=ax_colmiss, binwidth=0.02)
    ax_colmiss.set_xlim(0, 1)
    ax_colmiss.set_title('Bin MissRate')

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


threshold = 0.4
resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'

for sample_id in sampleid_list:

    fb_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.sprite'.format(data_dir, sample_id)
    fb_matrix_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.matrix'.format(save_dir, sample_id)
    chrom_data = chrom_info(resolution, threshold)
    sprite_matrix(fb_path, fb_matrix_path, chrom_data, resolution)

    plot_miss_path = '{0}Histplot_{1}_FB_0.4.missrate.pdf'.format(save_dir, sample_id)
    plot_filter5_path = '{0}Histplot_{1}_FB_0.4.missrate.filter5%.pdf'.format(save_dir, sample_id)
    filter_matrix_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.matrix'.format(save_dir, sample_id)

    data_matrix = pd.read_csv(fb_matrix_path, header=0, index_col=0, sep='\t')
    row_num, col_num = data_matrix.shape
    row_miss = (data_matrix == 0).sum(axis=1) / col_num
    col_miss = (data_matrix == 0).sum(axis=0) / row_num
    plot_missrate(row_miss, col_miss, plot_miss_path)

    # Filter Bin (All Nan)
    data_matrix = data_matrix.loc[:, col_miss < 1]
    row_miss = (data_matrix == 0).sum(axis=1) / col_num
    col_miss = (data_matrix == 0).sum(axis=0) / row_num
    # Filter Bin and Spot (5%)
    topbin_5p = col_miss.quantile(0.95)
    data_matrix = data_matrix.loc[:, col_miss < topbin_5p]
    topspot_5p = row_miss.quantile(0.95)
    data_matrix = data_matrix.loc[row_miss < topspot_5p, :]
    row_miss = (data_matrix == 0).sum(axis=1) / col_num
    col_miss = (data_matrix == 0).sum(axis=0) / row_num
    plot_missrate(row_miss, col_miss, plot_filter5_path)
    data_matrix.to_csv(filter_matrix_path, header=True, index=True, sep='\t')
    print(data_matrix.shape)