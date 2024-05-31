import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['pdf.fonttype'] = 42


# This script plots the contact frequency of long-distance (>40 MB) between compartments of the same type. 
# For easy comparison, we normalize (Max-Min) the contact frequency by tissue. 
# The color of the lines in the plot corresponds to the contact frequency.
# To avoid too many lines, we only show contacts with normalized (Max-Min) frequency above the threshold (>0.2).


def compartment_cmap(data, cate, vmin=None):
    contactnum = np.array(data.loc[data['cate'] == cate, 'contact_num'].values)
    cnum_max, cnum_min = max(contactnum), min(contactnum)
    if vmin:
        cnum_min = vmin
    contactnum_norm = (contactnum - cnum_min)/(cnum_max - cnum_min)
    if cate == 'A-A':
        com_cmap = mpl.colormaps.get_cmap('Reds')
    else:
        com_cmap = mpl.colormaps.get_cmap('Blues')
    com_cmap_list = com_cmap(contactnum_norm)
    data.loc[data['cate'] == cate, 'contact_norm'] = contactnum_norm
    return com_cmap_list, contactnum


resolution = 1_000_000
sample_id = 'MS0612-5'
tissue_id_list = [1, 17, 33, 3, 12, 16, 6, 18, 22, 30, 34]
chrom_list = ['2', '4', '5', '6', '8', '9', '11', '15', '17', '19']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Long_Compartment_Contact/'

# Prepare X label
chrlength_path = '{0}mm10_chrom_sizes.txt'.format(data_dir)
chrlength_data = pd.read_csv(chrlength_path, header=None, index_col=None, sep=',')
chrlength_data.columns = ['chrom', 'length']
chrlength_data['bin'] = chrlength_data['length'].apply(lambda x: int(x/resolution))
chrlength_data['chrom'] = chrlength_data['chrom'].apply(lambda x: x.strip('chr'))
chrlength_data = chrlength_data[chrlength_data['chrom'].isin(chrom_list)]
bin_num = chrlength_data['bin'].values
bin_accnum = [0]
chrom_interval = 20
for i, bin_i in enumerate(bin_num):
    bin_i = bin_accnum[i] + bin_i + chrom_interval
    bin_accnum.append(bin_i)
chrlength_data['x_startbin'] = bin_accnum[:-1]
chrlength_data['x_endbin'] = chrlength_data['x_startbin'] + chrlength_data['bin']
chrlength_data['x_midbin'] = chrlength_data['x_startbin'] + chrlength_data['x_endbin']
chrlength_data['x_midbin'] = chrlength_data['x_midbin'].apply(lambda x: int(x/2))

# Prepare Y label
tissue_yrange = 2
tissue_yinterval = 0.3
tissue_yloc = reversed([x * (tissue_yrange + tissue_yinterval) for x in range(len(tissue_id_list))])
y_data = pd.DataFrame({'tissue_yloc': tissue_yloc}, index=tissue_id_list)

# Inter-Contact
contact_path = '{0}Compartment_ContactNum_Inter.LongDistance.{1}.tsv'.format(data_dir, sample_id)
contact_data = pd.read_csv(contact_path, header=0, index_col=None, sep='\t')

fig = plt.figure(figsize=(1 * len(chrom_list), 0.7 * len(tissue_id_list) + 2))
ax = fig.add_subplot(111, aspect='auto')
x_ticks = chrlength_data['x_midbin'].values
x_label = chrom_list
y_ticks = y_data['tissue_yloc'].values
y_label = y_data['tissue_yloc'].index.values
ax.set_xticks(x_ticks, x_label, fontsize=1.2 * len(chrom_list), rotation=0)
ax.set_yticks(y_ticks, y_label, fontsize=1.2 * len(chrom_list), rotation=0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Chrom ID', fontsize=1.2 * len(chrom_list))
ax.set_ylabel('Label ID', fontsize=1.2 * len(chrom_list))
ax.set_title('Long Distance Inter-Compartment Contact', fontsize=1.2 * len(chrom_list))

for tissue_id in tissue_id_list:
    
    ab_linewidth = 0.10
    tissue_y_loc = y_data.loc[tissue_id, 'tissue_yloc']
    tissue_data = contact_data[contact_data['tissue'] == tissue_id].copy()
    # Tissue
    aa_cmap_list, aa_contactnum = compartment_cmap(tissue_data, 'A-A', vmin=0)
    bb_cmap_list, bb_contactnum = compartment_cmap(tissue_data, 'B-B', vmin=0)
    
    for chrom_id in chrom_list:

        chrom_data = tissue_data[tissue_data['chrom'] == int(chrom_id)].copy()
        x_startbin = chrlength_data.loc[chrlength_data['chrom'] == chrom_id, 'x_startbin'].values[0]
        chrom_data['start1_bin'] = chrom_data['start1_bin'] + x_startbin
        chrom_data['mid1_bin'] = chrom_data['mid1_bin'] + x_startbin
        chrom_data['end1_bin'] = chrom_data['end1_bin'] + x_startbin
        chrom_data['start2_bin'] = chrom_data['start2_bin'] + x_startbin
        chrom_data['mid2_bin'] = chrom_data['mid2_bin'] + x_startbin
        chrom_data['end2_bin'] = chrom_data['end2_bin'] + x_startbin

        # inter-compartment contact
        for row_i in range(chrom_data.shape[0]):

            data_i = chrom_data.iloc[row_i, :]
            start1_bin = data_i['start1_bin']
            mid1_bin = data_i['mid1_bin']
            end1_bin = data_i['end1_bin']
            start2_bin = data_i['start2_bin']
            mid2_bin = data_i['mid2_bin']
            end2_bin = data_i['end2_bin']
            bin_distance = data_i['bin_distance']
            compart_cate = data_i['cate']
            contact_num = data_i['contact_num']
            half_height = data_i['contact_norm']

            if half_height > ab_linewidth * 2:
                center = (mid1_bin + mid2_bin)/2
                if compart_cate == "A-A":
                    theta1 = 0
                    theta2 = 180
                    color_index = np.where(aa_contactnum == contact_num)[0][0]
                    color = aa_cmap_list[color_index]
                else:
                    theta1 = -180
                    theta2 = 0
                    color_index = np.where(bb_contactnum == contact_num)[0][0]
                    color = bb_cmap_list[color_index]
                arc=Arc(xy=(center, tissue_y_loc), width=bin_distance, height=2 * half_height, angle=0, 
                        theta1=theta1, theta2=theta2, color=color, lw=1.5, fill=False)
                ax.add_patch(arc)

        # A/B compartment
        compart_list = []

        for row_i in range(chrom_data.shape[0]):

            data_i = chrom_data.iloc[row_i, :]
            start1_bin = data_i['start1_bin']
            end1_bin = data_i['end1_bin']
            start2_bin = data_i['start2_bin']
            end2_bin = data_i['end2_bin']
            compart_cate = data_i['cate']
            compart1 = '{0}-{1}'.format(start1_bin, end1_bin)
            compart2 = '{0}-{1}'.format(start2_bin, end2_bin)

            if compart_cate == 'A-A':
                line_color = '#d33f3f'
            else:
                line_color = '#3b82db'
            if compart1 not in compart_list:
                compart_list.append(compart1)
                ax.fill_between(
                    [start1_bin, end1_bin], [tissue_y_loc-ab_linewidth, tissue_y_loc-ab_linewidth], 
                    [tissue_y_loc+ab_linewidth, tissue_y_loc+ab_linewidth], color=line_color)
            if compart2 not in compart_list:
                compart_list.append(compart2)
                ax.fill_between(
                    [start2_bin, end2_bin], [tissue_y_loc-ab_linewidth, tissue_y_loc-ab_linewidth], 
                    [tissue_y_loc+ab_linewidth, tissue_y_loc+ab_linewidth], color=line_color)

divider = make_axes_locatable(ax)
norm = mpl.colors.Normalize(vmin=0, vmax=1)
aa_cax = divider.append_axes("bottom", size="1%", pad=0.6)
aa_im = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colormaps.get_cmap('Reds'))
aa_cbar = fig.colorbar(aa_im, cax=aa_cax, orientation='horizontal', 
                       ticks=np.linspace(0, 1, 11))
bb_cax = divider.append_axes("bottom", size="1%", pad=0.25)
bb_im = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colormaps.get_cmap('Blues'))
bb_cbar = fig.colorbar(bb_im, cax=bb_cax, orientation='horizontal', 
                       ticks=np.linspace(0, 1, 11))

plot_path = '{0}Compartment_ContactNum_Inter.LongDistance.{1}.pdf'.format(data_dir, sample_id)
plt.savefig(plot_path)
