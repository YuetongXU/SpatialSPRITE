import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['pdf.fonttype'] = 42


# This script plots the frequency of contacts at different distances (3 MB, 10 MB, 100 MB) for each tissue.


resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Distance/'

for sample_id in sampleid_list:

    distance_path = '{0}Contact_Distance.{1}.tsv'.format(data_dir, sample_id)
    expected_data = pd.read_csv(distance_path, header=0, index_col=None, sep='\t')

    if sample_id == 'MS0612-5':
        figsize = (5, 5.5)
        vmax_list = [2.15e-2, 1.0e-2, 5e-3]
        vmin_list = [1.45e-2, 7e-3, 1e-3]
        point_size = 3
    elif sample_id == 'MS0612-3.odd70':
        figsize = (5, 5)
        vmax_list = [2.0e-2, 1.0e-2, 5e-3]
        vmin_list = [1.4e-2, 7e-3, 1e-3]
        point_size = 4
    distance_list = [3_000_000, 20_000_000, 100_000_000]
    title_list = ['3MB', '10MB', '100MB']

    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    plot_path = '{0}Heatmap_ContactDistance.{1}.pdf'.format(data_dir, sample_id)
    with PdfPages(plot_path) as pdf:
        fig, axs = plt.subplots(figsize=(figsize[0] * len(distance_list), figsize[1]), 
                                nrows=1, ncols=len(distance_list), 
                                gridspec_kw={"hspace": 0, "wspace": 0.4})
        for d, distance in enumerate(distance_list):
            expected_data1 = expected_data[expected_data['s_bp'] == distance]
            expected_data1 = expected_data1.loc[:, ['balanced.avg.smoothed', 'label']]
            expected_data1.set_index('label', inplace=True) 
            expected_dict = expected_data1.to_dict()['balanced.avg.smoothed']
            label_data['Contact_Frequency_{0}'.format(d)] = label_data['label'].apply(lambda x: expected_dict[x])
            ax_sp = axs[d]
            vmin, vmax = vmin_list[d], vmax_list[d]
            im = ax_sp.scatter(label_data['x'], label_data['y'], 
                               c=label_data['Contact_Frequency_{0}'.format(d)], 
                               s=point_size, marker='s', cmap='coolwarm', vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=axs[d], orientation='vertical', ticks=np.linspace(vmin, vmax, 6))
            axs[d].set_title(title_list[d] + ' Contact Frequency')
        pdf.savefig()