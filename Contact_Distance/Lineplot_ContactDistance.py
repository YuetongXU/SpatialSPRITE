import cooler
import cooltools
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.rcParams['pdf.fonttype'] = 42


# This script calculates the contact frequency for different contact distances at the genome level for each tissue.


resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Distance/'

for sample_id in sampleid_list:
    
    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    expected_data = pd.DataFrame()
    for label_i in label_set:
        cool_path = '{0}Data_Cool/{1}_C{2}.onechrom.cool'.format(data_dir, sample_id, label_i)
        clr = cooler.Cooler(cool_path)
        cooler.balance_cooler(clr, cis_only=True, store=True)

        view_df = pd.DataFrame({'chrom': clr.chromnames,
                                'start': 0,
                                'end': clr.chromsizes.values,
                                'name': clr.chromnames})

        cvd = cooltools.expected_cis(clr=clr, view_df=view_df, smooth=True, aggregate_smoothed=True, nproc=10)
        cvd = cvd[(cvd['dist'] >= 2) & (cvd['count.sum'] > 0)]
        cvd['s_bp'] = cvd['dist'] * resolution
        cvd['label'] = [label_i] * cvd.shape[0]
        expected_data = pd.concat([expected_data, cvd])
    expected_data = expected_data.dropna(subset=['balanced.avg.smoothed'])
    expected_data = expected_data.loc[:, ['s_bp', 'balanced.avg.smoothed', 'label']]
    distance_path = '{0}Contact_Distance.{1}.tsv'.format(data_dir, sample_id)
    expected_data.to_csv(distance_path, header=True, index=False, sep='\t')

    # Plot
    if sample_id == 'MS0612-5':
        vmin = 1.45e-2
        vmax = 2.15e-2
    elif sample_id == 'MS0612-3.odd70':
        vmin = 1.4e-2
        vmax = 2.0e-2
    distance = 3_000_000
    plt.rcParams.update({"font.size": 5})
    expected_data = pd.read_csv(distance_path, header=0, index_col=None, sep='\t')
    color_data = expected_data[expected_data['s_bp'] == distance]
    frequency_value = color_data['balanced.avg.smoothed'].values
    frequency_value = (frequency_value - vmin)/(vmax - vmin)
    com_cmap = mpl.colormaps.get_cmap('coolwarm')
    com_cmap_list = list(com_cmap(frequency_value))

    fig = plt.figure(figsize=(4, 4))
    gs = GridSpec(1, 20, figure=fig)

    ax_line = fig.add_subplot(gs[0, 0:19])
    sns.lineplot(expected_data, x='s_bp', y='balanced.avg.smoothed', hue='label', 
                 palette=com_cmap_list, legend=False, lw=3, ax=ax_line)
    ax_line.set_xscale('log')
    ax_line.set_yscale('log')
    ax_line.set_xlabel('Distance (bp)')
    ax_line.set_ylabel('Contact Frequency')

    ax_cbar = fig.add_subplot(gs[0, 19])
    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    im = mpl.cm.ScalarMappable(norm=norm, cmap='coolwarm')
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    
    plot_path = '{0}Lineplot_ContactDistance.{1}.pdf'.format(data_dir, sample_id)
    plt.savefig(plot_path)