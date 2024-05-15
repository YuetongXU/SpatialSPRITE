import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42


# The script calculated the proportion of the non-diffusive interaction frequency in the spot to the total interaction frequency ((sample-back)/sample). 
# Here, we used the mean interaction frequency within the spot in the background region as the diffusion interaction frequency.


sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'

plot_path = '{0}Histogram_ValidInformRatio.pdf'.format(data_dir)
fig, axs = plt.subplots(1, 2, figsize=(5 * len(sampleid_list), 4))

for i, sampleid in enumerate(sampleid_list):

    legnum_path = '{0}clusters_{1}.spot.legnum.tsv'.format(data_dir, sampleid)
    sample_spotid_path = '{0}clusters_{1}.sample.spotid'.format(data_dir, sampleid)
    legnum_data = pd.read_csv(legnum_path, header=0, index_col=None, sep='\t')
    sample_spotid_data = pd.read_csv(sample_spotid_path, header=None, index_col=None, sep='\t')
    sample_spotid = sample_spotid_data.values.flatten()

    sample_legnum = legnum_data[legnum_data['spot_id'].isin(sample_spotid)]
    back_legnum = legnum_data[~legnum_data['spot_id'].isin(sample_spotid)]
    back_legmean = back_legnum['legnum'].mean()
    sample_legnum['ratio'] = (sample_legnum['legnum'] - back_legmean) / sample_legnum['legnum']

    hist = sns.histplot(data=sample_legnum, x='ratio', stat='probability', bins=50, ax=axs[i])
    axs[i].set_title('Sample {0}'.format(i + 1), fontsize=15)
    axs[i].set_xlabel('Valid Information Ratio', fontsize=15)
    axs[i].set_ylabel('Frequency', fontsize=14)

plt.tight_layout()
plt.savefig(plot_path)
