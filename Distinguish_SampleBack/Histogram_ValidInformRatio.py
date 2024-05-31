import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42


# The script calculated the proportion of the non-diffusive interaction frequency in the spot to the total interaction frequency ((sample-back)/sample). 
# Here, we used the mean interaction frequency within the spot in the background region as the diffusion interaction frequency.


sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'

plot_path = '{0}Histogram_ValidInformRatio.pdf'.format(data_dir)
fig, axs = plt.subplots(1, 2, figsize=(5 * len(sampleid_list), 4))

for i, sampleid in enumerate(sampleid_list):

    endnum_path = '{0}clusters_{1}.spot.endnum.tsv'.format(data_dir, sampleid)
    sample_spotid_path = '{0}clusters_{1}.sample.spotid'.format(data_dir, sampleid)
    endnum_data = pd.read_csv(endnum_path, header=0, index_col=None, sep='\t')
    sample_spotid_data = pd.read_csv(sample_spotid_path, header=None, index_col=None, sep='\t')
    sample_spotid = sample_spotid_data.values.flatten()

    sample_endnum = endnum_data[endnum_data['spot_id'].isin(sample_spotid)]
    back_endnum = endnum_data[~endnum_data['spot_id'].isin(sample_spotid)]
    back_endmean = back_endnum['endnum'].mean()
    sample_endnum['ratio'] = (sample_endnum['endnum'] - back_endmean) / sample_endnum['endnum']

    hist = sns.histplot(data=sample_endnum, x='ratio', stat='probability', bins=50, ax=axs[i])
    axs[i].set_title('Sample {0}'.format(i + 1), fontsize=15)
    axs[i].set_xlabel('Valid Information Ratio', fontsize=15)
    axs[i].set_ylabel('Frequency', fontsize=14)

plt.tight_layout()
plt.savefig(plot_path)
