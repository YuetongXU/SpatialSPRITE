import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import resample  

matplotlib.rcParams['pdf.fonttype'] = 42


# The script uses a violinplot to show the distribution of contact frequency within and across compartments.
# 10 tissues were selected for analysis according to high, medium, and low compartmentalization strength.
# High compartmentalization strength: 17 (abdominal), 33 (abdominal).
# Medium compartmentalization strength: 3 (face), 12 (face), and 16 (lung).
# Low compartmentalization strength: 6 (brain), 18 (brain), 22 (brain), 30 (brain), 34 (brain).
# boostrapping is used to determine the difference in the distribution of contact frequency between different compartmentalization strengths


top_ratio = 0.05
sample_id = 'MS0612-5'
compartment_list = ['A', 'B']
contact_cate_list = ['Inter', 'Intra']
strength_list = ['high', 'high', 
                 'midden', 'midden', 'midden', 
                 'low', 'low', 'low', 'low', 'low']
tissue_id_list = [17, 33, 16, 3, 12, 18, 34, 6, 22, 30]
# The number of contacts in tissue
tissue_pairsnum = {17: 40367139, 33: 47212489, 
                   16: 34212374, 3: 137169077, 12: 18153921, 
                   18: 8180270, 34: 18722738, 6: 24548647, 22: 23970515, 30: 12355835}

data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartment_Contact/'
strength_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartmentalization_Strength/'

for contact_cate in contact_cate_list:

    compartcontact_path = '{0}Compartment_ContactNum_{1}.{2}.tsv'.format(data_dir, contact_cate, sample_id)
    compartcontact_data = pd.read_csv(compartcontact_path, header=0, index_col=None, sep='\t')
    compartcontact_data = compartcontact_data.loc[:, ['tissue', 'chrom', 'cate', 'contact_num']]

    # Color by Strength
    color_vmin = 1.0
    color_vmax = 2.0
    corner_extent = 5
    embryo_strength_path = '{0}Strength_{1}.tsv'.format(strength_dir, sample_id)
    embryo_strength = pd.read_csv(embryo_strength_path, header=0, index_col=0, sep='\t')
    embryo_strength = embryo_strength.loc[tissue_id_list, :]
    strength_e = embryo_strength.iloc[:, corner_extent+1].values
    strength_e = (strength_e - color_vmin)/(color_vmax - color_vmin)
    com_cmap = matplotlib.colormaps.get_cmap('inferno')
    com_cmap_list = com_cmap(strength_e)
    com_cmap_list = list(com_cmap_list[[0, 3, 7], :])

    # boostrapping params
    iter = 1000
    strength_cate_list = ['high', 'midden', 'low']
    strength_cate_dict = {'high': 'Abdomen', 'midden': 'Face', 'low': 'Brain'}

    # Plot
    for compartment_cate in compartment_list:

        if contact_cate == 'Inter':
            compartment_str = '{0}-{0}'.format(compartment_cate)
            violin_ylim = 3e-3
            boostrapping_ymax = 4.5e-4
            boostrapping_ymin = 2.0e-4
        else:
            compartment_str = compartment_cate
            violin_ylim = 8e-3
            boostrapping_ymax = 15e-4
            boostrapping_ymin = 3e-4

        # Prepare Data
        tissue_data = pd.DataFrame()
        for i, tissue_id in enumerate(tissue_id_list):
            tissue_strength = strength_list[i]
            contact_num = compartcontact_data.loc[
                (compartcontact_data['cate'] == compartment_str) & 
                (compartcontact_data['tissue'] == tissue_id), 'contact_num'].copy()
            if not contact_num.empty:
                contact_num = contact_num.to_frame()
                contact_sorted = np.sort(contact_num.values.flatten())
                threshold_idx = int(len(contact_sorted) * top_ratio)
                top_threshold = contact_sorted[-threshold_idx]
                contact_num = contact_num[contact_num['contact_num'] < top_threshold]
                if not contact_num.empty:
                    contact_num['tissue_id'] = [tissue_id] * contact_num.shape[0]
                    contact_num['cate'] = [compartment_cate] * contact_num.shape[0]
                    contact_num['strength'] = [tissue_strength] * contact_num.shape[0]
                    contact_num.columns = ['contact_num', 'tissue_id', 'cate', 'strength']
                    contact_num['contact_num'] = contact_num['contact_num']/tissue_pairsnum[tissue_id]
                    tissue_data = pd.concat([tissue_data, contact_num], ignore_index=True, axis=0)
        tissue_data['tissue_id'] = tissue_data['tissue_id'].apply(lambda x: str(x))

        # Violinplot
        fig = plt.figure(figsize=(8, 4))
        plot_data = tissue_data[tissue_data['cate'] == compartment_cate]
        sns.violinplot(x='strength', y='contact_num', data=plot_data, hue='strength', palette=com_cmap_list, inner='box')
        plt.title('{0}  {1}-Compartment'.format(compartment_cate, contact_cate))
        plt.xlabel('Tissue ID')
        plt.ylabel('Frequence')
        if compartment_cate == 'A':
            plt.ylim(0, violin_ylim)
        elif compartment_cate == 'B':
            plt.ylim(0, violin_ylim)
        plot_path = '{0}Violinplot_Compartment{1}_ContactFrequency_{2}.{3}.pdf'.format(data_dir, compartment_cate, contact_cate, sample_id)
        plt.savefig(plot_path)

        # boostrapping plot
        fig = plt.figure(figsize=(15, 7))
        ax_scatter = plt.subplot2grid((1, 3), (0, 0), rowspan=1, colspan=2)
        ax_distribution = plt.subplot2grid((1, 3), (0, 2), rowspan=1, colspan=1)

        plot_data = pd.DataFrame()
        for s, strength_i in enumerate(strength_cate_list):
            mean_list = []
            mean_data = pd.DataFrame()
            contact_num = tissue_data.loc[tissue_data['strength'] == strength_i, 'contact_num'].values
            for i in range(iter):
                samples_bootstrap = resample(contact_num, n_samples=int(len(contact_num)*0.7), replace=True)
                mean_list.append(np.mean(samples_bootstrap))
            mean_data = pd.DataFrame({'contact_num': mean_list, 'label': [strength_cate_dict[strength_i]] * iter})
            plot_data = pd.concat([plot_data, mean_data], ignore_index=True)
            ax_scatter.scatter(x=range(iter), y=mean_list, label=strength_cate_dict[strength_i], color=com_cmap_list[s])
        ax_scatter.set_xlabel('Repetitions')
        ax_scatter.set_ylabel('{0} {1}-Compartment Contact Frequence'.format(compartment_cate, contact_cate))
        ax_scatter.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
        ax_scatter.set_ylim(boostrapping_ymin, boostrapping_ymax)
        ax_scatter.legend()

        sns.kdeplot(data=plot_data, y="contact_num", hue='label', ax=ax_distribution, palette=com_cmap_list, fill=True, alpha=.85)
        ax_distribution.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
        ax_distribution.set_ylabel('{0} {1}-Compartment Contact Frequence'.format(compartment_cate, contact_cate))
        ax_distribution.set_ylim(boostrapping_ymin, boostrapping_ymax)

        plot_path = '{0}Bootstrapping_Compartment{1}_ContactFrequency_{2}.{3}.pdf'.format(data_dir, compartment_cate, contact_cate, sample_id)
        plt.savefig(plot_path)