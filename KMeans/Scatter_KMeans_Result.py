import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42


# The script annotated the categories of spots with different colors, with similar colors indicating categories belonging to the same organ.


sample1_color = {0: '#0072ff', 14: '#2a19e0', 17: '#6ba3e9', 33: '#0149e7',
                 1: '#1aa5b7',
                 2: '#e44f4f', 4: '#fe7777', 
                 3: '#bc58f1', 12: '#905cc8', 27: '#5504ab', 
                 5: '#938601', 
                 6: '#0ad376', 7: '#006603', 8: '#589101', 15: '#aed573', 19: '#a8f5d1', 21: '#37af58', 24: '#0dd710', 26: '#6dba53',
                 9: '#81e2e9', 25: '#7fd0f1', 31: '#62e9fa', 
                 10: '#e6c42d', 11: '#cd8519', 13: '#f5e903', 18: '#b6d810', 20: '#e07c00', 22: '#ebcb00', 29: '#d4e056', 30: '#bea400', 34: '#fe9800',
                 16: '#FFAEB9', 
                 23: '#ff0000', 
                 28: '#FFC904', 
                 32: '#f321e7', 
                 999: '#a69ca8'}
sample2_color = {0: '#ffd55c', 6: '#e77828', 9: '#fcd217', 14: '#DEE944', 15: '#F1D755', 
                 16: '#DBAD04', 20: '#D7EC3C', 21: '#FFA34D', 24: '#92611C', 31: '#DD870E',
                 8: '#80489F', 17: '#d4a3ee', 25: '#590F66', 28: '#B867E4',
                 2: '#C30490', 3: '#E85959', 12: '#ff3300', 30: '#7f1616',
                 13: '#CA3F3F', 19: '#FF73BA', 29: '#FFAEB9', 
                 5: '#76D6B6', 7: '#00B306', 10: '#449782', 11: '#0CAC91', 18: '#8A9442', 
                 1: '#677CBC', 4: '#167CCA', 22: '#a1afc9', 23: '#70B8F0', 26: '#307FBB', 27: '#063FEA', 
                 999: '#E3D696'}
sample1_label = {0: 'Abdomen', 14: 'Abdomen', 17: 'Abdomen', 33: 'Abdomen', 
                 1: 'Liver', 
                 2: 'Muscle', 4: 'Muscle', 
                 3: 'Face', 12: 'Face', 27: 'Face', 
                 5: 'Dorsal root ganglion', 
                 6: 'Cartilage primordium', 7: 'Cartilage primordium', 8: 'Cartilage primordium', 15: 'Cartilage primordium', 
                 19: 'Cartilage primordium', 21: 'Cartilage primordium', 24: 'Cartilage primordium', 26: 'Cartilage primordium', 
                 9: 'Connective tissue', 25: 'Connective tissue', 31: 'Connective tissue',
                 10: 'Brain', 11: 'Brain', 13: 'Brain', 18: 'Brain', 20: 'Brain ', 22: 'Brain', 29: 'Brain', 30: 'Brain', 34: 'Brain',
                 16: 'Blood vessel', 
                 23: 'Heart', 
                 28: 'Urogenital ridge',
                 32: 'Kidney', 
                 999: 'None'}
sample2_label = {0: 'Brain', 6: 'Brain', 9: 'Brain', 14: 'Brain', 15: 'Brain', 
                 16: 'Brain', 20: 'Brain', 21: 'Brain', 24: 'Brain', 31: 'Brain',
                 8: 'Face', 17: 'Face', 25: 'Face', 28: 'Face',
                 2: 'Neck', 3: 'Neck', 12: 'Neck', 30: 'Neck', 
                 13: 'Chest', 19: 'Chest', 29: 'Chest', 
                 5: 'Back', 7: 'Back', 10: 'Back', 11: 'Back', 18: 'Back', 
                 1: 'Abdomen', 4: 'Abdomen', 22: 'Abdomen', 23: 'Abdomen', 26: 'Abdomen', 27: 'Abdomen', 
                 999: 'None'}

k = 30
seed = 11
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
sampleid_list = ['MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'

for sampleid in sampleid_list:

    if sampleid == 'MS0612-5':
        seed = 11
        color_dict = sample1_color
        label_dict = sample1_label
        figsize = (3.5, 6)
        point_size = 2.5
    elif sampleid == 'MS0612-3.odd70':
        seed = 31
        color_dict = sample2_color
        label_dict = sample2_label
        figsize = (4, 6)
        point_size = 6

    label_path = '{0}KMeans_{1}.split.label'.format(data_dir, sampleid, seed)
    plot_path = '{0}Scatter_KMeans_Result.{1}.split.pdf'.format(data_dir, sampleid, seed)
    label_data = pd.read_csv(label_path, header=0, index_col=None, sep='\t')

    fig, ax_sp = plt.subplots(figsize=figsize)
    for label_i in color_dict.keys():
        if label_dict:
            label_id = label_dict[label_i]
        else:
            label_id = label_i
        ax_sp.scatter(label_data.loc[label_data['label'] == label_i, 'x'].values, 
                      label_data.loc[label_data['label'] == label_i, 'y'].values, 
                      c=color_dict[label_i], s=point_size, marker='s', label=label_id)
    de_lgnd = ax_sp.legend(prop={'size': 5}, ncol=2)
    for handle in de_lgnd.legend_handles:
        handle.set_sizes([10])
    plt.savefig(plot_path)




