import numpy as np
import pandas as pd


# The script uses the E1 value of the embryo to classify the chromosome segments into two classes: A/B Compartment.


resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
e1_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartmentalization_Strength/E1_Value/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartment_Contact/'

for sample_id in sampleid_list:

    e1_path = '{0}{1}.e1.value'.format(e1_dir, sample_id)
    embryo_e1_data = pd.read_csv(e1_path, header=0, index_col=None, sep='\t')
    chrom_id_list = sorted(list(set(embryo_e1_data['chrom'].values)))
    chrom_end_list = []
    for chrom_id in chrom_id_list:
        e1_idata = embryo_e1_data[embryo_e1_data['chrom'] == chrom_id]
        chrom_end = max(e1_idata['end'].values)
        chrom_end_list.append(chrom_end)

    border_index = np.where(np.diff(embryo_e1_data['embryo_e1']>0))[0]
    border_data = embryo_e1_data.iloc[border_index, :].copy()
    border_data.index = range(border_data.shape[0])

    chrom_list = []
    compartment_data = pd.DataFrame()
    for i, chrom_id in enumerate(chrom_id_list):
        border_idata = border_data[border_data['chrom'] == chrom_id]
        border_inum = border_idata.shape[0]
        for j in range(border_inum):
            border_j = border_idata.iloc[j, :]
            e1_j = border_j['embryo_e1']
            start_j = border_j['start']
            end_j = border_j['end']
            compartment_i = pd.DataFrame({'chrom': [chrom_id], 'start': [None], 'end': [end_j], 'cate': [None]})
            if j == 0:    
                compartment_i['start'] = 0
                if e1_j > 0:
                    compartment_i['cate'] = 'A'
                else:
                    compartment_i['cate'] = 'B'
            else:
                compartment_i['start'] = last_border
                if e1_j > 0:
                    compartment_i['cate'] = 'A'
                else:
                    compartment_i['cate'] = 'B'
            chrom_iend = chrom_end_list[i]
            if (j == border_inum - 1) and (end_j < chrom_iend):
                if e1_j > 0:
                    end_cate = 'B'
                else:
                    end_cate = 'A'
                compartment_iend = pd.DataFrame({'chrom': [chrom_id], 'start': [end_j], 
                                                'end': [chrom_end_list[i]], 'cate': [end_cate]})
                compartment_i = pd.concat([compartment_i, compartment_iend], ignore_index=True)
            last_border = end_j
            compartment_data = pd.concat([compartment_data, compartment_i], ignore_index=True)
    compartment_data['start_bin'] = compartment_data['start'].map(lambda x: int(int(x)/resolution))
    compartment_data['end_bin'] = compartment_data['end'].map(lambda x: int(int(x)/resolution))
    compartment_path = '{0}Compartment_AB_border.{1}.tsv'.format(save_dir, sample_id)
    compartment_data.to_csv(compartment_path, index=False, header=True, sep='\t')

