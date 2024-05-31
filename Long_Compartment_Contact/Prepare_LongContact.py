import pandas as pd


# The script extracts long-distance contacts (>40MB) in each tissue and 
# calculates the contact frequency between compartments of the same type.


add_tissue_id = 1
resolution = 1_000_000
sample_id = 'MS0612-5'
compartment_cate = ['A', 'B']
distance_threshold = 40_000_000
tissue_id_list = [17, 33, 3, 12, 16, 6, 18, 22, 30, 34]
chrom_id_list = ['2', '4', '5', '6', '8', '9', '11', '15', '17', '19']
pairs_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
border_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartment_Contact/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Long_Compartment_Contact/'

# Prepare Data
pairs_path = '{0}Tissue_Pairs/{1}_C{2}.pairs'.format(pairs_dir, sample_id, add_tissue_id)
for chrom_id in chrom_id_list:
    print('Sample ID: {0}\t\tTissue ID: {1}\t\tChrom ID: {2}'.format(sample_id, add_tissue_id, chrom_id))
    output_path = '{0}Tissue_Chrom_Pairs/{1}_C{2}_Chrom{3}.pairs'.format(pairs_dir, sample_id, add_tissue_id, chrom_id)
    with open(pairs_path, 'r') as pairs, open(output_path, 'w') as output:
        for line in pairs.readlines():
            if not line.startswith('#'):
                lines = line.split(' ')
                chrom_i = lines[1]
                if chrom_i == 'chr' + chrom_id:
                    output.write('\t'.join([lines[1], lines[2], lines[4]]) + '\n')

tissue_id_list = tissue_id_list + [add_tissue_id]
for tissue_id in tissue_id_list:
    for chrom_id in chrom_id_list:
        print('Tissue ID:', tissue_id, 'Chrom ID:', chrom_id)
        pairs_path = '{0}Tissue_Chrom_Pairs/{1}_C{2}_Chrom{3}.pairs'.format(pairs_dir, sample_id, tissue_id, chrom_id)
        long_path = '{0}Tissue_Chrom_Pairs_LongDistance/{1}_C{2}_Chrom{3}.LongDistance.pairs'.format(pairs_dir, sample_id, tissue_id, chrom_id)
        with open(pairs_path, 'r') as s, open(long_path, 'w') as p:
            for line in s.readlines():
                lines = line.split('\t')
                contact_start, contact_end = int(lines[1]), int(lines[2])
                if distance_threshold < contact_end - contact_start:
                    p.write(line)

# Calculate Contact Frequency
border_path = '{0}Compartment_AB_border.{1}.tsv'.format(border_dir, sample_id)
border_data = pd.read_csv(border_path, header=0, index_col=None, sep='\t')

inter_contact_data = pd.DataFrame()
for tissue_id in tissue_id_list:
    for chrom_id in chrom_id_list:

        pairs_path = '{0}Tissue_Chrom_Pairs_LongDistance/{1}_C{2}_Chrom{3}.LongDistance.pairs'.format(pairs_dir, sample_id, tissue_id, chrom_id)
        pairs_data = pd.read_csv(pairs_path, header=None, index_col=None, sep='\t')
        pairs_data.columns = ['chrom_id', 'start', 'end']
        pairs_data['start_bin'] = pairs_data['start'].apply(lambda x: int(int(x)/resolution))
        pairs_data['end_bin'] = pairs_data['end'].apply(lambda x: int(int(x)/resolution))

        for compartment in compartment_cate:      
            border_i = border_data[(border_data['chrom'] == 'chr' + chrom_id) &
                                   (border_data['cate'] == compartment)]
            compartment_num = border_i.shape[0]
            com1_index = list(range(compartment_num))
            com2_index = list(range(compartment_num))
            print('Tissue ID: {0}\tChrom ID: {1}\tCompartment: {2}'.format(tissue_id, chrom_id, compartment))

            for m in com1_index:
                compartment_m = border_i.iloc[m, :]
                start_bin_m, end_bin_m = compartment_m['start_bin'], compartment_m['end_bin']
                com2_index.remove(m)
                if len(com2_index) > 0:
                    for n in com2_index:
                        compartment_n = border_i.iloc[n, :]
                        start_bin_n, end_bin_n = compartment_n['start_bin'], compartment_n['end_bin']
                        pairs_n = pairs_data[(start_bin_m <= pairs_data['start_bin']) & 
                                             (pairs_data['start_bin'] <= end_bin_m) & 
                                             (start_bin_n <= pairs_data['end_bin']) & 
                                             (pairs_data['end_bin'] <= end_bin_n)]
                        inter_contact_num = pairs_n.shape[0]
                        if inter_contact_num > 0:
                            mid_bin_m = int((start_bin_m + end_bin_m) / 2)
                            mid_bin_n = int((start_bin_n + end_bin_n) / 2)
                            bin_distance = int(abs(mid_bin_n - mid_bin_m))
                            inter_contact_mn = pd.DataFrame({
                                'tissue': [tissue_id],
                                'chrom': [chrom_id], 
                                'start1_bin': [start_bin_m],
                                'mid1_bin': [mid_bin_m],
                                'end1_bin': [end_bin_m], 
                                'start2_bin': [start_bin_n], 
                                'mid2_bin': [mid_bin_n],
                                'end2_bin': [end_bin_n], 
                                'bin_distance': [bin_distance],
                                'cate': ['{0}-{1}'.format(compartment, compartment)],
                                'contact_num': [inter_contact_num]})
                            inter_contact_data = pd.concat([inter_contact_data, inter_contact_mn])

inter_contact_path = '{0}Compartment_ContactNum_Inter.LongDistance.{1}.tsv'.format(save_dir, sample_id)
inter_contact_data.to_csv(inter_contact_path, header=True, index=False, sep='\t')

