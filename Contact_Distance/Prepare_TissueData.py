import subprocess
import numpy as np
import pandas as pd


# The script converts SPRITE data into Pairs format. 
# To reduce the size of the Pairs file, we omit the Cluster ID information.


def sprtie2pairs(sprite_path: str, pairs_path: str, threshold: int=None) -> None:

    if not threshold:
        threshold = np.inf

    with open(sprite_path, 'r') as f, open(pairs_path, 'w') as h:
        for i, line in enumerate(f.readlines()):
            lines = line.strip().split('\t')
            sites = lines[1:]
            cluster_size = len(sites)
            if cluster_size < threshold:
                sitea_index = list(range(cluster_size))
                siteb_index = list(range(cluster_size))
                for a in sitea_index:
                    siteb_index.remove(a)
                    if len(siteb_index) > 0:
                        for b in siteb_index:
                            site_a_chr, site_a_loc = sites[a].split('_')[0:2]
                            site_b_chr, site_b_loc = sites[b].split('_')[0:2]
                            h.write('. chr{0} {1} chr{2} {3} . .\n'.format(site_a_chr, site_a_loc, site_b_chr, site_b_loc))

    subprocess.run("sort -t ' ' -k 2.4n -k 4.4n -k 3n -k 5n {0} -o {0}".format(pairs_path), shell=True)

    row1 = '## pairs format v1.0'
    row2 = '#columns: readID chr1 position1 chr2 position2 strand1 strand2'
    subprocess.run("sed -i '1i{0}\n' {1}".format(row2, pairs_path), shell=True)
    subprocess.run("sed -i '1i{0}\n' {1}".format(row1, pairs_path), shell=True)


sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
sprite_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'

for sample_id in sampleid_list:

    sprite_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.sprite'.format(sprite_dir, sample_id)
    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    for label_i in label_set:

        print('Sample ID: {0}\t Label ID: {1}'.format(sample_id, label_i))

        save_sprite = '{0}Tissue_SPRITE/{1}_C{2}.sprite'.format(save_dir, sample_id, label_i)
        save_pairs = '{0}Tissue_Pairs/{1}_C{2}.pairs'.format(save_dir, sample_id, label_i)

        label_i_data = label_data[label_data['label'] == label_i]
        spotid_i = label_i_data['spotid'].values
        with open(save_sprite, 'w') as ss, open(sprite_path, 'r') as sp:
            for line in sp.readlines():
                lines = line.strip().split('\t') 
                spot_id = lines[0].split('.')[4:-1]
                spot_id = '.'.join(spot_id)
                if spot_id in spotid_i:              
                    ss.write(line)

        print('SPRITE -> Pairs')
        sprtie2pairs(save_sprite, save_pairs)

