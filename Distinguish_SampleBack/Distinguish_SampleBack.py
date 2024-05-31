import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mpl.rcParams['pdf.fonttype'] = 42


# In this script, different thresholds are used to distinguish the categories of spots in the section, 
# in which the spot with endnum exceeding the threshold is the sample area, and the spot with endnum below the threshold is the background area.
# Finally, we used endnum=4000 as the threshold (10min), spots with endnum>4000 as the sample, and spots with endnum<=4000 as the background (which will be filtered out).
# Note that we manually eliminated spots (manually_elim_xy) with endnum>4000 but clearly in the background region.


def format_clean(data_path: str, save_path: str) -> None:

    with open(data_path, 'r') as d, open(save_path, 'w') as s:

        for line in d.readlines():
            lines = line.strip().split('\t')
            cluster_name = lines[0]
            sites = lines[1:]
            sites_clean = []

            for site_a in sites:
                site_a_chr = site_a.split(':')[0].split('chr')[1]
                site_a_site = site_a.split(':')[1].split('-')[0]
                site_a_strand = site_a.split(']')[0].split('[')[1]
                site_a_clean = '{0}_{1}_{2}'.format(site_a_chr, site_a_site, site_a_strand)
                sites_clean.append(site_a_clean)

            sites_clean = '\t'.join(sites_clean)
            s.write('{0}\t{1}\n'.format(cluster_name, sites_clean))


def stats_sprite(input_path: str, output_prefix: str) -> None:

    clusternum_dict, endnum_dict = {}, {}

    with open(input_path, 'r') as f:
        for line in f.readlines():
            lines = line.strip().split('\t')
            spot_id = '.'.join(lines[0].split('.')[4:-1])
            endnum = len(lines) - 1
            if spot_id not in clusternum_dict.keys():
                clusternum_dict[spot_id] = 1
                endnum_dict[spot_id] = endnum
            else:
                clusternum_dict[spot_id] += 1
                endnum_dict[spot_id] += endnum
                
    endnum_data = pd.DataFrame.from_dict(endnum_dict, orient='index', columns=['endnum'])
    endnum_data['spot_id'] = endnum_data.index.values
    endnum_data['x'] = endnum_data['spot_id'].apply(lambda x: 50 - int(x.split('.')[0].split('Bo')[-1]))
    endnum_data['y'] = endnum_data['spot_id'].apply(lambda x: 96 - int(x.split('.')[1].split('Bo')[-1]))
    endnum_data = endnum_data.loc[:, ['spot_id', 'x', 'y', 'endnum']]
    endnum_data.to_csv(output_prefix + '.endnum.tsv', header=True, index=False, sep='\t')

    clusternum_data = pd.DataFrame.from_dict(clusternum_dict, orient='index', columns=['clusternum'])
    clusternum_data['spot_id'] = clusternum_data.index.values
    clusternum_data['x'] = clusternum_data['spot_id'].apply(lambda x: 50 - int(x.split('.')[0].split('Bo')[-1]))
    clusternum_data['y'] = clusternum_data['spot_id'].apply(lambda x: 96 - int(x.split('.')[1].split('Bo')[-1]))
    clusternum_data = clusternum_data.loc[:, ['spot_id', 'x', 'y', 'clusternum']]
    clusternum_data.to_csv(output_prefix + '.clusternum.tsv', header=True, index=False, sep='\t')


def plot_sampleback(input_path: str, plot_path: str, threshold_list: list) -> None:

    spot_endnum = pd.read_csv(input_path, header=0, index_col=0, sep='\t')
    spot_endnum['spot_id'] = spot_endnum.index.values
    spot_endnum['x'] = spot_endnum['spot_id'].apply(lambda x: 50 - int(x.split('.')[0].split('Bo')[-1]))
    spot_endnum['y'] = spot_endnum['spot_id'].apply(lambda x: 96 - int(x.split('.')[1].split('Bo')[-1]))

    fig = plt.figure(figsize=(3.5 * len(threshold_list), 5))
    grid = plt.GridSpec(1, len(threshold_list), wspace=0.3, hspace=0.3)
    for i, threshold in enumerate(threshold_list):
        ax_de = plt.subplot(grid[0, i])
        ax_de.scatter(spot_endnum.loc[spot_endnum['endnum'] >= threshold, 'x'], 
                    spot_endnum.loc[spot_endnum['endnum'] >= threshold, 'y'], 
                    c='red', s=2)
        ax_de.scatter(spot_endnum.loc[spot_endnum['endnum'] < threshold, 'x'], 
                    spot_endnum.loc[spot_endnum['endnum'] < threshold, 'y'], 
                    c='blue', s=2)
        ax_de.set_title(threshold)
    plt.savefig(plot_path)
    plt.close()


def filter_sample(sprite_path: str, endnum_path: str, threshold: int, output_prefix: str, elim_spotid: list = []):
    
    endnum_data = pd.read_csv(endnum_path, header=0, index_col=None, sep='\t')
    sample_spotid = endnum_data[endnum_data['endnum'] > threshold]['spot_id'].values

    if elim_spotid:
        for elim_i in elim_spotid:
            elim_ix, elim_iy = elim_i[0], elim_i[1]
            elim_id = 'evenBo{0}.oddBo{1}'.format(50 - elim_ix, 96 - elim_iy)
            sample_spotid = np.delete(sample_spotid, np.where(sample_spotid == elim_id))
    sampleid_data = pd.DataFrame(sample_spotid)
    sampleid_data.to_csv(output_prefix + '.spotid', header=False, index=False, sep='\t')

    with open(sprite_path, 'r') as f, open(output_prefix + '.sprite', 'w') as o:
        for line in f.readlines():
            lines = line.strip().split('\t')
            spot_id = '.'.join(lines[0].split('.')[4:-1])
            if spot_id in sample_spotid:
                o.write(line)


sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
spotendnum_threshold_list = list(range(1000, 7000, 1000))
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'

for sampleid in sampleid_list:
    
    if sampleid == 'MS0612-5':
        spotendnum_threshold = 4000
        manually_elim_xy = [(5, 68), (5, 88), (6, 68), (8, 90), (8, 92), (12, 95), (23, 2), 
                            (23, 75), (23, 80), (23, 83), (26, 79), (26, 80), (26, 81), (31, 34), 
                            (32, 7), (32, 29), (32, 87), (35, 12), (35, 73), (36, 73), (37, 74), 
                            (39, 94), (41, 11), (42, 23), (46, 92), (47, 30), (48, 7)]
    elif sampleid == 'MS0612-3.odd70':
        spotendnum_threshold = 5000
        manually_elim_xy = []

    sprite_path = '{0}clusters_{1}'.format(data_dir, sampleid)
    clean_path = '{0}clusters_{1}.clean.sprite'.format(data_dir, sampleid)
    spot_prefix = '{0}clusters_{1}.spot'.format(save_dir, sampleid)
    plot_path = '{0}clusters_{1}.spot.endnum.pdf'.format(save_dir, sampleid)
    filter_prefix = '{0}clusters_{1}.sample'.format(data_dir)

    format_clean(sprite_path, clean_path)
    stats_sprite(clean_path, spot_prefix)
    plot_sampleback(spot_prefix + '.endnum.tsv', plot_path, spotendnum_threshold_list)
    filter_sample(clean_path, spot_prefix + '.endnum.tsv', spotendnum_threshold, filter_prefix, manually_elim_xy)
