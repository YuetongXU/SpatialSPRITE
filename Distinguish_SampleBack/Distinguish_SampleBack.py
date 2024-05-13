import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# In this script, different thresholds are used to distinguish the categories of spots on the slice, 
# in which the spot with legnum exceeding the threshold is the sample area, and the spot with legnum below the threshold is the background area.
# Finally, we chose LegNum=4000 as the threshold.


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


def sprite_matrix(input_path: str, output_path: str) -> None:

    spot_dict = {}

    with open(input_path, 'r') as f:
        for line in f.readlines():
            lines = line.strip().split('\t')
            spot_id = '.'.join(lines[0].split('.')[4:-1])
            legnum = len(lines) - 1
            if spot_id in spot_dict.keys():
                spot_dict[spot_id] = spot_dict[spot_id] + legnum
            else:
                spot_dict[spot_id] = legnum

    legnum_data = pd.DataFrame.from_dict(spot_dict, orient='index')
    legnum_data.columns = ['legnum']
    legnum_data.to_csv(output_path, header=True, index=True, sep='\t')


def plot_sampleback(input_path: str, plot_path: str, threshold_list: list) -> None:

    spot_legnum = pd.read_csv(input_path, header=0, index_col=0, sep='\t')
    spot_legnum['spot_id'] = spot_legnum.index.values
    spot_legnum['x'] = spot_legnum['spot_id'].apply(lambda x: 50 - int(x.split('.')[0].split('Bo')[-1]))
    spot_legnum['y'] = spot_legnum['spot_id'].apply(lambda x: 96 - int(x.split('.')[1].split('Bo')[-1]))

    with PdfPages(plot_path) as pdf:

        fig = plt.figure(figsize=(3.5 * len(threshold_list), 5))
        grid = plt.GridSpec(1, len(threshold_list), wspace=0.3, hspace=0.3)

        for i, threshold in enumerate(threshold_list):
            ax_de = plt.subplot(grid[0, i])
            ax_de.scatter(spot_legnum.loc[spot_legnum['legnum'] >= threshold, 'x'], 
                        spot_legnum.loc[spot_legnum['legnum'] >= threshold, 'y'], 
                        c='red', s=2)
            ax_de.scatter(spot_legnum.loc[spot_legnum['legnum'] < threshold, 'x'], 
                        spot_legnum.loc[spot_legnum['legnum'] < threshold, 'y'], 
                        c='blue', s=2)
            ax_de.set_title(threshold)

        pdf.savefig()


data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'
sprite_path = '{0}clusters_MS0612-5'.format(data_dir)
clean_path = '{0}clusters_MS0612-5.clean.sprite'.format(data_dir)
legnum_path = '{0}clusters_MS0612-5.spot.legnum.tsv'.format(save_dir)
plot_path = '{0}clusters_MS0612-5.spot.legnum.pdf'.format(save_dir)
spotlegnum_threshold_list = list(range(1000, 7000, 1000))

format_clean(sprite_path, clean_path)
sprite_matrix(clean_path, legnum_path)
plot_sampleback(legnum_path, plot_path, spotlegnum_threshold_list)


