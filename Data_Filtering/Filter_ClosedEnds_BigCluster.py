import numpy as np
import pandas as pd


# Since this experiment does not consider the interaction between chromosomes, the script first splits the multi-chromosome SPRITE into a single chromosome sprite format.
# Then filter the sprite data, containing: "cluster with closer end" and "large fragment cluster".

# Filter the ends that are close (30bp) within the cluster to avoid experimental errors caused by random primers.
# After ultrasound, there were still many chromosome fragments that were not sufficiently broken, and we remove the cluster containing too many ends. 
# Specifically, we divided the chromosomes into bins with a length of 1MB, and then counted the number of bins the cluster end is in. 
# If the number exceeded 40% of the total number of chromosomes, it was determined that the cluster represented a large fragment of the chromosome that had not been broken and would be filtered out.


def chrom_info(resolution: int, threshold: float):

    chrom_data = pd.read_csv('/home/xuyuetong/Reference/mm10_reference/mm10_chrom_sizes.txt', header=None, index_col=None)
    chrom_data.columns = ['chrID', 'chrLength']
    chrom_data['chrID'] = chrom_data['chrID'].apply(lambda x: x.strip('chr'))
    chrom_data['binNum'] = chrom_data['chrLength'].apply(lambda x: int(x / resolution) + 1)
    chrom_data['binThreshold'] = chrom_data['binNum'].apply(lambda x: int(x * threshold))
    return chrom_data


def sprite_intra(data_path: str, save_path: str) -> None:

    with open(data_path, 'r') as d, open(save_path, 'w') as s:

        for i, line in enumerate(d.readlines()):
            
            lines = line.strip().split('\t')
            cluster_id = lines[0]

            chr_list = []
            for site in lines[1:]:
                site_chr = site.split('_')[0]
                chr_list.append(site_chr)
            chr_set = set(chr_list)

            for chr_j in chr_set:
                chr_cluster = []
                cluster_id_chr = '{0}-{1}'.format(cluster_id, chr_j)
                for site in lines[1:]:
                    site_chr = site.split('_')[0]
                    if site_chr == chr_j:
                        chr_cluster.append(site)
                if len(chr_cluster) > 1:
                    s.write('\t'.join([cluster_id_chr] + chr_cluster) + '\n')
   

def filter_cl(input_path: str, output_path: str, threshold: int) -> None:

    with open(input_path, 'r') as d, open(output_path, 'w') as s:
        for i, line in enumerate(d.readlines()):
            lines = line.strip().split('\t')
            sites = lines[1:]
            endnum = len(sites)

            if endnum > 1:
                chr_list, start_list = [], []
                for site in sites:
                    site = site.split('_')
                    chr_list.append(site[0])
                    start_list.append(int(site[1]))

                line_data = pd.DataFrame({'chr':chr_list, 'start':start_list})
                line_data = line_data.sort_values(by=['chr', 'start'])
                line_data['start_s1'] = line_data['start'].shift(1)
                line_data['dis'] = line_data['start'] - line_data['start_s1']
                non_neighbor_data = line_data[(line_data['dis'] < 0) | (line_data['dis'] > threshold)].copy()
                non_neighbornum = non_neighbor_data.shape[0]

                if non_neighbornum > 1:
                    non_neighbor_data['end_id'] = non_neighbor_data['chr'] + '_' + non_neighbor_data['start'].apply(str)
                    writen_endid = non_neighbor_data['end_id'].values.tolist()
                    writen_line = '\t'.join([lines[0]] + writen_endid) + '\n'
                    s.write(writen_line)


def filter_bigfragment(ff_path: str, fb_path: str, resolution: int, threshold_list: list) -> None:

    with open(ff_path, 'r') as ff, open(fb_path, 'w') as fb:

        for line in ff.readlines():

            filter_pass = True
            lines = line.strip().split('\t')
            sites = lines[1:]

            sites_chr, sites_bin = np.array([]), np.array([])
            for site in sites:
                site = site.split('_')
                site_chr = site[0]
                site_loc = site[1]
                sites_chr = np.append(sites_chr, site_chr)
                sites_bin = np.append(sites_bin, int(int(site_loc)/resolution))
            
            sites_chrset = set(sites_chr)
            for chr_a in sites_chrset:
                sites_bin_a = sites_bin[np.where(sites_chr==chr_a)]
                binnum_a = len(set(sites_bin_a))
                if chr_a == 'X':
                    chr_a_index = 19
                elif chr_a == 'Y':
                    chr_a_index = 20
                else:
                    chr_a_index = int(chr_a) - 1
                chr_a_threshold = threshold_list[chr_a_index]
                if binnum_a > chr_a_threshold:
                    filter_pass = False
                    break

            if filter_pass:
                fb.write(line)

    return None


fb_threshold = 0.4
fcl_threshold = 30
resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'

for sample_id in sampleid_list:
        
    sprite_path = '{0}clusters_{1}.sample.sprite'.format(data_dir, sample_id)
    intra_path = '{0}clusters_{1}.sample.intra.sprite'.format(data_dir, sample_id)
    fcl_path = '{0}clusters_{1}.sample.intra.FRI.F1.sprite'.format(data_dir, sample_id)
    fb_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_{2}.sprite'.format(data_dir, sample_id, fb_threshold)

    # sprite -> intra-sprite
    sprite_intra(sprite_path, intra_path)

    # Filter closer end and single-end cluster
    filter_cl(intra_path, fcl_path, fcl_threshold)

    # Filter BigFragment
    chrom_size = chrom_info(resolution, fb_threshold)
    filter_bigfragment(fcl_path, fb_path, resolution, chrom_size['binThreshold'].values)





