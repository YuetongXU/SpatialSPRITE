import subprocess
import numpy as np


# The script converts SPRITE data into Pairs format. 
# To reduce the size of the Pairs file, we omit the Cluster ID information.
# Then, the script calls Juicer Tools Pre and hic2cool to convert the pairs into cool format.


def sprtie2pairs(sprite_path: str, pairs_path: str, threshold: int=None) -> None:

    if not threshold:
        threshold = np.inf

    with open(sprite_path, 'r') as f, open(pairs_path, 'w') as p:
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
                            p.write('. chr{0} {1} chr{2} {3} . .\n'.format(site_a_chr, site_a_loc, site_b_chr, site_b_loc))

    subprocess.run("sort -t ' ' -k 2.4n -k 4.4n -k 3n -k 5n {0} -o {0}".format(pairs_path), shell=True)

    row1 = '## pairs format v1.0'
    row2 = '#columns: readID chr1 position1 chr2 position2 strand1 strand2'
    subprocess.run("sed -i '1i{0}\n' {1}".format(row2, pairs_path), shell=True)
    subprocess.run("sed -i '1i{0}\n' {1}".format(row1, pairs_path), shell=True)


resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Heatmap/'

for sample_id in sampleid_list:

    data_prefix = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4'.format(data_dir, sample_id)
    save_prefix = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4'.format(save_dir, sample_id)
    sprite_path = data_prefix + '.sprite'
    pairs_path = data_prefix + '.pairs'
    hic_path = save_prefix + '.hic'
    cool_path = save_prefix + '.cool'

    sprtie2pairs(sprite_path, pairs_path)
    print('SPRITE to Pairs is Finish')
    subprocess.run("java -Xmx2g -jar /home/xuyuetong/Tools/JuicerTools/juicer_tools_1.22.01.jar pre {0} {1} mm10 -r {2} -j 30".format(
        pairs_path, hic_path, resolution), shell=True)
    print('Pairs to HiC is Finish')
    subprocess.run("hic2cool convert {0} {1} -p 30".format(hic_path, cool_path), shell=True)
    print('HiC to Cool is Finish')

