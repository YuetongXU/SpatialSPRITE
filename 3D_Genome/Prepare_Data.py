import subprocess


# The script utilizes the Hickit software to calculate the 3D structure of each chromosome in various tissues.


def write_header(file):
    file.write('## pairs format v1.0\n')
    file.write('#sorted: chr1-chr2-pos1-pos2\n')
    file.write('#shape: upper triangle\n')
    file.write('#chromosome: chr1 195471971\n')
    file.write('#chromosome: chr10 130694993\n')
    file.write('#chromosome: chr11 122082543\n')
    file.write('#chromosome: chr12 120129022\n')
    file.write('#chromosome: chr13 120421639\n')
    file.write('#chromosome: chr14 124902244\n')
    file.write('#chromosome: chr15 104043685\n')
    file.write('#chromosome: chr16 98207768\n')
    file.write('#chromosome: chr17 94987271\n')
    file.write('#chromosome: chr18 90702639\n')
    file.write('#chromosome: chr19 61431566\n')
    file.write('#chromosome: chr2 182113224\n')
    file.write('#chromosome: chr3 160039680\n')
    file.write('#chromosome: chr4 156508116\n')
    file.write('#chromosome: chr5 151834684\n')
    file.write('#chromosome: chr6 149736546\n')
    file.write('#chromosome: chr7 145441459\n')
    file.write('#chromosome: chr8 129401213\n')
    file.write('#chromosome: chr9 124595110\n')
    file.write('#chromosome: chrX 171031299\n')
    file.write('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n')


seed = 0
sample_id = 'MS0612-5'
tissue_id_list = [1, 17, 33, 3, 12, 16, 6, 18, 22, 30, 34]
chrom_id_list = ['2', '4', '5', '6', '8', '9', '11', '15', '17', '19']
pairs_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/Tissue_Chrom_Pairs/'
hickit_pairs_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/Tissue_Chrom_Pairs_3DGenome/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/3D_Genome/Data_Hickit/'

# Pairs -> Hickit Pairs
for tissue_id in tissue_id_list:
    for chr_id in chrom_id_list:
        print('Tissue ID:{0}\t\tChrom ID:{1}'.format(tissue_id, chr_id))
        paris_path = '{0}{1}_C{2}_Chrom{3}.pairs'.format(pairs_dir, sample_id, tissue_id, chr_id)
        hickit_pairs_path = '{0}{1}_C{2}_Chrom{3}.hickit.pairs'.format(hickit_pairs_dir, sample_id, tissue_id, chr_id)
        with open(paris_path, 'r') as s, open(hickit_pairs_path, 'w') as p:
            for line in s.readlines():
                lines = line.strip().split('\t')
                line = '\t'.join(['.', lines[0], lines[1], lines[0], lines[2], '.', '.']) + '\n'
                p.write(line)
        with open(hickit_pairs_path, 'r+') as p: 
            old = p.read()
            p.seek(0)
            write_header(p)
            p.write(old)

# Calculate 3D coordinate by Hickit
hickit_path = '/home/xuyuetong/Tools/hickit-0.1.1_x64-linux/hickit'
for tissue_id in tissue_id_list:
    for chr_id in chrom_id_list:
        print('Tissue ID:{0}\t\tChrom ID:{1}'.format(tissue_id, chr_id))
        hickit_pairs_path = '{0}{1}_C{2}_Chrom{3}.hickit.pairs'.format(hickit_pairs_dir, sample_id, tissue_id, chr_id)
        save_path = '{0}{1}_C{2}_Chrom{3}.seed{4}.3dg'.format(save_dir, sample_id, tissue_id, chr_id, seed)
        subprocess.call('{0} -s {1} -i {2} -r1m -c1 -r10m -c5 -b4m -b1m -O {3}'.format(
            hickit_path, seed, hickit_pairs_path, save_path), shell=True)