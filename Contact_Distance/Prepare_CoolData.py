import subprocess
import pandas as pd


# The script calls Juicer Tools Pre and hic2cool to convert the pairs into cool format.


resolution = 1_000_000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
lable_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'
pairs_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Data/Tissue_Pairs/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Contact_Distance/'

for sample_id in sampleid_list:

    label_path = '{0}KMeans_{1}.split.label'.format(lable_dir, sample_id)
    label_data = pd.read_csv(label_path, header=0, index_col=False, sep='\t')
    label_set = sorted(list(set(label_data['label'].values)))

    for label_i in label_set:
        
        print('Sample ID: {0}\t Label ID: {1}'.format(sample_id, label_i))

        pairs_path = '{0}{1}_C{2}.pairs'.format(pairs_dir, sample_id, label_i)
        onechrom_pairs_path = '{0}{1}_C{2}.onechrom.pairs'.format(pairs_dir, sample_id, label_i)
        save_hic = '{0}Data_HiC/{1}_C{2}.onechrom.hic'.format(save_dir, sample_id, label_i)
        save_cool = '{0}Data_Cool/{1}_C{2}.onechrom.cool'.format(save_dir, sample_id, label_i)
        
        with open(pairs_path, 'r') as p, open(onechrom_pairs_path, 'w') as o:
            for i, line in enumerate(p.readlines()):
                if line.startswith('. chr'):
                    lines = line.split(' ')
                    lines[1] = 'chr1'
                    lines[3] = 'chr1'
                    lines = ' '.join(lines)
                    o.write(lines)
                else:
                    o.write(line)

        print('Pairs -> HiC')
        subprocess.run("java -Xmx2g -jar /home/xuyuetong/Tools/JuicerTools/juicer_tools_1.22.01.jar pre {0} {1} mm10 -r {2} -j 30".format(
            onechrom_pairs_path, save_hic, resolution), shell=True)
        
        print('HiC -> Cool')
        subprocess.run("hic2cool convert {0} {1} -p 30".format(save_hic, save_cool), shell=True)
 