import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# This script shows the 3D structure of a chromosome using a scatter plot. 
# Different positions are given different colors to distinguish the head and tail of the chromosome.


seed = 0
sample_id = 'MS0612-5'
tissue_id_list = [1, 17, 33, 3, 12, 16, 6, 18, 22, 30, 34]
chrom_id_list = ['2', '4', '5', '6', '8', '9', '11', '15', '17', '19']
dir_path = '/home/xuyuetong/CRICK_Data/Github/Paper/3D_Genome/Data_Hickit/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/3D_Genome/Scatter_3D_Structure/'

for tissue_id in tissue_id_list:
    for chr_id in chrom_id_list:

        print('Tissue ID:{0}\t\tChrom ID:{1}'.format(tissue_id, chr_id))

        data_path = '{0}{1}_C{2}_Chrom{3}.seed{4}.3dg'.format(dir_path, sample_id, tissue_id, chr_id, seed)
        data_3d = pd.DataFrame()
        with open(data_path, 'r') as d:
            for row in d.readlines():
                if not row.startswith('#'):
                    row = row.strip().split('\t')
                    if (row[0] == 'chr' + chr_id) and (int(row[1]) > 0):
                        row = pd.DataFrame({
                            'chrID': [row[0].strip('chr')],
                            'BinID': [int(row[1])],
                            'x': [float(row[2])],
                            'y': [float(row[3])],
                            'z': [float(row[4])]})
                        data_3d = pd.concat([data_3d, row], axis=0)
        data_3d['chrID'] = data_3d['chrID'].apply(lambda x: int(x))

        # Plot
        plot_path = '{0}Sactter_3DGenome_{1}_C{2}_Chrom{3}.pdf'.format(save_dir, sample_id, tissue_id, chr_id)
        mp4_path = '{0}Sactter_3DGenome_{1}_C{2}_Chrom{3}.mp4'.format(save_dir, sample_id, tissue_id, chr_id)

        fig = plt.figure(figsize=(7, 7))
        ax = plt.axes(projection="3d")
        my_cmap = matplotlib.colormaps.get_cmap('gist_rainbow')
        sctt = ax.scatter3D(data_3d['x'], data_3d['y'], data_3d['z'], c=data_3d['BinID'], cmap=my_cmap,
                            s=400, alpha=1, edgecolor='black')
        
        my_cmap_list = my_cmap(np.linspace(0, 1, data_3d.shape[0]))
        for s in range(data_3d.shape[0]):
            line = ax.plot(data_3d['x'].values[s:s+2], data_3d['y'].values[s:s+2],
                           data_3d['z'].values[s:s+2], c=my_cmap_list[s], alpha=1)
            
        plt.title("Label {0}   1 MB   Chrom {1}".format(tissue_id, chr_id))

        def update(angle):
            angle_norm = (angle + 180) % 360 - 180
            roll = 0
            elev = 30
            azim = angle_norm
            # Update the axis view and title
            ax.view_init(elev, azim, roll)
            if angle == 240:
                plt.savefig(plot_path)

        ani = animation.FuncAnimation(fig, update, frames=range(0, 360 + 1, 10))
        writer = animation.FFMpegWriter()
        ani.save(mp4_path, writer=writer)
        plt.close()

    
        
