import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d


# 本实验的目的是观察肝脏和组织 AB Compartment 的 3D 结构差异。
# 隶属于 A Compartment 的染色体片段用红色表示，B Compartment 的用蓝色表示。
# 本脚本 E1 数据来源为 AB_Compartment_Tissue/Calculate_Saddle.py 脚本（n_eigs=30）。


seed = 0
sample_id = 'MS0612-5'
tissue_id_list = [1, 17, 33, 3, 12, 16, 6, 18, 22, 30, 34]
chrom_id_list = ['2', '4', '5', '6', '8', '9', '11', '15', '17', '19']
dir_path = '/home/xuyuetong/CRICK_Data/Github/Paper/3D_Genome/Data_Hickit/'
save_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/3D_Genome/Scatter_3D_Structure_AB/'
e1_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Compartmentalization_Strength/E1_Value/'
e1_path = '{0}{1}.e1.value'.format(e1_dir, sample_id)
e1_data = pd.read_csv(e1_path, header=0, index_col=None, sep='\t')

for tissue_id in tissue_id_list:
    for chr_id in chrom_id_list:

        print('Tissue ID:{0}\t\tChrom ID:{1}'.format(tissue_id, chr_id))

        e1_chr_data = e1_data[e1_data['chrom'] == 'chr' + chr_id].copy()
        e1_chr_data.dropna(subset=['embryo_e1'], inplace=True)
        e1_chr_data.loc[e1_chr_data['embryo_e1'] >= 0, 'color'] = '#d9544d'
        e1_chr_data.loc[e1_chr_data['embryo_e1'] < 0, 'color'] = '#0485d1'
        ab_data = e1_chr_data[['start', 'color']].copy()
        ab_data.rename(columns={'start': 'BinID'}, inplace=True)

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
        data_3d = pd.merge(data_3d, ab_data, on='BinID')

        # Plot
        plot_path = '{0}Sactter_3DGenome_{1}_C{2}_Chrom{3}.AB.pdf'.format(save_dir, sample_id, tissue_id, chr_id)
        mp4_path = '{0}Sactter_3DGenome_{1}_C{2}_Chrom{3}.AB.mp4'.format(save_dir, sample_id, tissue_id, chr_id)

        fig = plt.figure(figsize=(7, 7))
        ax = plt.axes(projection="3d")
        sctt = ax.scatter3D(data_3d['x'], data_3d['y'], data_3d['z'], c=data_3d['color'],
                            s=400, alpha=1, edgecolor='black')
        for s in range(data_3d.shape[0]):
            line = ax.plot(data_3d['x'].values[s:s+2], data_3d['y'].values[s:s+2],
                           data_3d['z'].values[s:s+2], c=data_3d['color'].values[s], alpha=0.5)
            
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
    
        
