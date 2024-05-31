import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


# We used the spot-bin matrix after dimensionality reduction as input and applied the K-means method (K=30) to cluster the spots. 
# Due to the instability of clustering results when using different seeds, 
# we selected results whose spatial positions and shapes were more consistent with organ structures. 
# After filtering, the parameters for sample 1 (MS0612-5) are K=30, seed=31; for sample 2 (MS0612-3.odd70), the parameters are K=30, seed=11. 
# Since some spots within same category were not spatially adjacent, for easier analysis, we divided them into multiple categories. 
# Therefore, the final number of categories exceeded 30. Sample 1 was divided into 35 categories, while sample 2 was divided into 32 categories.


def kmeans_plot(data: pd.DataFrame, plot_path: str) -> None:

    label_list = sorted(list(set(data['label'].values)))
    label_num = len(label_list)

    col_num = 6
    if label_num % col_num == 0:
        row_num = label_num // col_num
    else:
        row_num = label_num // col_num + 1
    f, ax = plt.subplots(figsize=(12, row_num * 3), ncols=col_num, nrows=row_num)

    for i, label_i in enumerate(label_list):
        col_idx = i % col_num 
        row_idx = i // col_num
        ax_sp = ax[row_idx, col_idx]
        ax_sp.scatter(data.loc[data['label'] != label_i, 'x'], 
                      data.loc[data['label'] != label_i, 'y'], 
                      c='black', s=0.6)
        ax_sp.scatter(data.loc[data['label'] == label_i, 'x'], 
                      data.loc[data['label'] == label_i, 'y'], 
                      c='red', s=0.6)
        ax_sp.set_title('Cluster {0}'.format(i))
        
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close


k = 30
dim_num = 200
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'

for sampleid in sampleid_list:

    if sampleid == 'MS0612-5':
        seed = 11
    elif sampleid == 'MS0612-3.odd70':
        seed = 31

    matrix_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix'.format(data_dir, sampleid)
    label_path = '{0}KMeans_{1}.label'.format(data_dir, sampleid)
    plot_label_path = '{0}KMeans_{1}.png'.format(data_dir, sampleid)
    label_split_path = '{0}KMeans_{1}.split.label'.format(data_dir, sampleid)
    plot_label_split_path = '{0}KMeans_{1}.split.png'.format(data_dir, sampleid)

    matrix_data = pd.read_csv(matrix_path, header=0, index_col=0, sep='\t')
    matrix_data = matrix_data.iloc[:, :dim_num]
    spot_id = matrix_data.index.values

    # KMeans
    cluster_model = KMeans(n_clusters=k, n_init=20, random_state=seed)
    cluster_model.fit(matrix_data)
    label= cluster_model.predict(matrix_data)
    label_data = pd.DataFrame({'spotid': spot_id, 'label': label})
    label_data['x'] = [50 - int(x.split('.')[0].split('Bo')[-1]) for x in spot_id]
    label_data['y'] = [96 - int(x.split('.')[1].split('Bo')[-1]) for x in spot_id]
    label_data.to_csv(label_path, header=True, index=False, sep='\t')
    kmeans_plot(label_data, plot_label_path)

    # Split
    if sampleid == 'MS0612-5':
        label_data.loc[(label_data['label'] == 3) & (label_data['y'] <= 40), 'label'] = 999
        label_data.loc[(label_data['label'] == 9) & (label_data['x'] <= 35), 'label'] = 30
        label_data.loc[(label_data['label'] == 15) & (label_data['x'] <= 20), 'label'] = 31
        label_data.loc[(label_data['label'] == 16) & (label_data['y'] <= 25), 'label'] = 32
        label_data.loc[(label_data['label'] == 23) & (label_data['y'] <= 30), 'label'] = 33
        label_data.loc[(label_data['label'] == 23) & (label_data['x'] <= 10), 'label'] = 999
        label_data.loc[(label_data['label'] == 25) & (label_data['y'] >= 50), 'label'] = 34
        label_data.loc[(label_data['label'] == 25) & (label_data['y'] >= 20), 'label'] = 999
        label_data.loc[(label_data['label'] == 33) & (label_data['x'] >= 20), 'label'] = 999
    elif sampleid == 'MS0612-3.odd70':
        label_data.loc[(label_data['label'] == 6) & (label_data['y'] <= 60), 'label'] = 30
        label_data.loc[(label_data['label'] == 30) & (label_data['x'] <= 25), 'label'] = 999
        label_data.loc[(label_data['label'] == 8) & (label_data['y'] >= 65), 'label'] = 999
        label_data.loc[(label_data['label'] == 24) & (label_data['x'] <= 20), 'label'] = 31
    label_data.to_csv(label_split_path, header=True, index=False, sep='\t')
    kmeans_plot(label_data, plot_label_split_path)
    



