import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from qnorm import quantile_normalize
from sklearn.decomposition import PCA
mpl.rcParams['pdf.fonttype'] = 42


# The script utilizes PCA for dimensionality reduction of the spot-bin matrix, 
# reducing the number of columns in the matrix to enhance the information content of the feature dimensions. 
# The script plots the accumulated explained variance ratio of the dimensions after PCA processing.


def clip_matrix(matrix_data: pd.DataFrame, threshold: float=0.005) -> pd.DataFrame:
    # Reduce the influence of extreme bin on distribution
    up_limit = matrix_data.quantile(1 - threshold, axis=1)
    matrix_data = matrix_data.clip(upper=up_limit, axis=0)
    return matrix_data


def pca_dim(data: pd.DataFrame, n_components: int, save_path: str=None) -> np.ndarray:
    
    pca = PCA(n_components=n_components)
    redim_data = pca.fit_transform(data)

    exp_array = pca.explained_variance_ratio_
    redim_data = pd.DataFrame(redim_data, index=data.index)
    if save_path:
        redim_data.to_csv(save_path, header=True, index=True, sep='\t')
    return exp_array


def pca_plot(exp_array: np.ndarray, plot_path: str, components_num: int, color=None) -> None:

    exp_sum = []
    for n in range(1, components_num + 1):
        exp_n = exp_array[:n]
        exp_sum.append(sum(exp_n))

    plt.figure(figsize=(6, 5))
    ax = plt.subplot(111)
    ax.scatter(range(1, components_num + 1), exp_sum, c=color)
    ax.set_xlabel('Dimension')
    ax.set_ylabel('Accumulated explained variance ratio')
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close


components_max = 1000
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'

for sampleid in sampleid_list:

    impu_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix'.format(data_dir, sampleid)
    impu_pca_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix'.format(data_dir, sampleid)
    impu_plot_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.pdf'.format(data_dir, sampleid)

    matrix_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.matrix'.format(data_dir, sampleid)
    plot_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.pca.pdf'.format(data_dir, sampleid)
    
    impu_data = pd.read_csv(impu_path, header=0, index_col=0, sep='\t')
    impu_exparray = pca_dim(impu_data, components_max, impu_pca_path)
    pca_plot(impu_exparray, impu_plot_path, components_max, color='#3170a7')

    data_matrix = pd.read_csv(matrix_path, header=0, index_col=0, sep='\t')
    clip_data = clip_matrix(data_matrix)
    zscore_data = stats.zscore(clip_data, axis=1, nan_policy='omit')
    norm_matrix = quantile_normalize(zscore_data, ncpus=10, axis=0)
    exp_array = pca_dim(norm_matrix, components_max)
    pca_plot(exp_array, plot_path, components_max, color='#f7c173')





