import math
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from qnorm import quantile_normalize


# This script performs z-score normalization, imputation, and qnorm normalization on the data. 
# First, it uses the clip_matrix() function to clip bins with extreme interaction counts within each spot (Top 0.5%). 
# Then, it normalizes using z-score, where bins with no interactions (end num=0) are treated as "None". 
# Next, it imputes missing values using the mean of interaction information from neighboring spots (N=2). 
# Finally, because the imputation process may alter the distribution of interaction data, the script applies qnorm for re-normalization of the data.


def clip_matrix(matrix_data: pd.DataFrame, threshold: float=0.005) -> pd.DataFrame:
    # Reduce the influence of extreme bin on distribution
    up_limit = matrix_data.quantile(1 - threshold, axis=1)
    matrix_data = matrix_data.clip(upper=up_limit, axis=0)
    return matrix_data


def normalized_zscore(matrix_data: pd.DataFrame) -> pd.DataFrame:
    from scipy import stats
    matrix_data.replace(0, np.nan, inplace=True)
    matrix_data = stats.zscore(matrix_data, axis=1, nan_policy='omit')
    return matrix_data


def ex_coordinate(data: pd.DataFrame) -> pd.DataFrame:

    coordinate_data = pd.DataFrame(index=data.index.values)
    coordinate_data['x'] = [50 - int(x.split('.')[0].split('Bo')[-1]) for x in data.index.values]
    coordinate_data['y'] = [96 - int(x.split('.')[1].split('Bo')[-1]) for x in data.index.values]
    
    return coordinate_data


def imputation(matrix_data: pd.DataFrame, neighbor_num: int) -> pd.DataFrame:

    coordinate_data = ex_coordinate(matrix_data)
    matrix_data = coordinate_data.merge(matrix_data, left_index=True, right_index=True)

    data_imputation = None
    spot_num = matrix_data.shape[0]
    
    for spot_i in range(spot_num):
        spot_i_series = matrix_data.iloc[spot_i, :]
        spot_i_x = spot_i_series['x']
        spot_i_y = spot_i_series['y']
        x_min, x_max = spot_i_x - neighbor_num, spot_i_x + neighbor_num
        y_min, y_max = spot_i_y - neighbor_num, spot_i_y + neighbor_num
        spot_i_neighbor = matrix_data[(x_min <= matrix_data['x']) & (matrix_data['x'] <= x_max) & 
                                      (y_min <= matrix_data['y']) & (matrix_data['y'] <= y_max)]
        spot_i_neighbor = spot_i_neighbor.iloc[:, 2:]
        spot_i_neighbor = spot_i_neighbor.mean().to_frame()
        spot_i_neighbor.columns = [spot_i_series.name]
        if spot_i == 0:
            data_imputation = spot_i_neighbor
        else:
            data_imputation = pd.merge(data_imputation, spot_i_neighbor, left_index=True, right_index=True, how='inner')

    return data_imputation.T


neighbor_num = 2
sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/KMeans/'

for sample_id in sampleid_list:

    matrix_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.matrix'.format(data_dir, sample_id)
    output_path = '{0}clusters_{1}.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix'.format(data_dir, sample_id)

    # Normalization zscore and Imputation
    data_matrix = pd.read_csv(matrix_path, header=0, index_col=0, sep='\t')
    clip_data = clip_matrix(data_matrix)
    zscore_data = normalized_zscore(clip_data)
    impu_data = imputation(zscore_data, neighbor_num)

    # Filter by missrate
    rowmiss = impu_data.isnull().sum(axis=1)
    filter_spotid = rowmiss[rowmiss == 0].index
    impu_data = impu_data.loc[filter_spotid, :]

    # Normalization qnorm
    norm_matrix = quantile_normalize(impu_data, ncpus=10, axis=0)
    norm_matrix.to_csv(output_path, header=True, index=True, sep='\t')

