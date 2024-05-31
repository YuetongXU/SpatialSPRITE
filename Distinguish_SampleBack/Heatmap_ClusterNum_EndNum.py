import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42


# The script shows the cluster num and end num of each spot in the section in the form of heatmap.


def add_right_cax(ax, pad, width, orientation='vertical'):
    '''
    在 ax 的右侧 (orientation='vertical') 或下方 (orientation='horizontal'）追加与 ax 等高或等宽的 cax。
    pad 是 cax 与 ax 的间距, width 是 cax 的宽度。
    '''
    axpos = ax.get_position()
    if orientation == 'vertical':
        caxpos = mpl.transforms.Bbox.from_extents(
            axpos.x1 + pad,
            axpos.y0,
            axpos.x1 + pad + width,
            axpos.y1)
    elif orientation == 'horizontal':
        caxpos = mpl.transforms.Bbox.from_extents(
            axpos.x0,
            axpos.y0 - pad - width,
            axpos.x1,
            axpos.y0 - pad)
    else:
        raise ValueError('The Value of Param: orientation is Error')
    cax = ax.figure.add_axes(caxpos)
    return cax


sampleid_list = ['MS0612-5', 'MS0612-3.odd70']
data_dir = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'

for sampleid in sampleid_list:
    
    if sampleid == 'MS0612-5':
        threshold = 0.15
        clustermin, clustermax = 150, 1500
        endmin, endmax = 100, 100000
    elif sampleid == 'MS0612-3.odd70':
        threshold = 0.25
        clustermin, clustermax = 100, 1200
        endmin, endmax = 100, 100000
    cdict1 = {
        'red': (
            (0.0, 0.0, 0.0),
            (threshold, 1.0, 1.0),
            (0.5, 1.0, 1.0),
            (1.0, 0.6, 1.0),
        ),
        'green': (
            (0.0, 0.0, 0.0),
            (threshold, 0.9, 0.9),
            (0.5, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ),
        'blue': (
            (0.0, 0.0, 0.7),
            (threshold, 1.0, 0.8),
            (0.5, 0.0, 0.0),
            (1.0, 0.0, 0.0)
        ),
    }
    cdict1 = {
        **cdict1,
        'alpha': (
            (0.0, 1.0, 1.0),
            # (0.25, 1.0, 1.0),
            (threshold, 0.5, 0.5),
            # (0.75, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        )
    }

    clusternum_path = '{0}clusters_{1}.spot.clusternum.tsv'.format(data_dir, sampleid)
    endnum_path = '{0}clusters_{1}.spot.endnum.tsv'.format(data_dir, sampleid)
    plot_path = '{0}Heatmap_ClusterNum_EndNum.{1}.pdf'.format(data_dir, sampleid)

    clusternum_data = pd.read_csv(clusternum_path, header=0, index_col=None, sep='\t')
    endnum_data = pd.read_csv(endnum_path, header=0, index_col=None, sep='\t')

    # Plot
    mpl.colormaps.register(LinearSegmentedColormap('BlueRed_{0}'.format(sampleid), cdict1))
    blue_red1 = LinearSegmentedColormap('BlueRed1_{0}'.format(sampleid), cdict1)
                                        
    with PdfPages(plot_path) as pdf:
        
        fig = plt.figure(figsize=(9, 6))

        # Cluster Num
        ax = fig.add_subplot(121)
        cax = add_right_cax(ax, pad=0.05, width=0.02, orientation='horizontal')
        norm = mpl.colors.Normalize(vmin=clustermin, vmax=clustermax)
        ax.scatter(clusternum_data['x'], clusternum_data['y'], c=clusternum_data['clusternum'], 
                s=10, norm=norm, cmap=blue_red1)
        ax.set_title('Cluster Num')
        im = mpl.cm.ScalarMappable(norm=norm, cmap=blue_red1)
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        
        # End Num   
        ax = fig.add_subplot(122)
        cax = add_right_cax(ax, pad=0.05, width=0.02, orientation='horizontal')
        norm = mpl.colors.LogNorm(vmin=endmin, vmax=endmax)
        ax.scatter(endnum_data['x'], endnum_data['y'], c=endnum_data['endnum'], 
                s=10, norm=norm, cmap='seismic')
        ax.set_title('End Num')
        im = mpl.cm.ScalarMappable(norm=norm, cmap='seismic')
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        pdf.savefig()
