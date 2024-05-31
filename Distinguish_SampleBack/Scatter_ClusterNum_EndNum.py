import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42


# This script displays the number of clusters and ends per spot in a scatter plot format. 
# The horizontal axis represents the Spot ID, and the vertical axis represents the number of clusters or ends. 
# The Spot IDs are sorted in ascending order by number.



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
dir_path = '/home/xuyuetong/CRICK_Data/Github/Paper/Distinguish_SampleBack/'
plot_path = '{0}Scatter_ClusterNum_EndNum.pdf'.format(dir_path)

fig, axs = plt.subplots(figsize=(9, 6), ncols=2, nrows=2, )
fig.subplots_adjust(wspace=0.5, hspace=0.4)


for i, sample_id in enumerate(sampleid_list):

    if sample_id == 'MS0612-5':
        threshold = 0.15
        titile = 'Sample 1'
    else:
        threshold = 0.25
        titile = 'Sample 2'
    cdict = {
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
    cdict = {
        **cdict,
        'alpha': (
            (0.0, 1.0, 1.0),
            # (0.25, 1.0, 1.0),
            (threshold, 0.5, 0.5),
            # (0.75, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        )
    }
    mpl.colormaps.register(LinearSegmentedColormap('BlueRed_{0}'.format(sample_id), cdict))
    blue_red = LinearSegmentedColormap('BlueRed_{0}'.format(sample_id), cdict)

    # Prepare Data
    clusternum_path = '{0}clusters_{1}.spot.clusternum.tsv'.format(dir_path, sample_id)
    spot_clusternum = pd.read_csv(clusternum_path, header=0, index_col=0, sep='\t')
    spot_clusternum.sort_values(by=['clusternum'], inplace=True, ascending=True)
    spot_clusternum['x_loc'] = range(1, spot_clusternum.shape[0] + 1)

    endnum_path = '{0}clusters_{1}.spot.endnum.tsv'.format(dir_path, sample_id)
    spot_endnum = pd.read_csv(endnum_path, header=0, index_col=0, sep='\t')
    spot_endnum.sort_values(by=['endnum'], inplace=True, ascending=True)
    spot_endnum['x_loc'] = range(1, spot_endnum.shape[0] + 1)

    # Plot Cluster Num
    ax_clusternum = axs[i, 0]
    norm1 = mpl.colors.Normalize(vmin=100, vmax=1500)
    ax_clusternum.scatter(spot_clusternum['x_loc'], spot_clusternum['clusternum'], 
                          c=spot_clusternum['clusternum'], cmap=blue_red, norm=norm1)
    ax_clusternum.set_title(titile)
    ax_clusternum.set_xlabel('Spot Num')
    ax_clusternum.set_ylabel('Cluster Num')
    ax_clusternum.set_yscale('log')
    # Plot Colorbar
    cax1 = add_right_cax(ax_clusternum, pad=0.01, width=0.01, orientation='vertical')
    im1 = mpl.cm.ScalarMappable(norm=norm1, cmap=blue_red)
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    
    # Plot end Num
    ax_endnum = axs[i, 1]
    norm2 = mpl.colors.LogNorm(vmin=1.6E2, vmax=1E5)
    ax_endnum.scatter(spot_endnum['x_loc'], spot_endnum['endnum'], 
                      c=spot_endnum['endnum'], cmap='seismic', norm=norm2)
    ax_endnum.set_title(titile)
    ax_endnum.set_xlabel('Spot Num')
    ax_endnum.set_ylabel('End Num')
    ax_endnum.set_yscale('log')
    # Plot Colorbar
    cax2 = add_right_cax(ax_endnum, pad=0.01, width=0.01, orientation='vertical')
    im2 = mpl.cm.ScalarMappable(norm=norm2, cmap='seismic')
    cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical')

plt.savefig(plot_path)