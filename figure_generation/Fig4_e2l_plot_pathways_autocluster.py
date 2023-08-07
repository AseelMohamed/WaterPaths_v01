# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 12:47:51 2020

Does the infrastructure pathway plotting for Descoberto and Santa Maria.

@author: dgold
@edited: lbl59
"""
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances_argmin
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib

SOL = 41
RDM = 481
NUM_WEEKS = 2080
NUM_REALIZATIONS = 1000
NUM_RDMS = 1000
FILEPATH = 'C:/Users/lbl59/Desktop/CAESB/Pathways/Pathways_s' + str(SOL) + '_RDM' + str(RDM)  + '.out'

def auto_cluster(cluster_input):
    '''
    Automatically determines the optimal number of clusters for a given dataset
    :param cluster_input: an array with each row a realization and each column
    '''
    # Reference: https://medium.com/analytics-vidhya/how-to-determine-the-optimal-k-for-k-means-708505d204eb 
    # Implement the silhouette score method to automatically
    # determine the optimal number of clusters

    silhouette = np.zeros(3, dtype=float)   # 2,3,4 clusters
    opt_num_clusters = 0
    curr_num_labels = 0

    # test 5 clusters, unlikely to be more 
    for num_clusters in range(2,len(silhouette)+2):
        k_mean_model = KMeans(init='k-means++', n_clusters = num_clusters, n_init=10).fit(cluster_input)
        labels = k_mean_model.labels_
        if len(np.unique(labels)) >= curr_num_labels and curr_num_labels < 3:
            curr_num_labels = len(np.unique(labels))

    opt_num_clusters = curr_num_labels
    print(opt_num_clusters)
        
    return opt_num_clusters


def get_clusters(sol, rdm, utility):
    '''
    Extracts the pathways from the output file and formats them for clustering
    :param sol: the solution number
    :param rdm: the realization number
    :param utility: the utility of interest

    :return: cluster_input: an array with each row a realization and each column
    '''
    # Get the pathways from the output file
    filepath =FILEPATH
    
    pathways_df = pd.read_csv(filepath, sep='\t')
    pathways_df['infra.'] = pathways_df['infra.'] - 6   # convert the infrastructure 6-12 ranking to a 0-6 ranking

    # set up mapping for each infrastructure to a number for each utility
    infra_dict = {'corumba1': 0, 'corumba2': 1, 'corumba3': 2, 'reservoir': 3, 
                  'paranoa1': 4, 'paranoa2': 5, 'paranoa3': 6}

    utility_dict = {'santa_maria': 0, 'descoberto': 1}

    infra_desc = list(map(infra_dict.get, ['corumba1', 'corumba2', 'corumba3', 'reservoir']))
    infra_sm = list(map(infra_dict.get, ['paranoa1', 'paranoa2', 'paranoa3']))
    
    # reformat for clustering (need an array with each row a realization and
    # each column a different infrastructure option. elements are const weeks)
    cluster_input = np.ones([NUM_REALIZATIONS,22])*NUM_WEEKS  

    # loop through each realization
    for real in range(0,NUM_REALIZATIONS):
    # extract the realization
        current_real = pathways_df[pathways_df['Realization']==real]
        # find the infrastructure option (ids 0-2 are already built, 6 is off)
        for inf in range(max(pathways_df['infra'])):   # Fed District infrastructure ranks from 6-12
            if not current_real.empty:
                for index, row in current_real.iterrows():
                    if row['infra.']==inf:
                        cluster_input[real, inf] = row['week']

    # post process to remove inf options never constructed and normalize weeks
    # to [0-1] by dividing by total weeks, 2080
    # cluster_input = cluster_input[:,[4,5,7,8,9,10,11,12]]/2080
    cluster_input = cluster_input/NUM_WEEKS
    print(cluster_input.shape)

    # extract columns for each utility
    if utility == 'Santa Maria':
        cluster_input = cluster_input[:, infra_sm]
    elif utility == 'Descoberto':
        cluster_input = cluster_input[:, infra_desc]
     
    return cluster_input 

def cluster_pathways(cluster_input, num_clusters):
    """
    Clusters infrastructure pathways be the week each option is constructed
    creates "representative pathways" for diagnostics and communication

    Parameters:
        solution: name of the solution to be plotted (should be folder name)
        compSol_full: name of the compromise solution to plot
        mode: baseline, worst regional robustness or best regional robustness
        utility: a string (all lowercase) of the name of the utility of interest
        num_clusters: number of clusters used (should be carefully assessed)

    returns:
        cluster_pathways: a 3-d list containing the pathways for each cluster
        cluster_medians: an array containing median construction weeks for each
        inf option in each cluster

    """
    k_means = KMeans(init='k-means++', n_clusters = num_clusters, n_init=10)

    k_means.fit(cluster_input)
    k_means_cluster_centers = k_means.cluster_centers_
    k_means_labels = pairwise_distances_argmin(cluster_input, k_means_cluster_centers)

    # assign each realization to a pathway, and calculate the median week
    # each infrstructure option is constructed in each cluster
    cluster_pathways = []
    cluster_medians = []

    for i in range(0, num_clusters):
        current_cluster =  cluster_input[k_means_labels==i,:]*NUM_WEEKS
        cluster_pathways.append(current_cluster)
        current_medians = np.zeros(len(current_cluster[0,:]))
        
        for j in range(0, len(current_cluster[0,:])):
            current_medians[j]= np.median(current_cluster[:,j])

        cluster_medians.append(current_medians)

    # sort clusters by average of medians to get heavy, mod and light clusters
    cluster_means = []
    for n in range(num_clusters):
        cluster_means.append(np.mean(cluster_medians[n]))

    sorted_indices = np.argsort(cluster_means)

    cluster_med = []
    for s in range(len(sorted_indices)):
        cluster_med.append(cluster_medians[sorted_indices[s]])

    cluster_medians = np.reshape(cluster_med, (len(sorted_indices),-1))
    '''
    cluster_medians = np.vstack((cluster_medians[sorted_indicies[2]],
                                     cluster_medians[sorted_indicies[1]],
                                     cluster_medians[sorted_indicies[0]]))
    '''
    return cluster_pathways, cluster_medians


def plot_single_pathway(cluster_medians, cluster_pathways, inf_options_idx,
                        c, cmap, ax, y_offset, plot_legend):
    """
    Makes a plot of an infrastructure Pathway

    Parameters:
        cluster_medians: an array with median weeks each option is built
        cluster_pathways: an array with every pathway in the cluster
        inf_options: an array with numbers representing each option (y-axis vals)
        should start at zero to represent the "baseline"
        c: color to plot pathway
        ax: axes object to plot pathway
        inf_names: a list of strings with the names of each pathway

    """
    # get array of the infrastructure options without baseline
    inf_options_idx_no_baseline=inf_options_idx[1:]

    sorted_inf = np.argsort(cluster_medians)
    print(sorted_inf)
    # plot heatmap of construction times
    cluster_pathways = np.rint(cluster_pathways/45)
    inf_im = np.zeros((45, np.shape(cluster_pathways)[1]+1))

    for k in range(1,np.shape(cluster_pathways)[1]+1) :
        for i in range(0,45):
            for j in range(0, len(cluster_pathways[:,k-1])):
                if cluster_pathways[j,k-1] == i:
                    inf_im[i,k] +=1
        #print(min(inf_im[:,k]))
        #print(max(inf_im[:,k]))

    ax.imshow((inf_im.T)/0.85, cmap=cmap, aspect='auto', alpha = 0.6)
    
    # sort by construction order
    #cluster_medians = np.rint(cluster_medians/45)
    #sorted_inf = np.argsort(cluster_medians)

    # plot pathways
    # create arrays to plot the pathway lines. 
    # To ensure pathways have corners, we need an array to have length 2*num_inf_options
    pathway_x = np.zeros(len(cluster_medians)*2+2)
    pathway_y = np.zeros(len(cluster_medians)*2+2)

    # to make corners, each inf option must be in the slot it is triggered, and
    # the one after
    cluster_medians = np.rint(cluster_medians/45)
    for i in range(0,len(cluster_medians)):
        for j in [1,2]:
            pathway_x[(i*2)+j] = cluster_medians[sorted_inf[i]]
            pathway_y[(i*2)+j+1] = inf_options_idx_no_baseline[sorted_inf[i]]
            

    # end case
    pathway_x[-1] = 45

    # plot the pathway line
    ax.plot(pathway_x, pathway_y+y_offset, color=c, linewidth=3,
            alpha = .9, zorder=1)

    ax.set_xlim([0,44])


def create_cluster_plots(util_meds, util_pathways,
                         n_clusters, util_inf, util_name,
                         cluster_colors_all, ax_num,
                         fig, gspec, fig_col, plot_legend):
    
    """
    creates a figure with three subplots, each representing a utility

    Parameters:
        
    """
    
    #fig.text(0.5, 0.01, 'Years', ha='center', va='center')
    #fig.text(0.01, 0.5, 'Infrastructure options number', ha='center', va='center', rotation='vertical')

    ax = fig.add_subplot(gspec[ax_num, fig_col])

    y_offsets = np.linspace(-0.15, 0.15, num=n_clusters)

    for i in np.arange(n_clusters):
        
        plot_single_pathway(util_meds[i], util_pathways[i], np.arange(0, len(util_inf)),
                              cluster_colors_all[i], 'bone_r', ax, 
                              y_offsets[i], plot_legend)
        
        if fig_col == 0:
            ax.set_ylabel(util_name, fontsize=14)
            ax.set_yticks(np.linspace(0, len(util_inf), num=len(util_inf)))
            ax.set_yticklabels(util_inf)

        else:
            #ax1.set_yticks(np.linspace(-0.65, 2.15, 5))
            ax.set_yticks(np.linspace(0, len(util_inf), len(util_inf)))
            ax.set_yticklabels(['']*len(util_inf))

 
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)
    #if fig_col == 3:
    #ax1.set_axis_off()
    if ax_num == len(util_name)-1:
        divider1 = make_axes_locatable(ax)
        cax1 = divider1.append_axes('bottom', size='5%', pad=0.3)
        cbar1 = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='bone_r'), 
                                shrink=0.89, fraction=0.1, pad=0.05,  \
                                orientation='horizontal', cax=cax1)
        cbar1.set_label(r'Constr. freq. $\rightarrow$', labelpad=5, fontsize=12)
    
    plt.tight_layout()

sm_inf = ['Baseline', 'Paranoa1', 'Paranoa2', 'Paranoa3']
desc_inf = ['Baseline', 'Corumba1', 'Corumba2', 'Corumba3', 'Reservoir']

util_names = ['Santa Maria', 'Descoberto']
util_infrastructure = {'Santa Maria':sm_inf, 'Descoberto': desc_inf}

cluster_colors_all = ['darkgoldenrod', 'orange', 'golden', 'navajowhite']

'''
Change these depending on the number of RDMs and solutions available
'''

for s in range(SOL,SOL+1):
    for i in range(RDM,RDM+1):
        rdm = i
        title_label = 'Pathways for soln ' + str(s) + ' in RDM ' + str(i)   

        fig_filepath = 'C:/Users/lbl59/Desktop/CAESB/Pathways/Figures/'
        fig_fileloc = fig_filepath + title_label + '.pdf'
        
        fig = plt.figure(figsize=(len(util_names),12), dpi=300)
        gspec = fig.add_gridspec(nrows=len(util_names), ncols=1, height_ratios =[1]*len(util_names), width_ratios=[1])
        
        heavy_inf_color = ""
        mid_inf_color = ""
        light_inf_color = ""

        for u in range(len(util_names)):
            utility_name = util_names[u]
            util_inf = util_infrastructure[utility_name]
            print(util_inf)

            cluster_input_util = get_clusters(s, rdm, utility_name)

            num_clusters_util = auto_cluster(cluster_input_util)

            color_util = []

            for n in range(num_clusters_util):
                #print(n)
                color_util.append(cluster_colors_all[n])
            '''
            for m in range(num_clusters_desc):
                color_desc.append(cluster_colors_all[m])
            '''


            util_cluster_pathways, util_cluster_meds = cluster_pathways(cluster_input_util, num_clusters_util)
            #desc_cluster_pathways, desc_cluster_meds = cluster_pathways(cluster_input_desc, num_clusters_desc)
            create_cluster_plots(util_cluster_meds, util_cluster_pathways, num_clusters_util, 
                                util_inf, utility_name, cluster_colors_all, u,
                                fig, gspec, 0, False)

        print(f"solution = {s} done")

        plt.suptitle(title_label, fontsize=12, y=1.00)
        #plt.tight_layout()
        plt.savefig(fig_fileloc)