#%%
import pandas as pd
from pandas.plotting import parallel_coordinates
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from copy import deepcopy
import seaborn as sns
from matplotlib.gridspec import GridSpec
sns.set(style='whitegrid')
from mpl_toolkits.axes_grid1 import make_axes_locatable

# load objectives
ref_set = pd.read_csv('formatted_objectives.csv')

# reverse the order of Rel, on the parallel axis plots we'll have down always
# equal "good"
ref_set['Desc Rel'] = 1- ref_set['Desc Rel']
ref_set['SM Rel'] = 1- ref_set['SM Rel']

#%%
fig = plt.figure(figsize=(12,8))
plt.scatter(1-ref_set['Desc Rel'], ref_set['Desc RF'], c=ref_set['Desc NPC'] + ref_set['SM NPC'], cmap = 'viridis',
            vmin = 0, vmax=6.5*10**8)
#plt.scatter(1-Desc_best['Desc Rel'], Desc_best['Desc RF'], c=Desc_best['Desc NPC']+ ref_set['SM NPC'], cmap='viridis', edgecolors='k', vmin = 0, vmax=6.5*10**8, linewidth=1.5)
plt.scatter(1-ref_set[ref_set['group']=='Economical']['Desc Rel'], ref_set[ref_set['group']=='Economical']['Desc RF'],
            c=ref_set[ref_set['group']=='Economical']['Desc NPC'] + ref_set[ref_set['group']=='Economical']['SM NPC'],
            cmap = 'viridis', vmin = 0, vmax=6.5*10**8, marker = '^', edgecolors='k', s=100)

plt.scatter(1-ref_set[ref_set['group']=='Conservative']['Desc Rel'], ref_set[ref_set['group']=='Conservative']['Desc RF'],
            c=ref_set[ref_set['group']=='Conservative']['Desc NPC'] + ref_set[ref_set['group']=='Conservative']['SM NPC'],
            cmap = 'viridis', vmin = 0, vmax=6.5*10**8, marker = 's', edgecolors='k', s=75)
plt.scatter(1-ref_set[ref_set['group']=='Baseline']['Desc Rel'], ref_set[ref_set['group']=='Baseline']['Desc RF'],
            c=ref_set[ref_set['group']=='Baseline']['Desc NPC'] + ref_set[ref_set['group']=='Baseline']['SM NPC'],
            cmap = 'viridis', vmin = 0, vmax=6.5*10**8, marker = 'v', edgecolors='w', s=150)
plt.scatter(1-ref_set[ref_set['group']=='Transfers']['Desc Rel'], ref_set[ref_set['group']=='Transfers']['Desc RF'],
            c=ref_set[ref_set['group']=='Transfers']['Desc NPC'] + ref_set[ref_set['group']=='Transfers']['SM NPC'],
            cmap = 'viridis', vmin = 0, vmax=6.5*10**8, marker = 'P', edgecolors='k', s=100)

plt.xlabel('Reliability', fontsize=16)
plt.ylabel('Restriction Frequency', fontsize=16)
plt.xlim([0,1])
plt.ylim([0,1])
plt.title('Descoberto Trade-offs', fontsize=16)
plt.colorbar(label='Inf NPC (R\$)')
plt.savefig('Paper_Descoberto_objectives.png', bbox_inches='tight', dpi=300)

#%%
# Normalize the objectives
# we want all of the axes to line up, so we'll transform all values 
# to have a range of 0-1. We'll do this by subtracting the minimum dividing 
# by the difference of the max and min
normalized_objectives = ref_set.copy(deep=True)

# set up bounds for each objective
REL_bounds = [1, 0] # no normalization needed
RF_bounds = [1, 0] # no normalization needed
NPC_bounds = [1100000000, 0]
PFC_bounds = [.35, 0]
WCC_bounds = [.5, 0]

# NPC
normalized_objectives['Desc NPC'] = (normalized_objectives['Desc NPC'] - \
        NPC_bounds[1])/(NPC_bounds[0]- NPC_bounds[1])
normalized_objectives['SM NPC'] = (normalized_objectives['SM NPC'] - \
        NPC_bounds[1])/(NPC_bounds[0]- NPC_bounds[1])

# PFC
normalized_objectives['Desc PFC'] = (normalized_objectives['Desc PFC'] - \
        PFC_bounds[1])/(PFC_bounds[0]- PFC_bounds[1])
normalized_objectives['SM PFC'] = (normalized_objectives['SM PFC'] - \
        PFC_bounds[1])/(PFC_bounds[0]- PFC_bounds[1])

# WCC
normalized_objectives['Desc WCC'] = (normalized_objectives['Desc WCC'] - \
        WCC_bounds[1])/(WCC_bounds[0]- WCC_bounds[1])
normalized_objectives['SM WCC'] = (normalized_objectives['SM WCC'] - \
        WCC_bounds[1])/(WCC_bounds[0]- WCC_bounds[1])

normalized_objectives = normalized_objectives.sort_values('group',  ascending=False)
Desc_normalized_objectives = normalized_objectives[['Desc NPC', 'Desc Rel', 'Desc RF', 'Desc WCC', 'Desc PFC', 'group']]
SM_normalized_objectives = normalized_objectives[['SM NPC', 'SM Rel', 'SM RF', 'SM WCC', 'SM PFC', 'group']]
#%%
# plot the utility's performance space
# First Descoberto
fig = plt.figure(figsize=(12,10))
ax1 = fig.add_subplot(211)
parallel_coordinates(Desc_normalized_objectives[Desc_normalized_objectives['group']=='A'], 'group', ax=ax1, color=['#CFCFCF'], linewidth=3, alpha=.25) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',
parallel_coordinates(Desc_normalized_objectives[Desc_normalized_objectives['group']=='B'], 'group', ax=ax1, color=['#5BBCD6'], linewidth=3, alpha=.9) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',
parallel_coordinates(Desc_normalized_objectives[Desc_normalized_objectives['group']=='C'], 'group', ax=ax1, color=['#F2AD00'], linewidth=3, alpha=.9) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',

ax1.set_ylim([0,1])
ax1.set_xticklabels(['Inf NPC', 'Rel', 'RF', 'WCC', 'PFC'], fontsize=14)
ax1.get_yaxis().set_ticks([])
ax1.set_ylabel('Objective Value\n$\longleftarrow$ Direction of preference', labelpad=10, fontsize =16)
ax1.get_legend().remove()
ax1.set_title('Descoberto', fontsize = 16)

# Now Santa Maria
ax2 = fig.add_subplot(212)
# Below is the green for no brush
parallel_coordinates(SM_normalized_objectives[SM_normalized_objectives['group']=='A'], 'group', ax=ax2, color=['#CFCFCF'], linewidth=3, alpha=.25) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',
parallel_coordinates(SM_normalized_objectives[SM_normalized_objectives['group']=='B'], 'group', ax=ax2, color=['#5BBCD6'], linewidth=3, alpha=.9) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',
parallel_coordinates(SM_normalized_objectives[SM_normalized_objectives['group']=='C'], 'group', ax=ax2, color=['#F2AD00'], linewidth=3, alpha=.9) #, color= ['#F98400','#00A08A', '#FF0000', '#5BBCD6',

ax2.set_ylim([0,1])
ax2.set_xticklabels(['INF NPC', 'Rel', 'RF' , 'WCC', 'PFC'], fontsize=14)
ax2.set_ylabel('Objective Value\n$\longleftarrow$ Direction of preference', labelpad=10, fontsize =16)
ax2.get_yaxis().set_ticks([])
ax2.legend(loc='upper right')
ax2.get_legend().remove()
ax2.set_title('Santa Maria', fontsize = 16)
plt.savefig('Paper_objectives.png', bbox_inches='tight', dpi=300)
