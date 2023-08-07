"""
Created on Fri Mar  6 15:53:05 2020

This script will create the visualization in Figure 3 of Gold et al. (2021)
"""


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import parallel_coordinates
import seaborn as sns
import matplotlib.cm as cm
#from plotting_functions_DG import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

#### plot decision variables on radial plots
dv_names = ['RT_D', 'RT_SM', 'RT_DD', 'RT_DSM', 'TT_D', 'TT_SM', 'AP_D', 'AP_SM',\
            'INF_D', 'INF_SM']

# load original solution (1-D array)
util = 'SM'
row_num = 3 
comp_name = "ModInf"
figname = 'C:/Users/lbl59/Desktop/CAESB/Figures/polar_' + util + '_' + comp_name + '.pdf'


dvs_directory = 'C:/Users/lbl59/Desktop/CAESB/'
dvs_filename = 'decision_vars.csv'

'''
Get dv filepaths then load the data
'''
dvs_filepath = dvs_directory + dvs_filename

dvs = np.loadtxt(dvs_filepath, delimiter=',', skiprows=1)
dvs_row = dvs[row_num, :]
# extract DVs for each utility (not including IP)
D_dvs = np.take(dvs_row,[0, 4, 6, 8])
SM_dvs = np.take(dvs_row,[1, 5, 7, 9])

# reverse direction of ROFs, since lower values mean increased use of the DV

# Descoberto RT
# Descoberto RT
D_dvs[0] = 1- D_dvs[0]

# Descoberto TT
D_dvs[1] = 1- D_dvs[1]

# Descoberto RC
D_dvs[2] = D_dvs[2]/.1

# Descoberto INF
D_dvs[3] = 1- D_dvs[3]

D_dvs = np.append(D_dvs, D_dvs[0])

# Santa Maria RT
# Santa Maria RT
SM_dvs[0] = 1- SM_dvs[0]

# Santa Maria TT
SM_dvs[1] = 1- SM_dvs[1]

# Santa Maria RC
SM_dvs[2] = SM_dvs[2]/.1

# Santa Maria INF
SM_dvs[3] = 1- SM_dvs[3]

SM_dvs = np.append(SM_dvs, SM_dvs[0])

c_best_cols = ['darkorange', 'purple', 'darkgreen']
c_worst_cols = ['gold', 'thistle', 'palegreen']
c_comp = ['darkorange', 'purple', 'darkgreen']


fig = plt.figure(figsize=(4,4))
plt.rcParams['axes.titlepad'] = 14
spec = fig.add_gridspec(ncols=1, nrows=1, figure=fig,
                        height_ratios = [1], \
                        width_ratios = [1])

SM_dvs_names = ['RT', 'TT', 'AP', 'INF']
ax3 = fig.add_subplot(spec[0,0], projection='polar')
theta3 = np.linspace(0,2*np.pi, len(SM_dvs))
ax3.plot(theta3, SM_dvs, color=c_comp[0], linewidth=4, alpha=0.6)
#ax3.plot(theta3, W_dvs_b, color=c_best, label='Best regional\nrobutness', linewidth=3)
lines, labels = ax3.set_thetagrids(range(0, 360, int(360/len(SM_dvs_names))), (SM_dvs_names))
ax3.set_yticklabels([])
ax3.set_title('Santa Maria', fontsize=12)

ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.suptitle(comp_name, size=14)
plt.tight_layout()
plt.savefig(figname, bbox_inches='tight')
plt.show()
