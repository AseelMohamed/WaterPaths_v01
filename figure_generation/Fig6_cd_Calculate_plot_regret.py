#%% Calculate Regret
''''
A script to calculate regret (type 2 as defined in Herman et al., 2015)

Regret is the deviation from the best solution in each state of the world

Following Herman et al., 2015, the 90th percentile deviation from the best solution in each SOW is used, though no
normalization is conducted
'''


import numpy as np

# First go through every solution and every SOW and store the best obj
num_sols = 91

#%% Calculate type I regret - deviatoin from base SOW
regrets_typeI = np.zeros((num_sols, 5))
for i in range(0,num_sols):
    base_RDM = np.loadtxt('RDM_by_sol/sol_' + str(i) + '.csv', delimiter=',')
    REL_regret = np.zeros(998)
    RF_regret = np.zeros(998)
    NPC_regret = np.zeros(998)
    PFC_regret = np.zeros(998)
    WCC_regret = np.zeros(998)

    current_sol = sol = np.loadtxt('RDM_by_sol/sol_' + str(i) + '.csv', delimiter=',')
    for j in range(0, 998):
        REL_regret[j] = base_RDM[481, 0] - current_sol[j,0]
        RF_regret[j] = current_sol[j, 1] - base_RDM[481, 1]
        NPC_regret[j] = current_sol[j, 2] - base_RDM[481, 2]
        PFC_regret[j] = current_sol[j, 3] - base_RDM[481, 3]
        WCC_regret[j] = current_sol[j, 4] - base_RDM[481, 4]

    regrets_typeI[i,:] = [np.percentile(REL_regret,90), np.percentile(RF_regret,90), np.percentile(NPC_regret, 90), np.percentile(PFC_regret,90), np.percentile(WCC_regret, 90)]

#%% Calculate type II regret - deviation from ideal
# initialize parameters to store best obj values
ideal_REL = np.ones(998)*0
ideal_RF = np.ones(998)*1
ideal_NPC = np.ones(998)*100000000000
ideal_PFC = np.ones(998)*1
ideal_WCC = np.ones(998)*1

# starting with this just for Descoberto
for i in range(0,num_sols):
    current_sol = np.loadtxt('RDM_by_sol/sol_' + str(i) + '.csv', delimiter=',')
    for j in range(0,998):
        if current_sol[j,0] > ideal_REL[j]:
            ideal_REL[j] = current_sol[j,0]
        if current_sol[j,1] < ideal_RF[j]:
            ideal_RF[j] = current_sol[j,1]
        if current_sol[j,2] < ideal_NPC[j]:
            ideal_NPC[j] = current_sol[j,2]
        if current_sol[j,3] < ideal_PFC[j]:
            ideal_PFC[j] = current_sol[j,3]
        if current_sol[j,4] < ideal_WCC[j]:
            ideal_WCC[j] = current_sol[j,4]

# initialize a matrix to store regret
regrets_typeII = np.zeros((num_sols, 5))

for i in range(0,num_sols):
    REL_regret = np.zeros(998)
    RF_regret = np.zeros(998)
    NPC_regret = np.zeros(998)
    PFC_regret = np.zeros(998)
    WCC_regret = np.zeros(998)

    current_sol = sol = np.loadtxt('RDM_by_sol/sol_' + str(i) + '.csv', delimiter=',')
    for j in range(0, 998):
        REL_regret[j] = ideal_REL[j] - current_sol[j,0]
        RF_regret[j] = current_sol[j, 1] - ideal_RF[j]
        NPC_regret[j] = current_sol[j, 2] - ideal_NPC[j]
        PFC_regret[j] = current_sol[j, 3] - ideal_PFC[j]
        WCC_regret[j] = current_sol[j, 4] - ideal_WCC[j]

    regrets_typeII[i,:] = [np.percentile(REL_regret,90), np.percentile(RF_regret,90), np.percentile(NPC_regret, 90), np.percentile(PFC_regret,90), np.percentile(WCC_regret, 90)]


#%% Plot regret type I
from matplotlib import pyplot as plt
best_sols = [3,4,24,25,26,30,39,41,60,76,82]

fig = plt.figure()

plt.scatter(regrets_typeI[:,0], regrets_typeI[:,1], c=regrets_typeI[:,2], cmap='viridis', vmin=0, vmax = 250000000)

plt.scatter(regrets_typeI[3,0], regrets_typeI[3,1], c = regrets_typeI[3,2], cmap = 'viridis', marker ='^',
            edgecolors='k', s=100, vmin=0, vmax = 250000000)
plt.scatter(regrets_typeI[41,0], regrets_typeI[41,1], c = regrets_typeI[41,2], cmap = 'viridis', marker ='s',
            edgecolors='k', s=100, vmin=0, vmax = 250000000)
plt.scatter(regrets_typeI[85,0], regrets_typeI[85,1], c = regrets_typeI[85,2], cmap = 'viridis', marker ='v',
            edgecolors='w', s=100, vmin=0, vmax = 250000000)
plt.scatter(regrets_typeI[20,0], regrets_typeI[20,1], c = regrets_typeI[20,2], cmap = 'viridis', marker ='P',
            edgecolors='w', s=150, vmin=0, vmax = 250000000)

plt.xlabel('Reliability Regret \n $\longleftarrow$ Preference ')
plt.ylabel('Restriction Frequency Regret \n $\longleftarrow$ Preference ')
plt.colorbar(label='Inf NPC Regret (R\$) \n $\longleftarrow$ Preference ')
plt.title('Type I Regret')

plt.tight_layout()
plt.savefig('Regret_typeI_tradeoffs.png', dpi=500)

#%% Plot regret type II
fig = plt.figure()

plt.scatter(regrets_typeII[:,0], regrets_typeII[:,1], c=regrets_typeII[:,2], cmap='viridis',
            vmin = 0, vmax=6.5*10**8)

plt.scatter(regrets_typeII[3,0], regrets_typeII[3,1], c = regrets_typeII[3,2], cmap = 'viridis', marker ='^',
            edgecolors='k', s=100,vmin = 0, vmax=6.5*10**8)
plt.scatter(regrets_typeII[41,0], regrets_typeII[41,1], c = regrets_typeII[41,2], cmap = 'viridis', marker ='s',
            edgecolors='k', s=100,vmin = 0, vmax=6.5*10**8)
plt.scatter(regrets_typeII[85,0], regrets_typeII[85,1], c = regrets_typeII[85,2], cmap = 'viridis', marker ='v',
            edgecolors='w', s=100,vmin = 0, vmax=6.5*10**8)
plt.scatter(regrets_typeII[20,0], regrets_typeII[20,1], c = regrets_typeII[20,2], cmap = 'viridis', marker ='P',
            edgecolors='w', s=150,vmin = 0, vmax=6.5*10**8)

plt.xlabel('Reliability Regret \n $\longleftarrow$ Preference ')
plt.ylabel('Restriction Frequency Regret \n $\longleftarrow$ Preference ')
plt.xlim([0,1])
plt.ylim([0,1])
plt.colorbar(label='Inf NPC Regret (R\$) \n $\longleftarrow$ Preference ')
plt.title('Type II Regret')

plt.tight_layout()
plt.savefig('Regret_typeII_tradeoffs.png', dpi=500)