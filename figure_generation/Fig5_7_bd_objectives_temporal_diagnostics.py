#%%
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

def calculate_sol_reliability(sol_num, RDM):
    fail_years = np.zeros([1000, 40])

    for rel in range(1000):
        cur_rel  = np.loadtxt('DiagnosticsPathways/RDM' + str(RDM) + '/utility_files/sol_' + str(sol_num) + '/Utilities_s' +
                              str(sol_num) + '_RDM' + str(RDM) + '_r' + str(rel) + '.csv', delimiter=',', skiprows=1)

        for year, i in enumerate(np.arange(0,2080, 52)):
            for w in range(52):
                if cur_rel[i+w, 0]/cur_rel[i+w, 1] < .2:
                    fail_years[rel, year] = 1

    # calculate reliability over time
    reliability = 1 - np.sum(fail_years, axis=0)/1000

    return fail_years, reliability

def calculate_annual_restriction_frequency(sol_num, RDM):
    restr_years = np.zeros([1000, 40])

    for rel in range(1000):
        cur_rel  = np.loadtxt('DiagnosticsPathways/RDM' + str(RDM) + '/policy_files/sol_' + str(sol_num) + '/Policies_s' +
                              str(sol_num) + '_RDM' + str(RDM) + '_r' + str(rel) + '.csv', delimiter=',', skiprows=1)

        for year, i in enumerate(np.arange(0,2080, 52)):
            for w in range(52):
                if cur_rel[i+w, 0] < 1:
                    restr_years[rel, year] = 1

    # calculate rf over time
    restriction_frequency =np.sum(restr_years, axis=0)/1000

    return restr_years, restriction_frequency

def calculate_ave_annual_transfers(sol_num, RDM):
    annual_transfers = np.zeros([1000, 40])

    for rel in range(1000):
        cur_rel = np.loadtxt('DiagnosticsPathways/RDM' + str(RDM) + '/policy_files/sol_' + str(sol_num) + '/Policies_s' +
                             str(sol_num) + '_RDM' + str(RDM) + '_r' + str(rel) + '.csv', delimiter=',', skiprows=1)

        for year, i in enumerate(np.arange(0, 2080, 52)):
            annual_transfers[rel, year] = np.sum(cur_rel[i:i+51, 2])

    # calculate average tranfser volume over time
    average_transfers = np.mean(annual_transfers, axis=0)

    return annual_transfers, average_transfers

def calculate_transfer_freq(sol_num, RDM):
    transfer_years = np.zeros([1000, 40])

    for rel in range(1000):
        cur_rel = np.loadtxt('DiagnosticsPathways/RDM' + str(RDM) + '/policy_files/sol_' + str(sol_num) + '/Policies_s' +
                             str(sol_num) + '_RDM' + str(RDM) + '_r' + str(rel) + '.csv', delimiter=',', skiprows=1)

        for year, i in enumerate(np.arange(0, 2080, 52)):
            for w in range(52):
                if cur_rel[i + w, 2] > 0:
                    transfer_years[rel, year] = 1

    # calculate transfer freq over time
    transfer_frequency =np.sum(transfer_years, axis=0)/1000

    return transfer_years, transfer_frequency

# Calculate Debt Service
def calculate_debt_service(sol_num, RDM):
    annual_debt_service = np.zeros([1000, 40])

    for rel in range(1000):
        cur_rel = np.loadtxt('DiagnosticsPathways/RDM' + str(RDM) + '/utility_files/sol_' + str(sol_num) + '/Utilities_s' +
                             str(sol_num) + '_RDM' + str(RDM) + '_r' + str(rel) + '.csv', delimiter=',', skiprows=1)

        for year, i in enumerate(np.arange(0, 2080, 52)):
            annual_debt_service[rel, year] = np.sum(cur_rel[i:i + 51, 14]) + np.sum(cur_rel[i:i + 51, 29])

    # calculate average debt service over time
    average_debt_service = np.mean(annual_debt_service, axis=0) / 1000

    return annual_debt_service, average_debt_service


def plotTemporalDiagnostics(sol_3_stats, sol_41_stats, sol_20_stats, sol_85_stats, RDM):
    '''
    Plots performance of four policies over time

    :param sol_3_stats:         list, performance of solution 3 - [rel, rf, trans, NPC]
    :param sol_41_stats:        list, performance of solution 41 - [rel, rf, trans, NPC]
    :param sol_20_stats:        list, performance of solution 20 - [rel, rf, trans, NPC]
    :param sol_85_stats:        list, performance of solution 85 - [rel, rf, trans, NPC]
    :param RDM:                 float, RDM number

    '''

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (8, 8))

    # plot reliability
    axes.flatten()[0].step(np.arange(40), sol_3_stats[0], linewidth=2)
    axes.flatten()[0].step(np.arange(40), sol_41_stats[0], linestyle='--', linewidth=2)
    axes.flatten()[0].step(np.arange(40), sol_20_stats[0], linestyle=':', linewidth=2)
    axes.flatten()[0].step(np.arange(40), sol_85_stats[0], linestyle='-.', linewidth=2)
    axes.flatten()[0].set_ylim([.5, 1.0])
    axes.flatten()[0].set_xlim([0,40])
    axes.flatten()[0].set_ylabel('Reliability')
    axes.flatten()[0].set_title('Reliability')

    # add a zoom in on reliability, make a 3x3 set of subplots and have the bars encompass two panels?

    axes.flatten()[1].step(np.arange(40), sol_3_stats[1], linewidth=2)
    axes.flatten()[1].step(np.arange(40), sol_41_stats[1], linestyle='--', linewidth=2)
    axes.flatten()[1].step(np.arange(40), sol_20_stats[1], linestyle=':', linewidth=2)
    axes.flatten()[1].step(np.arange(40), sol_85_stats[1], linestyle='-.', linewidth=2)
    axes.flatten()[1].set_ylim([-0.01, 1.01])
    axes.flatten()[1].set_xlim([0,40])
    axes.flatten()[1].legend(['Moderate Inf. Comp.', 'Heavy Inf. Comp.', 'Transfer Dominant',  'Baseline'], loc='right')
    axes.flatten()[1].set_ylabel('Restriction Frequency')
    axes.flatten()[1].set_title('Restriction Frequency')

    axes.flatten()[2].step(np.arange(40), sol_3_stats[2], linewidth=2)
    axes.flatten()[2].step(np.arange(40), sol_41_stats[2], linestyle='--', linewidth=2)
    axes.flatten()[2].step(np.arange(40), sol_20_stats[2], linestyle=':', linewidth=2)
    axes.flatten()[2].step(np.arange(40), sol_85_stats[2], linestyle='-.', linewidth=2)
    #axes.flatten()[2].set_ylim([0, 1.01])
    axes.flatten()[2].set_xlim([0,40])
    axes.flatten()[2].set_ylabel('Average Annual Transfer Volume ($hm^3$)')
    axes.flatten()[2].set_xlabel('Simulation Year')
    axes.flatten()[2].set_title('Transfers')

    axes.flatten()[3].step(np.arange(0, 40), sol_3_stats[3], linewidth=2)
    axes.flatten()[3].step(np.arange(0, 40), sol_41_stats[3], linestyle='--', linewidth=2, alpha =.7)
    axes.flatten()[3].step(np.arange(0, 40), sol_20_stats[3], linestyle=':', linewidth=2, alpha=.7)
    axes.flatten()[3].set_xlim([0,40])
    axes.flatten()[3].set_ylim([0,80000])
    axes.flatten()[3].set_ylabel('Annual Debt Service ($R)')
    axes.flatten()[3].set_title('Debt Service Payments')
    axes.flatten()[3].set_xlabel('Simulation Year')

    plt.tight_layout()
    #plt.show()
    plt.savefig('Obj_Diagnostics_' + str(RDM) + '.pdf', format='pdf')




#%% RDM 481
RDM = 481
fail_years_3_481, sol_3_rel_481 = calculate_sol_reliability(3, RDM)
fail_years_41_481, sol_41_rel_481 = calculate_sol_reliability(41, RDM)
fail_years_20_481, sol_20_rel_481 = calculate_sol_reliability(20, RDM)
fail_years_85_481, sol_85_rel_481 = calculate_sol_reliability(85, RDM)
print('Rel Done')

annual_transfers_3_481, sol_3_trans_481 = calculate_ave_annual_transfers(3, RDM)
annual_transfers_41_481, sol_41_trans_481 = calculate_ave_annual_transfers(41, RDM)
annual_transfers_20_481, sol_20_trans_481 = calculate_ave_annual_transfers(20, RDM)
annual_transfers_85_481, sol_85_trans_481 = calculate_ave_annual_transfers(85, RDM)
print('Trans Done')

rest_years_3_481, sol_3_rf_481 = calculate_annual_restriction_frequency(3, RDM)
rest_years_41_481, sol_41_rf_481 = calculate_annual_restriction_frequency(41, RDM)
rest_years_20_481, sol_20_rf_481 = calculate_annual_restriction_frequency(20, RDM)
rest_years_85_481, sol_85_rf_481 = calculate_annual_restriction_frequency(85, RDM)
print('Rest Done')

annuLdebt_3_481, sol_3_ds_481 = calculate_debt_service(3, RDM)
annuLdebt_41_481, sol_41_ds_481 = calculate_debt_service(41, RDM)
annuLdebt_20_481, sol_20_ds_481 = calculate_debt_service(20, RDM)
annuLdebt_85_481, sol_85_ds_481 = calculate_debt_service(85, RDM)
print('NPC Done')

plotTemporalDiagnostics([sol_3_rel_481, sol_3_rf_481, sol_3_trans_481, sol_3_ds_481], [sol_41_rel_481, sol_41_rf_481,
                        sol_41_trans_481, sol_41_ds_481], [sol_20_rel_481, sol_20_rf_481, sol_20_trans_481,
                                                           sol_20_ds_481], [sol_85_rel_481, sol_85_rf_481,
                                                            sol_85_trans_481, sol_85_ds_481], 481)

print('Plotting Done')
