//
// Created by bernardoct on 5/1/17.
//

#ifndef TRIANGLEMODEL_INSURANCESTORAGETOROF_H
#define TRIANGLEMODEL_INSURANCESTORAGETOROF_H

#include <string>
#include "Base/DroughtMitigationPolicy.h"
#include "../ContinuityModels/ContinuityModelROF.h"
#include "../ContinuityModels/ContinuityModelRealization.h"

class InsuranceStorageToROF : public DroughtMitigationPolicy,
                              public ContinuityModelROF {
private:
    vector<double> rof_triggers;
    const unsigned long total_simulation_time;
    const double insurance_premium;
    vector<double> payout_multiplier;
    vector<double> insurance_price;
    const vector<double> &fixed_payouts;
    vector<double> utilities_revenue_update;
    vector<double> utilities_revenue_last_year;
    vector<DroughtMitigationPolicy *> drought_mitigation_policies;

    // --- Static lookup table pricing ---
    // Four infrastructure states:
    //   0 = no infra built for either WSS
    //   1 = infra built for Descoberto only
    //   2 = infra built for TortoSM only
    //   3 = infra built for both WSSes
    //
    // Each table row: [rof_threshold, fail_freq_period_0, fail_freq_period_1, ...]
    //   fail_freq is the COMBINED expected trigger-weeks/year from BOTH WSSes
    //   (single utility perspective). Period width = 40 / n_time_periods years.
    vector<vector<vector<double>>> rof_freq_tables; // [4 states][n_rows][1 + n_time_periods]
    vector<string> rof_freq_table_files;            // stored for copy constructor
    int n_time_periods  = 0;                        // inferred from CSV column count
    int years_per_period = 0;                       // = 40 / n_time_periods
    vector<double> initial_treatment_capacity;      // per wss_ids order, set at addSystemComponents
    vector<double> initial_storage_capacity;        // per wss_ids order, set at addSystemComponents

    static vector<vector<double>> loadROFFreqTable(const string &filepath);
    int getInfraState() const;
    double lookupFailFreq(int state, int time_period_idx, double rof_trigger) const;

public:

    InsuranceStorageToROF(const int id, vector<WaterSource *> &water_sources,
                              const Graph &water_sources_graph,
                              const vector<vector<int>> &water_sources_to_wss,
                              vector<WaterSupplySystems *> &wss,
                              vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
                              vector<MinEnvFlowControl *> min_env_flow_controls,
                              vector<vector<double>>& wss_rdm,
                              vector<vector<double>>& water_sources_rdm,
                              vector<vector<double>>& policy_rdm, vector<double> &rof_triggers,
                              const double insurance_premium, const vector<double> &fixed_payouts,
                              unsigned long total_simulation_time,
                              const vector<string> &rof_freq_table_files);

    InsuranceStorageToROF(InsuranceStorageToROF &insurance);

    ~InsuranceStorageToROF() override;

    void priceInsurance(int week);

    void applyPolicy(int week) override;

    void addSystemComponents(vector<WaterSupplySystems *> wss,
                                 vector<WaterSource *> water_sources,
                                 vector<MinEnvFlowControl *> min_env_flow_controls) override;

    void setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                        vector<double> &water_sources_rdm, vector<double> &policy_rdm) override;

    vector<double> calculateShortTermROFTable(int week);

    void updateOnlineInfrastructure(int week) override;
};


#endif //TRIANGLEMODEL_INSURANCESTORAGETOROF_H
