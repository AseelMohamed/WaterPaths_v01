//
// Created by bernardoct on 5/1/17.
//

#ifndef TRIANGLEMODEL_INSURANCESTORAGETOROF_H
#define TRIANGLEMODEL_INSURANCESTORAGETOROF_H


#include "Base/DroughtMitigationPolicy.h"
#include "../ContinuityModels/ContinuityModelROF.h"
#include "../ContinuityModels/ContinuityModelRealization.h"

class InsuranceStorageToROF : public DroughtMitigationPolicy,
                              public ContinuityModelROF {
private:
    vector<double> rof_triggers;// = vector<double>(40, 0.);
    const unsigned long total_simulation_time;
    const double insurance_premium;
    vector<double> payout_multiplier;
    vector<double> insurance_price;
    const vector<double> &fixed_payouts;
    vector<double> utilities_revenue_update;
    vector<double> utilities_revenue_last_year;
    vector<DroughtMitigationPolicy *> drought_mitigation_policies;

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
                              unsigned long total_simulation_time);

    InsuranceStorageToROF(InsuranceStorageToROF &insurance);

    ~InsuranceStorageToROF() override;

    void priceInsurance(int week);

    void applyPolicy(int week) override;

    void addSystemComponents(vector<WaterSupplySystems *> wss,
                                 vector<WaterSource *> water_sources,
                                 vector<MinEnvFlowControl *> min_env_flow_controls) override;

    void setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                        vector<double> &water_sources_rdm, vector<double> &policy_rdm) override;

    vector<double> calculateShortTermROFTable(int week, const vector<WaterSupplySystems *> &wss, const int &n_wss);

    void updateOnlineInfrastructure(int week) override;
};


#endif //TRIANGLEMODEL_INSURANCESTORAGETOROF_H
