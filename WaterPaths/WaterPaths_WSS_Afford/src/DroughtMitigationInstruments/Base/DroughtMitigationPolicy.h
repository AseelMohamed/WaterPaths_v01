//
// Created by bernardo on 2/6/17.
//

#ifndef TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H
#define TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H

#include "../../SystemComponents/Utility/Utility.h"
#include "../../SystemComponents/Utility/WaterSupplySystems.h"
#include "../../Utils/Constants.h"
#include "../../Utils/Graph/Graph.h"
#include "../../Utils/Matrices.h"
#include "../../Controls/Base/MinEnvFlowControl.h"

class DroughtMitigationPolicy {
private:
    vector<Matrix2D<int>> *storage_to_rof_table_;

protected:
    DroughtMitigationPolicy(const DroughtMitigationPolicy &drought_mitigation_policy);

    vector<int> wss_ids;
    vector<WaterSupplySystems *> realization_wss;
    vector<WaterSource *> realization_water_sources;
    vector<MinEnvFlowControl *> realization_min_env_flow_controls;
    vector<vector<double>> *rdm_factors_all;
    double *rdm_factors_realization;
    bool use_imported_tables;

    double getRofFromRealizationTable(int wss_id, int week, int tier);

public:
    const int id;
    const int type;

    DroughtMitigationPolicy(const int id, const int type);

    virtual void applyPolicy(int week)= 0;

    virtual void addSystemComponents(vector<WaterSupplySystems *> wss,
                                         vector<WaterSource *> water_sources,
                                         vector<MinEnvFlowControl *> min_env_flow_controls)= 0;

    const vector<int> &getWSS_ids() const;

    bool operator<(const DroughtMitigationPolicy *other);

    bool operator>(const DroughtMitigationPolicy *other);

    virtual ~DroughtMitigationPolicy();

    void setStorage_to_rof_table_(vector<Matrix2D<int>> &storage_to_rof_table_, int use_imported_tables);

    virtual void setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                                vector<double> &water_sources_rdm, vector<double> &policy_rdm)= 0;

};


#endif //TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H
