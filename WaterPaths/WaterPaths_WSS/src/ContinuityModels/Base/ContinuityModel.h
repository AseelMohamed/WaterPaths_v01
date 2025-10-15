//
// Created by bernardo on 1/12/17.
//

#ifndef TRIANGLEMODEL_CONTINUITYMODEL_H
#define TRIANGLEMODEL_CONTINUITYMODEL_H

#include "../../SystemComponents/WaterSources/Base/WaterSource.h"
#include "../../Utils/Constants.h"
#include "../../SystemComponents/Utility/Utility.h"
#include "../../SystemComponents/Utility/WaterSupplySystems.h"
#include "../../Utils/Graph/Graph.h"
#include "../../Controls/Base/MinEnvFlowControl.h"
#include <vector>

using namespace Constants;
using namespace std;

class ContinuityModel {

protected:

    vector<WaterSource *> continuity_water_sources;
    // vector<Utility *> continuity_wss;
    vector<WaterSupplySystems *> continuity_wss;
    vector<MinEnvFlowControl *> min_env_flow_controls;
    Graph water_sources_graph;
    vector<vector<int> > water_sources_to_wss;
    vector<vector<int> > water_sources_online_to_wss;
    vector<vector<int> > wss_to_water_sources;
    vector<int> downstream_sources;
    const vector<int> sources_topological_order;
    vector<double> water_sources_capacities;
    vector<double> wss_capacities;
    vector<vector<double>> demands;
    vector<double> wss_rdm;
    vector<double> water_sources_rdm;
    const int n_wss;
    const int n_sources;
    int delta_realization_weeks[NUMBER_REALIZATIONS_ROF + 1];
//    int delta_realization_weeks[NUMBER_REALIZATIONS_ROF];

public:
    const unsigned long realization_id;

    ContinuityModel(vector<WaterSource *> &water_sources, vector<WaterSupplySystems *> &wss,
                    vector<MinEnvFlowControl *> &min_env_flow_controls,
                    const Graph &water_sources_graph,
                    const vector<vector<int>> &water_sources_to_wss,
                    vector<double> &wss_rdm,
                    vector<double> &water_sources_rdm,
                    unsigned long realization_id);

    ContinuityModel(ContinuityModel &continuity_model);

    void continuityStep(
            int week, int rof_realization = -1, bool
    apply_demand_buffer = false);

    virtual ~ContinuityModel();

    void setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                        vector<double> &water_sources_rdm);

    vector<int> getOnlineDownstreamSources();

    const vector<WaterSource *> &getContinuity_water_sources() const;

    const vector<WaterSupplySystems *> &getContinuity_wss() const;
};


#endif //TRIANGLEMODEL_REGION_H
