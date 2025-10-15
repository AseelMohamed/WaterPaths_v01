//
// Created by bernardo on 1/26/17.
//

#ifndef TRIANGLEMODEL_CONTINUITYMODELROF_H
#define TRIANGLEMODEL_CONTINUITYMODELROF_H



#include "Base/ContinuityModel.h"
#include "../Utils/Matrices.h"


class ContinuityModelROF : public ContinuityModel {
private:
    vector<int> online_downstream_sources;
    bool *storage_wout_downstream;
    const int n_topo_sources;
    const int use_precomputed_rof_tables;

protected:
    int beginning_tier = 0;
    vector<WaterSource *> realization_water_sources;
    vector<WaterSupplySystems *> realization_wss;
    vector<Matrix2D<int>> ut_storage_to_rof_table;

    vector<vector<double>> table_storage_shift;
    vector<double> wss_base_storage_capacity;
    vector<double> wss_base_delta_capacity_table;
    vector<double> current_and_base_storage_capacity_ratio;
    vector<double> current_storage_table_shift;

public:
    ContinuityModelROF(vector<WaterSource *> water_sources, const Graph &water_sources_graph,
                       const vector<vector<int>> &water_sources_to_wss, vector<WaterSupplySystems *> wss,
                       vector<MinEnvFlowControl *> min_env_flow_controls, vector<double>& wss_rdm,
                       vector<double>& water_sources_rdm, unsigned long total_weeks_simulation,
                       const int use_precomputed_rof_tables, const unsigned long realization_id);

//    ContinuityModelROF(ContinuityModelROF &continuity_model_rof);

    vector<double> calculateShortTermROF(int week, int import_export_rof_tables);

    vector<double> calculateShortTermROFFullCalcs(int week);

    vector<double> calculateShortTermROFTable(int week);

    vector<double> calculateLongTermROF(int week);

    // WSS-level methods - return vector<vector<double>> where outer vector is wss, inner vector is WSS within each utility
    vector<vector<double>> calculateShortTermROF_WSS(int week, int import_export_rof_tables);

    vector<vector<double>> calculateShortTermROFFullCalcs_WSS(int week);

    vector<vector<double>> calculateShortTermROFTable_WSS(int week);

    vector<vector<double>> calculateLongTermROF_WSS(int week);

    void resetWSSAndReservoirs(int rof_type);

    void connectRealizationWaterSources(const vector<WaterSource *> &realization_water_sources);

    void connectRealizationWSS(const vector<WaterSupplySystems *> &realization_wss);

    virtual void updateOnlineInfrastructure(int week);

    virtual ~ContinuityModelROF();

    void updateStorageToROFTable(double storage_percent_decrement,
                                 int week,
                                 const double *to_full_toposort, int rof_realization_number);

    vector<Matrix2D<int>> &getUt_storage_to_rof_table();

    void shiftStorages(double *available_volumes_shifted, const double
    *delta_storage);

    void printROFTable(const string &folder);

    void setROFTablesAndShifts(const vector<Matrix2D<int>> &storage_to_rof_table,
                               const vector<vector<double>> &table_storage_shift);

    void calculateEmptyVolumes(vector<WaterSource *> &realization_water_sources, double *to_full);

};


#endif //TRIANGLEMODEL_CONTINUITYMODELROF_H
