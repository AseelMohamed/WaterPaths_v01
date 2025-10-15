#ifndef FDB_WSS_H
#define FDB_WSS_H


#include <vector>
#include <string>
#include <memory>
#include "../WaterSources/Reservoir.h"
#include "../../Utils/Constants.h"
#include "../../Controls/WwtpDischargeRule.h"
#include "InfrastructureManager.h"

using namespace std;

class Utility; // forward declaration

class WaterSupplySystems {
public:
    const int system_id;     // Unique system identifier
    const int utility_id;    // Owner utility identifier  
    const double demand_buffer;
    const int number_of_week_demands;
    const char *name;


    // Constructor & Destructor
    WaterSupplySystems(const string& name,
                      int system_id,
                      int utility_id, 
                      Utility* owner_utility,
                      const WwtpDischargeRule& wwtp_rule);

    WaterSupplySystems(const string& name,
                      int system_id,
                      int utility_id,
                      Utility* owner_utility,
                      vector<vector<double>> &demands_all_realizations,
                      int number_of_week_demands,
                      WwtpDischargeRule wwtp_discharge_rule,
                      double demand_buffer,
                      vector<vector<int>> water_source_to_wtp,
                      vector<double> wss_owned_wtp_capacities);

    // Constructor for when there is infrastructure to be built
    WaterSupplySystems(const string& name,
                      int system_id,
                      int utility_id,
                      Utility* owner_utility,
                      vector<vector<double>> &demands_all_realizations,
                      int number_of_week_demands,
                      WwtpDischargeRule wwtp_discharge_rule,
                      double demand_buffer,
                      vector<vector<int>> water_source_to_wtp,
                      vector<double> wss_owned_wtp_capacities,
                      const vector<int> &rof_infra_construction_order,
                      const vector<int> &demand_infra_construction_order,
                      const vector<double> &infra_construction_triggers);

    ~WaterSupplySystems();
    
    // Water source management
    void clearWaterSources();
    void addWaterSource(WaterSource* water_source); //Checked
    void setWaterSourceOnline(unsigned int source_id, int week);
    bool hasTreatmentConnected(int ws);
    void checkErrorsAddWaterSourceOnline(WaterSource* water_source); //Checked
//     void copyWaterSourceConnections(const WaterSupplySystems& other);

    // Operational calculations
    void updateTreatmentAndNumberOfStorageSources();
    void updateTotalAvailableVolume(); //Checked
    void calculateWastewater_releases(int week, double* discharges); //Checked
//     void resetOperationalVariables();
    
    // Realization management
    void setRealization(unsigned long r, vector<double> &rdm_factors); //Checked
    
    // Demand management
    void splitDemands(int week, 
                     vector<vector<double>>& demands,
                     bool apply_demand_buffer = false); //Checked

    static vector<double> calculateWeeklyPeakingFactor(
                        vector<double> *demands); //Checked

    // Infrastructure mapping
    void unrollWaterSourceToWtpVector(
            const vector<vector<int>>& water_source_to_wtp,
            const vector<double>& wss_owned_wtp_capacities);
    void reconnectInfrastructureManager();
    int infrastructureConstructionHandler(double long_term_rof, int week); //checked

    // Demand splitting algorithms
    static bool idealDemandSplitUnconstrained(
            double* split_demands,
            const double* available_treated_flow_rate,
            double total_demand,
            const double* storage,
            double total_storage,
            int n_storage_sources);

    static bool idealDemandSplitConstrained(
            double* split_demands,
            bool* over_allocated,
            bool* has_spare_capacity,
            const double* available_treated_flow_rate,
            double total_demand,
            const double* storage,
            double total_storage,
            int n_storage_sources);

    // State accessors
    double getTotal_available_volume() const;
    double getTotal_stored_volume() const;
    double getTotal_storage_capacity() const;
    double getTotal_treatment_capacity() const; //checked
    double getWaste_water_discharge() const;
    double getUnfulfilled_demand() const;
    double getNet_stream_inflow() const;
    double getShort_term_risk_of_failure() const;
    double getLong_term_risk_of_failure() const;
    double getRisk_of_failure() const; //checked
    double getStorageToCapacityRatio() const; //checked
    double getTotal_storage_treatment_capacity() const;
    double getUnrestrictedDemand(int week = -1) const; //checked
    double getDemand_multiplier() const; //checked  
    double getDemand_offset() const; //checked
    void resetTotal_storage_capacity(); //checked
    
    // Setter methods for delegation
    void setRisk_of_failure(double risk_of_failure); //Checked
    
    // Metadata accessors
    int getSystemId() const; //Checked
//     int getUtilityId() const { return utility_id; }

    Utility* getOwner() const { return owner; }
    const vector<WaterSource*>& getWater_sources() const; //Checked

    // Operational drought response
    void setDemand_multiplier(double demand_multiplier); //checked
    void setDemand_offset(double demand_offset, double offset_rate_per_volume); //checked
    double getRestrictedDemand() const; //checked
    

private:
    Utility* owner;
    vector<int> priority_draw_water_source;
    vector<int> non_priority_draw_water_source;
    vector<double> weekly_peaking_factor;
    double short_term_risk_of_failure = 0;
    double long_term_risk_of_failure = 0;
    double total_storage_capacity = 0;
    double total_available_volume = 0;
    double total_stored_volume = 0;
    double total_treatment_capacity = 0;
    double total_storage_treatment_capacity = 0;
    double waste_water_discharge = 0;
    double unfulfilled_demand = 0.0;
    double net_stream_inflow = 0.0;
    double *available_treated_flow_rate = new double[0];
    bool used_for_realization = true;
    unsigned short n_storage_sources = 0;
    vector<WaterSource*> water_sources;
    WwtpDischargeRule wwtp_discharge_rule;
    vector<vector<double>> &demands_all_realizations;
    vector<double> demand_series_realization;
    vector<double> wss_owned_wtp_capacities;
    vector<int> water_source_to_wtp;
    InfrastructureManager infrastructure_construction_manager;


    /// Drought mitigation variables
    double demand_multiplier = 1;
    double demand_offset = 0;
    double offset_rate_per_volume = 0;
    double restricted_demand = 0;
    double unrestricted_demand = 0;
    int n_sources = 0;
    double max_capacity = 0;
    
    bool hasTreatmentCapacity() const;
};

#endif // FDB_WSS_H