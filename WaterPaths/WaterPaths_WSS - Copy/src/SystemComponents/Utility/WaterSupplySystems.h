#ifndef FDB_WSS_H
#define FDB_WSS_H


#include <vector>
#include <string>
#include <memory>
#include <map>
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
                      vector<vector<double>>& demands_all_realizations,
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

    // Copy constructor for proper deep copying
    WaterSupplySystems(const WaterSupplySystems& other);

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
    
    // Infrastructure decision-making (WSS-level triggers based on individual ROF)
    int infrastructureConstructionHandler(double long_term_rof, int week);
    
    // Force infrastructure construction when another WSS triggers shared infrastructure
    void forceInfrastructureConstruction(int week, const vector<int>& new_infra_triggered);

    // Infrastructure management accessors
    const vector<int>& getRof_infra_construction_order() const;
    
    // Method to set infrastructure parameters after WSS creation
    void setInfrastructureParameters(const vector<int>& rof_infra_construction_order,
                                    const vector<int>& demand_infra_construction_order, 
                                    const vector<double>& infra_construction_triggers);

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
    double getDemand_offset() const;

    void resetTotal_storage_capacity(); //checked
    
    // Setter methods for delegation
    void setRisk_of_failure(double risk_of_failure); //Checked
    void setTotal_available_volume(double volume); // For synchronization
    
    // Metadata accessors
    int getSystemId() const; //Checked
//     int getUtilityId() const { return utility_id; }

    Utility* getOwner() const { return owner; }
    const InfrastructureManager& getInfrastructure_construction_manager() const { return infrastructure_construction_manager; }
    const vector<WaterSource*>& getWater_sources() const; //Checked
    bool isUsedForRealization() const { return used_for_realization; }
    void setUsedForRealization(bool value) { used_for_realization = value; }
    unsigned long getCurrentRealizationId() const { return current_realization_id; }

    // Per-source failure thresholds for reliability calculation
    void setSourceFailureThreshold(int source_id, double threshold);
    double getSourceFailureThreshold(int source_id) const;
    const std::map<int, double>& getSourceFailureThresholds() const { return source_failure_thresholds; }

    // Operational drought response
    void setDemand_multiplier(double demand_multiplier); //checked
    void setDemand_offset(double demand_offset, double offset_rate_per_volume); //checked
    double getRestrictedDemand() const; //checked
    double getOffset_rate_per_volume() const; //NEW - for transfer cost calculations
    
    // Financial data storage (calculated by Utility, stored in WSS)
    void setWssGrossRevenue(double revenue);
    void setWssDroughtMitigationCost(double cost);
    void setWssContingencyFundShare(double share);
    void setWssDebtServiceShare(double share);
    void setWssInfrastructureNPC(double npc);
    void setWssResidentialPrice(double price);
    
    double getWssGrossRevenue() const;
    double getWssDroughtMitigationCost() const;
    double getWssContingencyFundShare() const;
    double getWssDebtServiceShare() const;
    double getWssInfrastructureNPC() const;
    double getWssResidentialPrice() const;
    
    void setAverageMonthlyIncome(double income);
    double getAverageMonthlyIncome() const;
    void setInitialHouseholds(double households);
    double getInitialHouseholds() const;
    
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
    double infra_discount_rate = 0.0;  // Realization-specific discount rate for bond NPC calculations
    unsigned long current_realization_id = 0;  // Track which realization this WSS copy belongs to
    double bond_term_multiplier = 1.0;  // Realization-specific bond term multiplier
    double bond_interest_rate_multiplier = 1.0;  // Realization-specific bond interest rate multiplier
    vector<double> available_treated_flow_rate;
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
    
    /// Financial data storage (set by Utility, no calculation logic in WSS)
    double wss_gross_revenue = 0.0;
    double wss_drought_mitigation_cost = 0.0;
    double wss_contingency_fund_share = 0.0;
    double wss_debt_service_share = 0.0;
    double wss_infrastructure_npc = 0.0;
    double wss_residential_price = 0.0;  // Current residential first-tier weekly price
    double average_monthly_income = 0.0;  // Average monthly income for affordability calculation
    double initial_households = 0.0;  // Initial households for affordability scaling
    
    // Per-source failure thresholds (source_id -> threshold ratio)
    // If a source is not in this map, it uses STORAGE_CAPACITY_RATIO_FAIL default
    std::map<int, double> source_failure_thresholds;
    
    bool hasTreatmentCapacity() const;
};

#endif // FDB_WSS_H