//
// Created by AI Assistant on 10/13/25.
//

#ifndef TRIANGLEMODEL_WSSDATACOLLECTOR_H
#define TRIANGLEMODEL_WSSDATACOLLECTOR_H

#include "Base/DataCollector.h"
#include "../SystemComponents/Utility/WaterSupplySystems.h"
#include "../SystemComponents/Utility/Utility.h"

class WSSDataCollector : public DataCollector {
private:
    vector<double> combined_storage;
    vector<double> storage_capacity;
    vector<double> storage_to_capacity_ratio;
    vector<double> unrestricted_demand;
    vector<double> restricted_demand;
    vector<double> demand_multiplier;
    vector<double> demand_offset;
    vector<double> waste_water_discharge;
    vector<double> unfulfilled_demand;
    vector<double> net_stream_inflow;
    vector<double> total_treatment_capacity;
    vector<double> total_storage_treatment_capacity;
    vector<int> water_sources_count;
    vector<double> short_term_rof;
    vector<double> long_term_rof;
    
    // Affordability index (water price / average income)
    vector<double> weekly_affordability_index;
    vector<double> aggregated_current_price;
    vector<double> residential_current_price;
    
    // Financial data vectors (WSS-level)
    vector<double> gross_revenues;
    vector<double> drought_mitigation_cost;
    vector<double> contingency_fund_contribution;
    vector<double> debt_service_payments;
    vector<double> net_present_infrastructure_cost;
    vector<double> contingency_fund_size;
    double average_monthly_income = 0.0;
    double initial_households = 0.0;
    
    // Infrastructure pathways built by this WSS
    vector<vector<int>> pathways;
    
    // Per-source failure flag: 1 if ANY source in WSS failed that week, 0 otherwise
    vector<int> weekly_failure_flag;
    
    // Per-source failure tracking: source_id -> weekly failure flag (1=failed, 0=ok)
    map<int, vector<int>> per_source_failure_flag;
    // Per-source names for CSV headers (source_id -> name)
    map<int, string> source_names;
    
    const WaterSupplySystems *wss;
    const Utility *owner;  // Store owner directly to avoid accessing deleted WSS
    int owner_id;  // Store owner ID to avoid dangling pointer access

public:

    explicit WSSDataCollector(const WaterSupplySystems *wss, unsigned long realization);

    WSSDataCollector &operator=(const WSSDataCollector &wss_data_collector);

    string printTabularString(int week) override;

    string printCompactString(int week) override;

    void collect_data() override;

    string printTabularStringHeaderLine1() override;

    string printTabularStringHeaderLine2() override;

    string printCompactStringHeader() override;

    void checkForNans() const;

    // Getters - Operational data
    const vector<double> &getCombined_storage() const;
    const vector<double> &getStorage_capacity() const;
    const vector<double> &getStorage_to_capacity_ratio() const;
    const vector<double> &getUnrestricted_demand() const;
    const vector<double> &getRestricted_demand() const;
    const vector<double> &getDemand_multiplier() const;
    const vector<double> &getDemand_offset() const;
    const vector<double> &getWaste_water_discharge() const;
    const vector<double> &getUnfulfilled_demand() const;
    const vector<double> &getNet_stream_inflow() const;
    const vector<double> &getTotal_treatment_capacity() const;
    const vector<double> &getTotal_storage_capacity() const;
    const vector<int> &getWater_sources_count() const;
    const WaterSupplySystems *getWss() const;
    const vector<double> &getLong_term_rof() const;
    
    // Getters - Financial data
    const vector<double> &getGross_revenues() const;
    const vector<double> &getDrought_mitigation_cost() const;
    const vector<double> &getContingency_fund_contribution() const;
    const vector<double> &getDebt_service_payments() const;
    const vector<double> &getNet_present_infrastructure_cost() const;
    const vector<double> &getContingency_fund_size() const;
    const vector<double> &getWeekly_affordability_index() const;
    const vector<double> &getResidential_current_price() const;
    double getAverage_monthly_income() const;
    double getInitial_households() const;
    
    const Utility *getOwner() const;  // Direct access to owner (safe after WSS deletion)
    int getOwnerId() const;  // Get stored owner ID (safe even after WSS deletion)
    const vector<vector<int>> &getPathways() const;
    const vector<int> &getWeekly_failure_flag() const;
    const map<int, vector<int>> &getPer_source_failure_flag() const;
    const map<int, string> &getSource_names() const;
};

#endif //TRIANGLEMODEL_WSSDATACOLLECTOR_H