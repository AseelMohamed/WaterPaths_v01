//
// Created by AI Assistant on 10/13/25.
//

#include <iomanip>
#include <sstream>
#include <algorithm>
#include "WSSDataCollector.h"

WSSDataCollector::WSSDataCollector(const WaterSupplySystems *wss, unsigned long realization)
        : DataCollector(wss->getSystemId(), wss->name, realization, UTILITY, 11 * COLUMN_WIDTH),
          wss(wss),
          owner(wss->getOwner()),
          owner_id(wss->getOwner() ? wss->getOwner()->id : -1) {  // Store owner ID immediately while WSS is still valid
}

int WSSDataCollector::getOwnerId() const {
    return owner_id;
}

string WSSDataCollector::printTabularString(int week) {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << combined_storage[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << storage_capacity[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << storage_to_capacity_ratio[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << net_stream_inflow[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << unfulfilled_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << restricted_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << unrestricted_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << demand_multiplier[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << waste_water_discharge[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << total_treatment_capacity[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << total_storage_treatment_capacity[week];

    return outStream.str();
}

string WSSDataCollector::printCompactString(int week) {

    stringstream outStream;

    outStream << combined_storage[week]
              << ","
              << storage_capacity[week]
              << ","
              << net_stream_inflow[week]
              << ","
              << short_term_rof[week]  
              << ","
              << long_term_rof[week]  
              << ","
              << restricted_demand[week]
              << ","
              << unrestricted_demand[week]
              << ","
              << unfulfilled_demand[week]
              << ","
              << waste_water_discharge[week]
              << ","
              << total_treatment_capacity[week]
              << ",";

    return outStream.str();
}

string WSSDataCollector::printTabularStringHeaderLine1() {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << "Stored"
              << setw(COLUMN_WIDTH) << "Storage"
              << setw(COLUMN_WIDTH) << "Storage"
              << setw(COLUMN_WIDTH) << "Net"
              << setw(COLUMN_WIDTH) << "Rest."
              << setw(COLUMN_WIDTH) << "Unrest."
              << setw(COLUMN_WIDTH) << "Unfulf."
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "W. Water"
              << setw(COLUMN_WIDTH) << "Treatment"
              << setw(COLUMN_WIDTH) << "Storage Treat.";

    return outStream.str();
}

string WSSDataCollector::printTabularStringHeaderLine2() {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << "Volume"
              << setw(COLUMN_WIDTH) << "Capacity"
              << setw(COLUMN_WIDTH) << "Ratio"
              << setw(COLUMN_WIDTH) << "Inflow"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Multiplier"
              << setw(COLUMN_WIDTH) << "Discharge"
              << setw(COLUMN_WIDTH) << "Capacity"
              << setw(COLUMN_WIDTH) << "Capacity";

    return outStream.str();
}

string WSSDataCollector::printCompactStringHeader() {
    stringstream outStream;

    outStream << id << "st_vol" << ","
              << id << "capacity" << ","
              << id << "net_inf" << ","
              << id << "st_rof" << ","
              << id << "lt_rof" << ","
              << id << "rest_demand" << ","
              << id << "unrest_demand" << ","
              << id << "unfulf_demand" << ","
              << id << "wastewater" << ","
              << id << "treat_capacity" << ",";

    return outStream.str();
}

void WSSDataCollector::collect_data() {
    // Collect available volume to match Original model behavior
    // Available volume includes flow-through sources (like intakes)
    double available_vol = wss->getTotal_available_volume();
    combined_storage.push_back(available_vol);
    storage_capacity.push_back(wss->getTotal_storage_capacity());
    storage_to_capacity_ratio.push_back(wss->getStorageToCapacityRatio());
    unrestricted_demand.push_back(wss->getUnrestrictedDemand());
    restricted_demand.push_back(wss->getRestrictedDemand());
    demand_multiplier.push_back(wss->getDemand_multiplier());
    demand_offset.push_back(wss->getDemand_offset());
    waste_water_discharge.push_back(wss->getWaste_water_discharge());
    unfulfilled_demand.push_back(wss->getUnfulfilled_demand());
    net_stream_inflow.push_back(wss->getNet_stream_inflow());
    total_treatment_capacity.push_back(wss->getTotal_treatment_capacity());
    total_storage_treatment_capacity.push_back(0.0); // Placeholder - implement getTotal_storage_treatment_capacity() if needed
    water_sources_count.push_back((int)wss->getWater_sources().size());
    short_term_rof.push_back(wss->getShort_term_risk_of_failure());
    long_term_rof.push_back(wss->getLong_term_risk_of_failure());

    // Calculate affordability index: weekly_water_price / weekly_average_income
    // Income is stored as monthly, so convert to weekly by dividing by WEEKS_IN_MONTH
    double affordability_index = 0.0;
    double average_monthly_income = wss->getAverageMonthlyIncome();
    
    if (wss->getWssGrossRevenue() > 0.0 && wss->getRestrictedDemand() > 0.0 && average_monthly_income > 0.0) {
        // Weekly water price = weekly revenue / weekly demand (both in same time unit)
        double weekly_water_price = wss->getWssGrossRevenue() / wss->getRestrictedDemand();
        // Convert monthly income to weekly income
        double weekly_average_income = average_monthly_income / WEEKS_IN_MONTH;
        // Affordability index = weekly water price / weekly income
        affordability_index = weekly_water_price / weekly_average_income;
    }
    weekly_affordability_index.push_back(affordability_index);
    
    double wss_npc = wss->getWssInfrastructureNPC();
    
    // Validate NPC value - catch memory corruption early
    if (std::isnan(wss_npc) || std::isinf(wss_npc) || wss_npc < -1e100 || wss_npc > 1e100) {
        char error[512];
        sprintf(error, "Invalid WSS NPC value %.2e for WSS %d (system_id=%d), week %zu, realization %lu. "
                "This indicates memory corruption or uninitialized access. WSS pointer: %p",
                wss_npc, id, wss->system_id, combined_storage.size(), realization, (void*)wss);
        throw std::runtime_error(error);
    }
    net_present_infrastructure_cost.push_back(wss_npc);
    
    // Collect infrastructure pathways from THIS WSS's infrastructure manager
    // This WSS data collector points to the REALIZATION WSS where infrastructure is actually built
    vector<vector<int>> infrastructure_pathways = 
        const_cast<InfrastructureManager&>(wss->getInfrastructure_construction_manager())
            .getAllAndClearInfraBuilt();
    
    for (const auto& pathway : infrastructure_pathways) {
        pathways.push_back(pathway);
    }

    // checkForNans();
}

void WSSDataCollector::checkForNans() const {
    string error = "nan collecting data for WSS " + to_string(id) + " in week " + to_string(combined_storage.size()) + ", realization " + to_string(realization);
    
    if (std::isnan(combined_storage.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(storage_capacity.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(storage_to_capacity_ratio.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(unrestricted_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(restricted_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(demand_multiplier.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(waste_water_discharge.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(unfulfilled_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(net_stream_inflow.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(total_treatment_capacity.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(total_storage_treatment_capacity.back()))
        throw_with_nested(runtime_error(error.c_str()));
}

// Getters implementation
const vector<double> &WSSDataCollector::getCombined_storage() const {
    return combined_storage;
}

const vector<double> &WSSDataCollector::getStorage_capacity() const {
    return storage_capacity;
}

const vector<double> &WSSDataCollector::getStorage_to_capacity_ratio() const {
    return storage_to_capacity_ratio;
}

const vector<double> &WSSDataCollector::getUnrestricted_demand() const {
    return unrestricted_demand;
}

const vector<double> &WSSDataCollector::getRestricted_demand() const {
    return restricted_demand;
}

const vector<double> &WSSDataCollector::getDemand_multiplier() const {
    return demand_multiplier;
}

const vector<double> &WSSDataCollector::getDemand_offset() const {
    return demand_offset;
}

const vector<double> &WSSDataCollector::getWaste_water_discharge() const {
    return waste_water_discharge;
}

const vector<double> &WSSDataCollector::getUnfulfilled_demand() const {
    return unfulfilled_demand;
}

const vector<double> &WSSDataCollector::getNet_stream_inflow() const {
    return net_stream_inflow;
}

const vector<double> &WSSDataCollector::getTotal_treatment_capacity() const {
    return total_treatment_capacity;
}

// const vector<double> &WSSDataCollector::getTotal_storage_treatment_capacity() const {
//     return total_storage_treatment_capacity;
// }

const vector<int> &WSSDataCollector::getWater_sources_count() const {
    return water_sources_count;
}

const WaterSupplySystems *WSSDataCollector::getWss() const {
    return wss;
}

// Financial data getters
const vector<double> &WSSDataCollector::getGross_revenues() const {
    return gross_revenues;
}

const vector<double> &WSSDataCollector::getDrought_mitigation_cost() const {
    return drought_mitigation_cost;
}

const vector<double> &WSSDataCollector::getContingency_fund_contribution() const {
    return contingency_fund_contribution;
}

const vector<double> &WSSDataCollector::getDebt_service_payments() const {
    return debt_service_payments;
}

const vector<double> &WSSDataCollector::getNet_present_infrastructure_cost() const {
    return net_present_infrastructure_cost;
}

const vector<double> &WSSDataCollector::getContingency_fund_size() const {
    return contingency_fund_size;
}

const vector<double> &WSSDataCollector::getWeekly_affordability_index() const {
    return weekly_affordability_index;
}

const Utility *WSSDataCollector::getOwner() const {
    return owner;
}

const vector<vector<int>> &WSSDataCollector::getPathways() const {
    return pathways;
}