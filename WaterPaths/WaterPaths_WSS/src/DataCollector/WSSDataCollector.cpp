//
// Created by AI Assistant on 10/13/25.
//

#include <iomanip>
#include <sstream>
#include "WSSDataCollector.h"

WSSDataCollector::WSSDataCollector(const WaterSupplySystems *wss, unsigned long realization)
        : DataCollector(wss->getSystemId(), wss->name, realization, UTILITY, 11 * COLUMN_WIDTH),
          wss(wss) {
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
              << storage_to_capacity_ratio[week]
              << ","
              << net_stream_inflow[week]
              << ","
              << restricted_demand[week]
              << ","
              << unrestricted_demand[week]
              << ","
              << unfulfilled_demand[week]
              << ","
              << demand_multiplier[week]
              << ","
              << waste_water_discharge[week]
              << ","
              << total_treatment_capacity[week]
              << ","
              << total_storage_treatment_capacity[week]
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
              << id << "st_capacity" << ","
              << id << "st_ratio" << ","
              << id << "net_inf" << ","
              << id << "rest_demand" << ","
              << id << "unrest_demand" << ","
              << id << "unfulf_demand" << ","
              << id << "demand_mult" << ","
              << id << "wastewater" << ","
              << id << "treat_capacity" << ","
              << id << "st_treat_capacity" << ",";

    return outStream.str();
}

void WSSDataCollector::collect_data() {
    combined_storage.push_back(wss->getTotal_available_volume());
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
    // total_storage_treatment_capacity.push_back(wss->getTotal_storage_treatment_capacity());
    water_sources_count.push_back((int)wss->getWater_sources().size());

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