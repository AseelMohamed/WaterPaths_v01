#include "WaterSupplySystems.h"
#include "Utility.h"
#include "../../Utils/Utils.h"
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstdio>

/**
 * Basic constructor for the WaterSupplySystem class.
 * @param name Water Supply System name (e.g. Raleigh_WSS_1)
 * @param system_id Numeric ID assigned to that water supply system.
 * @param utility_id Numeric ID of the parent utility that owns this WSS.
 * @param owner_utility Pointer to the parent utility object.
 * @param wwtp_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 */
WaterSupplySystems::WaterSupplySystems(
        const string& name,
        int system_id,
        int utility_id,
        Utility* owner_utility,
        const WwtpDischargeRule& wwtp_rule) :
        system_id(system_id),
        utility_id(utility_id),
        demand_buffer(0.0),
        number_of_week_demands(0),
        name(name.c_str()),
        owner(owner_utility),
        wwtp_discharge_rule(const_cast<WwtpDischargeRule&>(wwtp_rule)),
        demands_all_realizations(*(new vector<vector<double>>())) {
    
    // Initialize default values for basic constructor
    total_storage_capacity = 0.0;
    total_available_volume = 0.0;
    total_stored_volume = 0.0;
    total_treatment_capacity = 0.0;
    n_sources = 0;
    max_capacity = 0.0;
    
    // Initialize infrastructure construction manager (basic constructor with no infrastructure)
    infrastructure_construction_manager = InfrastructureManager(
        system_id, 
        vector<double>(),  // empty construction triggers
        vector<vector<int>>(),  // empty if_built_remove
        0.0,  // no discount rate 
        0.0,  // no bond term
        0.0,  // no bond interest rate
        vector<int>(),  // empty rof construction order
        vector<int>()   // empty demand construction order
    );
}

/**
 * Main constructor for the WaterSupplySystem class.
 * @param name Water Supply System name (e.g. Raleigh_WSS_1)
 * @param system_id Numeric ID assigned to that water supply system.
 * @param utility_id Numeric ID of the parent utility that owns this WSS.
 * @param owner_utility Pointer to the parent utility object.
 * @param demands_all_realizations Vector containing WSS's demand series for all realizations.
 * @param number_of_week_demands Length of weeks in demand series.
 * @param wwtp_discharge_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 * @param demand_buffer Buffer factor for demand calculations.
 * @param water_source_to_wtp Mapping of water sources to water treatment plants.
 * @param wss_owned_wtp_capacities Treatment capacities owned by this WSS for each WTP.
 */
WaterSupplySystems::WaterSupplySystems(
        const string& name,
        int system_id,
        int utility_id,
        Utility* owner_utility,
        vector<vector<double>>& demands_all_realizations,
        int number_of_week_demands,
        WwtpDischargeRule wwtp_discharge_rule,
        double demand_buffer,
        vector<vector<int>> water_source_to_wtp,
        vector<double> wss_owned_wtp_capacities) :
        system_id(system_id),
        utility_id(utility_id),
        demand_buffer(demand_buffer),
        number_of_week_demands(number_of_week_demands),
        name(name.c_str()),
        owner(owner_utility),
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations) {
    
    // Initialize default values following Utility pattern
    total_storage_capacity = 0.0;
    total_available_volume = 0.0;
    total_stored_volume = 0.0;
    total_treatment_capacity = 0.0;
    n_sources = 0;
    max_capacity = 0.0;
    
    // Copy WTP capacities
    this->wss_owned_wtp_capacities = wss_owned_wtp_capacities;
    
    // Initialize infrastructure construction manager (with empty construction orders for basic constructor)
    infrastructure_construction_manager = InfrastructureManager(
        system_id, 
        vector<double>(),  // empty construction triggers
        vector<vector<int>>(),  // empty if_built_remove
        0.0,  // no discount rate 
        0.0,  // no bond term
        0.0,  // no bond interest rate
        vector<int>(),  // empty rof construction order
        vector<int>()   // empty demand construction order
    );
    
    // Setup water source to WTP mapping, following Utility pattern
    unrollWaterSourceToWtpVector(water_source_to_wtp,
                                wss_owned_wtp_capacities);
    
    // printf("DEBUG: WSS Constructor - about to connect water sources vectors\n");
    // printf("DEBUG: water_sources vector address = %p, size = %zu\n", (void*)&water_sources, water_sources.size());
    
    // Connect water sources vectors to the infrastructure manager
    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
            
    // printf("DEBUG: WSS Constructor - connected water sources vectors\n");
}

/**
 * Constructor for when there is infrastructure to be built.
 * @param name Water Supply System name (e.g. Raleigh_WSS_1)
 * @param system_id Numeric ID assigned to that water supply system.
 * @param utility_id Numeric ID of the parent utility that owns this WSS.
 * @param owner_utility Pointer to the parent utility object.
 * @param demands_all_realizations Vector containing WSS's demand series for all realizations.
 * @param number_of_week_demands Length of weeks in demand series.
 * @param wwtp_discharge_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 * @param demand_buffer Buffer factor for demand calculations.
 * @param water_source_to_wtp Mapping of water sources to water treatment plants.
 * @param wss_owned_wtp_capacities Treatment capacities owned by this WSS for each WTP.
 * @param rof_infra_construction_order Order of infrastructure construction based on risk of failure.
 * @param demand_infra_construction_order Order of infrastructure construction based on demand.
 * @param infra_construction_triggers Triggers for infrastructure construction.
 */
WaterSupplySystems::WaterSupplySystems(
        const string& name,
        int system_id,
        int utility_id,
        Utility* owner_utility,
        vector<vector<double>>& demands_all_realizations,
        int number_of_week_demands,
        WwtpDischargeRule wwtp_discharge_rule,
        double demand_buffer,
        vector<vector<int>> water_source_to_wtp,
        vector<double> wss_owned_wtp_capacities,
        const vector<int> &rof_infra_construction_order,
        const vector<int> &demand_infra_construction_order,
        const vector<double> &infra_construction_triggers) :
        system_id(system_id),
        utility_id(utility_id),
        demand_buffer(demand_buffer),
        number_of_week_demands(number_of_week_demands),
        name(name.c_str()),
        owner(owner_utility),
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations) {
    
    // Initialize default values
    total_storage_capacity = 0.0;
    total_available_volume = 0.0;
    total_stored_volume = 0.0;
    total_treatment_capacity = 0.0;
    n_sources = 0;
    max_capacity = 0.0;
    
    // Copy WTP capacities
    this->wss_owned_wtp_capacities = wss_owned_wtp_capacities;
    
    // Initialize infrastructure construction manager with infrastructure parameters (no financial params)
    infrastructure_construction_manager = InfrastructureManager(
        system_id, 
        infra_construction_triggers,
        vector<vector<int>>(),  // empty infra_if_built_remove (WSS doesn't handle complex financial rules)
        0.0,  // no discount rate (handled by utility)
        0.0,  // no bond term (handled by utility)
        0.0,  // no bond interest rate (handled by utility)
        rof_infra_construction_order,
        demand_infra_construction_order
    );
    
    // Setup water source to WTP mapping
    unrollWaterSourceToWtpVector(water_source_to_wtp,
                                wss_owned_wtp_capacities);
    
    // Connect water sources vectors to the infrastructure manager
    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Validation checks for infrastructure orders
    if (rof_infra_construction_order.empty() &&
        demand_infra_construction_order.empty())
        throw std::invalid_argument("At least one infrastructure construction "
                                    "order vector must have at least "
                                    "one water source ID. If there's "
                                    "no infrastructure to be built, "
                                    "use other constructor "
                                    "instead.");
    
    if (demands_all_realizations.empty()) {
        char error[256];
        sprintf(error, "Empty demand vectors passed to WSS %d", system_id);
        throw std::invalid_argument(error);
    }
}

WaterSupplySystems::~WaterSupplySystems() {
    water_sources.clear();
}

///////////   =============================================== ///////////

void WaterSupplySystems::unrollWaterSourceToWtpVector(
        const vector<vector<int>> &water_source_to_wtp,
        const vector<double> &wss_owned_wtp_capacities) {

    if (water_source_to_wtp.size() != wss_owned_wtp_capacities.size()) {
        char error[512];
        sprintf(error, "WSS %d has %zu WTPs but %zu water sources (or "
                       "groups of) assigned to WTPs.", system_id,
                wss_owned_wtp_capacities.size(),
                water_source_to_wtp.size());
        throw invalid_argument(error);
    }

    for (int i = 0; i < water_source_to_wtp.size(); ++i) {
        for (int ws : water_source_to_wtp[i]) {
            if (ws >= this->water_source_to_wtp.size()) {
                this->water_source_to_wtp.resize(ws + 1, NON_INITIALIZED);
            }
            this->water_source_to_wtp[ws] = i;
        }
    }
}

void WaterSupplySystems::reconnectInfrastructureManager() {
    // printf("DEBUG: reconnectInfrastructureManager called\n");
    // printf("DEBUG: Reconnecting with water_sources vector address = %p, size = %zu\n", 
        //    (void*)&water_sources, water_sources.size());
    
    // Reconnect the infrastructure manager to the current water sources vectors
    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // printf("DEBUG: Reconnection complete\n");
}

void WaterSupplySystems::updateTreatmentAndNumberOfStorageSources() {
    n_storage_sources = non_priority_draw_water_source.size();
    delete[] available_treated_flow_rate;
    available_treated_flow_rate = new double[non_priority_draw_water_source.size()];
    for (int i = 0; i < n_storage_sources; ++i) {
        auto ws = water_sources[non_priority_draw_water_source[i]];
        available_treated_flow_rate[i] = wss_owned_wtp_capacities[water_source_to_wtp[ws->id]];
        total_storage_treatment_capacity += available_treated_flow_rate[i];
    }

    total_treatment_capacity = accumulate(wss_owned_wtp_capacities.begin(),
                                          wss_owned_wtp_capacities.end(),
                                          0.);
}

void WaterSupplySystems::updateTotalAvailableVolume() {
    total_available_volume = 0.0;
    total_stored_volume = 0.0;
    net_stream_inflow = 0.0;

    for (int ws : priority_draw_water_source) {
        total_available_volume +=
                max(1.0e-6,
                    water_sources[ws]->getAvailableAllocatedVolume(system_id));
        net_stream_inflow += water_sources[ws]->getAllocatedInflow(system_id);
    }

    for (int i = 0; i < non_priority_draw_water_source.size(); ++i) {
        auto ws = water_sources[non_priority_draw_water_source[i]];
        double stored_volume = max(1.0e-6,
                                   ws->getAvailableAllocatedVolume(system_id));
        total_available_volume += stored_volume;
        total_stored_volume += stored_volume;
        net_stream_inflow += ws->getAllocatedInflow(system_id);
        available_treated_flow_rate[i] = wss_owned_wtp_capacities[water_source_to_wtp[ws->id]];
    }
}

void WaterSupplySystems::clearWaterSources() {
    water_sources.clear();
}

/**
 * Connects a reservoir to the WSS.
 * @param water_source
 */
void WaterSupplySystems::addWaterSource(WaterSource* water_source) {
    // printf("DEBUG: addWaterSource called for water_source id=%d\n", water_source->id);
    // printf("DEBUG: current water_sources vector address = %p, size = %zu\n", (void*)&water_sources, water_sources.size());
    
    checkErrorsAddWaterSourceOnline(water_source);

    // Add water sources with their IDs matching the water sources vector
    // indexes.
    if (water_source->id > (int) water_sources.size() - 1) {
        // printf("DEBUG: Resizing water_sources from %zu to %d\n", water_sources.size(), water_source->id + 1);
        water_sources.resize((unsigned int) water_source->id + 1);
    }

    // Add water source
    water_sources[water_source->id] = water_source;
    // printf("DEBUG: Added water source to vector, new size = %zu\n", water_sources.size());

    // Add water source to infrastructure construction manager.
    infrastructure_construction_manager.addWaterSource(water_source);

    // If watersource is online and the WSS owns some of its installed
    // treatment capacity, make it online.
    if (water_source->isOnline() && 
        water_source->id < water_source_to_wtp.size() &&
        water_source_to_wtp[water_source->id] != NON_INITIALIZED &&
        water_source_to_wtp[water_source->id] < wss_owned_wtp_capacities.size() &&
        wss_owned_wtp_capacities[water_source_to_wtp[water_source->id]] > 0) {
        // printf("DEBUG: Water source is online and has WTP capacity, calling addWaterSourceToOnlineLists\n");
        infrastructure_construction_manager.addWaterSourceToOnlineLists(
                water_source->id, total_storage_capacity,
                total_available_volume,
                total_stored_volume);
    }

    n_sources++;
    max_capacity += water_source->getAllocatedCapacity(system_id);

    updateTreatmentAndNumberOfStorageSources();
}

void WaterSupplySystems::checkErrorsAddWaterSourceOnline(WaterSource* water_source) {
    for (WaterSource *ws : water_sources) {
        if ((ws != nullptr) && ws->id == water_source->id) {
            cout << "Water source ID: " << water_source->id << endl <<
                 "WSS ID: " << system_id << endl;
            throw invalid_argument("Attempt to add water source with "
                                   "duplicate ID to WSS.");
        }
    }
}


bool WaterSupplySystems::idealDemandSplitUnconstrained(
        double* split_demands,
        const double* available_treated_flow_rate,
        double total_demand,
        const double* storage,
        double total_storage,
        int n_storage_sources) {
    
    bool treatment_capacity_violated = false;
    for (int i = 0; i < n_storage_sources; ++i) {
        split_demands[i] = total_demand * storage[i] / total_storage;
        if (split_demands[i] - 1e-9 > available_treated_flow_rate[i]) {
            treatment_capacity_violated = true;
        }
    }
    return treatment_capacity_violated;
}

bool WaterSupplySystems::idealDemandSplitConstrained(
        double* split_demands,
        bool* over_allocated,
        bool* has_spare_capacity,
        const double* available_treated_flow_rate,
        double total_demand,
        const double* storage,
        double total_storage,
        int n_storage_sources) {
    
    // Consider only storage of sources that are still not at provision capacity.
    total_storage = 0;
    for (int j = 0; j < n_storage_sources; ++j) {
        if (has_spare_capacity[j]) total_storage += storage[j];
    }

    // Split demands not fulfilled by sources at provision capacity across
    // sources with spare capacity while checking for over allocation.
    bool treatment_capacity_violated = false;
    for (int i = 0; i < n_storage_sources; ++i) {
        if (has_spare_capacity[i]) {
            split_demands[i] = total_demand * storage[i] / total_storage;
        }
        over_allocated[i] =
                split_demands[i] - 1e-9 > available_treated_flow_rate[i];
        has_spare_capacity[i] =
                split_demands[i] + 1e-9 < available_treated_flow_rate[i];
        if (over_allocated[i]) treatment_capacity_violated = true;
    }
    return treatment_capacity_violated;
}

void WaterSupplySystems::splitDemands(
    int week, vector<vector<double>> &demands,
    bool apply_demand_buffer) {
    auto wss_owned_wtp_capacities = this->wss_owned_wtp_capacities;
    unrestricted_demand = demand_series_realization[week] +
                          apply_demand_buffer * demand_buffer *
                          weekly_peaking_factor[Utils::weekOfTheYear(week)];
    restricted_demand = unrestricted_demand * demand_multiplier - demand_offset;
    unfulfilled_demand = max(max(restricted_demand - total_available_volume,
                                 restricted_demand - total_treatment_capacity),
                             0.);
    restricted_demand -= unfulfilled_demand;
    double demand_non_priority_sources = restricted_demand;
    double total_serviced_demand = 0;

    // Allocates demand to intakes and reuse based on allocated volume to
    // this watersupply system.
    for (int &ws : priority_draw_water_source) {
        double max_source_output = min(
                water_sources[ws]->getAvailableAllocatedVolume(system_id),
                wss_owned_wtp_capacities[water_source_to_wtp[ws]]);
        double source_demand =
                min(demand_non_priority_sources,
                    max_source_output);
    // Write demand contribution into the owning utility's column
    demands[ws][this->utility_id] = source_demand;
        demand_non_priority_sources -= source_demand;
        total_serviced_demand += source_demand;
        wss_owned_wtp_capacities[water_source_to_wtp[ws]] -= source_demand;
    }

    vector<double> storages(n_storage_sources);
    double total_available_flow_rate = 0;
    for (int i = 0; i < n_storage_sources; ++i) {
        auto ws = water_sources[non_priority_draw_water_source[i]];
        storages[i] = ws->getAvailableAllocatedVolume(system_id);
        available_treated_flow_rate[i] = min(
                storages[i],
                wss_owned_wtp_capacities[water_source_to_wtp[ws->id]]
        );
        total_available_flow_rate += available_treated_flow_rate[i];
    }

    bool treatment_capacity_violation = false;
    if (demand_non_priority_sources > total_available_flow_rate) {
        // If the WSS's demand is greater than the sum of treatment
        // capacities of all water sources, all WTPs will be fully used.
        for (int i = 0; i < n_storage_sources; ++i) {
        // Write demand contribution into the owning utility's column
        demands[non_priority_draw_water_source[i]][utility_id] =
            available_treated_flow_rate[i];
        }
        treatment_capacity_violation = true;
    } else if (demand_non_priority_sources > 0) {
        // If a given WTP cannot fulfill its ideal demand but there is spare
        // treatment capacity available in other WTPs, use it.

        // Create auxiliary variables and check which sources are over allocated
        // and which have spare capacity.
        bool* has_spare_flow_rate_array = new bool[n_storage_sources];
        bool* over_allocated_array = new bool[n_storage_sources];
        
        // Initialize arrays
        for (int i = 0; i < n_storage_sources; ++i) {
            has_spare_flow_rate_array[i] = true;
            over_allocated_array[i] = false;
        }
        vector<double> split_demands(n_storage_sources);

        treatment_capacity_violation = idealDemandSplitUnconstrained(
                split_demands.data(),
                available_treated_flow_rate,
                demand_non_priority_sources,
                storages.data(),
                total_stored_volume,
                n_storage_sources);

        if (treatment_capacity_violation) {
            // Check which sources are over allocated or have spare capacity.
            for (int i = 0; i < n_storage_sources; ++i) {
                over_allocated_array[i] = split_demands[i] - 1e-9 >
                                    available_treated_flow_rate[i];
                has_spare_flow_rate_array[i] =
                        split_demands[i] + 1e-9 <
                        available_treated_flow_rate[i];
            }

            // Redistribute demands across water sources that may still have
            // spare capacity.
            while (treatment_capacity_violation) {
                double remainder_demand = demand_non_priority_sources;
                for (int i = 0; i < n_storage_sources; ++i) {
                    if (over_allocated_array[i]) {
                        split_demands[i] = available_treated_flow_rate[i];
                    }
                    if (!has_spare_flow_rate_array[i]) {
                        remainder_demand -= split_demands[i];
                    }
                }
                treatment_capacity_violation = idealDemandSplitConstrained(
                        split_demands.data(),
                        over_allocated_array,
                        has_spare_flow_rate_array,
                        available_treated_flow_rate,
                        remainder_demand,
                        storages.data(),
                        total_stored_volume,
                        non_priority_draw_water_source.size());
            }
        }

        for (int j = 0; j < n_storage_sources; ++j) {
            // Write demand contribution into the owning utility's column
            demands[non_priority_draw_water_source[j]][utility_id] = split_demands[j];
        }
        
        // Clean up dynamically allocated arrays
        delete[] has_spare_flow_rate_array;
        delete[] over_allocated_array;
    }
}

void WaterSupplySystems::setWaterSourceOnline(unsigned int source_id, int week) {
    infrastructure_construction_manager.setWaterSourceOnline(
            source_id, week, wss_owned_wtp_capacities, water_source_to_wtp,
            total_storage_capacity, total_available_volume,
            total_stored_volume);

    updateTreatmentAndNumberOfStorageSources();
}

/**
 * Check if new infrastructure is to be triggered based on long-term risk of failure and, if so, handle
 * the beginning of construction, issue corresponding bonds and update debt.
 * @param long_term_rof
 * @param week
 * @return
 */
int WaterSupplySystems::infrastructureConstructionHandler(double long_term_rof, int week) {
    double past_year_average_demand = 0;
    if (week >= (int) WEEKS_IN_YEAR) {
        //     past_year_average_demand =
        //            std::accumulate(demand_series_realization.begin() + week - (int) WEEKS_IN_YEAR,
        //                            demand_series_realization.begin() + week, 0.0) / WEEKS_IN_YEAR;

        for (int w = week - (int) WEEKS_IN_YEAR; w < week; ++w) {
            past_year_average_demand += demand_series_realization.at(w);
        }
    }

    long_term_risk_of_failure = long_term_rof;

    // Check if new infrastructure is to be triggered and, if so, trigger it.
    int new_infra_triggered = infrastructure_construction_manager.infrastructureConstructionHandler(
            long_term_rof, week,
            past_year_average_demand,
            wss_owned_wtp_capacities,
            water_source_to_wtp,
            total_storage_capacity,
            total_available_volume,
            total_stored_volume);

    // Issue and add bond of triggered water source to list of outstanding bonds, and update total new
    // infrastructure NPV.
    // issueBond(new_infra_triggered, week);

    updateTreatmentAndNumberOfStorageSources();

    return new_infra_triggered;
}

void WaterSupplySystems::calculateWastewater_releases(int week, double* discharges) {
    double discharge;
    waste_water_discharge = 0;

    for (int &id_disch : wwtp_discharge_rule.discharge_to_source_ids) {
        discharge = restricted_demand * wwtp_discharge_rule
                .get_dependent_variable(id_disch, Utils::weekOfTheYear(week));
        discharges[id_disch] += discharge;

        waste_water_discharge += discharge;
    }
}

void WaterSupplySystems::setDemand_offset(double demand_offset, double offset_rate_per_volume) {
    this->demand_offset = demand_offset;
    this->offset_rate_per_volume = offset_rate_per_volume;
}

/**
 * Get time series corresponding to realization index and eliminate reference to
 * comprehensive demand data set.
 * @param r
 */
void WaterSupplySystems::setRealization(unsigned long r, vector<double> &rdm_factors) {
    unsigned long n_weeks = demands_all_realizations.at(r).size();
    demand_series_realization = vector<double>(n_weeks);

    // Apply demand multiplier and copy demands pertaining to current realization.
    double delta_demand =
            demands_all_realizations.at(r)[0] * (1. - rdm_factors.at(0));
    for (unsigned long w = 0; w < n_weeks; ++w) {
        demand_series_realization[w] = demands_all_realizations.at(r)[w] *
                                       rdm_factors.at(0)
                                       + delta_demand;
    }

    // Set peaking demand factor.
    weekly_peaking_factor = calculateWeeklyPeakingFactor
            (&demands_all_realizations.at(r));

}

vector<double> WaterSupplySystems::calculateWeeklyPeakingFactor(vector<double> *demands) {
    unsigned long n_weeks = (unsigned long) WEEKS_IN_YEAR + 1;
    int n_years = (int) (demands->size() / WEEKS_IN_YEAR - 1);
    vector<double> year_averages(n_weeks,
                                 0.0);

    double year_average_demand;
    for (int y = 0; y < n_years; ++y) {
        year_average_demand = accumulate(
                demands->begin() + y * WEEKS_IN_YEAR,
                demands->begin() + (y + 1) * WEEKS_IN_YEAR,
                0.0) /
                              ((int) ((y + 1) * WEEKS_IN_YEAR) -
                               (int) (y * WEEKS_IN_YEAR));
        for (int w = 0; w < n_weeks; ++w) {
            year_averages[w] += (*demands)[y * WEEKS_IN_YEAR + w] /
                                year_average_demand / n_years;
        }
    }

    return year_averages;
}


//========================= GETTERS AND SETTERS =============================//

int WaterSupplySystems::getSystemId() const { 
    return system_id; 
}

bool WaterSupplySystems::hasTreatmentConnected(int ws) {
    return wss_owned_wtp_capacities[water_source_to_wtp[ws]] > 0.;
}

double WaterSupplySystems::getStorageToCapacityRatio() const { 
    return total_storage_capacity > 0 ? total_stored_volume / total_storage_capacity : 0.0; 
}

double WaterSupplySystems::getTotal_available_volume() const {
     return total_available_volume; 
}

double WaterSupplySystems::getTotal_stored_volume() const { 
    return total_stored_volume; 
}

double WaterSupplySystems::getTotal_storage_capacity() const { 
    return total_storage_capacity; 
}

double WaterSupplySystems::getRisk_of_failure() const { 
    return short_term_risk_of_failure; 
}

void WaterSupplySystems::setRisk_of_failure(double risk_of_failure) {
    this->short_term_risk_of_failure = risk_of_failure;
}

double WaterSupplySystems::getTotal_treatment_capacity() const { 
    return total_treatment_capacity; 
}

void WaterSupplySystems::setDemand_multiplier(double demand_multiplier) {
    WaterSupplySystems::demand_multiplier = demand_multiplier;
}

double WaterSupplySystems::getUnrestrictedDemand(int week) const {
    if (week == -1) {
        return unrestricted_demand;
    } else {
        return demand_series_realization[week];
    }
}

double WaterSupplySystems::getRestrictedDemand() const {
    return restricted_demand;
}

double WaterSupplySystems::getDemand_multiplier() const {
    return demand_multiplier;
}

double WaterSupplySystems::getLong_term_risk_of_failure() const { 
    return long_term_risk_of_failure; 
}

const vector<WaterSource*>& WaterSupplySystems::getWater_sources() const { 
    return water_sources; 
}

double WaterSupplySystems::getWaste_water_discharge() const { 
    return waste_water_discharge; 
}

void WaterSupplySystems::resetTotal_storage_capacity() {
    WaterSupplySystems::total_storage_capacity = 0;
}

double WaterSupplySystems::getUnfulfilled_demand() const { 
    return unfulfilled_demand; 
}

double WaterSupplySystems::getNet_stream_inflow() const { 
    return net_stream_inflow; 
}

// double WaterSupplySystems::getShort_term_risk_of_failure() const {
//     return short_term_risk_of_failure;
// }

// double WaterSupplySystems::getTotal_storage_treatment_capacity() const {
//     return total_storage_treatment_capacity;
// }

double WaterSupplySystems::getDemand_offset() const {
    return demand_offset;
}