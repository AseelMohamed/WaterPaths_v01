#include "WaterSupplySystems.h"
#include "Utility.h"
#include "../../Utils/Utils.h"
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <iostream>
#include <cstdio>

/**
 * Basic constructor for the WaterSupplySystem class.
 * @param name Water Supply System name (e.g. Raleigh_WSS_1)
 * @param system_id Numeric ID assigned to that water supply system.
 * @param utility_id Numeric ID of the parent utility that owns this WSS.
 * @param owner_utility Pointer to the parent utility object.
 * @param demands_all_realizations Demands vector (externally provided)
 * @param wwtp_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 */
WaterSupplySystems::WaterSupplySystems(
        const string& name,
        int system_id,
        int utility_id,
        Utility* owner_utility,
        vector<vector<double>>& demands_all_realizations,
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
    wss_infrastructure_npc = 0.0;  // Explicit initialization
    
    // Initialize demand_series_realization as empty vector to prevent crashes
    // This will be properly set when setRealization() is called
    demand_series_realization = vector<double>();
    
    // Initialize infrastructure construction manager (basic constructor with no infrastructure)
    infrastructure_construction_manager = InfrastructureManager(
        system_id,  // Use system_id so treatment capacity lookup is correct for this WSS
        vector<double>(),  // empty construction triggers
        vector<vector<int>>(),  
        0.0,  
        0.0,  
        0.0,  
        vector<int>(),  
        vector<int>()   
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
    wss_infrastructure_npc = 0.0;  // Explicit initialization
    
    // Copy WTP capacities
    this->wss_owned_wtp_capacities = wss_owned_wtp_capacities;
    
    // Initialize infrastructure construction manager (with empty construction orders for basic constructor)
    infrastructure_construction_manager = InfrastructureManager(
        system_id,  // Use system_id so treatment capacity lookup is correct
        vector<double>(),  // empty construction triggers
        vector<vector<int>>(),  
        0.0,  
        0.0,  
        0.0,  
        vector<int>(), 
        vector<int>()   
    );
    
    // Setup water source to WTP mapping, following Utility pattern
    unrollWaterSourceToWtpVector(water_source_to_wtp,
                                wss_owned_wtp_capacities);
 
    // NOTE: Do not connect infrastructure manager here with empty water_sources vector.
    // The connection will happen automatically when water sources are added via addWaterSource()
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
        demands_all_realizations(demands_all_realizations) {  // Reference to externally-owned demands
    
    // Initialize default values
    total_storage_capacity = 0.0;
    total_available_volume = 0.0;
    total_stored_volume = 0.0;
    total_treatment_capacity = 0.0;
    n_sources = 0;
    max_capacity = 0.0;
    wss_infrastructure_npc = 0.0;  // Explicit initialization
    
    // Copy WTP capacities
    this->wss_owned_wtp_capacities = wss_owned_wtp_capacities;
    
    // Initialize infrastructure construction manager with infrastructure parameters (no financial params)
    infrastructure_construction_manager = InfrastructureManager(
        system_id,  // Use system_id so treatment capacity lookup is correct
        infra_construction_triggers,
        vector<vector<int>>(),  
        0.0,  
        0.0,  
        0.0,  
        rof_infra_construction_order,
        demand_infra_construction_order
    );
    
    // Setup water source to WTP mapping
    unrollWaterSourceToWtpVector(water_source_to_wtp,
                                wss_owned_wtp_capacities);
    
    // NOTE: Do not connect infrastructure manager here with empty water_sources vector.
    // The connection will happen automatically when water sources are added via addWaterSource()
    
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

/**
 * Copy constructor - Performs deep copy of WaterSupplySystems object.
 * This is critical because Simulation.cpp creates copies for realization and ROF models.
 * Without proper copy constructor, C++ does shallow copy which causes double-delete crashes.
 *
 * @param other The WaterSupplySystems object to copy from
 */
WaterSupplySystems::WaterSupplySystems(const WaterSupplySystems& other) :
        system_id(other.system_id),
        utility_id(other.utility_id),
        demand_buffer(other.demand_buffer),
        number_of_week_demands(other.number_of_week_demands),
        name(other.name),  // Shallow copy of const char* is OK (points to string literal)
        owner(other.owner),  // Shallow copy of Utility pointer is OK (not owned by WSS)
        wwtp_discharge_rule(const_cast<WwtpDischargeRule&>(other.wwtp_discharge_rule)),  // Copy via copy constructor
        demands_all_realizations(const_cast<vector<vector<double>>&>(other.demands_all_realizations)) {  // Reference to same externally-owned vector

    
    // Validate that the copied owner is not null
    if (owner == nullptr) {
        throw std::invalid_argument("WaterSupplySystems copy constructor: copied owner cannot be null");
    }
    
    // Copy all simple member variables
    short_term_risk_of_failure = other.short_term_risk_of_failure;
    long_term_risk_of_failure = other.long_term_risk_of_failure;
    total_storage_capacity = other.total_storage_capacity;
    total_available_volume = other.total_available_volume;
    total_stored_volume = other.total_stored_volume;
    total_treatment_capacity = other.total_treatment_capacity;
    total_storage_treatment_capacity = other.total_storage_treatment_capacity;
    waste_water_discharge = other.waste_water_discharge;
    unfulfilled_demand = other.unfulfilled_demand;
    net_stream_inflow = other.net_stream_inflow;
    used_for_realization = other.used_for_realization;
    n_storage_sources = other.n_storage_sources;
    demand_multiplier = other.demand_multiplier;
    demand_offset = other.demand_offset;
    offset_rate_per_volume = other.offset_rate_per_volume;
    restricted_demand = other.restricted_demand;
    unrestricted_demand = other.unrestricted_demand;
    n_sources = other.n_sources;
    max_capacity = other.max_capacity;
    // Copy discount rate and bond multipliers from original - critical for realization isolation
    infra_discount_rate = other.infra_discount_rate;
    bond_term_multiplier = other.bond_term_multiplier;
    bond_interest_rate_multiplier = other.bond_interest_rate_multiplier;
    // Copy average monthly income for affordability calculation
    average_monthly_income = other.average_monthly_income;
    // NOTE: current_realization_id is NOT copied - it will be set correctly by setRealization() 
    // after this copy is created, avoiding race condition from shared original WSS
    // NOTE: wss_infrastructure_npc deliberately NOT copied - realization WSS should get NPC from ROF WSS via references
    
    // Deep copy all vectors (these are value types, so std::vector handles deep copy)
    priority_draw_water_source = other.priority_draw_water_source;
    non_priority_draw_water_source = other.non_priority_draw_water_source;
    weekly_peaking_factor = other.weekly_peaking_factor;
    
    // This prevents copying uninitialized or partially initialized vectors during parallel execution
    demand_series_realization = vector<double>();
    
    wss_owned_wtp_capacities = other.wss_owned_wtp_capacities;
    water_source_to_wtp = other.water_source_to_wtp;
    
    // Initialize empty water_sources vector
    // The actual water source pointers will be set later via addWaterSource() calls
    // in ContinuityModel constructor. DO NOT copy pointers directly as they're owned
    // by ContinuityModel and will be different copies for realization vs ROF models.
    water_sources.clear();
    water_sources.resize(other.water_sources.size(), nullptr);
    
    // Deep copy the infrastructure manager
    // Note: The id in infrastructure_construction_manager should be utility_id, not system_id
    infrastructure_construction_manager = other.infrastructure_construction_manager;
    
    // Deep copy the vector for available_treated_flow_rate
    available_treated_flow_rate = other.available_treated_flow_rate;
}

WaterSupplySystems::~WaterSupplySystems() {
    water_sources.clear();
    // No manual cleanup needed - std::vector handles its own memory
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
    // Infrastructure manager connection now happens automatically in addWaterSource()
    // when the first water source is added, similar to how Original model works.
    // No explicit reconnection needed here.
}

void WaterSupplySystems::updateTreatmentAndNumberOfStorageSources() {
    // Count actual storage sources by checking which water sources exist and are storage sources
    vector<int> actual_storage_sources;
    for (int i = 0; i < non_priority_draw_water_source.size(); ++i) {
        int ws_id = non_priority_draw_water_source[i];
        if (ws_id < water_sources.size() && water_sources[ws_id] != nullptr) {
            actual_storage_sources.push_back(ws_id);
        }
    }
    
    n_storage_sources = actual_storage_sources.size();
    available_treated_flow_rate.assign(n_storage_sources, 0.0);
    for (int i = 0; i < n_storage_sources; ++i) {
        auto ws = water_sources[actual_storage_sources[i]];
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
    double priority_vol = 0.0;
    double non_priority_vol = 0.0;

    for (int ws : priority_draw_water_source) {
        if (ws < water_sources.size() && water_sources[ws] != nullptr) {
            double vol = max(1.0e-6, water_sources[ws]->getAvailableAllocatedVolume(system_id));
            total_available_volume += vol;
            priority_vol += vol;
            net_stream_inflow += water_sources[ws]->getAllocatedInflow(system_id);
        }
    }

    for (int i = 0; i < non_priority_draw_water_source.size(); ++i) {
        int ws_id = non_priority_draw_water_source[i];
        if (ws_id < water_sources.size() && water_sources[ws_id] != nullptr) {
            auto ws = water_sources[ws_id];
            double stored_volume = max(1.0e-6, ws->getAvailableAllocatedVolume(system_id));
            total_available_volume += stored_volume;
            total_stored_volume += stored_volume;
            non_priority_vol += stored_volume;
            net_stream_inflow += ws->getAllocatedInflow(system_id);
        }
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
    checkErrorsAddWaterSourceOnline(water_source);

    // Add water sources with their IDs matching the water sources vector
    // indexes.
    if (water_source->id > (int) water_sources.size() - 1) {
        water_sources.resize((unsigned int) water_source->id + 1);
    }

    // Add water source
    water_sources[water_source->id] = water_source;

    // Connect infrastructure manager to water sources vectors when first water source is added
    // Check if this is the first water source by counting non-null entries
    int non_null_count = 0;
    for (auto* ws : water_sources) {
        if (ws != nullptr) non_null_count++;
    }
    
    if (non_null_count == 1) {  // This is the first water source
        infrastructure_construction_manager.connectWaterSourcesVectorsToWSS(
                water_sources,
                priority_draw_water_source,
                non_priority_draw_water_source);
    }

    // Add water source to infrastructure construction manager.
    infrastructure_construction_manager.addWaterSource(water_source);

    // If watersource is online and the WSS owns some of its installed
    // treatment capacity, make it online.
    if (water_source->isOnline() && 
        water_source->id < water_source_to_wtp.size() &&
        water_source_to_wtp[water_source->id] != NON_INITIALIZED &&
        water_source_to_wtp[water_source->id] < wss_owned_wtp_capacities.size() &&
        wss_owned_wtp_capacities[water_source_to_wtp[water_source->id]] > 0) {
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
            printf("Water source ID: %d\nWSS ID: %d\n", water_source->id, system_id);
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
    
    // Bounds check for demand_series_realization access
    if (week < 0 || week >= (int)demand_series_realization.size()) {
        char error_msg[512];
        sprintf(error_msg, "WSS %d: Week %d out of bounds for demand_series_realization (size=%zu)", 
                system_id, week, demand_series_realization.size());
        printf("ERROR: %s\n", error_msg);
        throw std::out_of_range(error_msg);
    }
    
    // Bounds check for weekly_peaking_factor access
    int week_of_year = Utils::weekOfTheYear(week);
    if (week_of_year < 0 || week_of_year >= (int)weekly_peaking_factor.size()) {
        char error_msg[512];
        sprintf(error_msg, "WSS %d: WeekOfYear %d (from week %d) out of bounds for weekly_peaking_factor (size=%zu)", 
                system_id, week_of_year, week, weekly_peaking_factor.size());
        printf("ERROR: %s\n", error_msg);
        throw std::out_of_range(error_msg);
    }
    
    unrestricted_demand = demand_series_realization[week] +
                          apply_demand_buffer * demand_buffer *
                          weekly_peaking_factor[week_of_year];
           
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
                available_treated_flow_rate.data(),
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
                        available_treated_flow_rate.data(),
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
    
    // Note: Financial calculations are handled at utility level, not WSS level
    // This follows the Original model architecture where Utility aggregates all WSS data
    // for financial calculations in Utility::updateFinancialCalculations()
}

void WaterSupplySystems::setWaterSourceOnline(unsigned int source_id, int week) {
    infrastructure_construction_manager.setWaterSourceOnline(
            source_id, week, wss_owned_wtp_capacities, water_source_to_wtp,
            total_storage_capacity, total_available_volume,
            total_stored_volume);

    updateTreatmentAndNumberOfStorageSources();
}

/**
 * Calculate wastewater releases for the Water Supply System.
 * @param week Week for which to calculate discharges
 * @param discharges Array to store discharge values for each source
 */
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

double WaterSupplySystems::getOffset_rate_per_volume() const {
    return offset_rate_per_volume;
}

/**
 * Get time series corresponding to realization index and eliminate reference to
 * comprehensive demand data set.
 * @param r
 */
void WaterSupplySystems::setRealization(unsigned long r, vector<double> &rdm_factors) {
    // Store which realization this WSS copy belongs to (for bond tracking)
    current_realization_id = r;
    
    // BOUNDS CHECK: Verify realization index is valid
    if (r >= demands_all_realizations.size()) {
        char error_msg[512];
        sprintf(error_msg, "ERROR [WSS %d]: setRealization called with realization index %lu, "
                "but demands_all_realizations only has %zu realizations available. "
                "This typically means ROF model is trying to access demand data it shouldn't need.",
                system_id, r, demands_all_realizations.size());
        cerr << error_msg << endl;
        
        // For ROF models, we don't actually need to set realization-specific demands
        // since they run with synthetic short-term scenarios. Just skip this operation.
        if (!used_for_realization) {
            return;
        }
        
        throw std::out_of_range(error_msg);
    }
    
    // Store realization-specific discount rate and bond multipliers for infrastructure NPC calculations
    // rdm_factors[1] = bond_interest_rate_multiplier
    // rdm_factors[2] = bond_term_multiplier  
    // rdm_factors[3] = infra_discount_rate_multiplier
    if (owner && rdm_factors.size() > 3) {
        double base_rate = owner->getBaseInfraDiscountRate();
        double rdm_multiplier = rdm_factors.at(3);
        infra_discount_rate = base_rate * rdm_multiplier;
        bond_interest_rate_multiplier = rdm_factors.at(1);
        bond_term_multiplier = rdm_factors.at(2);
    } else {
        // Fallback: if no owner or insufficient RDM factors, use default
        char error_msg[512];
        sprintf(error_msg, "WARNING: WSS %d setRealization(%lu): Cannot set discount rate (owner=%p, rdm_size=%zu)",
                system_id, r, (void*)owner, rdm_factors.size());
        fprintf(stderr, "%s\n", error_msg);
        infra_discount_rate = 0.05;  // Fallback default rate
    }
    
    // Simple, clean implementation like the original
    unsigned long n_weeks = demands_all_realizations.at(r).size();
    
    // Validate n_weeks is reasonable before allocating
    if (n_weeks == 0 || n_weeks > 100000) {
        char error_msg[512];
        sprintf(error_msg, "ERROR: WSS %d setRealization got invalid n_weeks=%lu from demands_all_realizations[%lu]",
                system_id, n_weeks, r);
        throw std::runtime_error(error_msg);
    }
    
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

void WaterSupplySystems::setTotal_available_volume(double volume) {
    total_available_volume = volume;
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
        // Check if vector is empty (not yet initialized by setRealization)
        size_t vec_size = demand_series_realization.size();
        if (vec_size == 0) {
            // Vector not yet initialized - return 0 demand rather than crashing
            // This can happen if ROF calculation accesses WSS before setRealization completes
            printf("WARNING: WSS %d getUnrestrictedDemand called before setRealization initialized demand vector (week=%d)\n",
                   system_id, week);
            return 0.0;
        }
        
        // Check for corrupted vector (impossibly large size indicates memory corruption)
        if (vec_size > 100000) {
            char error_msg[512];
            sprintf(error_msg, "FATAL: WSS %d (used_for_realization=%d) has corrupted demand_series_realization vector (size=%zu). "
                    "Week=%d. demands_all_realizations.size()=%zu. This indicates the vector was not properly initialized by setRealization().",
                    system_id, used_for_realization, vec_size, week, demands_all_realizations.size());
            printf("ERROR: %s\n", error_msg);
            throw std::runtime_error(error_msg);
        }
        
        // Bounds check for demand_series_realization access
        if (week < 0 || week >= (int)vec_size) {
            char error_msg[512];
            sprintf(error_msg, "WSS %d: Week %d out of bounds for demand_series_realization (size=%zu) in getUnrestrictedDemand", 
                    system_id, week, vec_size);
            printf("ERROR: %s\n", error_msg);
            throw std::out_of_range(error_msg);
        }
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

double WaterSupplySystems::getShort_term_risk_of_failure() const {
    return short_term_risk_of_failure;
}

// double WaterSupplySystems::getTotal_storage_treatment_capacity() const {
//     return total_storage_treatment_capacity;
// }

double WaterSupplySystems::getDemand_offset() const {
    return demand_offset;
}

/**
 * WSS-level infrastructure construction handler.
 * Each WSS makes infrastructure decisions based on its individual ROF and demand conditions.
 * This is the correct location for infrastructure decision-making, not at the Utility level.
 * 
 * @param long_term_rof WSS-specific risk of failure
 * @param week Current simulation week
 * @return Infrastructure ID if new infrastructure is triggered, NON_INITIALIZED otherwise
 */
int WaterSupplySystems::infrastructureConstructionHandler(double long_term_rof, int week) {
    // Calculate WSS-specific past year average demand
    double past_year_average_demand = 0;
    if (week >= (int) WEEKS_IN_YEAR) {
        for (int w = week - (int) WEEKS_IN_YEAR; w < week; ++w) {
            past_year_average_demand += demand_series_realization.at(w);
        }
    }

    long_term_risk_of_failure = long_term_rof;

    // WSS makes infrastructure decision based on its specific conditions
    int new_infra_triggered = infrastructure_construction_manager.infrastructureConstructionHandler(
            long_term_rof, week,
            past_year_average_demand,
            wss_owned_wtp_capacities,
            water_source_to_wtp,
            total_storage_capacity,
            total_available_volume,
            total_stored_volume,
            system_id);  // Pass WSS ID for pathway tracking

    // Handle bond issuance immediately here since WSS has the infrastructure sources
    if (new_infra_triggered != NON_INITIALIZED) {
        // THREAD-SAFE: Only issue bonds for realization models, not ROF models
        // ROF models should not issue actual bonds as they're just calculating risks
        if (used_for_realization && owner) {
            // Find the infrastructure source in this WSS
            WaterSource* target_source = nullptr;
            for (const auto& source : water_sources) {
                if (source && source->id == new_infra_triggered) {
                    target_source = source;
                    break;
                }
            }
            
            if (target_source) {
                // Check if bond has already been issued for this source in this WSS
                try {
                    Bond &bond = target_source->getBond(system_id);
                    if (!bond.isIssued()) {
                        // Issue bond directly on the infrastructure source using utility's method
                        // but searching only in this WSS (which has the source)
                        // Pass the realization-specific discount rate and bond multipliers
                        owner->issueBond(new_infra_triggered, week, this, infra_discount_rate,
                                        bond_term_multiplier, bond_interest_rate_multiplier);
                    } 
                } catch (const std::exception& e) {
                    printf("ERROR [WSS] Failed to access bond for infrastructure %d in WSS %d: %s\n", 
                           new_infra_triggered, system_id, e.what());
                }
            } else {
                printf("ERROR [WSS] Infrastructure %d triggered but not found in WSS %d\n", 
                       new_infra_triggered, system_id);
            }
        }
    }

    updateTreatmentAndNumberOfStorageSources();
    return new_infra_triggered;
}

/**
 * Get the ROF-based infrastructure construction order for this WSS.
 * Used by ContinuityModel to coordinate shared infrastructure decisions.
 * 
 * @return Reference to the ROF infrastructure construction order vector
 */
const vector<int>& WaterSupplySystems::getRof_infra_construction_order() const {
    return infrastructure_construction_manager.getRof_infra_construction_order();
}

/**
 * Set infrastructure parameters for this WSS after creation.
 * This recreates the infrastructure manager with the provided parameters.
 * 
 * @param rof_infra_construction_order Infrastructure construction order based on ROF
 * @param demand_infra_construction_order Infrastructure construction order based on demand 
 * @param infra_construction_triggers Infrastructure trigger thresholds
 */
void WaterSupplySystems::setInfrastructureParameters(const vector<int>& rof_infra_construction_order,
                                                     const vector<int>& demand_infra_construction_order, 
                                                     const vector<double>& infra_construction_triggers) {
    // Recreate the infrastructure manager with the provided parameters
    infrastructure_construction_manager = InfrastructureManager(
        system_id,  // Use system_id so treatment capacity lookup is correct
        infra_construction_triggers,
        vector<vector<int>>(),  
        0.0,  
        0.0,  
        0.0, 
        rof_infra_construction_order,
        demand_infra_construction_order
    );
    
    // Reconnect water sources vectors to the new infrastructure manager
    infrastructure_construction_manager.connectWaterSourcesVectorsToWSS(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Re-add all existing water sources to the new infrastructure manager
    for (WaterSource* water_source : water_sources) {
        if (water_source != nullptr) {
            infrastructure_construction_manager.addWaterSource(water_source);
        }
    }
}

// ============================================================================
// Financial Data Storage Methods (Set by Utility, no calculation logic in WSS)
// ============================================================================

void WaterSupplySystems::setWssGrossRevenue(double revenue) {
    wss_gross_revenue = revenue;
}

void WaterSupplySystems::setWssDroughtMitigationCost(double cost) {
    wss_drought_mitigation_cost = cost;
}

void WaterSupplySystems::setWssContingencyFundShare(double share) {
    wss_contingency_fund_share = share;
}

void WaterSupplySystems::setWssDebtServiceShare(double share) {
    wss_debt_service_share = share;
}

void WaterSupplySystems::setWssInfrastructureNPC(double npc) {
    // Validate input to catch corruption at source
    if (std::isnan(npc) || std::isinf(npc) || npc < -1e100 || npc > 1e100) {
        char error[512];
        sprintf(error, "Attempting to add invalid NPC value %.2e to WSS %d (system_id=%d). "
                "Previous value: %.2e. This indicates memory corruption or calculation error.",
                npc, utility_id, system_id, wss_infrastructure_npc);
        throw std::runtime_error(error);
    }
    // ACCUMULATE instead of overwrite to avoid race conditions
    wss_infrastructure_npc += npc;
}

double WaterSupplySystems::getWssGrossRevenue() const {
    return wss_gross_revenue;
}

double WaterSupplySystems::getWssDroughtMitigationCost() const {
    return wss_drought_mitigation_cost;
}

double WaterSupplySystems::getWssContingencyFundShare() const {
    return wss_contingency_fund_share;
}

double WaterSupplySystems::getWssDebtServiceShare() const {
    return wss_debt_service_share;
}

double WaterSupplySystems::getWssInfrastructureNPC() const {
    return wss_infrastructure_npc;
}

void WaterSupplySystems::setAverageMonthlyIncome(double income) {
    average_monthly_income = income;
}

double WaterSupplySystems::getAverageMonthlyIncome() const {
    return average_monthly_income;
}