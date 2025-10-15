//
// Created by bernardo on 1/13/17.
//

#include <algorithm>
#include "Utility.h"
#include "../../Utils/Utils.h"
#include "WaterSupplySystems.h"

/**
 * Main constructor for the Utility class.
 * @param name Utility name (e.g. Raleigh_water)
 * @param id Numeric ID assigned to that utility.
 * @param demands_all_realizations Text file containing utility's demand series.
 * @param number_of_week_demands Length of weeks in demand series.
 * @param typesMonthlyDemandFraction Table of size 12 (months in year) by
 * number of consumer tiers with the fraction of the total demand consumed by
 * each tier in each month of the year. The last column must be the fraction
 * of the demand treated as sewage. The summation of all number in a row but
 * the last one, therefore, must sum to 1.
 * @param typesMonthlyWaterPrice Monthly water price for each tier. The last
 * column is the price charged for sewage treatment.
 * @param wwtp_discharge_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 */

void Utility::setTotal_storage_capacity() {
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->setTotal_storage_capacity();
    }
}

void Utility::connectWaterSources(std::vector<WaterSource*>& water_sources_ref,
                                  std::vector<int>& priority_sources,
                                  std::vector<int>& non_priority_sources) {
    // Connect to infrastructure manager for infrastructure utilities
    if (!infrastructure_construction_manager.getRof_infra_construction_order().empty() ||
        !infrastructure_construction_manager.getDemand_infra_construction_order().empty()) {
        infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
                water_sources_ref, priority_sources, non_priority_sources);
    }
    
    // Connect to water supply system for storage capacity calculation
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->connectWaterSources(water_sources_ref,
                                                     priority_sources,
                                                     non_priority_sources);
        
        // Force recalculation of storage capacity
        water_supply_systems[0]->setTotal_storage_capacity();
    }
}

Utility::Utility(
        const char *name, int id,
        vector<vector<double>> &demands_all_realizations,
        int number_of_week_demands,
        const double percent_contingency_fund_contribution,
        const vector<vector<double>> &typesMonthlyDemandFraction,
        const vector<vector<double>> &typesMonthlyWaterPrice,
        WwtpDischargeRule wwtp_discharge_rule,
        double demand_buffer,
        vector<vector<int>> water_source_to_wtp,
        vector<double> utility_owned_wtp_capacities) :
        // Removed operational variables - now in WaterSupplySystems
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations),
        infra_discount_rate(NON_INITIALIZED),
        bond_term_multiplier(NON_INITIALIZED),
        bond_interest_rate_multiplier(NON_INITIALIZED),
        id(id),
        number_of_week_demands(number_of_week_demands),
        name(name),
        percent_contingency_fund_contribution(
                percent_contingency_fund_contribution),
        demand_buffer(demand_buffer),
        utility_owned_wtp_capacities(utility_owned_wtp_capacities) {
    
    // Initialize water supply systems (for now, create one default system)
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, id, id, this, demands_all_realizations, number_of_week_demands,
        typesMonthlyDemandFraction, wwtp_discharge_rule, demand_buffer,
        water_source_to_wtp, utility_owned_wtp_capacities));
    
    calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                      typesMonthlyWaterPrice);
    
    // Note: unrollWaterSourceToWtpVector is already called in WSS constructor
}

/*
// MOVED TO WaterSupplySystems CLASS
void Utility::unrollWaterSourceToWtpVector(
        const vector<vector<int>> &water_source_to_wtp,
        const vector<double> &utility_owned_wtp_capacities) {

    if (water_source_to_wtp.size() != utility_owned_wtp_capacities.size()) {
        string error = "Utility " + to_string(id) + " has " + to_string(utility_owned_wtp_capacities.size()) + " WTPs but " + to_string(water_source_to_wtp.size()) + " water sources (or groups of) assigned to WTPs.";
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
*/

/**
 * Constructor for when there is infrastructure to be built.
 * @param name Utility name (e.g. Raleigh_water)
 * @param id Numeric id assigned to that utility.
 * @param demands_all_realizations Text file containing utility's demand series.
 * @param number_of_week_demands Length of weeks in demand series.
 * @param percent_contingency_fund_contribution
 * @param typesMonthlyDemandFraction Table of size 12 (months in year) by
 * number of consumer tiers with the fraction of the total demand consumed by
 * each tier in each month of the year. The last column must be the fraction
 * of the demand treated as sewage. The summation of all number in a row but
 * the last one, therefore, must sum to 1.
 * @param typesMonthlyWaterPrice Monthly water price for each tier. The last
 * column is the price charged for waste water treatment.
 * @param wwtp_discharge_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 * @param rof_infra_construction_order
 * @param infra_discount_rate
 * @param infra_if_built_remove if infra option in position 0 of a row is
 * built, remove infra options of IDs in remaining positions of the same row.
 */
Utility::Utility(const char *name, int id,
                 vector<vector<double>> &demands_all_realizations,
                 int number_of_week_demands,
                 const double percent_contingency_fund_contribution,
                 const vector<vector<double>> &typesMonthlyDemandFraction,
                 const vector<vector<double>> &typesMonthlyWaterPrice,
                 WwtpDischargeRule wwtp_discharge_rule,
                 double demand_buffer,
                 vector<vector<int>> water_source_to_wtp,
                 vector<double> utility_owned_wtp_capacities,
                 const vector<int> &rof_infra_construction_order,
                 const vector<int> &demand_infra_construction_order,
                 const vector<double> &infra_construction_triggers,
                 double infra_discount_rate,
                 const vector<vector<int>> &infra_if_built_remove,
                 double bond_term, double bond_interest_rate) :
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations),
        infra_discount_rate(infra_discount_rate),
        bond_term_multiplier(bond_term),
        bond_interest_rate_multiplier(bond_interest_rate),
        id(id),
        number_of_week_demands(number_of_week_demands),
        name(name),
        percent_contingency_fund_contribution(
                percent_contingency_fund_contribution),
        demand_buffer(demand_buffer),
        utility_owned_wtp_capacities(utility_owned_wtp_capacities) {
    infrastructure_construction_manager =
            InfrastructureManager(id, infra_construction_triggers,
                                  infra_if_built_remove,
                                  infra_discount_rate, bond_term,
                                  bond_interest_rate,
                                  rof_infra_construction_order,
                                  demand_infra_construction_order);

    // Initialize water supply systems
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, id, id, this, wwtp_discharge_rule));
        
    // Delegate to water supply systems
    water_supply_systems[0]->unrollWaterSourceToWtpVector(water_source_to_wtp,
                                                       utility_owned_wtp_capacities);

    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Also connect water sources to the water supply system for storage capacity calculation
    water_supply_systems[0]->connectWaterSources(water_sources,
                                                 priority_draw_water_source,
                                                 non_priority_draw_water_source);

    if (rof_infra_construction_order.empty() &&
        demand_infra_construction_order.empty())
        throw std::invalid_argument("At least one infrastructure construction "
                                    "order vector  must have at least "
                                    "one water source ID. If there's "
                                    "not infrastructure to be build, "
                                    "use other constructor "
                                    "instead.");
    if (infra_discount_rate <= 0)
        throw std::invalid_argument("Infrastructure discount rate must be "
                                    "greater than 0.");

    calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                      typesMonthlyWaterPrice);
}


/**
 * Constructor for when there is infrastructure to be built.
 * @param name Utility name (e.g. Raleigh_water)
 * @param id Numeric id assigned to that utility.
 * @param demands_all_realizations Text file containing utility's demand series.
 * @param number_of_week_demands Length of weeks in demand series.
 * @param percent_contingency_fund_contribution
 * @param typesMonthlyDemandFraction Table of size 12 (months in year) by
 * number of consumer tiers with the fraction of the total demand consumed by
 * each tier in each month of the year. The last column must be the fraction
 * of the demand treated as sewage. The summation of all number in a row but
 * the last one, therefore, must sum to 1.
 * @param typesMonthlyWaterPrice Monthly water price for each tier. The last
 * column is the price charged for waste water treatment.
 * @param wwtp_discharge_rule 53 weeks long time series according to which
 * fractions of sewage is discharged in different water sources (normally one
 * for each WWTP).
 * @param rof_infra_construction_order
 * @param infra_discount_rate
 */
Utility::Utility(const char *name, int id,
                 vector<vector<double>> &demands_all_realizations,
                 int number_of_week_demands,
                 const double percent_contingency_fund_contribution,
                 const vector<vector<double>> &typesMonthlyDemandFraction,
                 const vector<vector<double>> &typesMonthlyWaterPrice,
                 WwtpDischargeRule wwtp_discharge_rule,
                 double demand_buffer,
                 vector<vector<int>> water_source_to_wtp,
                 vector<double> utility_owned_wtp_capacities,
                 const vector<int> &rof_infra_construction_order,
                 const vector<int> &demand_infra_construction_order,
                 const vector<double> &infra_construction_triggers,
                 double infra_discount_rate, double bond_term,
                 double bond_interest_rate) :
        // Removed operational variables - now in WaterSupplySystems
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations),
        infra_discount_rate(infra_discount_rate),
        bond_term_multiplier(bond_term),
        bond_interest_rate_multiplier(bond_interest_rate),
        id(id),
        number_of_week_demands(number_of_week_demands),
        name(name),
        percent_contingency_fund_contribution(
                percent_contingency_fund_contribution),
        demand_buffer(demand_buffer),
        utility_owned_wtp_capacities(utility_owned_wtp_capacities) {

    // Initialize water supply systems (create one default system)
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, id, id, this, wwtp_discharge_rule));
        
    // Delegate to water supply systems
    water_supply_systems[0]->unrollWaterSourceToWtpVector(water_source_to_wtp,
                                                          utility_owned_wtp_capacities);

    infrastructure_construction_manager = InfrastructureManager(id,
                                                                infra_construction_triggers,
                                                                vector<vector<int>>(),
                                                                infra_discount_rate,
                                                                bond_term,
                                                                bond_interest_rate,
                                                                rof_infra_construction_order,
                                                                demand_infra_construction_order);

    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Also connect water sources to the water supply system for storage capacity calculation
    water_supply_systems[0]->connectWaterSources(water_sources,
                                                 priority_draw_water_source,
                                                 non_priority_draw_water_source);

    if (rof_infra_construction_order.empty() &&
        demand_infra_construction_order.empty())
        throw std::invalid_argument("At least one infrastructure construction "
                                    "order vector must have at least "
                                    "one water source ID. If there's "
                                    "not infrastructure to be build, "
                                    "use other constructor "
                                    "instead.");
    if (infra_discount_rate <= 0)
        throw std::invalid_argument("Infrastructure discount rate must be "
                                    "greater than 0.");

    if (demands_all_realizations.empty()) {
        string error = "Empty demand vectors passed to utility " + to_string(id);
        throw std::invalid_argument(error);
    }

    calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                      typesMonthlyWaterPrice);
}

Utility::Utility(Utility &utility) :
        weekly_average_volumetric_price(
                utility.weekly_average_volumetric_price),
        // Removed operational variables - now in WaterSupplySystems
        wwtp_discharge_rule(utility.wwtp_discharge_rule),
        demands_all_realizations(utility.demands_all_realizations),
        demand_series_realization(utility.demand_series_realization),
        infra_discount_rate(utility.infra_discount_rate),
        bond_term_multiplier(utility.bond_term_multiplier),
        bond_interest_rate_multiplier(utility.bond_interest_rate_multiplier),
        id(utility.id),
        number_of_week_demands(utility.number_of_week_demands),
        name(utility.name),
        percent_contingency_fund_contribution(
                utility.percent_contingency_fund_contribution),
        demand_buffer(utility.demand_buffer),
        infrastructure_construction_manager(
                utility.infrastructure_construction_manager),
        water_source_to_wtp(
                utility.water_source_to_wtp),
        utility_owned_wtp_capacities(utility.utility_owned_wtp_capacities) {

    // Initialize water supply systems like in other constructors
    if (!utility.water_supply_systems.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            utility.name, utility.id, utility.id, this, utility.wwtp_discharge_rule));
        
        // Copy the water source connections from the original utility's water supply system
        if (!utility.water_supply_systems.empty()) {
            auto& original_wss = utility.water_supply_systems[0];
            auto& new_wss = water_supply_systems[0];
            
            // Copy water sources and connection lists from original
            new_wss->copyWaterSourceConnections(*original_wss);
        }
    }

    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Note: Don't clear water_sources here as WaterSupplySystems needs them
}

Utility::~Utility() {
    water_sources.clear();
}

Utility &Utility::operator=(const Utility &utility) {
    demand_series_realization = vector<double>(
            (unsigned long) utility.number_of_week_demands);

    // Initialize water supply systems like in other constructors
    water_supply_systems.clear();
    if (!utility.water_supply_systems.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            utility.name, utility.id, utility.id, this, utility.wwtp_discharge_rule));
        
        // Copy the water source connections from the original utility's water supply system
        if (!utility.water_supply_systems.empty()) {
            auto& original_wss = utility.water_supply_systems[0];
            auto& new_wss = water_supply_systems[0];
            
            // Copy water sources and connection lists from original
            new_wss->copyWaterSourceConnections(*original_wss);
        }
    }

    infrastructure_construction_manager.connectWaterSourcesVectorsToUtilities(
            water_sources,
            priority_draw_water_source,
            non_priority_draw_water_source);
    
    // Note: Don't clear water_sources here as WaterSupplySystems needs them

    return *this;
}

bool Utility::operator<(const Utility *other) const {
    return id < other->id;
}

bool Utility::operator>(const Utility *other) const {
    return id > other->id;
}

bool Utility::compById(Utility *a, Utility *b) {
    return a->id < b->id;
}

/*
// MOVED TO WaterSupplySystems CLASS
void Utility::updateTreatmentAndNumberOfStorageSources() {
    n_storage_sources = non_priority_draw_water_source.size();
    delete[] available_treated_flow_rate;
    available_treated_flow_rate = new double[non_priority_draw_water_source.size()];
    for (int i = 0; i < n_storage_sources; ++i) {
        auto ws = water_sources[non_priority_draw_water_source[i]];
        available_treated_flow_rate[i] = utility_owned_wtp_capacities[water_source_to_wtp[ws->id]];
        total_storage_treatment_capacity += available_treated_flow_rate[i];
    }

    total_treatment_capacity = accumulate(utility_owned_wtp_capacities.begin(),
                                          utility_owned_wtp_capacities.end(),
                                          0.);

    //TODO: IMPLEMENT QP HERE
//    P_x = new double[n_storage_sources];
//    A_x = new double[n_storage_sources];
}
*/

/**
 * Calculates average water price from consumer types and respective prices.
 * @param typesMonthlyDemandFraction
 * @param typesMonthlyWaterPrice
 */
void Utility::calculateWeeklyAverageWaterPrices(
        const vector<vector<double>> &typesMonthlyDemandFraction,
        const vector<vector<double>> &typesMonthlyWaterPrice) {
    priceCalculationErrorChecking(typesMonthlyDemandFraction,
                                  typesMonthlyWaterPrice);

    weekly_average_volumetric_price = vector<double>((int) WEEKS_IN_YEAR + 1,
                                                     0.);
    double monthly_average_price[NUMBER_OF_MONTHS] = {};
    int n_tiers = static_cast<int>(typesMonthlyWaterPrice.at(0).size());

    // Calculate monthly average prices across consumer types.
    for (int m = 0; m < NUMBER_OF_MONTHS; ++m) {
        for (int t = 0; t < n_tiers; ++t) {
            monthly_average_price[m] += typesMonthlyDemandFraction[m][t] *
                                        typesMonthlyWaterPrice[m][t];
            if (monthly_average_price[m] < 1e-6) {
                string error = "Utility " + to_string(id) + " has $0.00 water price.";
                throw runtime_error(error);
            }
        }
    }
    // Create weekly price table from monthly prices.
    bool issued_high_tariff_warning = false;
    for (int w = 0; w < (int) (WEEKS_IN_YEAR + 1); ++w) {
        weekly_average_volumetric_price[w] =
                monthly_average_price[(int) (w / WEEKS_IN_MONTH)] /
                WEEKS_IN_MONTH;
    }
}

/**
 * Checks price calculation input matrices for errors.
 * @param typesMonthlyDemandFraction
 * @param typesMonthlyWaterPrice
 */
void Utility::priceCalculationErrorChecking(
        const vector<vector<double>> &typesMonthlyDemandFraction,
        const vector<vector<double>> &typesMonthlyWaterPrice) {
    if (typesMonthlyDemandFraction.size() != NUMBER_OF_MONTHS)
        throw invalid_argument(
                "There must be 12 total_demand fractions per tier.");
    if (typesMonthlyWaterPrice.size() != NUMBER_OF_MONTHS) {
        throw invalid_argument("There must be 12 water prices per tier.");
    }
    if ((&typesMonthlyWaterPrice)[0].size() !=
        (&typesMonthlyDemandFraction)[0].size()) {
        throw invalid_argument("There must be Demand fractions and water "
                               "prices for the same number of tiers.");
    }
}

// MOVED TO WaterSupplySystems CLASS - updateTotalAvailableVolume() method commented out

// MOVED TO WaterSupplySystems CLASS - clearWaterSources() method commented out

/**
 * Connects a reservoir to the utility.
 * @param water_source
 */
void Utility::addWaterSource(WaterSource *water_source) {
    // Check if the utility has any WTPs if water sources are being added.
    if (utility_owned_wtp_capacities.empty() && water_source != nullptr) {
        throw std::invalid_argument("Utility " + std::to_string(id) +
                                    " has no WTPs, but water source " +
                                    std::to_string(water_source->id) + " is being added.");
    }

    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->addWaterSource(water_source);
    }
    
    // Handle infrastructure management at utility level
    infrastructure_construction_manager.addWaterSource(water_source);
    
    // If water source is online and the utility owns some of its installed
    // treatment capacity, add it to online lists.
    if (water_source->isOnline()) {
        // Check if we have valid mapping and capacity
        if (water_source->id < water_source_to_wtp.size()) {
            int wtp_index = water_source_to_wtp[water_source->id];
            if (wtp_index >= 0 && wtp_index < utility_owned_wtp_capacities.size() &&
                utility_owned_wtp_capacities[wtp_index] > 0) {
                
                double total_storage_capacity = getTotal_storage_capacity();
                double total_available_volume = getTotal_available_volume();
                double total_stored_volume = getTotal_stored_volume();
                
                infrastructure_construction_manager.addWaterSourceToOnlineLists(
                        water_source->id, total_storage_capacity,
                        total_available_volume, total_stored_volume);
            }
        }
    }
    
    // Update utility-level tracking
    n_sources++;
    max_capacity += water_source->getAllocatedCapacity(id);
}

void Utility::addWaterSupplySystem(const std::string& name, int system_id, int utility_id,
                                   const std::vector<std::vector<double>>& demands_all_realizations,
                                   int number_of_week_demands,
                                   const std::vector<std::vector<double>>& typesMonthlyDemandFraction,
                                   const WwtpDischargeRule& wwtp_discharge_rule,
                                   double demand_buffer,
                                   const std::vector<std::vector<int>>& water_source_to_wtp,
                                   const std::vector<double>& utility_owned_wtp_capacities) {
    // Use the simpler constructor like in lines 162 and 248, then configure the system
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, system_id, utility_id, this, wwtp_discharge_rule));
    
    // Configure the newly created system with additional parameters
    auto& new_system = water_supply_systems.back();
    new_system->unrollWaterSourceToWtpVector(water_source_to_wtp, utility_owned_wtp_capacities);
}

void Utility::checkErrorsAddWaterSourceOnline(WaterSource *water_source) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->checkErrorsAddWaterSourceOnline(water_source);
    }
}

void Utility::updateTotalAvailableVolume() {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->updateTotalAvailableVolume();
    }
}

void Utility::clearWaterSources() {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->clearWaterSources();
    }
}

void Utility::setWaterSourceOnline(unsigned int source_id, int week) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->setWaterSourceOnline(source_id, week);
    }
}

void Utility::calculateWastewater_releases(int week, double *discharges) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->calculateWastewater_releases(week, discharges);
    }
}

void Utility::updateTreatmentAndNumberOfStorageSources() {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->updateTreatmentAndNumberOfStorageSources();
    }
}

bool Utility::hasTreatmentConnected(int ws) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->hasTreatmentConnected(ws);
    }
    return false;
}

void Utility::unrollWaterSourceToWtpVector(
        const vector<vector<int>> &water_source_to_wtp,
        const vector<double>& utility_owned_wtp_capacities) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->unrollWaterSourceToWtpVector(water_source_to_wtp, utility_owned_wtp_capacities);
    }
}

vector<double> Utility::calculateWeeklyPeakingFactor(vector<double> *demands) {
    // Delegate to water supply systems - convert pointer to reference
    return WaterSupplySystems::calculateWeeklyPeakingFactor(*demands);
}

void Utility::splitDemands(
        int week, vector<vector<double>> &demands, bool apply_demand_buffer) {
    // Delegate to water supply systems (for now, use first system)
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->splitDemands(week, demands, apply_demand_buffer);
    }
}

//#pragma GCC optimize("O3")
bool Utility::idealDemandSplitUnconstrained(double *split_demands,
                                            const double *available_treated_flow_rate,
                                            double total_demand,
                                            const double *storage,
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

//#pragma GCC optimize("O3")
bool Utility::idealDemandSplitConstrained(double *split_demands,
                                          bool *over_allocated,
                                          bool *has_spare_capacity,
                                          const double *available_treated_flow_rate,
                                          double total_demand,
                                          const double *storage,
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

/**
 * Splits demands among sources. Demand is allocated so that river intakes
 * and reuse are first used to their capacity before requesting water from
 * allocations in reservoirs.
 * @param week
 */
/**
 * Update contingency fund based on regular contribution, restrictions, and
 * transfers. This function works for both sources and receivers of
 * transfers, and the transfer water prices are different than regular prices
 * for both sources and receivers. It also stores the cost of drought
 * mitigation.
 * @param unrestricted_demand
 * @param demand_multiplier
 * @param demand_offset
 * @return contingency fund contribution or draw.
 */
//#pragma GCC optimize("O3")
void Utility::updateContingencyFundAndDebtService(
        double unrestricted_demand, double demand_multiplier,
        double demand_offset, double unfulfilled_demand, int week) {
    int week_of_year = Utils::weekOfTheYear(week);
    double unrestricted_price = weekly_average_volumetric_price[week_of_year];
    double current_price;

    // Clear yearly updated data collecting variables.
    if (week_of_year == 0) {
        insurance_purchase = 0.;
    } else if (week_of_year == 1) {
        infra_net_present_cost = 0.;
        current_debt_payment = 0.;
    }

    // Set current water price, contingent on restrictions being enacted.
    if (restricted_price == NON_INITIALIZED)
        current_price = unrestricted_price;
    else
        current_price = restricted_price;

    if (current_price < unrestricted_price)
        // throw logic_error("Prices under surcharge cannot be smaller than "
        //                   "prices w/o restrictions enacted.");
        current_price = unrestricted_price; //This command is temporal, 
        // intended to solve the zero water price

    // calculate fund contributions if there were no shortage.
    double projected_fund_contribution = percent_contingency_fund_contribution *
                                         unrestricted_demand *
                                         unrestricted_price;

    // Calculate actual gross revenue.
    gross_revenue = restricted_demand * current_price;

    // Calculate losses due to restrictions and transfers.
    double lost_demand_vol_sales =
            (unrestricted_demand * (1 - demand_multiplier) +
             unfulfilled_demand);
    double revenue_losses = lost_demand_vol_sales * unrestricted_price;
    double transfer_costs = demand_offset * (offset_rate_per_volume -
                                             unrestricted_price);
    double recouped_loss_price_surcharge =
            restricted_demand * (current_price - unrestricted_price);

    // contingency fund cannot get negative.
    contingency_fund = max(contingency_fund + projected_fund_contribution -
                           revenue_losses - transfer_costs +
                           recouped_loss_price_surcharge,
                           0.0);


//    if (demand_multiplier < 1.0 && demand_offset != 0 && week > 285) {
//        int i = 0;
//    }
//    if (week > 1028) {
//        int i = 0;
//    }

    // Update variables for data collection and next iteration.
    drought_mitigation_cost = max(revenue_losses + transfer_costs -
                                  insurance_payout -
                                  recouped_loss_price_surcharge,
                                  0.0);
    fund_contribution =
            projected_fund_contribution - revenue_losses - transfer_costs +
            recouped_loss_price_surcharge;


    resetDroughtMitigationVariables();

    // Calculate current debt payment to be made on that week (if first
    // week of year), if any.
    current_debt_payment = updateCurrent_debt_payment(week);
}

void Utility::resetDroughtMitigationVariables() {
    restricted_price = NON_INITIALIZED;
    offset_rate_per_volume = NONE;
    this->demand_offset = NONE;
}


/**
 * Calculates total debt payments to be made in a week, if that's the first week
 * of the year.
 * @param week
 * @param debt_payment_streams
 * @return
 */
double Utility::updateCurrent_debt_payment(int week) {
    double updated_debt_payment = 0;

    // Checks if it's the first week of the year, when outstanding debt
    // payments should be made.
    for (Bond *bond : issued_bonds) {
        updated_debt_payment += bond->getDebtService(week);
    }

    return updated_debt_payment;
}

void Utility::issueBond(int new_infra_triggered, int week) {
    if (new_infra_triggered != NON_INITIALIZED) {
        Bond &bond = water_sources.at((unsigned long) new_infra_triggered)
                ->getBond(id);
        if (!bond.isIssued()) {
            double construction_time = water_sources
                    .at((unsigned long) new_infra_triggered)->construction_time;
            bond.issueBond(week, (int) construction_time, bond_term_multiplier,
                           bond_interest_rate_multiplier);
            issued_bonds.push_back(&bond);
            infra_net_present_cost += bond.getNetPresentValueAtIssuance(
                    infra_discount_rate, week);
        }
    }
}

void Utility::forceInfrastructureConstruction(int week,
                                              const vector<int> &new_infra_triggered) {
    // Build all triggered infrastructure
    infrastructure_construction_manager.forceInfrastructureConstruction(week,
                                                                        new_infra_triggered);

    // Issue bonds for triggered infrastructure
    auto under_construction = infrastructure_construction_manager.getUnder_construction();
    for (int ws : new_infra_triggered) {
        if (under_construction.size() > ws &&
            under_construction.at((unsigned long) ws)) {
            issueBond(ws, week);
        }
    }
}

/**
 * Check if new infrastructure is to be triggered based on long-term risk of failure and, if so, handle
 * the beginning of construction, issue corresponding bonds and update debt.
 * @param long_term_rof
 * @param week
 * @return
 */
int Utility::infrastructureConstructionHandler(double long_term_rof, int week) {
    double past_year_average_demand = 0;
    if (week >= (int) WEEKS_IN_YEAR) {
        //     past_year_average_demand =
        //            std::accumulate(demand_series_realization.begin() + week - (int) WEEKS_IN_YEAR,
        //                            demand_series_realization.begin() + week, 0.0) / WEEKS_IN_YEAR;

        for (int w = week - (int) WEEKS_IN_YEAR; w < week; ++w) {
            past_year_average_demand += demand_series_realization.at(w);
        }
    }

    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->setLong_term_risk_of_failure(long_term_rof);
    }

    // Get current values from water supply systems
    double current_storage_capacity = getTotal_storage_capacity();
    double current_available_volume = getTotal_available_volume();
    double current_stored_volume = getTotal_stored_volume();
    
    // Check if new infrastructure is to be triggered and, if so, trigger it.
    int new_infra_triggered = infrastructure_construction_manager.infrastructureConstructionHandler(
            long_term_rof, week,
            past_year_average_demand,
            utility_owned_wtp_capacities,
            water_source_to_wtp,
            current_storage_capacity,
            current_available_volume,
            current_stored_volume);

    // Issue and add bond of triggered water source to list of outstanding bonds, and update total new
    // infrastructure NPV.
    issueBond(new_infra_triggered, week);

    updateTreatmentAndNumberOfStorageSources();

    return new_infra_triggered;
}

void Utility::addInsurancePayout(double payout_value) {
    contingency_fund += payout_value;
    insurance_payout = payout_value;
}

void Utility::purchaseInsurance(double insurance_price) {
    contingency_fund -= insurance_price;
    insurance_purchase = insurance_price;
}

void
Utility::setDemand_offset(double demand_offset, double offset_rate_per_volume) {
    this->demand_offset = demand_offset;
    this->offset_rate_per_volume = offset_rate_per_volume;
}

/**
 * Get time series corresponding to realization index and eliminate reference to
 * comprehensive demand data set.
 * @param r
 */
void Utility::setRealization(unsigned long r, vector<double> &rdm_factors) {
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

    bond_interest_rate_multiplier = rdm_factors.at(1);
    bond_term_multiplier = rdm_factors.at(2);
    infra_discount_rate *= rdm_factors.at(3);

    // Set peaking demand factor.
    weekly_peaking_factor = calculateWeeklyPeakingFactor
            (&demands_all_realizations.at(r));

    price_rdm_multiplier = rdm_factors.at(4);
    for (double &awp : weekly_average_volumetric_price) {
        awp *= price_rdm_multiplier;
    }
}

//========================= GETTERS AND SETTERS =============================//

double Utility::getStorageToCapacityRatio() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getStorageToCapacityRatio();
    }
    return 0.0;
}

double Utility::getTotal_available_volume() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getTotal_available_volume();
    }
    return 0.0;
}

double Utility::getTotal_stored_volume() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getTotal_stored_volume();
    }
    return 0.0;
}

double Utility::getTotal_storage_capacity() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getTotal_storage_capacity();
    }
    return 0.0;
}

double Utility::getRisk_of_failure() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getRisk_of_failure();
    }
    return 0.0;
}

void Utility::setRisk_of_failure(double risk_of_failure) {
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->setRisk_of_failure(risk_of_failure);
    }
}

double Utility::getTotal_treatment_capacity() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getTotal_treatment_capacity();
    }
    return 0.0;
}

void Utility::setDemand_multiplier(double demand_multiplier) {
    Utility::demand_multiplier = demand_multiplier;
}

double Utility::getContingency_fund() const {
    return contingency_fund;
}

double Utility::getUnrestrictedDemand() const {
    return unrestricted_demand;
}

double Utility::getRestrictedDemand() const {
    return restricted_demand;
}

double Utility::getGrossRevenue() const {
    return gross_revenue;
}

double Utility::getDemand_multiplier() const {
    return demand_multiplier;
}

double Utility::getUnrestrictedDemand(int week) const {
    return demand_series_realization[week];
}

double Utility::getInfrastructure_net_present_cost() const {
    return infra_net_present_cost;
}

double Utility::getCurrent_debt_payment() const {
    return current_debt_payment;
}

double Utility::getCurrent_contingency_fund_contribution() const {
    return fund_contribution;
}

double Utility::getDrought_mitigation_cost() const {
    return drought_mitigation_cost;
}

double Utility::getInsurance_payout() const {
    return insurance_payout;
}

double Utility::getInsurance_purchase() const {
    return insurance_purchase;
}

const vector<int> &Utility::getRof_infrastructure_construction_order()
const {
    return infrastructure_construction_manager.getRof_infra_construction_order();
}

const vector<int> &Utility::getDemand_infra_construction_order() const {
    return infrastructure_construction_manager.getDemand_infra_construction_order();
}

vector<int> Utility::getInfrastructure_built() const {
    return infrastructure_construction_manager.getInfra_built_last_week();
}

double Utility::waterPrice(int week) {
    return weekly_average_volumetric_price[week];
}

void Utility::setRestricted_price(double restricted_price) {
    Utility::restricted_price = restricted_price * price_rdm_multiplier;
}

void Utility::setNoFinaicalCalculations() {
    used_for_realization = false;
}

double Utility::getLong_term_risk_of_failure() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getLong_term_risk_of_failure();
    }
    return NONE;
}

const vector<WaterSource *> &Utility::getWater_sources() const {
    return water_sources;
}

double Utility::getWaste_water_discharge() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getWaste_water_discharge();
    }
    return NONE;
}

void Utility::resetTotal_storage_capacity() {
    if (!water_supply_systems.empty()) {
        water_supply_systems[0]->setTotal_storage_capacity();
    }
}

double Utility::getUnfulfilled_demand() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getUnfulfilled_demand();
    }
    return NONE;
}

double Utility::getNet_stream_inflow() const {
    if (!water_supply_systems.empty()) {
        return water_supply_systems[0]->getNet_stream_inflow();
    }
    return NONE;
}

const InfrastructureManager &
Utility::getInfrastructure_construction_manager() const {
    return infrastructure_construction_manager;
}

double Utility::getDemand_offset() const {
    return demand_offset;
}

double Utility::getInfraDiscountRate() const {
    return infra_discount_rate;
}
