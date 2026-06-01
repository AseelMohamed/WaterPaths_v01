//
// Created by bernardo on 1/13/17.
//

#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include "Utility.h"
#include "WaterSupplySystems.h"
#include "../../Utils/Utils.h"

// THREAD-SAFE: Define static member for globally tracking issued bonds
std::set<std::string> Utility::globally_issued_bonds;
std::map<std::string, double> Utility::globally_issued_bond_npcs;

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
        vector<double> utility_owned_wtp_capacities,
        double infra_discount_rate,
        double bond_term,
        double bond_interest_rate) :
        // Removed operational variables - now in WaterSupplySystems
        wwtp_discharge_rule(wwtp_discharge_rule),
        demands_all_realizations(demands_all_realizations),
        infra_discount_rate(infra_discount_rate),  // Set from constructor parameter
        base_infra_discount_rate(infra_discount_rate),  // Store base rate for RDM multiplication
        bond_term_multiplier(bond_term),
        bond_interest_rate_multiplier(bond_interest_rate),
        base_bond_term_multiplier(bond_term),
        base_bond_interest_rate_multiplier(bond_interest_rate),
        id(id),
        number_of_week_demands(number_of_week_demands),
        name(name),
        percent_contingency_fund_contribution(
                percent_contingency_fund_contribution),
        demand_buffer(demand_buffer),
        utility_owned_wtp_capacities(utility_owned_wtp_capacities) {
    
    // Only initialize water supply systems if water_source_to_wtp is not empty
    // If empty, it means WSS will be added manually later
    if (!water_source_to_wtp.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            name, id, id, this, demands_all_realizations, number_of_week_demands,
            wwtp_discharge_rule, demand_buffer,
            water_source_to_wtp, utility_owned_wtp_capacities));
    }
    
    // Only calculate prices if financial parameters are provided
    // If empty, WSS-specific prices will be set later via addWaterSupplySystem
    if (!typesMonthlyDemandFraction.empty() && !typesMonthlyWaterPrice.empty()) {
        calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                          typesMonthlyWaterPrice);
    }
    
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
        base_infra_discount_rate(infra_discount_rate),
        base_bond_term_multiplier(bond_term),
        base_bond_interest_rate_multiplier(bond_interest_rate),
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

    // Initialize water supply systems with INFRASTRUCTURE PARAMETERS
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, id, id, this, demands_all_realizations, number_of_week_demands,
        wwtp_discharge_rule, demand_buffer, water_source_to_wtp, utility_owned_wtp_capacities,
        rof_infra_construction_order, demand_infra_construction_order,
        infra_construction_triggers));
        
    // unrollWaterSourceToWtpVector is already called in the WSS constructor above

    // InfrastructureManager will work through utility to coordinate with WaterSupplySystems
    // Water sources are managed by WaterSupplySystems, not directly by Utility

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

    // Only calculate prices if financial parameters are provided
    // If empty, WSS-specific prices will be set later via addWaterSupplySystem
    if (!typesMonthlyDemandFraction.empty() && !typesMonthlyWaterPrice.empty()) {
        calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                          typesMonthlyWaterPrice);
    }
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
        base_infra_discount_rate(infra_discount_rate),
        base_bond_term_multiplier(bond_term),
        base_bond_interest_rate_multiplier(bond_interest_rate),
        id(id),
        number_of_week_demands(number_of_week_demands),
        name(name),
        percent_contingency_fund_contribution(
                percent_contingency_fund_contribution),
        demand_buffer(demand_buffer),
        utility_owned_wtp_capacities(utility_owned_wtp_capacities) {

    // Initialize water supply systems (create one default system with infrastructure parameters)
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, id, id, this, demands_all_realizations, number_of_week_demands,
        wwtp_discharge_rule, demand_buffer, water_source_to_wtp, 
        utility_owned_wtp_capacities, rof_infra_construction_order,
        demand_infra_construction_order, infra_construction_triggers));
        
    // Note: unrollWaterSourceToWtpVector is now called within the WSS constructor

    infrastructure_construction_manager = InfrastructureManager(id,
                                                                infra_construction_triggers,
                                                                vector<vector<int>>(),
                                                                infra_discount_rate,
                                                                bond_term,
                                                                bond_interest_rate,
                                                                rof_infra_construction_order,
                                                                demand_infra_construction_order);

    // InfrastructureManager will coordinate through utility to work with WaterSupplySystems
    // Water sources are managed by WaterSupplySystems, not directly by Utility

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

    // Only calculate prices if financial parameters are provided
    // If empty, WSS-specific prices will be set later via addWaterSupplySystem
    if (!typesMonthlyDemandFraction.empty() && !typesMonthlyWaterPrice.empty()) {
        calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                          typesMonthlyWaterPrice);
    }
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
            utility.name, utility.id, utility.id, this, demands_all_realizations, utility.wwtp_discharge_rule));
        
        // Copy the water source connections from the original utility's water supply system
        if (!utility.water_supply_systems.empty()) {
            auto& original_wss = utility.water_supply_systems[0];
            auto& new_wss = water_supply_systems[0];
        }
    }

    // InfrastructureManager coordinates through utility to work with WaterSupplySystems
    // Water sources are managed by WaterSupplySystems, not directly by Utility
}

Utility::~Utility() {

}

Utility &Utility::operator=(const Utility &utility) {
    demand_series_realization = vector<double>(
            (unsigned long) utility.number_of_week_demands);

    weekly_average_volumetric_price = utility.weekly_average_volumetric_price;

    // Initialize water supply systems like in other constructors
    water_supply_systems.clear();
    if (!utility.water_supply_systems.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            utility.name, utility.id, utility.id, this, demands_all_realizations, utility.wwtp_discharge_rule));
        
        // Copy the water source connections from the original utility's water supply system
        if (!utility.water_supply_systems.empty()) {
            auto& original_wss = utility.water_supply_systems[0];
            auto& new_wss = water_supply_systems[0];
        }
    }

    // InfrastructureManager coordinates through utility to work with WaterSupplySystems
    // Water sources are managed by WaterSupplySystems, not directly by Utility

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

    base_weekly_average_volumetric_price = vector<double>((int) WEEKS_IN_YEAR + 1,
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
    for (int w = 0; w < (int) (WEEKS_IN_YEAR + 1); ++w) {
        int month_index = min((int) (w / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
        base_weekly_average_volumetric_price[w] =
                monthly_average_price[month_index] /
                WEEKS_IN_MONTH;
    }

    weekly_average_volumetric_price = base_weekly_average_volumetric_price;
    for (double &awp : weekly_average_volumetric_price) {
        awp *= price_rdm_multiplier;
    }
}

/**
 * Overloaded version: Calculate weekly average water prices for a specific WSS and store in output parameter.
 * @param typesMonthlyDemandFraction Demand class fractions by month
 * @param typesMonthlyWaterPrice Water prices by month and tier
 * @param output_weekly_prices Output vector to store calculated weekly prices
 */
void Utility::calculateWeeklyAverageWaterPrices(
        const vector<vector<double>> &typesMonthlyDemandFraction,
        const vector<vector<double>> &typesMonthlyWaterPrice,
        vector<double>& output_weekly_prices) {
    output_weekly_prices = vector<double>((int) WEEKS_IN_YEAR + 1, 0.);
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
    for (int w = 0; w < (int) (WEEKS_IN_YEAR + 1); ++w) {
        int month_index = min((int) (w / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
        output_weekly_prices[w] = monthly_average_price[month_index] / WEEKS_IN_MONTH;
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
    if (typesMonthlyWaterPrice[0].size() !=
        typesMonthlyDemandFraction[0].size()) {
        throw invalid_argument("There must be Demand fractions and water "
                               "prices for the same number of tiers.");
    }
}

/////////////////////////////////////////////////
// --------------- CHECK THIS --------------- //
///////////////////////////////////////////////

WaterSupplySystems& Utility::systemById(int system_id) {
    for (auto& wss : water_supply_systems) {
        if (wss->getSystemId() == system_id) {
            return *wss;
        }
    }
    throw std::out_of_range("Water supply system with ID " 
        + std::to_string(system_id) + " not found in utility " 
        + this->name);
}

const WaterSupplySystems& Utility::systemById(int system_id) const {
    for (const auto& wss : water_supply_systems) {
        if (wss->getSystemId() == system_id) {
            return *wss;
        }
    }
    throw std::out_of_range("Water supply system with ID " 
        + std::to_string(system_id) 
        + " not found in utility " 
        + this->name);
}

WaterSupplySystems& Utility::systemForSource(int source_id) {
    for (auto& wss : water_supply_systems) {
        const auto& sources = wss->getWater_sources();
        for (const auto& source : sources) {
            if (source && source->id == source_id) {
                return *wss;
            }
        }
    }
    throw std::invalid_argument("Source not found in any water supply system");
}


void Utility::addWaterSupplySystem(const std::string& name, int system_id, int utility_id,
                                   std::vector<std::vector<double>>& demands_all_realizations,
                                   int number_of_week_demands,
                                   const WwtpDischargeRule& wwtp_discharge_rule,
                                   double demand_buffer,
                                   const std::vector<std::vector<int>>& water_source_to_wtp,
                                   const std::vector<double>& utility_owned_wtp_capacities,
                                   double wss_contingency_percent,
                                   const std::vector<std::vector<double>>& wss_demand_class_fractions,
                                   const std::vector<std::vector<double>>& wss_water_prices_param) {
    // Use the main constructor with all parameters
    // Create a mutable copy of const references for the constructor
    auto water_source_to_wtp_copy = water_source_to_wtp;
    auto utility_owned_wtp_capacities_copy = utility_owned_wtp_capacities;
    WwtpDischargeRule wwtp_discharge_rule_copy(const_cast<WwtpDischargeRule&>(wwtp_discharge_rule));
    
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, system_id, utility_id, this, demands_all_realizations, 
        number_of_week_demands, wwtp_discharge_rule_copy, demand_buffer,
        water_source_to_wtp_copy, utility_owned_wtp_capacities_copy));
    
    // Reconnect the infrastructure manager after the WSS is in its final location
    // This ensures the pointer is correct after any potential moves/copies
    auto& wss = water_supply_systems.back();
    wss->reconnectInfrastructureManager();
    
    // Store WSS-specific financial parameters
    wss_contingency_percentages[system_id] = wss_contingency_percent;
    wss_demand_fractions[system_id] = wss_demand_class_fractions;
    wss_water_prices[system_id] = wss_water_prices_param;
    
    // Calculate weekly average water prices for this WSS
    priceCalculationErrorChecking(wss_demand_class_fractions, wss_water_prices_param);
    vector<double> wss_weekly_prices;
    calculateWeeklyAverageWaterPrices(wss_demand_class_fractions, wss_water_prices_param, wss_weekly_prices);
    base_wss_weekly_average_prices[system_id] = wss_weekly_prices;
    for (double &price : wss_weekly_prices) {
        price *= price_rdm_multiplier;
    }
    wss_weekly_average_prices[system_id] = wss_weekly_prices;
}

void Utility::clearWaterSupplySystems() {
    water_supply_systems.clear();
}

/////////////////////////////////////////////////
////////////////////////////////////////////////

void Utility::updateContingencyFundAndDebtService(
        double utility_unrestricted_demand, double demand_multiplier,
        double demand_offset, double utility_unfulfilled_demand, int week) {
    
    int week_of_year = Utils::weekOfTheYear(week);
    
    // Add bounds checking for weekly_average_volumetric_price
    if (week_of_year < 0 || week_of_year >= (int)weekly_average_volumetric_price.size()) {
        char error_msg[512];
        sprintf(error_msg, "Utility %d: week_of_year %d (from week %d) out of bounds for weekly_average_volumetric_price (size=%zu)", 
                id, week_of_year, week, weekly_average_volumetric_price.size());
        printf("ERROR: %s\n", error_msg);
        throw std::out_of_range(error_msg);
    }
    
    this->utility_unrestricted_demand = utility_unrestricted_demand;
    this->demand_multiplier = demand_multiplier;
    this->demand_offset = demand_offset;
    double unrestricted_price = weekly_average_volumetric_price[week_of_year];
    double current_price;

    // Clear yearly updated data collecting variables.
    if (week_of_year == 0) {
        insurance_purchase = 0.;
    } else if (week_of_year == 1) {
        // NOTE: infra_net_present_cost should NOT be reset - it accumulates throughout simulation
        // Only reset current_debt_payment which is a per-week value
        current_debt_payment = 0.;
    }

    // Set current water price, contingent on restrictions being enacted.
    if (restricted_price == NON_INITIALIZED)
        current_price = unrestricted_price;
    else
        current_price = restricted_price;

    // Validation: ensure restricted price is not lower than unrestricted price
    // In some edge cases during transition periods, this might occur due to timing
    if (current_price < unrestricted_price) {
        current_price = unrestricted_price;
        restricted_price = NON_INITIALIZED;  // Reset to prevent repeated warnings
    }

    // calculate fund contributions if there were no shortage.
    double projected_fund_contribution = percent_contingency_fund_contribution *
                          this->utility_unrestricted_demand *
                          unrestricted_price;

    // Calculate actual gross revenue using the aggregated restricted demand
    gross_revenue = this->utility_restricted_demand * current_price;    
    // Calculate losses due to restrictions and transfers.
    double lost_demand_vol_sales =
            (utility_unrestricted_demand * (1 - demand_multiplier) +
             utility_unfulfilled_demand);
    double revenue_losses = lost_demand_vol_sales * unrestricted_price;
    double transfer_costs = demand_offset * (offset_rate_per_volume -
                                             unrestricted_price);
    double recouped_loss_price_surcharge =
            utility_restricted_demand * (current_price - unrestricted_price);

    

    // Update variables for data collection and next iteration.
    drought_mitigation_cost = max(revenue_losses + transfer_costs -
                                  insurance_payout -
                                  recouped_loss_price_surcharge-
                                  contingency_fund - projected_fund_contribution,
                                  0.0);
    
    // contingency fund cannot get negative.
    contingency_fund = max(contingency_fund + projected_fund_contribution -
                           revenue_losses - transfer_costs +
                           recouped_loss_price_surcharge +
                           insurance_payout,
                           0.0);

    fund_contribution = projected_fund_contribution - 
                          revenue_losses - transfer_costs +
                          recouped_loss_price_surcharge +
                          insurance_payout;


    resetDroughtMitigationVariables();

    // NOTE: current_debt_payment calculation moved to updateUtilityFinancialCalculations()
    // to avoid duplicate calls that would exhaust bond payments prematurely.
    // The updateCurrent_debt_payment(week) is not called here!
    
}

/**
 * REFACTORED: Calculate financial data at WSS level first, then aggregate to utility level.
 * This method:
 * 1. Loops through each WSS to calculate WSS-specific revenues and costs
 * 2. Aggregates totals for utility-level debt and contingency fund calculations
 */
void Utility::updateUtilityFinancialCalculations(int week, const std::vector<WaterSupplySystems*>& realization_wss) {
    // Early exit if this is not for realization (e.g., ROF model)
    if (!used_for_realization) {
        return;
    }
    
    int week_of_year = Utils::weekOfTheYear(week);
    
    // ========================================================================
    // STEP 1: Calculate financial data FOR EACH WSS using WSS-specific prices
    // ========================================================================
    
    double total_unrestricted_demand = 0.0;
    double total_restricted_demand = 0.0;
    double total_unfulfilled_demand = 0.0;
    double total_wss_gross_revenue = 0.0;
    double total_wss_drought_cost = 0.0;
    double avg_demand_multiplier = 0.0;
    double total_demand_offset = 0.0;
    // Track aggregated losses for contingency fund (instead of recalculating with averages)
    double total_revenue_losses = 0.0;
    double total_transfer_costs = 0.0;
    double total_recouped = 0.0;
    double total_transfer_volume = 0.0;
    double weighted_transfer_price_sum = 0.0;

    // Use passed realization WSS (thread-safe - each realization has its own copies)
    size_t num_wss = realization_wss.size();
    
    // Sanity check - ensure we have WSS
    if (num_wss == 0) {
        printf("WARNING [Financial] Utility %d has zero WSS!\n", id);
        return;
    }
    
    for (size_t i = 0; i < num_wss; ++i) {
        if (i >= realization_wss.size()) {
            char error[256];
            sprintf(error, "WSS index %zu out of bounds (size=%zu) for utility %d",
                    i, realization_wss.size(), id);
            throw std::out_of_range(error);
        }
        
        WaterSupplySystems* wss = realization_wss[i];
        
        if (wss == nullptr) {
            printf("ERROR [Financial] WSS %zu is null!\n", i);
            continue;
        }
        
        int system_id = wss->getSystemId();
        
        // Get WSS-specific price (lookup from stored parameters)
        double wss_unrestricted_price = 0.0;
        double wss_current_price = 0.0;
        
        if (wss_weekly_average_prices.find(system_id) != wss_weekly_average_prices.end()) {
            auto& wss_prices = wss_weekly_average_prices[system_id];
            if (week_of_year >= 0 && week_of_year < (int)wss_prices.size()) {
                wss_unrestricted_price = wss_prices[week_of_year];

                // Look up restricted price if a surcharge is currently active for this WSS
                double fetched_restricted_price = NON_INITIALIZED;
                #pragma omp critical(utility_restriction_prices)
                {
                    auto it = wss_restricted_prices.find(system_id);
                    if (it != wss_restricted_prices.end()) {
                        fetched_restricted_price = it->second;
                    }
                }

                // Use restricted price when active, otherwise fall back to unrestricted price
                wss_current_price = (fetched_restricted_price == NON_INITIALIZED)
                        ? wss_unrestricted_price : fetched_restricted_price;

                // Ensure restricted price never falls below unrestricted price
                if (wss_current_price < wss_unrestricted_price) {
                    wss_current_price = wss_unrestricted_price;
                }
            } else {
                printf("WARNING: week_of_year %d out of bounds for WSS %d prices\n", week_of_year, system_id);
                continue;
            }
        } else {
            printf("WARNING: No prices stored for WSS %d\n", system_id);
            continue;
        }

        // Calculate current residential first-tier weekly price for affordability
        double wss_residential_base_price = 0.0;
        auto wss_monthly_price_it = wss_water_prices.find(system_id);
        if (wss_monthly_price_it != wss_water_prices.end()) {
            const auto& monthly_prices = wss_monthly_price_it->second;
            int month_index = min((int) (week_of_year / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
            if (month_index >= 0 && month_index < (int)monthly_prices.size() &&
                !monthly_prices[month_index].empty()) {
                wss_residential_base_price = monthly_prices[month_index][0] / WEEKS_IN_MONTH * price_rdm_multiplier;
            }
        }

        double wss_residential_current_price = wss_residential_base_price;
        #pragma omp critical(utility_restriction_prices)
        {
            auto wss_res_restricted_it = wss_restricted_residential_prices.find(system_id);
            if (wss_res_restricted_it != wss_restricted_residential_prices.end()) {
                wss_residential_current_price = wss_res_restricted_it->second;
            }
        }

        if (wss_residential_current_price < wss_residential_base_price) {
            wss_residential_current_price = wss_residential_base_price;
        }
        wss->setWssResidentialPrice(wss_residential_current_price);
        
        // Get WSS operational data
        double wss_unrestricted = wss->getUnrestrictedDemand(-1);
        double wss_restricted = wss->getRestrictedDemand();
        double wss_unfulfilled = wss->getUnfulfilled_demand();
        double wss_demand_mult = wss->getDemand_multiplier();
        double wss_offset = wss->getDemand_offset();
        double wss_offset_rate = wss->getOffset_rate_per_volume();
        double wss_transfer_volume = max(wss_offset, 0.0);
        
        // Revenue is based on what customers are billed (unrestricted × multiplier), NOT on
        // wss_restricted (= unrestricted × mult − demand_offset).  demand_offset only controls
        // which physical source supplies the water; it does not change what customers pay.
        // Using wss_restricted would inflate the sender's revenue and deflate the receiver's
        // revenue during intra-utility transfers, which is incorrect.
        double customer_demand = wss_unrestricted * wss_demand_mult;

        // Calculate WSS-specific gross revenue using WSS-specific price
        double wss_gross_revenue = customer_demand * wss_current_price;
        
        // Calculate WSS-specific drought mitigation costs using WSS-specific price
        double lost_demand_vol_sales = (wss_unrestricted * (1 - wss_demand_mult) + wss_unfulfilled);
        double revenue_losses = lost_demand_vol_sales * wss_unrestricted_price;
        double transfer_costs = wss_transfer_volume * (wss_offset_rate - wss_unrestricted_price);
        // Surcharge revenue is also based on the customer-facing demand volume.
        double recouped_loss_price_surcharge = customer_demand * (wss_current_price - wss_unrestricted_price);
        double wss_drought_cost = max(revenue_losses + transfer_costs - recouped_loss_price_surcharge, 0.0);

        double wss_contingency_pct = 0.0;
        if (wss_contingency_percentages.find(system_id) != wss_contingency_percentages.end()) {
            wss_contingency_pct = wss_contingency_percentages[system_id];
        }
        double wss_projected_fund_contribution = wss_contingency_pct * wss_unrestricted * wss_unrestricted_price;

        
        // Store WSS-level financial data (for data collection and objective calculations)
        wss->setWssGrossRevenue(wss_gross_revenue);
        // wss->setWssDroughtMitigationCost(wss_drought_cost);
        
        // Aggregate for utility-level calculations
        total_unrestricted_demand += wss_unrestricted;
        total_restricted_demand += wss_restricted;
        total_unfulfilled_demand += wss_unfulfilled;
        total_wss_gross_revenue += wss_gross_revenue;
        // total_wss_drought_cost += wss_drought_cost;
        avg_demand_multiplier += wss_demand_mult;
        total_demand_offset += wss_offset;
        // Aggregate actual losses (don't recalculate with averages later)
        total_revenue_losses += revenue_losses;
        total_transfer_costs += transfer_costs;
        total_recouped += recouped_loss_price_surcharge;
        total_transfer_volume += wss_transfer_volume;
        weighted_transfer_price_sum += wss_transfer_volume * wss_offset_rate;
    }

    
    
    // Calculate average demand multiplier
    if (num_wss > 0) {
        avg_demand_multiplier /= num_wss;
    } else {
        avg_demand_multiplier = 1.0;
    }
    
    // ========================================================================
    // STEP 2: Calculate utility-level financials (debt, contingency fund)
    // ========================================================================
    
    // Store aggregated values at utility level
    utility_unrestricted_demand = total_unrestricted_demand;
    utility_restricted_demand = total_restricted_demand;
    gross_revenue = total_wss_gross_revenue;
    
    // Clear yearly updated data collecting variables
    if (week_of_year == 0) {
        insurance_purchase = 0.0;
    } else if (week_of_year == 1) {
        // NOTE: infra_net_present_cost should NOT be reset - it accumulates throughout simulation
        // Only reset current_debt_payment which is a per-week value
        current_debt_payment = 0.0;
    }
    
    // Calculate WSS-specific contingency fund contributions and aggregate
    double total_projected_fund_contribution = 0.0;
    
    for (size_t i = 0; i < num_wss; ++i) {
        WaterSupplySystems* wss = realization_wss[i];
        if (wss == nullptr) continue;
        
        int system_id = wss->getSystemId();
        double wss_unrestricted = wss->getUnrestrictedDemand(-1);
        
        // Get WSS-specific contingency percentage and price
        double wss_contingency_pct = 0.0;
        double wss_price = 0.0;
        
        if (wss_contingency_percentages.find(system_id) != wss_contingency_percentages.end()) {
            wss_contingency_pct = wss_contingency_percentages[system_id];
        }
        
        if (wss_weekly_average_prices.find(system_id) != wss_weekly_average_prices.end()) {
            auto& wss_prices = wss_weekly_average_prices[system_id];
            if (week_of_year >= 0 && week_of_year < (int)wss_prices.size()) {
                wss_price = wss_prices[week_of_year];
            }
        }
        
        // Calculate WSS-specific fund contribution
        double wss_fund_contribution = wss_contingency_pct * wss_unrestricted * wss_price;
        total_projected_fund_contribution += wss_fund_contribution;
        
        // Store WSS-specific contingency fund share (for data collection)
        wss->setWssContingencyFundShare(wss_fund_contribution);
    }
    
    // Use aggregated losses from WSS calculations (already computed above)
    // Do NOT recalculate with averaged demand_multiplier - that causes incorrect results
    // The total_revenue_losses, total_transfer_costs, and total_recouped were summed from each WSS
    
    // Update drought mitigation cost at UTILITY level using aggregated variables
    // DMC_i = RL_i + TC_i - INP_i - RV_{sur,i} - CF_{i-1} - CF_{IN,i}
    double contingency_fund_before = contingency_fund;
    drought_mitigation_cost = max(total_revenue_losses + total_transfer_costs -
                          insurance_payout - total_recouped -
                          contingency_fund_before - total_projected_fund_contribution, 0.0);

    // Update contingency fund (utility-wide) AFTER drought mitigation cost
    contingency_fund = max(contingency_fund + total_projected_fund_contribution -
                          total_revenue_losses - total_transfer_costs + total_recouped, 0.0);
    
    // Fund contribution calculation
    fund_contribution = total_projected_fund_contribution - total_revenue_losses - 
                       total_transfer_costs + total_recouped;

    // // Update drought mitigation cost (aggregated from WSS, subtract insurance at utility level)
    // drought_mitigation_cost = max(total_wss_drought_cost - insurance_payout, 0.0);
    
    // Calculate current debt payment
    current_debt_payment = updateCurrent_debt_payment(week, realization_wss);
    
}

void Utility::resetDroughtMitigationVariables() {
    #pragma omp critical(utility_restriction_prices)
    {
        restricted_price = NON_INITIALIZED;
        offset_rate_per_volume = NONE;
        this->demand_offset = NONE;
        insurance_payout = 0.0;  // Reset payout so it doesn't persist across weeks
        wss_restricted_prices.clear();
        wss_restricted_residential_prices.clear();
    }
}

/**
 * Calculates total debt payments to be made in a week, if that's the first week
 * of the year.
 * @param week
 * @param current_wss The water supply systems for the current model (realization or ROF)
 * @return Total debt payment for this week
 */
double Utility::updateCurrent_debt_payment(int week, const std::vector<WaterSupplySystems*>& current_wss) {
    double updated_debt_payment = 0;

    // THREAD-SAFE FIX: Instead of using global issued_bonds list (which contains bonds from ALL models),
    // iterate through bonds from the CURRENT model's water sources only.
    // This prevents accessing deleted bonds when parallel threads destroy their ContinuityModel instances.
    for (const auto& wss : current_wss) {
        if (wss == nullptr) continue;
        
        const auto& sources = wss->getWater_sources();
        for (const auto& source : sources) {
            if (source == nullptr) continue;
            
            // Access bonds from the water source using getter method
            const auto& source_bonds = source->getBonds();
            for (Bond* bond : source_bonds) {
                if (bond != nullptr) {
                    double debt_service = bond->getDebtService(week);
                    updated_debt_payment += debt_service;
                }
            }
        }
    }
    
    return updated_debt_payment;
}


void Utility::issueBond(int new_infra_triggered, int week) {
    // THREAD-SAFE: Use critical section to prevent race conditions on bond issuance
    #pragma omp critical(bond_issuance)
    {   
        if (new_infra_triggered != NON_INITIALIZED) {
            // Try to use the non-owning references first (points to realization WSS with actual water sources)
            // Fall back to owned WSS if refs are empty
            vector<WaterSupplySystems*> wss_to_search;
            
            if (!water_supply_systems_refs.empty()) {
                wss_to_search = water_supply_systems_refs;
            } else {
                for (const auto& wss : water_supply_systems) {
                    wss_to_search.push_back(wss.get());
                }
            }
            
            // Find the water source in the WSS and remember the owning WSS so we
            // can request the bond corresponding to that system (not the utility)
            WaterSource* target_source = nullptr;
            WaterSupplySystems* owner_wss = nullptr;
            for (const auto& wss : wss_to_search) {
                if (wss == nullptr) continue;
                const auto& sources = wss->getWater_sources();
                for (const auto& source : sources) {
                    if (source && source->id == new_infra_triggered) {
                        target_source = source;
                        owner_wss = wss;
                        break;
                    }
                }
                if (target_source) break;
            }

            if (target_source && owner_wss) {
                // In the WSS design multiple systems live under the same Utility
                // (same utility id). The water source holds a vector of bonds
                // indexed by system id (WSS id). Use the owner's system id when
                // available; fall back to utility id for legacy cases.
                int bond_index = owner_wss->getSystemId();
                
                // Create unique key for global tracking: "utility_id:source_id:wss_id"
                char bond_key[64];
                sprintf(bond_key, "%d:%d:%d", id, new_infra_triggered, bond_index);
                std::string bond_key_str(bond_key);
                
                // Check if this bond has already been issued globally
                bool should_skip = (globally_issued_bonds.find(bond_key_str) != globally_issued_bonds.end());
                
                if (should_skip) {                    
                    // THREAD-SAFE: Even though bond was already issued, this realization needs its cost recorded
                    if (used_for_realization && current_realization_id != (unsigned long)NON_INITIALIZED) {
                        auto npc_it = globally_issued_bond_npcs.find(bond_key_str);
                        if (npc_it != globally_issued_bond_npcs.end()) {
                            double npc = npc_it->second;
                            realization_infra_costs[current_realization_id] += npc;
                        }
                    }
                } else {
                    try {
                        Bond &bond = target_source->getBond(bond_index);
                        if (!bond.isIssued()) {
                        double construction_time = target_source->construction_time;
                        
                        bond.issueBond(week, (int) construction_time, bond_term_multiplier,
                                       bond_interest_rate_multiplier);
                        
                        // Add bond to issued_bonds so debt service can be calculated
                        issued_bonds.push_back(&bond);
                        
                        // Calculate NPC first
                        double npc = bond.getNetPresentValueAtIssuance(infra_discount_rate, week);
                        
                        // Mark this bond as globally issued and store its NPC
                        globally_issued_bonds.insert(bond_key_str);
                        globally_issued_bond_npcs[bond_key_str] = npc;
                        
                        // THREAD-SAFE: Only accumulate costs for realization models
                        if (used_for_realization) {
                            infra_net_present_cost += npc;
                            
                            // THREAD-SAFE: Also track per-realization cost
                            if (current_realization_id != (unsigned long)NON_INITIALIZED) {
                                realization_infra_costs[current_realization_id] += npc;
                            }
                        }
                        
                        // NEW: Track infrastructure NPC at WSS level
                        if (owner_wss) {
                            double current_wss_npc = owner_wss->getWssInfrastructureNPC();
                            owner_wss->setWssInfrastructureNPC(current_wss_npc + npc);
                        }
                        } else {
                            // Still mark it globally to prevent future attempts
                            globally_issued_bonds.insert(bond_key_str);
                        }
                    } catch (const std::exception& e) {
                        printf("ERROR [WSS] Failed to get bond for source %d with index %d: %s\n", 
                               new_infra_triggered, bond_index, e.what());
                    }
                }
            } else {
                printf("ERROR [WSS] Could not find water source %d in any WSS\n", new_infra_triggered);
            }
        } 
    }
}

// Overloaded version to issue bond for infrastructure in specific WSS
void Utility::issueBond(int new_infra_triggered, int week, WaterSupplySystems* target_wss, double discount_rate,
                        double bond_term_mult, double bond_interest_rate_mult) {
    // THREAD-SAFE: Use critical section to prevent race conditions on bond issuance
    #pragma omp critical(bond_issuance)
    {
        if (new_infra_triggered != NON_INITIALIZED && target_wss != nullptr) {
            // Create unique key for global tracking: "utility_id:source_id:wss_id:realization_id"
            // Including realization_id allows each realization to issue bonds independently
            char bond_key[64];
            sprintf(bond_key, "%d:%d:%d:%lu", id, new_infra_triggered, target_wss->getSystemId(), 
                    target_wss->getCurrentRealizationId());
            std::string bond_key_str(bond_key);
            
            // Check if this bond has already been issued globally
            bool should_skip = (globally_issued_bonds.find(bond_key_str) != globally_issued_bonds.end());
            
            if (should_skip) {
                // THREAD-SAFE: Even though bond was already issued, this realization needs its cost recorded
                if (used_for_realization && current_realization_id != (unsigned long)NON_INITIALIZED) {
                    auto npc_it = globally_issued_bond_npcs.find(bond_key_str);
                    if (npc_it != globally_issued_bond_npcs.end()) {
                        double npc = npc_it->second;
                        realization_infra_costs[current_realization_id] += npc;
                        
                        // NEW: Track infrastructure NPC at WSS level
                        if (target_wss) {
                            double current_wss_npc = target_wss->getWssInfrastructureNPC();
                            target_wss->setWssInfrastructureNPC(current_wss_npc + npc);
                        }
                    }
                }
            } else {
                // Look for the infrastructure source only in the target WSS
                const auto& water_sources = target_wss->getWater_sources();
                WaterSource* target_source = nullptr;
            
            for (const auto& source : water_sources) {
                if (source && source->id == new_infra_triggered) {
                    target_source = source;
                    break;
                }
            }
            
            if (target_source) {
                // Issue bond from thread-local source but don't add pointer to issued_bonds
                // This avoids the problem of persistent storage while still tracking costs
                try {
                    Bond &bond = target_source->getBond(target_wss->getSystemId());
                    if (!bond.isIssued()) {
                        double construction_time = target_source->construction_time;
                        
                        // Use passed multipliers (realization-specific from WSS) instead of member variables (shared)
                        bond.issueBond(week, (int) construction_time, bond_term_mult,
                                       bond_interest_rate_mult);
                        
                        // Add bond to issued_bonds so debt service can be calculated
                        issued_bonds.push_back(&bond);
                        
                        // Calculate NPC using realization-specific discount rate passed from WSS
                        double npc = bond.getNetPresentValueAtIssuance(discount_rate, week);
                        // Mark this bond as issued for this realization and store its NPC
                        globally_issued_bonds.insert(bond_key_str);
                        globally_issued_bond_npcs[bond_key_str] = npc;
                        // But note: this will have race conditions if used across realizations
                        if (used_for_realization) {
                            double old_util_npc = infra_net_present_cost;
                            infra_net_present_cost += npc;
                            // THREAD-SAFE: Track per-realization cost using WSS's realization ID (not utility's)
                            // Use target_wss->getCurrentRealizationId() because utility's current_realization_id
                            // is shared and gets overwritten by multiple threads
                            unsigned long wss_realization_id = target_wss->getCurrentRealizationId();
                            if (wss_realization_id != (unsigned long)NON_INITIALIZED) {
                                double old_real_npc = realization_infra_costs[wss_realization_id];
                                realization_infra_costs[wss_realization_id] += npc;                            }
                        }
                        // Track infrastructure NPC at WSS level for objective calculations
                        if (target_wss) {
                            target_wss->setWssInfrastructureNPC(npc);
                        }
                    } else {
                        // Still mark it globally to prevent future attempts
                        globally_issued_bonds.insert(bond_key_str);
                    }
                } catch (const std::exception& e) {
                    printf("ERROR [WSS] Failed to get bond for source %d in WSS %d: %s\n", 
                           new_infra_triggered, target_wss->getSystemId(), e.what());
                }
            }
            }  
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

void Utility::purchaseInsurance(double insurance_price) {
    contingency_fund -= insurance_price;
    insurance_purchase = insurance_price;
}

void Utility::addInsurancePayout(double payout_value) {
    contingency_fund += payout_value;
    insurance_payout = payout_value;
}

void Utility::setDemand_offset(double demand_offset, double offset_rate_per_volume) {
    this->demand_offset = demand_offset;
    this->offset_rate_per_volume = offset_rate_per_volume;
}

/**
 * Get time series corresponding to realization index and eliminate reference to
 * comprehensive demand data set.
 * @param r
 */
void Utility::setRealization(unsigned long r, vector<double> &rdm_factors) {
    // THREAD-SAFE: DO NOT clear global bond tracking here!
    // globally_issued_bonds and globally_issued_bond_npcs are shared across ALL realizations
    // to prevent duplicate bond issuance. Clearing them here causes race conditions where
    // different threads see inconsistent state.
    //
    // These globals should only be cleared ONCE before running all realizations (in Problem),
    // not per-realization. Each realization will check if a bond was globally issued and either:
    // 1) Issue it for the first time and mark it globally
    // 2) Reuse the NPC from the global map if already issued
    //
    // issued_bonds can still be cleared per-utility since each thread gets its own utility objects
    issued_bonds.clear();
    
    // Prefer operating on non-owning WSS references that point to realization-specific
    // WaterSupplySystems copies (populated by Simulation::createContinuityModels).
    // Fall back to owned water_supply_systems if refs are not set.
    if (!water_supply_systems_refs.empty()) {
        // Call setRealization on each WSS - they will check their own used_for_realization flag
        for (auto* wss : water_supply_systems_refs) {
            if (wss) wss->setRealization(r, rdm_factors);
        }

        // Aggregate demand series realization from referenced WSS
        utility_demand_series_realization = vector<double>(number_of_week_demands, 0.0);
        for (int week = 0; week < number_of_week_demands; ++week) {
            for (const auto* wss : water_supply_systems_refs) {
                if (wss) utility_demand_series_realization[week] += wss->getUnrestrictedDemand(week);
            }
        }
    } else {
        // Set realization for owned water supply systems
        // Each WSS will check its own used_for_realization flag
        for (auto& wss : water_supply_systems) {
            wss->setRealization(r, rdm_factors);
        }

        // Then aggregate demand series realization from owned WSS at utility level
        utility_demand_series_realization = vector<double>(number_of_week_demands, 0.0);
        for (int week = 0; week < number_of_week_demands; ++week) {
            for (const auto& wss : water_supply_systems) {
                utility_demand_series_realization[week] += wss->getUnrestrictedDemand(week);
            }
        }
    }
    
    // THREAD-SAFE: Track current realization for per-realization cost tracking
    current_realization_id = r;
    
    // Initialize this realization's cost tracking if not already present
    if (realization_infra_costs.find(r) == realization_infra_costs.end()) {
        realization_infra_costs[r] = 0.0;
    }
    
    // Reset all financial state variables to prevent accumulation across realizations
    // Since utilities are shared across threads, each realization must start with clean state
    gross_revenue = 0.0;
    current_debt_payment = 0.0;
    fund_contribution = 0.0;
    contingency_fund = 0.0;
    drought_mitigation_cost = 0.0;
    insurance_payout = 0.0;
    insurance_purchase = 0.0;
    infra_net_present_cost = 0.0;
    
    bond_interest_rate_multiplier = rdm_factors.at(1);
    bond_term_multiplier = rdm_factors.at(2);
    
    // THREAD-SAFE FIX: Apply RDM multiplier to BASE discount rate (not accumulated)
    // This prevents accumulation when setRealization is called multiple times on shared utility
    infra_discount_rate = base_infra_discount_rate * rdm_factors.at(3);
    price_rdm_multiplier = rdm_factors.at(4);

    if (!base_weekly_average_volumetric_price.empty()) {
        weekly_average_volumetric_price = base_weekly_average_volumetric_price;
        for (double &awp : weekly_average_volumetric_price) {
            awp *= price_rdm_multiplier;
        }
    }

    if (!base_wss_weekly_average_prices.empty()) {
        wss_weekly_average_prices.clear();
        for (const auto &pair : base_wss_weekly_average_prices) {
            vector<double> scaled_prices = pair.second;
            for (double &price : scaled_prices) {
                price *= price_rdm_multiplier;
            }
            wss_weekly_average_prices[pair.first] = std::move(scaled_prices);
        }
    }

}

void Utility::clearGlobalBondTracking() {
    // Clear global bond tracking before running realizations
    // This must be called ONCE before the parallel region starts, not per-realization
    globally_issued_bonds.clear();
    globally_issued_bond_npcs.clear();
}

void Utility::resetFinancialState() {
    // Reset all financial state variables to zero
    // Called at the start of each realization in ordered fashion to prevent race conditions
    gross_revenue = 0.0;
    current_debt_payment = 0.0;
    fund_contribution = 0.0;
    contingency_fund = 0.0;
    drought_mitigation_cost = 0.0;
    insurance_payout = 0.0;
    insurance_purchase = 0.0;
    infra_net_present_cost = 0.0;
    utility_unrestricted_demand = 0.0;
    utility_restricted_demand = 0.0;
}

//========================= GETTERS AND SETTERS =============================//

double Utility::getTotal_treatment_capacity() const {
    double utility_total_treatment_capacity = 0.0;
        for (const auto& wss : water_supply_systems) {
            utility_total_treatment_capacity += wss->getTotal_treatment_capacity();
    }
    return utility_total_treatment_capacity;
}

double Utility::getContingency_fund() const {
    return contingency_fund;
}

double Utility::getUnrestrictedDemand(int week) const {
    if (week == -1) {
        // no-parameter logic
        double utility_unrestricted_demand = 0.0;
            for (const auto& wss : water_supply_systems) {
            double wss_unrestricted_demand = wss->getUnrestrictedDemand();
            utility_unrestricted_demand += wss_unrestricted_demand;
        }
        return utility_unrestricted_demand;

    } else {
        // week-specific logic - aggregate from all water supply systems
        return utility_demand_series_realization[week];
    }
}

double Utility::getRestrictedDemand() const {
    double utility_restricted_demand = 0.0;
    for (const auto& wss : water_supply_systems) {
        double wss_restricted_demand = wss->getRestrictedDemand();
        utility_restricted_demand += wss_restricted_demand;
    }
    return utility_restricted_demand;
}

double Utility::getGrossRevenue() const {
    return gross_revenue;
}

double Utility::getInfrastructure_net_present_cost() const {
    // Return the main accumulated infrastructure net present cost
    // This value is accumulated whenever bonds are issued (in used_for_realization mode)
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
    // Return empty - pathway collection not supported in WSS model
    return vector<int>();
}

vector<vector<int>> Utility::getAllAndClearInfraBuilt() const {
    return const_cast<InfrastructureManager&>(infrastructure_construction_manager).getAllAndClearInfraBuilt();
}

double Utility::waterPrice(int week) const {
    int week_of_year = Utils::weekOfTheYear(week);
    
    // If utility-level price exists, use it (legacy single-WSS utilities)
    if (!weekly_average_volumetric_price.empty() && week_of_year >= 0 && 
        week_of_year < (int)weekly_average_volumetric_price.size()) {
        return weekly_average_volumetric_price[week_of_year];
    }
    
    // Otherwise, calculate average price across all WSS (multi-WSS utilities)
    if (!wss_weekly_average_prices.empty()) {
        double total_price = 0.0;
        int count = 0;
        
        for (const auto& wss_price_pair : wss_weekly_average_prices) {
            const auto& wss_prices = wss_price_pair.second;
            if (week_of_year >= 0 && week_of_year < (int)wss_prices.size()) {
                total_price += wss_prices[week_of_year];
                count++;
            }
        }
        
        if (count > 0) {
            return total_price / count;
        }
    }
    
    // Fallback: return a default price or throw an error
    char error_msg[256];
    sprintf(error_msg, "Utility %d: No water price available for week %d (week_of_year %d)", 
            id, week, week_of_year);
    throw std::runtime_error(error_msg);
}

double Utility::getCurrentWaterPrice(int week) const {
    int week_of_year = Utils::weekOfTheYear(week);

    // If WSS prices exist, sum current prices across all WSS in this utility
    if (!wss_weekly_average_prices.empty()) {
        double total_current_price = 0.0;

        for (const auto& wss_price_pair : wss_weekly_average_prices) {
            int system_id = wss_price_pair.first;
            const auto& wss_prices = wss_price_pair.second;

            if (week_of_year < 0 || week_of_year >= (int)wss_prices.size()) {
                continue;
            }

            double wss_unrestricted_price = wss_prices[week_of_year];
            double wss_restricted_price = NON_INITIALIZED;
            double local_restricted_price = NON_INITIALIZED;
            #pragma omp critical(utility_restriction_prices)
            {
                auto wss_restricted_it = wss_restricted_prices.find(system_id);
                if (wss_restricted_it != wss_restricted_prices.end()) {
                    wss_restricted_price = wss_restricted_it->second;
                }
                local_restricted_price = restricted_price;
            }

            double wss_current_price = (wss_restricted_price == NON_INITIALIZED)
                    ? ((local_restricted_price == NON_INITIALIZED)
                       ? wss_unrestricted_price : local_restricted_price)
                    : wss_restricted_price;

            if (wss_current_price < wss_unrestricted_price) {
                wss_current_price = wss_unrestricted_price;
            }

            total_current_price += wss_current_price;
        }

        return total_current_price;
    }

    // Fallback: single utility price behavior
    double unrestricted_price = waterPrice(week);
    double local_restricted_price = NON_INITIALIZED;
    #pragma omp critical(utility_restriction_prices)
    {
        local_restricted_price = restricted_price;
    }
    double current_price = (local_restricted_price == NON_INITIALIZED)
            ? unrestricted_price : local_restricted_price;

    if (current_price < unrestricted_price) {
        current_price = unrestricted_price;
    }

    return current_price;
}

void Utility::setRestricted_price(double restricted_price) {
    #pragma omp critical(utility_restriction_prices)
    {
        Utility::restricted_price = restricted_price * price_rdm_multiplier;
    }
}

void Utility::setNoFinancialCalculations() {
    used_for_realization = false;
    // Clear any bonds that may have been issued during ROF table generation
    // ROF models should not persist bonds into realization runs
    issued_bonds.clear();
}

double Utility::getUnfulfilled_demand() const {
    double utility_unfulfilled_demand = 0.0;
        for (const auto& wss : water_supply_systems) {
            utility_unfulfilled_demand += wss->getUnfulfilled_demand();
    }
    return utility_unfulfilled_demand;
}

const InfrastructureManager &Utility::getInfrastructure_construction_manager() const {
    return infrastructure_construction_manager;
}

double Utility::getInfraDiscountRate() const {
    return infra_discount_rate;
}

double Utility::getBaseInfraDiscountRate() const {
    return base_infra_discount_rate;
}

double Utility::getTotal_storage_capacity() const {
    double utility_total_storage_capacity = 0.0;
        for (const auto& wss : water_supply_systems) {
            utility_total_storage_capacity += wss->getTotal_storage_capacity();
    }
    return utility_total_storage_capacity;
}

double Utility::getTotal_available_volume() const {
    double utility_total_available_volume = 0.0;
    
    // Force WSS to update their volumes before aggregating
        for (auto& wss : water_supply_systems) {
            const_cast<WaterSupplySystems*>(wss.get())->updateTotalAvailableVolume();
            utility_total_available_volume += wss->getTotal_available_volume();
    }
    
    return utility_total_available_volume;
}

double Utility::getTotal_stored_volume() const {
    double utility_total_stored_volume = 0.0;
    
    // Force WSS to update their volumes before aggregating
        for (auto& wss : water_supply_systems) {
            // Force update of volumes from water sources
            const_cast<WaterSupplySystems*>(wss.get())->updateTotalAvailableVolume();
        double wss_stored = wss->getTotal_stored_volume();
        utility_total_stored_volume += wss_stored;
    }
    
    return utility_total_stored_volume;
}

void Utility::updateTotalAvailableVolume() {
    // Delegate to all water supply systems to update their available volumes
    for (auto& wss : water_supply_systems) {
        wss->updateTotalAvailableVolume();
    }
}

const vector<unique_ptr<WaterSupplySystems>> &Utility::getWaterSupplySystems() const {
    return water_supply_systems;
}

void Utility::updateWSSReferences(const vector<WaterSupplySystems*>& new_wss) {
    // Update the non-owning references to point to the WSS objects that have water sources
    water_supply_systems_refs.clear();
    
    for (auto* wss : new_wss) {
        if (wss && wss->getOwner() == this) {
            water_supply_systems_refs.push_back(wss);
        }
    }
}

const vector<WaterSupplySystems*>& Utility::getWSSReferences() const {
    return water_supply_systems_refs;
}

bool Utility::isUsedForRealization() const {
    return used_for_realization;
}

void Utility::setRestrictedPriceForWss(int system_id, double restricted_price) {
    #pragma omp critical(utility_restriction_prices)
    {
        wss_restricted_prices[system_id] = restricted_price * price_rdm_multiplier;
    }
}

double Utility::getRestrictedPriceForWss(int system_id) const {
    auto it = wss_restricted_prices.find(system_id);
    if (it != wss_restricted_prices.end()) {
        return it->second;
    }
    return NON_INITIALIZED;
}

void Utility::setRestrictedResidentialPriceForWss(int system_id, double restricted_price) {
    #pragma omp critical(utility_restriction_prices)
    {
        wss_restricted_residential_prices[system_id] = restricted_price * price_rdm_multiplier;
    }
}

double Utility::getRestrictedResidentialPriceForWss(int system_id) const {
    auto it = wss_restricted_residential_prices.find(system_id);
    if (it != wss_restricted_residential_prices.end()) {
        return it->second;
    }
    return NON_INITIALIZED;
}

double Utility::calculateRestrictedWeeklyPriceForWss(int system_id, int stage,
                                                    int week_of_year,
                                                    const vector<vector<double>> &priceMultipliers) const {
    auto demand_it = wss_demand_fractions.find(system_id);
    auto price_it = wss_water_prices.find(system_id);
    if (demand_it == wss_demand_fractions.end() || price_it == wss_water_prices.end()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Missing demand fractions or water prices for WSS %d", id, system_id);
        throw std::runtime_error(error_msg);
    }

    const auto &fractions = demand_it->second;
    const auto &prices = price_it->second;

    if (stage < 0 || stage >= (int)priceMultipliers.size()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Restriction stage %d out of bounds for price multipliers (size=%zu)",
                id, stage, priceMultipliers.size());
        throw std::out_of_range(error_msg);
    }

    int month_index = min((int) (week_of_year / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
    double monthly_avg_price = 0.0;
    int n_tiers = static_cast<int>(prices.at(0).size());

    for (int t = 0; t < n_tiers; ++t) {
        monthly_avg_price += fractions[month_index][t] * prices[month_index][t] * priceMultipliers[stage][t];
    }

    return monthly_avg_price / WEEKS_IN_MONTH;
}

double Utility::calculateRestrictedWeeklyResidentialPriceForWss(int system_id, int stage,
                                                                int week_of_year,
                                                                const vector<vector<double>> &priceMultipliers) const {
    auto price_it = wss_water_prices.find(system_id);
    if (price_it == wss_water_prices.end()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Missing water prices for WSS %d", id, system_id);
        throw std::runtime_error(error_msg);
    }

    if (stage < 0 || stage >= (int)priceMultipliers.size()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Restriction stage %d out of bounds for price multipliers (size=%zu)",
                id, stage, priceMultipliers.size());
        throw std::out_of_range(error_msg);
    }

    const auto &prices = price_it->second;
    int month_index = min((int) (week_of_year / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
    if (month_index < 0 || month_index >= (int)prices.size() || prices[month_index].empty()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Missing residential price for WSS %d in month %d", id, system_id, month_index);
        throw std::runtime_error(error_msg);
    }

    if (priceMultipliers[stage].empty()) {
        char error_msg[256];
        sprintf(error_msg, "Utility %d: Missing residential price multiplier for stage %d", id, stage);
        throw std::runtime_error(error_msg);
    }

    double residential_price = prices[month_index][0] * priceMultipliers[stage][0];
    return residential_price / WEEKS_IN_MONTH;
}