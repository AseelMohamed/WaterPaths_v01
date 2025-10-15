//
// Created by bernardo on 1/13/17.
//

#include <algorithm>
#include "Utility.h"
#include "WaterSupplySystems.h"
#include "../../Utils/Utils.h"

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
    
    // Only initialize water supply systems if water_source_to_wtp is not empty
    // If empty, it means WSS will be added manually later
    if (!water_source_to_wtp.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            name, id, id, this, demands_all_realizations, number_of_week_demands,
            wwtp_discharge_rule, demand_buffer,
            water_source_to_wtp, utility_owned_wtp_capacities));
    }
    
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

    // Initialize water supply systems like in other constructors
    water_supply_systems.clear();
    if (!utility.water_supply_systems.empty()) {
        water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
            utility.name, utility.id, utility.id, this, utility.wwtp_discharge_rule));
        
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
                                   const std::vector<double>& utility_owned_wtp_capacities) {
    // Use the main constructor with all parameters
    // Create a mutable copy of const references for the constructor
    auto water_source_to_wtp_copy = water_source_to_wtp;
    auto utility_owned_wtp_capacities_copy = utility_owned_wtp_capacities;
    WwtpDischargeRule wwtp_discharge_rule_copy(const_cast<WwtpDischargeRule&>(wwtp_discharge_rule));
    
    water_supply_systems.emplace_back(std::make_unique<WaterSupplySystems>(
        name, system_id, utility_id, this, demands_all_realizations, 
        number_of_week_demands, wwtp_discharge_rule_copy, demand_buffer,
        water_source_to_wtp_copy, utility_owned_wtp_capacities_copy));
    
    // CRITICAL: Reconnect the infrastructure manager after the WSS is in its final location
    // This ensures the pointer is correct after any potential moves/copies
    auto& wss = water_supply_systems.back();
    wss->reconnectInfrastructureManager();
}

void Utility::clearWaterSupplySystems() {
    water_supply_systems.clear();
}

/////////////////////////////////////////////////
////////////////////////////////////////////////

void Utility::updateContingencyFundAndDebtService(
        double utility_restricted_demand, double demand_multiplier,
        double demand_offset, double utility_unfulfilled_demand, int week) {
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
        throw logic_error("Prices under surcharge cannot be smaller than "
                          "prices w/o restrictions enacted.");

    // calculate fund contributions if there were no shortage.
    double projected_fund_contribution = percent_contingency_fund_contribution *
                                         utility_restricted_demand *
                                         unrestricted_price;

    // Calculate actual gross revenue.
    gross_revenue = utility_restricted_demand * current_price;

    // Calculate losses due to restrictions and transfers.
    double lost_demand_vol_sales =
            (utility_restricted_demand * (1 - demand_multiplier) +
             utility_unfulfilled_demand);
    double revenue_losses = lost_demand_vol_sales * unrestricted_price;
    double transfer_costs = demand_offset * (offset_rate_per_volume -
                                             unrestricted_price);
    double recouped_loss_price_surcharge =
            utility_restricted_demand * (current_price - unrestricted_price);

    // contingency fund cannot get negative.
    contingency_fund = max(contingency_fund + projected_fund_contribution -
                           revenue_losses - transfer_costs +
                           recouped_loss_price_surcharge,
                           0.0);

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

/////////////////////////////////////////////////
// --------------- CHECK THIS --------------- //
///////////////////////////////////////////////

void Utility::issueBond(int new_infra_triggered, int week) {
    if (new_infra_triggered != NON_INITIALIZED) {
        // Find the water source across all water supply systems using the new_infra_triggered ID
        WaterSource* water_source = nullptr;
        for (const auto& wss : water_supply_systems) {
            const auto& sources = wss->getWater_sources();
            for (WaterSource* source : sources) {
                if (source && source->id == new_infra_triggered) {
                    water_source = source;
                    break;
                }
            }
            if (water_source) break;
        }
        
        if (water_source) {
            Bond &bond = water_source->getBond(id);
            if (!bond.isIssued()) {
                double construction_time = water_source->getConstruction_time();
                bond.issueBond(week, (int) construction_time, bond_term_multiplier,
                               bond_interest_rate_multiplier);
                issued_bonds.push_back(&bond);
                infra_net_present_cost += bond.getNetPresentValueAtIssuance(
                        infra_discount_rate, week);
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
    // Set realization for all water supply systems
    for (auto& wss : water_supply_systems) {
        wss->setRealization(r, rdm_factors);
    }

    // Then aggregate demand series realization from all WSS at utility level
    utility_demand_series_realization = vector<double>(number_of_week_demands, 0.0);
    for (int week = 0; week < number_of_week_demands; ++week) {
        for (const auto& wss : water_supply_systems) {
            utility_demand_series_realization[week] += wss->getUnrestrictedDemand(week);
        }
    }

    // Set financial/policy realization factors at the utility level
    bond_interest_rate_multiplier = rdm_factors.at(1);
    bond_term_multiplier = rdm_factors.at(2);
    infra_discount_rate *= rdm_factors.at(3);
    price_rdm_multiplier = rdm_factors.at(4);

    for (double &awp : weekly_average_volumetric_price) {
        awp *= price_rdm_multiplier;
    }

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
        utility_unrestricted_demand += wss->getUnrestrictedDemand();
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
        utility_restricted_demand += wss->getRestrictedDemand();
    }
    return utility_restricted_demand;
}

double Utility::getGrossRevenue() const {
    return gross_revenue;
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

void Utility::setNoFinancialCalculations() {
    used_for_realization = false;
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

double Utility::getTotal_storage_capacity() const {
    double utility_total_storage_capacity = 0.0;
    for (const auto& wss : water_supply_systems) {
        utility_total_storage_capacity += wss->getTotal_storage_capacity();
    }
    return utility_total_storage_capacity;
}

const vector<unique_ptr<WaterSupplySystems>> &Utility::getWaterSupplySystems() const {
    return water_supply_systems;
}