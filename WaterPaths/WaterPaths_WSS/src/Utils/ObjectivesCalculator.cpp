//
// Created by bernardoct on 8/25/17.
//

#include <numeric>
#include <algorithm>
#include "ObjectivesCalculator.h"
#include "Utils.h"

double ObjectivesCalculator::calculateReliabilityObjective(
        const vector<UtilitiesDataCollector *> &utility_collector,
        vector<unsigned long> realizations) {
    unsigned long n_realizations = utility_collector.size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }

    unsigned long n_weeks = utility_collector[realizations[0]]->getCombined_storage().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);

    vector<vector<int>> realizations_year_reliabilities(
            utility_collector.size(), vector<int>(n_years, NON_INITIALIZED));
    vector<int> year_reliabilities(n_years, 0);

    // Creates a table with years that failed in each realization.
    for (const unsigned long &r : realizations) {
        for (unsigned long y = 0; y < n_years; ++y) {
            for (int w = (int) round(y * WEEKS_IN_YEAR);
                 w <
                 (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR));
                 ++w) {
                if (utility_collector[r]->getCombined_storage()[w] /
                    utility_collector[r]->getCapacity()[w] <
                    STORAGE_CAPACITY_RATIO_FAIL) {
                    realizations_year_reliabilities[r][y] = FAILURE;
                }
            }
        }
    }

    // Creates a vector with the number of realizations that failed for each year.
    for (unsigned long y = 0; y < n_years; ++y) {
        for (const unsigned long &r : realizations) {
            if (realizations_year_reliabilities[r][y] == FAILURE)
                year_reliabilities[y]++;
        }
    }

    double check_non_zero = accumulate(year_reliabilities.rbegin(),
                                       year_reliabilities.rend(),
                                       0.0);

    // Returns year with most realization failures, divided by the number of realizations (reliability objective).
    if (check_non_zero > 0) {
        double obj_value =
                1. - (double) *max_element(year_reliabilities.begin(),
                                           year_reliabilities.end()) /
                     n_realizations;

        if (std::isinf(obj_value)) {
            string error_inf = "Infinite reliability.";
            throw logic_error(error_inf.c_str());
        } else {
            return obj_value;
        }
    } else {
        return 1.;
    }
}

double ObjectivesCalculator::calculateRestrictionFrequencyObjective(
        const vector<RestrictionsDataCollector *> &restriction_data,
        vector<unsigned long> realizations) {

    unsigned long n_realizations = restriction_data.size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }

    // Check if there were restriction policies in place.
    if (!restriction_data.empty()) {
        unsigned long n_weeks = restriction_data[realizations[0]]->getRestriction_multipliers().size();
        unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);

        double restriction_frequency = 0;

        // Counts how many years across all realizations had restrictions.
        for (const unsigned long &r : realizations) {
            for (unsigned long y = 0; y < n_years; ++y) {
                for (int w = (int) round(y * WEEKS_IN_YEAR);
                     w < (int) min((int) n_weeks,
                                   (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                    if (restriction_data[r]->getRestriction_multipliers()[w] !=
                        1.0) {
                        restriction_frequency++;
                        break;
                    }
                }
            }
        }

        double obj_value = restriction_frequency / (n_realizations * n_years);

        if (std::isinf(obj_value)) {
            string error_inf = "Infinite restriction frequency.";
            throw logic_error(error_inf.c_str());
        } else {
            return obj_value;
        }
    } else
        return NONE;
}

double ObjectivesCalculator::calculateNetPresentCostInfrastructureObjective(
        const vector<UtilitiesDataCollector *> &utility_data,
        vector<unsigned long> realizations) {

    unsigned long n_realizations = utility_data.size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }

    double infrastructure_npc = 0;
    // WSS FIX: Infrastructure cost is cumulative and stored weekly
    // We only need the FINAL value (last week), not the sum of all weeks
    // The cost accumulates as bonds are issued, so the last value = total NPC
    for (const unsigned long &r : realizations) {
        if (!utility_data[r]->getNet_present_infrastructure_cost().empty()) {
            infrastructure_npc += utility_data[r]->getNet_present_infrastructure_cost().back();
        }
    }
    return infrastructure_npc / n_realizations;
}

double ObjectivesCalculator::calculatePeakFinancialCostsObjective(
        const vector<UtilitiesDataCollector *> &utility_data,
        vector<unsigned long> realizations) {

    unsigned long n_realizations = utility_data.size();
    
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }

    unsigned long n_weeks = utility_data[realizations[0]]->getGross_revenues().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    double discount_rate = utility_data[0]->getUtility()->getInfraDiscountRate();

    double realizations_year_debt_payment = 0;
    double realizations_year_cont_fund_contribution = 0;
    double realizations_year_gross_revenue = 1e-6;
    double realizations_year_insurance_contract_cost = 0;
    vector<double> year_financial_costs;
    vector<double> realization_financial_costs(utility_data.size(), 0);

    // Creates a table with years that failed in each realization.
    int y;
    for (const unsigned long &r : realizations) {
        year_financial_costs.assign(n_years, 0.0);
        y = 0;
        
        // Check vector sizes and sample values for first realization
        if (r == 0) {
            auto &debt = utility_data[r]->getDebt_service_payments();
            auto &cont = utility_data[r]->getContingency_fund_contribution();
            auto &rev = utility_data[r]->getGross_revenues();
            auto &ins = utility_data[r]->getInsurance_contract_cost();
        }
        
        for (unsigned long w = 0; w < n_weeks; ++w) {
            // accumulate year's info by summing weekly amounts.
            realizations_year_debt_payment +=
                    utility_data[r]->getDebt_service_payments()[w];
            realizations_year_cont_fund_contribution +=
                    utility_data[r]->getContingency_fund_contribution()[w];
            realizations_year_gross_revenue +=
                    utility_data[r]->getGross_revenues()[w];
            realizations_year_insurance_contract_cost +=
                    utility_data[r]->getInsurance_contract_cost()[w];

            // if last week of the year, close the books and calculate
            // financial cost for the year.
            if (Utils::isFirstWeekOfTheYear(w + 1)) {
                year_financial_costs[y] +=
                        (realizations_year_debt_payment +
                         realizations_year_cont_fund_contribution +
                         realizations_year_insurance_contract_cost) /
//                                realizations_year_gross_revenue;
                        (realizations_year_gross_revenue *
                         (1. + pow(1. + discount_rate, y)));
                // update year count.
                y++;

                // reset accounts.
                realizations_year_debt_payment = 0;
                realizations_year_cont_fund_contribution = 0;
                realizations_year_gross_revenue = 1e-6;
                realizations_year_insurance_contract_cost = 0;
            }
        }
        // store highest year cost as the cost financial cost of the realization.
        realization_financial_costs[r] =
                *max_element(year_financial_costs.begin(),
                             year_financial_costs.end());
        if (realization_financial_costs[r] > 1e10) {
            printf("Absurdly high financial cost in realization %lu.\n", r);
        }
    }

    double obj_value = accumulate(realization_financial_costs.begin(),
                                  realization_financial_costs.end(),
                                  0.0) / n_realizations;

    if (std::isinf(obj_value)) {
        string error_inf = "Infinite peak financial cost.";
        throw logic_error(error_inf.c_str());
    } else {
        return obj_value;
    }
}

double ObjectivesCalculator::calculateWorseCaseCostsObjective(
        const vector<UtilitiesDataCollector *> &utility_data,
        vector<unsigned long> realizations) {

    unsigned long n_realizations = utility_data.size();
    
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }

    unsigned long n_weeks = utility_data[realizations[0]]->getGross_revenues().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    double discount_rate = utility_data[0]->getUtility()->getInfraDiscountRate();

    vector<double> worse_year_financial_costs;
    vector<double> year_financial_costs;
    double year_drought_mitigation_cost = 0;
    double year_gross_revenue = 1e-6;

    // Creates a table with years that failed in each realization.
    int y;
    for (const unsigned long &r : realizations) {
        y = 0;
        year_financial_costs.assign(n_years, 0);
        
        for (unsigned long w = 0; w < n_weeks; ++w) {
            // accumulate year's info by summing weekly amounts.
            year_drought_mitigation_cost +=
                    utility_data[r]->getDrought_mitigation_cost()[w];
            year_gross_revenue += utility_data[r]->getGross_revenues()[w];

            // if last week of the year, close the books and calculate financial cost for the year.
            if (Utils::isFirstWeekOfTheYear(w + 1)) {
                year_financial_costs[y] =
                        max(year_drought_mitigation_cost
                            - utility_data[r]->getContingency_fund_size()[w],
//                            0.0) / year_gross_revenue;
                            0.0) / (year_gross_revenue * (1. + pow(1. + discount_rate, y)));

                year_gross_revenue = 1e-6;
                year_drought_mitigation_cost = 0;
                // update year count.
                y++;
            }
        }
        // store highest year cost as the drought mitigation cost of the realization.
        worse_year_financial_costs.push_back(*max_element(
                year_financial_costs.begin(),
                year_financial_costs.end()));
    }

    // sort costs to get the worse 1 percentile.
    sort(worse_year_financial_costs.begin(),
         worse_year_financial_costs.end());

    double obj_value = worse_year_financial_costs.at(
            (unsigned long) floor(WORSE_CASE_COST_PERCENTILE * n_realizations));

    if (std::isinf(obj_value)) {
        string error_inf = "Infinite worse case cost.";
        throw logic_error(error_inf.c_str());
    } else {
        return obj_value;
    }
}

/**
 * Calculate infrastructure NPC at WSS level, then aggregate across all WSS and realizations.
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param utility_data Utility data collectors (used for discount rate - not dereferenced, just passed for safety)
 * @param realizations Optional vector of realization indices to include
 * @return Average infrastructure NPC across all WSS and realizations
 */
double ObjectivesCalculator::calculateNetPresentCostInfrastructureObjective_WSS(
        const vector<vector<WSSDataCollector *>>& wss_data,
        const vector<UtilitiesDataCollector *>& utility_data,
        vector<unsigned long> realizations) {
    
    if (wss_data.empty() || wss_data[0].empty()) {
        return 0.0;
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    double total_infrastructure_npc = 0.0;
    
    // Loop through each WSS
    for (const auto& wss_realization_data : wss_data) {
        // Loop through selected realizations for this WSS
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                // Infrastructure NPC accumulates weekly, so take the final value
                const auto& npc_vector = wss_realization_data[r]->getNet_present_infrastructure_cost();
                if (!npc_vector.empty()) {
                    total_infrastructure_npc += npc_vector.back();
                }
            }
        }
    }
    
    // Return average across all WSS and realizations
    return total_infrastructure_npc / n_realizations;
}

/**
 * Calculate worst case cost at WSS level, then aggregate across all WSS.
 * For each realization, sum worst case costs across all WSS, then take the 99th percentile.
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param utility_data Utility data collectors (used to safely get discount rate)
 * @param realizations Optional vector of realization indices to include
 * @return 99th percentile worst case cost across realizations
 */
double ObjectivesCalculator::calculateWorseCaseCostsObjective_WSS(
        const vector<vector<WSSDataCollector *>>& wss_data,
        const vector<UtilitiesDataCollector *>& utility_data,
        vector<unsigned long> realizations) {
    
    if (wss_data.empty() || wss_data[0].empty()) {
        return 0.0;
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    unsigned long n_weeks = wss_data[0][realizations[0]]->getGross_revenues().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    
    // Get discount rate safely from utility data collector
    double discount_rate = 0.05; // Default fallback
    if (!utility_data.empty() && utility_data[realizations[0]] != nullptr &&
        utility_data[realizations[0]]->getUtility() != nullptr) {
        discount_rate = utility_data[realizations[0]]->getUtility()->getInfraDiscountRate();
    }
    
    vector<double> worse_year_financial_costs;
    
    // For each realization, calculate worst case cost
    for (const unsigned long &r : realizations) {
        vector<double> year_financial_costs(n_years, 0.0);
        
        // Aggregate financial data across all WSS for this realization
        for (int y = 0; y < (int)n_years; ++y) {
            double year_drought_mitigation_cost = 0.0;
            double year_gross_revenue = 1e-6;
            double year_end_contingency_fund = 0.0;
            
            // Sum weekly values for each WSS across the year
            for (const auto& wss_realization_data : wss_data) {
                if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                    for (int w = (int) round(y * WEEKS_IN_YEAR);
                         w < (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                        
                        year_drought_mitigation_cost += 
                            wss_realization_data[r]->getDrought_mitigation_cost()[w];
                        year_gross_revenue += 
                            wss_realization_data[r]->getGross_revenues()[w];
                        
                        // Contingency fund is utility-wide, but stored in each WSS
                        // Use the value from the last week of the year
                        if (w == (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR)) - 1) {
                            year_end_contingency_fund = 
                                wss_realization_data[r]->getContingency_fund_size()[w];
                        }
                    }
                }
            }
            
            // Calculate year cost: max(drought_cost - contingency_fund, 0) / (revenue * discount_factor)
            year_financial_costs[y] = 
                max(year_drought_mitigation_cost - year_end_contingency_fund, 0.0) / 
                (year_gross_revenue * (1.0 + pow(1.0 + discount_rate, y)));
        }
        
        // Store the worst (maximum) year cost for this realization
        worse_year_financial_costs.push_back(*max_element(
            year_financial_costs.begin(), year_financial_costs.end()));
    }
    
    // Sort costs to get the 99th percentile
    sort(worse_year_financial_costs.begin(), worse_year_financial_costs.end());
    
    double obj_value = worse_year_financial_costs.at(
        (unsigned long) floor(WORSE_CASE_COST_PERCENTILE * n_realizations));
    
    if (std::isinf(obj_value)) {
        string error_inf = "Infinite worse case cost (WSS-level).";
        throw logic_error(error_inf.c_str());
    } else {
        return obj_value;
    }
}

/**
 * Calculate reliability at WSS level, then take MINIMUM reliability across all WSS.
 * This ensures that utility reliability = reliability of its weakest WSS.
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param realizations Optional vector of realization indices to include
 * @return Minimum reliability across all WSS (utility is only as reliable as weakest WSS)
 */
double ObjectivesCalculator::calculateReliabilityObjective_WSS(
        const vector<vector<WSSDataCollector *>>& wss_data,
        vector<unsigned long> realizations) {
    
    if (wss_data.empty() || wss_data[0].empty()) {
        return 1.0; // No WSS means perfect reliability (no failure)
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    unsigned long n_weeks = wss_data[0][realizations[0]]->getCombined_storage().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    
    vector<double> wss_reliabilities; // Store reliability for each WSS
    
    // Calculate reliability for EACH WSS independently
    for (const auto& wss_realization_data : wss_data) {
        vector<vector<int>> realizations_year_reliabilities(
                n_realizations, vector<int>(n_years, NON_INITIALIZED));
        vector<int> year_reliabilities(n_years, 0);
        
        // Check failures for this WSS across all realizations
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                for (unsigned long y = 0; y < n_years; ++y) {
                    for (int w = (int) round(y * WEEKS_IN_YEAR);
                         w < (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                        
                        double storage = wss_realization_data[r]->getCombined_storage()[w];
                        double capacity = wss_realization_data[r]->getStorage_capacity()[w];
                        
                        if (capacity > 0 && storage / capacity < STORAGE_CAPACITY_RATIO_FAIL) {
                            realizations_year_reliabilities[r][y] = FAILURE;
                            break; // Year already failed, no need to check more weeks
                        }
                    }
                }
            }
        }
        
        // Count failures per year for this WSS
        for (unsigned long y = 0; y < n_years; ++y) {
            for (const unsigned long &r : realizations) {
                if (r < realizations_year_reliabilities.size() && 
                    realizations_year_reliabilities[r][y] == FAILURE) {
                    year_reliabilities[y]++;
                }
            }
        }
        
        // Calculate reliability for this WSS
        double check_non_zero = accumulate(year_reliabilities.begin(),
                                           year_reliabilities.end(), 0.0);
        
        double wss_reliability;
        if (check_non_zero > 0) {
            int max_failures = *max_element(year_reliabilities.begin(),
                                           year_reliabilities.end());
            wss_reliability = 1.0 - (double)max_failures / n_realizations;
        } else {
            wss_reliability = 1.0; // No failures
        }
        
        wss_reliabilities.push_back(wss_reliability);
    }
    
    // Utility reliability = MINIMUM reliability among all its WSS
    // (utility is only as reliable as its weakest system)
    double utility_reliability = *min_element(wss_reliabilities.begin(),
                                              wss_reliabilities.end());
    
    if (std::isinf(utility_reliability)) {
        string error_inf = "Infinite reliability (WSS-level).";
        throw logic_error(error_inf.c_str());
    }
    
    return utility_reliability;
}
