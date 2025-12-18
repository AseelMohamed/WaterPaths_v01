//
// Created by bernardoct on 8/25/17.
//

#include <numeric>
#include <algorithm>
#include <limits>
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
    double obj_value =
            1. - (double) *max_element(year_reliabilities.begin(),
                                       year_reliabilities.end()) /
                 n_realizations;

    if (std::isinf(obj_value)) {
        string error_inf = "Infinite reliability.";
        throw logic_error(error_inf.c_str());
    }
    
    return obj_value;
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
    // NPC is cumulative, so take only the final value (last week)
    for (const unsigned long &r : realizations) {
        const auto& npc_vector = utility_data[r]->getNet_present_infrastructure_cost();
        if (!npc_vector.empty()) {
            double realization_npc = npc_vector.back();
            infrastructure_npc += realization_npc;
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

    // Use minimum n_weeks across all realizations to prevent reading uninitialized memory
    unsigned long n_weeks = std::numeric_limits<unsigned long>::max();
    for (const auto &r : realizations) {
        n_weeks = std::min(n_weeks,
                           (unsigned long)utility_data[r]->getGross_revenues().size());
    }
    if (n_weeks == std::numeric_limits<unsigned long>::max()) {
        throw std::runtime_error("No weeks available for peak financial cost objective.");
    }

#ifndef NDEBUG
    // In debug builds, ensure consistency across all realizations
    unsigned long first_weeks = utility_data[realizations[0]]->getGross_revenues().size();
    for (const auto &r : realizations) {
        const auto &g = utility_data[r]->getGross_revenues();
        if (g.size() != first_weeks) {
            fprintf(stderr, "Warning: Inconsistent number of weeks between realizations in calculatePeakFinancialCostsObjective. "
                    "Realization %lu has %lu weeks, realization %lu has %lu weeks.\n",
                    realizations[0], first_weeks, r, (unsigned long)g.size());
        }
    }
#endif

    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    double discount_rate = utility_data[0]->getInfraDiscountRate();
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

    // Use minimum n_weeks across all realizations to prevent reading uninitialized memory
    unsigned long n_weeks = std::numeric_limits<unsigned long>::max();
    for (const auto &r : realizations) {
        n_weeks = std::min(n_weeks,
                           (unsigned long)utility_data[r]->getGross_revenues().size());
    }
    if (n_weeks == std::numeric_limits<unsigned long>::max()) {
        throw std::runtime_error("No weeks available for worse case cost objective.");
    }

#ifndef NDEBUG
    // In debug builds, ensure consistency across all realizations
    unsigned long first_weeks = utility_data[realizations[0]]->getGross_revenues().size();
    for (const auto &r : realizations) {
        const auto &g = utility_data[r]->getGross_revenues();
        if (g.size() != first_weeks) {
            fprintf(stderr, "Warning: Inconsistent number of weeks between realizations in calculateWorseCaseCostsObjective. "
                    "Realization %lu has %lu weeks, realization %lu has %lu weeks.\n",
                    realizations[0], first_weeks, r, (unsigned long)g.size());
        }
    }
#endif

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
        double max_cost = *max_element(
            year_financial_costs.begin(), 
            year_financial_costs.end());

        worse_year_financial_costs.push_back(max_cost);
    }

    // sort costs to get the worse 1 percentile.
    sort(worse_year_financial_costs.begin(),
         worse_year_financial_costs.end());

    unsigned long percentile_index = (unsigned long) floor(WORSE_CASE_COST_PERCENTILE * n_realizations);
    double obj_value = worse_year_financial_costs.at(percentile_index);

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
    
    // Loop through realizations first, then sum across WSS within each realization
    for (const unsigned long &r : realizations) {
        double realization_npc = 0.0;
        
        // Sum NPC across all WSS for this realization
        int wss_idx = 0;
        for (const auto& wss_realization_data : wss_data) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                // NPC is cumulative, so take only the final value (last week)
                const auto& npc_vector = wss_realization_data[r]->getNet_present_infrastructure_cost();
                if (!npc_vector.empty()) {
                    double wss_npc = npc_vector.back();
                    realization_npc += wss_npc;
                }
            }
            wss_idx++;
        }
        total_infrastructure_npc += realization_npc;
    }
    // Return average across realizations
    return total_infrastructure_npc / n_realizations;
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
    
    if (wss_data.empty()) {
        throw std::runtime_error("ERROR: wss_data is empty in calculateReliabilityObjective_WSS");
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    if (realizations.empty()) {
        throw std::runtime_error("ERROR: realizations vector is empty in calculateReliabilityObjective_WSS");
    }
    
    // Check if first realization data exists
    if (realizations[0] >= wss_data[0].size() || wss_data[0][realizations[0]] == nullptr) {
        char error[512];
        sprintf(error, "ERROR: First realization %lu is invalid or nullptr in wss_data[0] (size=%zu)",
                realizations[0], wss_data[0].size());
        throw std::runtime_error(error);
    }
    
    unsigned long n_weeks = wss_data[0][realizations[0]]->getCombined_storage().size();
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    
    vector<double> wss_reliabilities; // Store reliability for each WSS
    int wss_skipped_no_data = 0;
    int wss_skipped_no_storage = 0;
    int wss_processed = 0;
    
    // Calculate reliability for EACH WSS independently
    int wss_idx = 0;
    for (const auto& wss_realization_data : wss_data) {
        // Check if this WSS has any storage capacity > 0 (actual storage, not flow-through)
        bool has_storage_data = false;
        bool has_any_data = false;
        
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                const auto& storage_vec = wss_realization_data[r]->getCombined_storage();
                const auto& capacity_vec = wss_realization_data[r]->getStorage_capacity();
                
                if (!storage_vec.empty() && !capacity_vec.empty()) {
                    has_any_data = true;
                    
                    // Check if ANY capacity value is > 0 (actual storage exists)
                    // Only need to find one non-zero capacity to know this WSS has storage
                    for (const auto& cap : capacity_vec) {
                        if (cap > 0.0) {
                            has_storage_data = true;
                            break;
                        }
                    }
                    if (has_storage_data) break;
                }
            }
        }
        
        // Skip WSS with no storage capacity (e.g., flow-through intakes)
        // We only skip if we found data but all capacities are zero
        // If we found NO data at all, something is wrong - don't skip silently
        if (!has_any_data) {
            wss_skipped_no_data++;
            continue;
        }
        
        if (!has_storage_data) {
            // This WSS has data but no storage capacity - it's a flow-through system
            wss_skipped_no_storage++;
            continue;
        }
        
        wss_processed++;
        
        vector<vector<int>> realizations_year_reliabilities(
                n_realizations, vector<int>(n_years, NON_INITIALIZED));
        vector<int> year_reliabilities(n_years, 0);
        
        // Check failures for this WSS across all realizations
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                const auto& storage_vec = wss_realization_data[r]->getCombined_storage();
                const auto& capacity_vec = wss_realization_data[r]->getStorage_capacity();
                
                // Skip if empty
                if (storage_vec.empty() || capacity_vec.empty()) {
                    continue;
                }
                
                for (unsigned long y = 0; y < n_years; ++y) {
                    for (int w = (int) round(y * WEEKS_IN_YEAR);
                         w < (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                        
                        // Bounds check before accessing vectors
                        if (w >= (int)storage_vec.size() || w >= (int)capacity_vec.size()) {
                            char error[512];
                            sprintf(error, "ERROR: Week %d out of bounds for WSS data (storage size=%zu, capacity size=%zu) "
                                    "in realization %lu, year %lu. Expected %lu weeks.",
                                    w, storage_vec.size(), capacity_vec.size(), r, y, n_weeks);
                            throw std::out_of_range(error);
                        }
                        
                        double storage = storage_vec[w];
                        double capacity = capacity_vec[w];
                        
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
        
        int max_failures = *max_element(year_reliabilities.begin(),
                                       year_reliabilities.end());
        
        // Validate max_failures doesn't exceed n_realizations
        if (max_failures > (int)n_realizations) {
            char error[512];
            sprintf(error, "ERROR: max_failures (%d) exceeds n_realizations (%lu). "
                    "This indicates a counting error in reliability calculation.",
                    max_failures, (unsigned long)n_realizations);
            throw std::runtime_error(error);
        }
        
        double wss_reliability = 1.0 - (double)max_failures / n_realizations;
        wss_reliabilities.push_back(wss_reliability);
        wss_idx++;
    }
    
    // Utility is only as reliable as its weakest system (among those with storage)
    // Check if wss_reliabilities is empty before calling min_element
    if (wss_reliabilities.empty()) {
        // No WSS with storage found - this can happen when:
        // 1. All WSS are flow-through (no storage capacity)
        // 2. Data collection was skipped (e.g., during ROF table export mode)
        // 3. WSS data collectors are not properly initialized
        // No WSS with storage found - return worst-case reliability
        // This can happen when data collection was skipped or all WSS are flow-through
        #ifdef PARALLEL
        printf("WARNING: wss_reliabilities is EMPTY - returning 0.0 reliability!\n");
        #endif
        return 0.0;
    }
    
    double utility_reliability = *min_element(wss_reliabilities.begin(),
                                              wss_reliabilities.end());
    
    if (std::isinf(utility_reliability)) {
        string error_inf = "Infinite reliability (WSS-level).";
        throw logic_error(error_inf.c_str());
    }
    
    return utility_reliability;
}
