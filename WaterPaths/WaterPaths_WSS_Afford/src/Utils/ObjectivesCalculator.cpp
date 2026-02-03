//
// Created by bernardoct on 8/25/17.
//

#include <numeric>
#include <algorithm>
#include <limits>
#include <fstream>
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
    // double realizations_year_cont_fund_contribution = 0;
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
            // realizations_year_cont_fund_contribution +=
            //         utility_data[r]->getContingency_fund_contribution()[w];
            realizations_year_gross_revenue +=
                    utility_data[r]->getGross_revenues()[w];
            realizations_year_insurance_contract_cost +=
                    utility_data[r]->getInsurance_contract_cost()[w];

            // if last week of the year, close the books and calculate
            // financial cost for the year.
            if (Utils::isFirstWeekOfTheYear(w + 1)) {
                year_financial_costs[y] +=
                        (realizations_year_debt_payment +
                         realizations_year_insurance_contract_cost) /
                        (realizations_year_gross_revenue *
                         (1. + pow(1. + discount_rate, y)));
                // update year count.
                y++;

                // reset accounts.
                realizations_year_debt_payment = 0;
                // realizations_year_cont_fund_contribution = 0;
                realizations_year_gross_revenue = 1e-6;
                realizations_year_insurance_contract_cost = 0;
            }
        }
        // store highest year cost as the cost financial cost of the realization.
        realization_financial_costs[r] =
                *max_element(year_financial_costs.begin(),
                             year_financial_costs.end());
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

// #ifndef NDEBUG
//     // In debug builds, ensure consistency across all realizations
//     unsigned long first_weeks = utility_data[realizations[0]]->getGross_revenues().size();
//     for (const auto &r : realizations) {
//         const auto &g = utility_data[r]->getGross_revenues();
//         if (g.size() != first_weeks) {
//             fprintf(stderr, "Warning: Inconsistent number of weeks between realizations in calculateWorseCaseCostsObjective. "
//                     "Realization %lu has %lu weeks, realization %lu has %lu weeks.\n",
//                     realizations[0], first_weeks, r, (unsigned long)g.size());
//         }
//     }
// #endif

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
                        max(year_drought_mitigation_cost,0.0) / 
                        (year_gross_revenue * (1. + pow(1. + discount_rate, y)));

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
    // CSV output disabled (reliability)
    // const char* reliability_csv_path = "output/reliability_yearly_wss.csv";
    // std::ofstream reliability_out(reliability_csv_path, std::ios::trunc);
    // if (reliability_out.good()) {
    //     reliability_out << "wss_id,year,reliability\n";
    // } else {
    //     fprintf(stderr, "WARNING: Could not open %s for writing.\n", reliability_csv_path);
    // }
    
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

        // CSV output disabled (reliability)
        // for (unsigned long y = 0; y < n_years; ++y) {
        //     double year_reliability = 1.0 - (double)year_reliabilities[y] / n_realizations;
        //     printf("Reliability (year %lu, WSS %d): %.10f\n", y, wss_idx, year_reliability);
        //     if (reliability_out.good()) {
        //         reliability_out << wss_idx << "," << y << "," << year_reliability << "\n";
        //     }
        // }

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

/**
 * Calculate affordability index using year-based approach for all WSS.
 * Affordability index = water_price / average_income
 * For each realization: find worst weekly affordability per year.
 * For each year: average the worst affordability across realizations.
 * For each WSS: average across all years in the simulation period.
 * The objective is the maximum affordability across all WSS (worst case).
 * 
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param realizations Optional vector of realization indices to include
 * @return Year-based affordability index (worst case across WSS)
 */
double ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS(
        const vector<vector<WSSDataCollector *>>& wss_data,
        vector<unsigned long> realizations) {
    
    if (wss_data.empty()) {
        throw std::runtime_error("ERROR: wss_data is empty in calculateAffordabilityIndexObjective_WSS");
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    if (realizations.empty()) {
        throw std::runtime_error("ERROR: realizations vector is empty in calculateAffordabilityIndexObjective_WSS");
    }
    
    // Get number of weeks and years from first WSS, first realization
    unsigned long n_weeks = 0;
    for (const auto& wss_realization_data : wss_data) {
        if (!wss_realization_data.empty() && realizations[0] < wss_realization_data.size() && 
            wss_realization_data[realizations[0]] != nullptr) {
            const auto& affordability_vec = wss_realization_data[realizations[0]]->getWeekly_affordability_index();
            if (!affordability_vec.empty()) {
                n_weeks = affordability_vec.size();
                break;
            }
        }
    }
    
    if (n_weeks == 0) {
        #ifdef PARALLEL
        printf("WARNING: No affordability data found - returning 0.0 affordability!\n");
        #endif
        return 0.0;
    }
    
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    const double affordability_scale = 2.44 / 1e6; // Convert R$/hm3 to R$/m3 and per-capita to household
    // CSV output disabled (affordability)
    // const char* affordability_csv_path = "output/affordability_yearly_wss.csv";
    // std::ofstream affordability_out(affordability_csv_path, std::ios::trunc);
    // if (affordability_out.good()) {
    //     affordability_out << "wss_id,year,affordability_index\n";
    // }
    
    vector<double> wss_affordability_avg; // Store average affordability for each WSS
    
    // Calculate average affordability for EACH WSS independently
    int wss_idx = 0;
    for (const auto& wss_realization_data : wss_data) {
        // For each year, store the average affordability across realizations
        vector<double> year_avg_affordability(n_years, 0.0);
        vector<int> year_realization_count(n_years, 0);
        
        // For each realization, find worst weekly affordability per year
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                const auto& affordability_vec = wss_realization_data[r]->getWeekly_affordability_index();
                
                if (!affordability_vec.empty()) {
                    // Find worst affordability for each year in this realization
                    for (unsigned long y = 0; y < n_years; ++y) {
                        double max_affordability_in_year = 0.0;
                        
                        for (int w = (int) round(y * WEEKS_IN_YEAR);
                             w < (int) min((int) affordability_vec.size(), 
                                          (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                            max_affordability_in_year = max(max_affordability_in_year, affordability_vec[w]);
                        }
                        
                        // Accumulate for averaging across realizations
                        year_avg_affordability[y] += max_affordability_in_year;
                        year_realization_count[y]++;
                    }
                }
            }
        }
        
        // Calculate average across realizations for each year, then average across years for this WSS
        double wss_avg_affordability = 0.0;
        int valid_years = 0;
        
        for (unsigned long y = 0; y < n_years; ++y) {
            if (year_realization_count[y] > 0) {
                double year_affordability = (year_avg_affordability[y] / year_realization_count[y]) * affordability_scale;
                // if (affordability_out.good()) {
                //     affordability_out << wss_idx << "," << y << "," << year_affordability << "\n";
                // }
                wss_avg_affordability += year_affordability;
                valid_years++;
            }
        }
        
        if (valid_years > 0) {
            wss_avg_affordability /= valid_years;
            wss_affordability_avg.push_back(wss_avg_affordability);
        }
        wss_idx++;
    }
    
    // Return the maximum (worst case) affordability across all WSS
    if (wss_affordability_avg.empty()) {
        #ifdef PARALLEL
        printf("WARNING: wss_affordability_avg is EMPTY - returning 0.0 affordability!\n");
        #endif
        return 0.0;
    }
    
    double worst_affordability = *max_element(wss_affordability_avg.begin(),
                                             wss_affordability_avg.end());
    
    if (std::isinf(worst_affordability)) {
        string error_inf = "Infinite affordability index (WSS-level).";
        throw logic_error(error_inf.c_str());
    }
    
    return worst_affordability;  // Values already scaled by 2.44/1e6.
}

/**
 * Configurable reliability calculation with choice of aggregation method.
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param realizations Vector of realization indices to include
 * @param aggregation_method 0=MIN (worst case), 1=AVERAGE
 * @return Reliability based on chosen aggregation method
 */
double ObjectivesCalculator::calculateReliabilityObjective_WSS_Configurable(
        const vector<vector<WSSDataCollector *>>& wss_data,
        vector<unsigned long> realizations,
        int aggregation_method) {
    
    if (wss_data.empty()) {
        throw std::runtime_error("ERROR: wss_data is empty in calculateReliabilityObjective_WSS_Configurable");
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    if (realizations.empty()) {
        throw std::runtime_error("ERROR: realizations vector is empty");
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
    // CSV output disabled (reliability configurable)
    // const char* reliability_csv_path = "output/reliability_yearly_wss.csv";
    // std::ofstream reliability_out(reliability_csv_path, std::ios::trunc);
    // if (reliability_out.good()) {
    //     reliability_out << "wss_id,year,reliability\n";
    // } else {
    //     fprintf(stderr, "WARNING: Could not open %s for writing.\n", reliability_csv_path);
    // }
    
    vector<double> wss_reliabilities; // Store reliability for each WSS
    int wss_skipped_no_data = 0;
    int wss_skipped_no_storage = 0;
    int wss_processed = 0;
    
    // Calculate reliability for EACH WSS independently
    for (const auto& wss_realization_data : wss_data) {
        // Check if this WSS has any storage capacity > 0
        bool has_storage_data = false;
        bool has_any_data = false;
        
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                const auto& storage_vec = wss_realization_data[r]->getCombined_storage();
                const auto& capacity_vec = wss_realization_data[r]->getStorage_capacity();
                
                if (!storage_vec.empty() && !capacity_vec.empty()) {
                    has_any_data = true;
                    
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
        
        if (!has_any_data) {
            wss_skipped_no_data++;
            continue;
        }
        
        if (!has_storage_data) {
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
                
                if (storage_vec.empty() || capacity_vec.empty()) {
                    continue;
                }
                
                for (unsigned long y = 0; y < n_years; ++y) {
                    for (int w = (int) round(y * WEEKS_IN_YEAR);
                         w < (int) min((int) n_weeks, (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                        
                        if (w >= (int)storage_vec.size() || w >= (int)capacity_vec.size()) {
                            char error[512];
                            sprintf(error, "ERROR: Week %d out of bounds for WSS data", w);
                            throw std::out_of_range(error);
                        }
                        
                        double storage = storage_vec[w];
                        double capacity = capacity_vec[w];
                        
                        if (capacity > 0 && storage / capacity < STORAGE_CAPACITY_RATIO_FAIL) {
                            realizations_year_reliabilities[r][y] = FAILURE;
                            break;
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

        // CSV output disabled (reliability configurable)
        // if (reliability_out.good()) {
        //     for (unsigned long y = 0; y < n_years; ++y) {
        //         double year_reliability = 1.0 - (double)year_reliabilities[y] / n_realizations;
        //         reliability_out << wss_processed << "," << y << "," << year_reliability << "\n";
        //     }
        // }
        
        int max_failures = *max_element(year_reliabilities.begin(),
                                       year_reliabilities.end());
        
        if (max_failures > (int)n_realizations) {
            char error[512];
            sprintf(error, "ERROR: max_failures (%d) exceeds n_realizations (%lu)",
                    max_failures, (unsigned long)n_realizations);
            throw std::runtime_error(error);
        }
        
        double wss_reliability = 1.0 - (double)max_failures / n_realizations;
        wss_reliabilities.push_back(wss_reliability);
    }
    
    if (wss_reliabilities.empty()) {
        #ifdef PARALLEL
        printf("WARNING: wss_reliabilities is EMPTY - returning 0.0 reliability!\n");
        #endif
        return 0.0;
    }
    
    double utility_reliability;
    
    if (aggregation_method == 1) { // AVERAGE
        utility_reliability = accumulate(wss_reliabilities.begin(), 
                                        wss_reliabilities.end(), 0.0) / wss_reliabilities.size();
    } else { // MIN (worst case) - default
        utility_reliability = *min_element(wss_reliabilities.begin(),
                                          wss_reliabilities.end());
    }
    
    if (std::isinf(utility_reliability)) {
        string error_inf = "Infinite reliability (WSS-level configurable).";
        throw logic_error(error_inf.c_str());
    }
    
    return utility_reliability;
}

/**
 * Configurable affordability calculation with choice of aggregation method.
 * @param wss_data Vector of vectors: outer index = WSS ID, inner index = realization
 * @param realizations Vector of realization indices to include
 * @param aggregation_method 0=MAX (worst case), 1=AVERAGE
 * @return Affordability based on chosen aggregation method
 */
double ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_Configurable(
        const vector<vector<WSSDataCollector *>>& wss_data,
        vector<unsigned long> realizations,
        int aggregation_method) {
    
    if (wss_data.empty()) {
        throw std::runtime_error("ERROR: wss_data is empty in calculateAffordabilityIndexObjective_WSS_Configurable");
    }
    
    unsigned long n_realizations = wss_data[0].size();
    if (realizations.empty()) {
        realizations = vector<unsigned long>(n_realizations);
        iota(realizations.begin(), realizations.end(), 0);
    } else {
        n_realizations = realizations.size();
    }
    
    if (realizations.empty()) {
        throw std::runtime_error("ERROR: realizations vector is empty");
    }
    
    // Get number of weeks and years from first WSS, first realization
    unsigned long n_weeks = 0;
    for (const auto& wss_realization_data : wss_data) {
        if (!wss_realization_data.empty() && realizations[0] < wss_realization_data.size() && 
            wss_realization_data[realizations[0]] != nullptr) {
            const auto& affordability_vec = wss_realization_data[realizations[0]]->getWeekly_affordability_index();
            if (!affordability_vec.empty()) {
                n_weeks = affordability_vec.size();
                break;
            }
        }
    }
    
    if (n_weeks == 0) {
        #ifdef PARALLEL
        printf("WARNING: No affordability data found - returning 0.0 affordability!\n");
        #endif
        return 0.0;
    }
    
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    const double affordability_scale = 2.44 / 1e6; // Convert R$/hm3 to R$/m3 and per-capita to household
    // CSV output disabled (affordability configurable)
    // const char* affordability_csv_path = "output/affordability_yearly_wss.csv";
    // std::ofstream affordability_out(affordability_csv_path, std::ios::trunc);
    // if (affordability_out.good()) {
    //     affordability_out << "wss_id,year,affordability_index\n";
    // }
    
    vector<double> wss_affordability_avg; // Store average affordability for each WSS
    
    // Calculate average affordability for EACH WSS independently
    int wss_idx = 0;
    for (const auto& wss_realization_data : wss_data) {
        // For each year, store the average affordability across realizations
        vector<double> year_avg_affordability(n_years, 0.0);
        vector<int> year_realization_count(n_years, 0);
        
        // For each realization, find worst weekly affordability per year
        for (const unsigned long &r : realizations) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                const auto& affordability_vec = wss_realization_data[r]->getWeekly_affordability_index();
                
                if (!affordability_vec.empty()) {
                    // Find worst affordability for each year in this realization
                    for (unsigned long y = 0; y < n_years; ++y) {
                        double max_affordability_in_year = 0.0;
                        
                        for (int w = (int) round(y * WEEKS_IN_YEAR);
                             w < (int) min((int) affordability_vec.size(), 
                                          (int) round((y + 1) * WEEKS_IN_YEAR)); ++w) {
                            max_affordability_in_year = max(max_affordability_in_year, affordability_vec[w]);
                        }
                        
                        // Accumulate for averaging across realizations
                        year_avg_affordability[y] += max_affordability_in_year;
                        year_realization_count[y]++;
                    }
                }
            }
        }
        
        // Calculate average across realizations for each year, then average across years for this WSS
        double wss_avg_affordability = 0.0;
        int valid_years = 0;
        
        for (unsigned long y = 0; y < n_years; ++y) {
            if (year_realization_count[y] > 0) {
                double year_affordability = (year_avg_affordability[y] / year_realization_count[y]) * affordability_scale;
                // if (affordability_out.good()) {
                //     affordability_out << wss_idx << "," << y << "," << year_affordability << "\n";
                // }
                wss_avg_affordability += year_affordability;
                valid_years++;
            }
        }
        
        if (valid_years > 0) {
            wss_avg_affordability /= valid_years;
            wss_affordability_avg.push_back(wss_avg_affordability);
        }
        wss_idx++;
    }
    
    if (wss_affordability_avg.empty()) {
        #ifdef PARALLEL
        printf("WARNING: wss_affordability_avg is EMPTY - returning 0.0 affordability!\n");
        #endif
        return 0.0;
    }
    
    double result_affordability;
    
    if (aggregation_method == 1) { // AVERAGE
        result_affordability = accumulate(wss_affordability_avg.begin(), 
                                         wss_affordability_avg.end(), 0.0) / wss_affordability_avg.size();
    } else { // MAX (worst case) - default
        result_affordability = *max_element(wss_affordability_avg.begin(),
                                           wss_affordability_avg.end());
    }
    
    if (std::isinf(result_affordability)) {
        string error_inf = "Infinite affordability index (WSS-level configurable).";
        throw logic_error(error_inf.c_str());
    }
    
    return result_affordability;  // Values already scaled by 2.44/1e6.
}
