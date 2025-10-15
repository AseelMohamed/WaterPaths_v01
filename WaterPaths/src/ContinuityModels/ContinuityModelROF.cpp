//
// Created by bernardo on 1/26/17.
//

#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "ContinuityModelROF.h"
#include "../Utils/Utils.h"

ContinuityModelROF::ContinuityModelROF(vector<WaterSource *> water_sources,
                                       const Graph &water_sources_graph,
                                       const vector<vector<int>> &water_sources_to_wss,
                                       vector<WaterSupplySystems *> wss,
                                       vector<MinEnvFlowControl *> min_env_flow_controls,
                                       vector<double> &wss_rdm,
                                       vector<double> &water_sources_rdm,
                                       unsigned long total_weeks_simulation,
                                       const int use_precomputed_rof_tables,
                                       const unsigned long realization_id)
        : ContinuityModel(water_sources, wss, min_env_flow_controls,
                          water_sources_graph, water_sources_to_wss,
                          wss_rdm,
                          water_sources_rdm,
                          realization_id),
          n_topo_sources((int) sources_topological_order.size()),
          use_precomputed_rof_tables(use_precomputed_rof_tables) {
    // update wss' total stored volume
    for (WaterSupplySystems *u : this->continuity_wss) {
        u->updateTotalAvailableVolume();
        u->getOwner()->setNoFinancialCalculations();
    }

    for (int u = 0; u < n_wss; ++u) {
        if (use_precomputed_rof_tables != IMPORT_ROF_TABLES) {
            ut_storage_to_rof_table.emplace_back(
                    total_weeks_simulation,
                    (unsigned long) NUMBER_REALIZATIONS_ROF);
        }
    }

    // Record which sources have no downstream sources.
    storage_wout_downstream = new bool[sources_topological_order.size()];
    for (int ws : sources_topological_order)
        storage_wout_downstream[ws] = downstream_sources[ws] != NON_INITIALIZED;

    // Get next online downstream source for each source.
    online_downstream_sources = getOnlineDownstreamSources();

    // Calculate wss' base delta storage corresponding to one table
    // tier and status-quo base storage capacity.
    if (use_precomputed_rof_tables == IMPORT_ROF_TABLES) {
        for (int u = 0; u < n_wss; ++u) {
            wss_base_storage_capacity.push_back(
                    continuity_wss[u]->getTotal_storage_capacity() *
                    BASE_STORAGE_CAPACITY_MULTIPLIER);
            wss_base_delta_capacity_table.push_back(
                    wss_base_storage_capacity[u] /
                    NO_OF_INSURANCE_STORAGE_TIERS);
        }

        current_and_base_storage_capacity_ratio =
                vector<double>((unsigned long) n_wss);
        current_storage_table_shift =
                vector<double>((unsigned long) n_wss);
    }
}

ContinuityModelROF::~ContinuityModelROF() {
    // Prevent base class destructor from deleting shared objects
    // by nullifying pointers to shared resources
    for (auto& ws : continuity_water_sources) {
        ws = nullptr;  // Null out pointers to prevent deletion
    }
    for (auto& u : continuity_wss) {
        u = nullptr;   // Null out pointers to prevent deletion  
    }
    for (auto& mef : min_env_flow_controls) {
        mef = nullptr; // Null out pointers to prevent deletion
    }
    
    // CRITICAL: Also null out realization_wss to prevent double deletion
    // These WSS objects are owned by the realization model
    for (auto& wss : realization_wss) {
        wss = nullptr; // Null out pointers to prevent deletion
    }
    
    delete[] storage_wout_downstream;
}

/**
 * Runs one the full rof calculations for realization #realization_id for a
 * given week.
 * @param week for which rof is to be calculated.
 */
vector<double> ContinuityModelROF::calculateLongTermROF(int week) {
    // Validate input week
    if (week < 0) {
        printf("ERROR: Negative week passed to calculateLongTermROF: week=%d\n", week);
        throw std::runtime_error("Negative week in calculateLongTermROF");
    }
    
    // vector where risks of failure will be stored.
    vector<double> risk_of_failure((unsigned long) n_wss, 0.0);
    vector<double> year_failure((unsigned long) n_wss, 0.0);

    // checks if new infrastructure became available and, if so, set the
    // corresponding realization
    // infrastructure online.
    updateOnlineInfrastructure(week);

    // perform a continuity simulation for NUMBER_REALIZATIONS_ROF (50) yearly
    // realization.
    for (int yr = 0; yr < NUMBER_REALIZATIONS_ROF; ++yr) {
        // reset current reservoirs' and wss' storage and combined
        // storage, respectively, in the corresponding realization simulation.
        resetWSSAndReservoirs(LONG_TERM_ROF);

        for (int w = 0; w < WEEKS_ROF_LONG_TERM; ++w) {
            // one week continuity time-step.
            int actual_week = w + week;
            if (actual_week < 0) {
                printf("ERROR: Negative week detected in long-term ROF: w=%d, week=%d, actual_week=%d\n", 
                       w, week, actual_week);
                throw std::runtime_error("Negative week in long-term ROF calculation");
            }
            continuityStep(actual_week, yr, APPLY_DEMAND_BUFFER);

            // check total available storage for each utility and, if smaller
            // than the fail ration, increase the number of failed years of
            // that utility by 1 (FAILURE).
            for (int u = 0; u < n_wss; ++u) {
                auto storage_condition =
                        continuity_wss[u]->getStorageToCapacityRatio() <=
                        STORAGE_CAPACITY_RATIO_FAIL;
                auto treatment_condition =
                        continuity_wss[u]->getUnrestrictedDemand() >
                        0.9 * continuity_wss[u]->getTotal_treatment_capacity();
                if (storage_condition || treatment_condition) {
                    year_failure[u] = FAILURE;
                }
            }
        }

        // Count failures and reset failures counter.
        for (int uu = 0; uu < n_wss; ++uu) {
            risk_of_failure[uu] += year_failure[uu];
            year_failure[uu] = NON_FAILURE;
        }
    }

    // Finish ROF calculations
    for (int i = 0; i < n_wss; ++i) {
        risk_of_failure[i] /= NUMBER_REALIZATIONS_ROF;
    }

    return risk_of_failure;
}

vector<double> ContinuityModelROF::calculateShortTermROF(int week,
                                                         int import_export_rof_tables) {
    // Validate input week
    if (week < 0) {
        printf("ERROR: Negative week passed to calculateShortTermROF: week=%d\n", week);
        throw std::runtime_error("Negative week in calculateShortTermROF");
    }
    
    vector<double> risk_of_failure;
    if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        return ContinuityModelROF::calculateShortTermROFTable(week);
    } else {
        return ContinuityModelROF::calculateShortTermROFFullCalcs(week);
    }
}

/**
 * Runs one the full rof calculations for realization #realization_id for a
 * given week.
 * @param week for which rof is to be calculated.
 */
#pragma GCC optimize("O0")
vector<double> ContinuityModelROF::calculateShortTermROFTable(int week) {
    // vector where risks of failure will be stored.
    auto n_wss = realization_wss.size();
    vector<double> risk_of_failure(n_wss, 0.0);
    double m;
    for (int u = 0; u < n_wss; ++u) {
        // Get current stored volume for utility u.
        double utility_storage =
                realization_wss[u]->getTotal_stored_volume();
        // Ratio of current and status-quo utility storage capacities
        //        double m = current_and_base_storage_capacity_ratio[u];
        m = realization_wss[u]->getTotal_storage_capacity() /
                wss_base_storage_capacity[u];
        // Calculate base table tier that contains the desired ROF by
        // shifting the table around based on new infrastructure -- the
        // shift is made by the part (m - 1) * STORAGE_CAPACITY_RATIO_FAIL *
        // wss_base_storage_capacity[u] - current_storage_table_shift[u]
        double storage_convert = utility_storage +
                                 STORAGE_CAPACITY_RATIO_FAIL *
                                         wss_base_storage_capacity[u] *
                                 (1. - m) + current_storage_table_shift[u];
        double tier = (storage_convert * NO_OF_INSURANCE_STORAGE_TIERS /
                wss_base_storage_capacity[u]);
        // Mean ROF between the two tiers of the ROF table where
        // current storage is located.

        //TODO: Make this more efficient if profiling shows the need.
        for (int r = 0; r < NUMBER_REALIZATIONS_ROF; ++r) {
            if (realization_wss[u]->getUnrestrictedDemand(week) > 0.9 * realization_wss[u]->getTotal_treatment_capacity()) {
                risk_of_failure[u] = 1.;
            } else {
                auto x = ut_storage_to_rof_table[u](week, r);
                if (tier < (double) NO_OF_INSURANCE_STORAGE_TIERS - x) {
                    risk_of_failure[u] += 0.5 / NUMBER_REALIZATIONS_ROF;
                }
                if (tier < (double) NO_OF_INSURANCE_STORAGE_TIERS - x + 1) {
                    risk_of_failure[u] += 0.5 / NUMBER_REALIZATIONS_ROF;
                }
            }
        }
    }

    return risk_of_failure;
}

/**
 * Runs one the full rof calculations for realization #realization_id for a
 * given week.
 * @param week for which rof is to be calculated.
 */
vector<double> ContinuityModelROF::calculateShortTermROFFullCalcs(int week) {
    // vector where risks of failure will be stored.
    vector<double> risk_of_failure((unsigned long) n_wss, 0.0);
    vector<double> year_failure((unsigned long) n_wss, 0.0);
    double *to_full = new double[n_sources];

    // Empty volumes are later used to update ROF tables.
    calculateEmptyVolumes(realization_water_sources, to_full);

    // checks if new infrastructure became available and, if so, set the
    // corresponding realization infrastructure online.
    updateOnlineInfrastructure(week);

    // perform a continuity simulation for NUMBER_REALIZATIONS_ROF (50)
    // yearly realization.
    for (int yr = 0; yr < NUMBER_REALIZATIONS_ROF; ++yr) {

        // reset current reservoirs' and wss' storage and combined
        // storage, respectively, in the corresponding realization simulation.
        resetWSSAndReservoirs(SHORT_TERM_ROF);

        for (int w = 0; w < WEEKS_ROF_SHORT_TERM; ++w) {
            // one week continuity time-step.
            int actual_week = w + week;
            if (actual_week < 0) {
                printf("ERROR: Negative week detected in short-term ROF: w=%d, week=%d, actual_week=%d\n", 
                       w, week, actual_week);
                throw std::runtime_error("Negative week in short-term ROF calculation");
            }
            continuityStep(actual_week, yr, !APPLY_DEMAND_BUFFER);

            // check total available storage for each utility and, if smaller
            // than the fail ration, increase the number of failed years of
            // that utility by 1 (FAILURE).
            for (int u = 0; u < n_wss; ++u) {
//                if (continuity_wss[u]->getStorageToCapacityRatio() <=
//                    STORAGE_CAPACITY_RATIO_FAIL) {
//                    year_failure[u] = FAILURE;
//                }
                auto storage_condition =
                        continuity_wss[u]->getStorageToCapacityRatio() <=
                        STORAGE_CAPACITY_RATIO_FAIL;
                auto treatment_condition =
                        continuity_wss[u]->getUnrestrictedDemand() >
                        0.9 * continuity_wss[u]->getTotal_treatment_capacity();
                if (storage_condition || treatment_condition) {
                    year_failure[u] = FAILURE;
                }
            }

            // calculated week of storage-rof table
            updateStorageToROFTable(INSURANCE_SHIFT_STORAGE_CURVES_THRESHOLD,
                                    week, to_full, yr);
        }

        // Count failures and reset failures counter.
        for (int uu = 0; uu < n_wss; ++uu) {
            risk_of_failure[uu] += year_failure[uu];
            year_failure[uu] = NON_FAILURE;
        }
    }

    // Finish ROF calculations
    for (int u = 0; u < n_wss; ++u) {
        risk_of_failure[u] /= NUMBER_REALIZATIONS_ROF;
        if (std::isnan(risk_of_failure[u])) {
            string error_m = "nan rof imported tables. Realization " +
                             to_string(realization_id) + ", week " +
                             to_string(week) + ", utility " + to_string(u);
            printf("%s", error_m.c_str());
            throw_with_nested(logic_error(error_m.c_str()));
        }
    }

    delete[] to_full;
    return risk_of_failure;
}

/**
 * Updates approximate ROF table based on continuity realization ran for
 * simulation based ROF calculations.
 * @param storage_percent_decrement
 * @param week_of_the_year
 * @param to_full empty volume of all reservoir in ID order.
 */
#pragma GCC optimize("O3")
void ContinuityModelROF::updateStorageToROFTable(
        double storage_percent_decrement, int week,
        const double *to_full, int rof_realization_number) {
    double *available_volumes = new double[n_sources];
    for (int ws = 0; ws < n_sources; ++ws) {
        available_volumes[ws] =
                continuity_water_sources[ws]->getAvailableSupplyVolume();
    }

    // loops over the percent storage levels to populate table. The loop
    // begins from one level above the level  where at least one failure was
    // observed in the last iteration. This saves a lot of computational time.
    for (int s = beginning_tier; s <= NO_OF_INSURANCE_STORAGE_TIERS; ++s) {
        // calculate delta storage for all reservoirs and array that will
        // receive the shifted storage curves.
        double percent_decrement_storage_level =
                (double) s * storage_percent_decrement;
        double *delta_storage = new double[n_sources];
        double *available_volumes_shifted = new double[n_sources];
        memcpy(available_volumes_shifted, available_volumes,
               sizeof(double) * n_sources);

        // calculate the difference between the simulated available water and
        // the one for the table calculated above based on the percent
        // decrement.
        for (int wss = 0; wss < n_sources; ++wss) {
            delta_storage[wss] = to_full[wss] - water_sources_capacities[wss] *
                                                percent_decrement_storage_level;
        }

        // Shift storages.
        shiftStorages(available_volumes_shifted, delta_storage);

        // Checks for wss failures.
        int count_fails = 0;
        // printf("DEBUG: Starting WSS failure checks, n_wss = %d, realization_wss.size() = %lu\n", 
            //    n_wss, realization_wss.size());
        for (int u = 0; u < n_wss; ++u) {
            // printf("DEBUG: Processing WSS %d\n", u);
            double utility_storage = 0;
            // Calculate combined stored volume for each utility based on
            // shifted storages.
            // printf("DEBUG: water_sources_online_to_wss[%d].size() = %lu\n", u, water_sources_online_to_wss[u].size());
            for (int ws : water_sources_online_to_wss[u]) {
                // printf("DEBUG: Processing water source %d for WSS %d\n", ws, u);
                bool has_treatment = realization_wss[u]->hasTreatmentConnected(ws);
                utility_storage += available_volumes_shifted[ws] *
                                   continuity_water_sources[ws]->getSupplyAllocatedFraction(
                                           u) *
                                   (has_treatment && realization_water_sources[ws]->isOnline());
            }
            // printf("DEBUG: Completed WSS %d\n", u);

            // Register failure in the table for each utility meeting
            // failure criteria. The treatment capacity criterion is INTENTIONALLY
            // missing, as it will be accounted for when using the tables--the
            // tables can theoretically only account for storage risk.
            if (utility_storage / wss_capacities[u] < STORAGE_CAPACITY_RATIO_FAIL) {
                auto p = ut_storage_to_rof_table[u].getPointerToElement(
                        week,
                        rof_realization_number
                        );
                if (*p == (double) NONE) {
                    *p = s;
                } else {
                    *p = min(*p, s);
                }
                count_fails++;
            }
        }
        delete[] delta_storage;
        delete[] available_volumes_shifted;
    }
    delete[] available_volumes;
}

//FIXME: MAKE THIS MORE EFFICIENT. THIS METHOD IS THE MOST EXPENSIVE ONE IN THE CODE.
#pragma GCC optimize("O3")
void ContinuityModelROF::shiftStorages(
        double *available_volumes_shifted,
        const double *delta_storage) {
    // Add deltas to all sources following the topological order, so that
    // upstream is calculated before downstream.
    for (int ws : sources_topological_order) {
        // calculate initial estimate for shifted
        available_volumes_shifted[ws] += delta_storage[ws];

        double available_volume_to_full = water_sources_capacities[ws] -
                                          available_volumes_shifted[ws];

        // if not full, retrieve spill to downstream source.
        if (available_volume_to_full > 0) {
            // Calculate spilled water. Since the curves are shifted as the
            // weeks of the rof realizations are calculated, the minimum
            // environmental outflows below will be the ones at the time when
            // the storage is being shifted.
            double spillage =
                    continuity_water_sources[ws]->getTotal_outflow() -
                    continuity_water_sources[ws]
                            ->getMin_environmental_outflow();

            double spillage_retrieved = min(available_volume_to_full, spillage);

            available_volumes_shifted[ws] += spillage_retrieved;

            if (online_downstream_sources[ws] > 0)
                available_volumes_shifted[online_downstream_sources[ws]] -=
                        spillage_retrieved;
        } else {
            double spillage = -available_volume_to_full;
            if (online_downstream_sources[ws] > 0) {
                available_volumes_shifted[ws] -= spillage;
                available_volumes_shifted[online_downstream_sources[ws]] +=
                        spillage;
            }
        }
    }
}


/**
 * Prints a binary file with the rof_table for a given realization in a
 * given week.
 * @param week
 */
void ContinuityModelROF::printROFTable(const string &folder) {
    for (int u = 0; u < n_wss; ++u) {

        string file_name =
                folder + "tables_r" + to_string(realization_id) + "_u" +
                to_string(u) + ".csv";
        ofstream output_file(file_name);

        auto num_weeks = ut_storage_to_rof_table[u].get_i();
        for (int w = 0; w < num_weeks; ++w) {
            auto data = ut_storage_to_rof_table[u].getPointerToElement(w, 0);
            std::ostringstream week_table;
            week_table << std::fixed;
            week_table << std::setprecision(0);
            for (int t = 0; t < NO_OF_INSURANCE_STORAGE_TIERS; ++t) {
                week_table << to_string(data[t]) + ",";
            }

            string line = week_table.str();
            line.pop_back();
            output_file << line;
            if (w < num_weeks - 1)
                output_file << endl;
        }

        output_file.close();
    }
}

/**
 * reset reservoirs' and wss' storage and last release, and
 * combined storage, respectively, they currently have in the
 * corresponding realization simulation.
 */
void ContinuityModelROF::resetWSSAndReservoirs(int rof_type) {
    // update water sources info. If short-term rof, return to current
    // storage; if long-term, make them full.
    if (rof_type == SHORT_TERM_ROF)
        for (int i = 0; i < n_sources; ++i) {   // Current available volume
            continuity_water_sources[i]->setAvailableAllocatedVolumes
                    (realization_water_sources[i]
                             ->getAvailable_allocated_volumes(),
                     realization_water_sources[i]->getAvailableVolume());
            continuity_water_sources[i]->setOutflow_previous_week(
                    realization_water_sources[i]->getTotal_outflow());
        }
    else
        for (int i = 0; i < n_sources; ++i) {   // Full capacity
            continuity_water_sources[i]->setFull();
            continuity_water_sources[i]->setOutflow_previous_week(
                    realization_water_sources[i]->getTotal_outflow());
        }

    // update wss combined storage.
    for (WaterSupplySystems *u : continuity_wss) {
        u->updateTotalAvailableVolume();
    }
}

/**
 * Pass to the rof continuity model the locations of the wss
 * of the realization it calculated rofs for.
 * @param realization_water_sources
 */
void ContinuityModelROF::connectRealizationWaterSources(
        const vector<WaterSource *> &realization_water_sources) {
    ContinuityModelROF::realization_water_sources =
            realization_water_sources;
}

/**
 * Pass to the rof continuity model the locations of the wss
 * of the realization it calculated rofs for.
 * @param realization_wss
 */
void ContinuityModelROF::connectRealizationWSS(
        const vector<WaterSupplySystems *> &realization_wss) {
    ContinuityModelROF::realization_wss = realization_wss;
}

/**
 * Checks if new infrastructure became online.
 */
#pragma GCC optimize("O0")
void ContinuityModelROF::updateOnlineInfrastructure(int week) {
    for (unsigned long ws = 0; ws < (unsigned long) n_sources; ++ws) {
        // Check if any infrastructure option is online in the
        // realization model and not in the ROF model.
        if (realization_water_sources.at(ws)->isOnline() &&
            !continuity_water_sources.at(ws)->isOnline()) {
            // If so, set it online in the ROF calculation model.
            for (int uu : wss_to_water_sources[ws]) {
                auto u = (unsigned long) uu;
                if (continuity_water_sources[ws]->source_type < NON_STRUCTURAL_SOURCES) {
                    water_sources_online_to_wss.at(u).push_back((int) ws);
                }
                continuity_wss.at(u)->setWaterSourceOnline((int) ws, week);

                // Update the shift in storage to be used to calculate the
                // tier in precomputed ROF tables corresponding to the
                // current storage of a given utility.
                if (use_precomputed_rof_tables == IMPORT_ROF_TABLES) {
                    //TODO: SHIFT NÃO ATUALIZADO QUANDO ETA É CONSTRUIDA EM UM RESERVATÓRIO EXISTENTE MAS DESCONECTADO OU SEM CAPACIDADE DE TRATAMENTO. PROVAVELMENTE OUTRAS COISAS TAMBÉM NÃO ESTÃO SENDO ATUALIZADAS.
                    current_storage_table_shift.at(u) += table_storage_shift.at(u).at(ws);
                }
            }

            // Update water source capacities in case a reservoir expansion
            // was built.
            water_sources_capacities.at(ws) =
                    continuity_water_sources.at(ws)->getSupplyCapacity();
        }
    }

    // Update wss' storage capacities and their ratios to status-quo
    // capacities in case new infrastructure has been built.
    if (Utils::isFirstWeekOfTheYear(week) || week == 0) {
        for (unsigned long u = 0; u < (unsigned long) n_wss; ++u) {
            wss_capacities.at(u) =
                    continuity_wss.at(u)->getTotal_storage_capacity();
        }

        if (use_precomputed_rof_tables == IMPORT_ROF_TABLES) {
            for (unsigned long u = 0; u < (unsigned long) n_wss; ++u) {
                current_and_base_storage_capacity_ratio.at(u) =
                        wss_capacities.at(u) /
                        wss_base_storage_capacity.at(u);
            }
        }
    }

    // Update list of downstream sources of each source.
    online_downstream_sources = getOnlineDownstreamSources();
}

void ContinuityModelROF::setROFTablesAndShifts(
        const vector<Matrix2D<int>> &storage_to_rof_table,
        const vector<vector<double>> &table_storage_shift) {
    this->ut_storage_to_rof_table = storage_to_rof_table;
    this->table_storage_shift = table_storage_shift;
}


vector<Matrix2D<int>> &ContinuityModelROF::getUt_storage_to_rof_table() {
    return ut_storage_to_rof_table;
}

/**
 * Calculate empty volume in storage-based water sources. This information is later used for updating the ROF tables.
 * @param realization_water_sources
 * @param to_full
 */
void ContinuityModelROF::calculateEmptyVolumes(
        vector<WaterSource *> &realization_water_sources, double *to_full) {
    for (int ws = 0; ws < n_sources; ++ws) {
        if (realization_water_sources[ws]->isOnline()) {
            to_full[ws] = realization_water_sources[ws]->getSupplyCapacity() -
                          realization_water_sources[ws]->getAvailableSupplyVolume();
        } else {
            to_full[ws] = 0;
        }
    }
}
