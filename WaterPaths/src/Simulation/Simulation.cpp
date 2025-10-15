//
// Created by bernardo on 1/25/17.
//

#include "Simulation.h"
#include "../Utils/Utils.h"
#include <algorithm>
#include <omp.h>
#include <set>
#include <cstdio>

#ifdef  PARALLEL
#include <mpi.h>
#endif

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


Simulation::Simulation(
        vector<WaterSource *> &water_sources, Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_utilities,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &utilities_rdm,
        vector<vector<double>> &water_sources_rdm,
        vector<vector<double>> &policies_rdm,
        const unsigned long total_simulation_time,
        vector<unsigned long> &realizations_to_run) :
        total_simulation_time(total_simulation_time),
        realizations_to_run(realizations_to_run),
        import_export_rof_tables(DO_NOT_EXPORT_OR_IMPORT_ROF_TABLES),
        n_realizations(realizations_to_run.size()),
        water_sources(water_sources),
        water_sources_graph(water_sources_graph),
        water_sources_to_utilities(water_sources_to_utilities),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        utilities_rdm(utilities_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm) {
    setupSimulation(
            water_sources, water_sources_graph,
            water_sources_to_utilities, utilities, drought_mitigation_policies,
            min_env_flow_controls,
            utilities_rdm, water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

Simulation::Simulation(
        vector<WaterSource *> &water_sources, Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_utilities,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &utilities_rdm,
        vector<vector<double>> &water_sources_rdm,
        vector<vector<double>> &policies_rdm,
        const unsigned long total_simulation_time,
        vector<unsigned long> &realizations_to_run,
        vector<vector<Matrix2D<int>>> &precomputed_rof_tables,
        vector<vector<double>> &table_storage_shift,
        string &rof_tables_folder) :
        total_simulation_time(total_simulation_time),
        realizations_to_run(realizations_to_run),
        import_export_rof_tables(IMPORT_ROF_TABLES),
        n_realizations(realizations_to_run.size()),
        water_sources(water_sources),
        water_sources_graph(water_sources_graph),
        water_sources_to_utilities(water_sources_to_utilities),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        utilities_rdm(utilities_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm),
        precomputed_rof_tables(&precomputed_rof_tables),
        table_storage_shift(&table_storage_shift) {
    setRof_tables_folder(rof_tables_folder);

    setupSimulation(
            water_sources, water_sources_graph,
            water_sources_to_utilities, utilities, drought_mitigation_policies,
            min_env_flow_controls,
            utilities_rdm, water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

Simulation::Simulation(
        vector<WaterSource *> &water_sources,
        Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_utilities,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &utilities_rdm,
        vector<vector<double>> &water_sources_rdm,
        vector<vector<double>> &policies_rdm,
        const unsigned long total_simulation_time,
        vector<unsigned long> &realizations_to_run,
        string &rof_tables_folder) :
        total_simulation_time(total_simulation_time),
        realizations_to_run(realizations_to_run),
        import_export_rof_tables(EXPORT_ROF_TABLES),
        n_realizations(realizations_to_run.size()),
        water_sources(water_sources),
        water_sources_graph(water_sources_graph),
        water_sources_to_utilities(water_sources_to_utilities),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        utilities_rdm(utilities_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm) {
    setRof_tables_folder(rof_tables_folder);

    setupSimulation(
            water_sources,
            water_sources_graph,
            water_sources_to_utilities,
            utilities,
            drought_mitigation_policies,
            min_env_flow_controls,
            utilities_rdm,
            water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

void Simulation::setupSimulation(vector<WaterSource *> &water_sources,
                                 Graph &water_sources_graph,
                                 const vector<vector<int>> &water_sources_to_utilities,
                                 vector<Utility *> &utilities,
                                 const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
                                 vector<MinEnvFlowControl *> &min_env_flow_controls,
                                 vector<vector<double>> &utilities_rdm,
                                 vector<vector<double>> &water_sources_rdm,
                                 vector<vector<double>> &policies_rdm,
                                 vector<unsigned long> &realizations_to_run) {
    // Sort water sources and utilities by their IDs.
    //FIXME: THERE IS A STUPID MISTAKE HERE IN THE SORT FUNCTION THAT IS PREVENTING IT FROM WORKING UNDER WINDOWS AND LINUX.
    std::sort(water_sources.begin(), water_sources.end(), WaterSource::compare);
    std::sort(utilities.begin(), utilities.end(), Utility::compById);

    // Check if IDs are sequential.
    for (int ws = 1; ws < (int) water_sources.size(); ++ws) {
        if (water_sources[ws]->id != water_sources[ws - 1]->id + 1) {
            cout << "The IDs of water sources " << water_sources[ws]->id << " "
                                                                            "and "
                 << water_sources[ws - 1]->id << " do not follow a "
                                                 "unit progression." << endl;
            throw_with_nested(
                    invalid_argument("Improper water source ID sequencing"));
        }
    }

    for (int u = 1; u < (int) utilities.size(); ++u) {
        if (utilities[u]->id != utilities[u - 1]->id + 1) {
            cout << "The IDs of utilities " << utilities[u]->id << " "
                                                                   "and "
                 << utilities[u - 1]->id << " do not follow a "
                                            "unit progression." << endl;
            throw_with_nested(
                    invalid_argument("Improper utility ID sequencing"));
        }
    }

    // Check if sources listed in construction order array of a utility are
    // listed as belonging to that utility
    for (int u = 0; u < (int) utilities.size(); ++u) {
        // Create a vector with rof and demand triggered infrastructure for
        // utility u.
        vector<int> demand_rof_infra_order =
                utilities[u]->getRof_infrastructure_construction_order();
        demand_rof_infra_order.insert(
                demand_rof_infra_order.begin(),
                utilities[u]->getDemand_infra_construction_order().begin(),
                utilities[u]->getDemand_infra_construction_order().end());
        // Iterate over demand and rof combined infrastructure vector
        // looking for sources declared as to be constructed that were not
        // declared as belonging to utility u.
        for (int ws :
                demand_rof_infra_order)
            if (std::find(water_sources_to_utilities[u].begin(),
                          water_sources_to_utilities[u].end(),
                          ws)
                == water_sources_to_utilities[u].end()) {
                cout << "Water source #" << ws << " is listed in the "
                                                  "construction order for utility "
                     << utilities[u]->id
                     << " (" << utilities[u]->name << ")  but is  not  "
                                                      "present in  utility's list of water sources."
                     << endl;
                throw invalid_argument("Utility's construction order and "
                                       "owned sources mismatch.");
            }

        for (int ws : water_sources_to_utilities[u])
            if (find_if(water_sources.begin(),
                        water_sources.end(),
                        [&ws](
                                const WaterSource *
                                obj) { return obj->id == ws; }) ==
                water_sources.end()) {
                cout << "Water source #" << ws << " not present in "
                                                  "comprehensive water sources vector."
                     << endl;
                throw invalid_argument("Water sources declared to belong to"
                                       " a utility is not present "
                                       "in vector of water sources.");
            }
    }

    // Water source connections are handled by the WaterSupplySystems constructors
    // and their infrastructure managers. No additional setup needed here.

    // Creates the data collector for the simulation.
    master_data_collector = new MasterDataCollector(realizations_to_run);
}

Simulation::~Simulation() = default;

/**
 * Assignment constructor
 * @param simulation
 * @return
 * @todo implement assignment constructor.
 */
Simulation &Simulation::operator=(const Simulation &simulation) {
    this->n_realizations = simulation.n_realizations;
    return *this;
}

void Simulation::createContinuityModels(unsigned long realization,
                                        ContinuityModelRealization *&realization_model,
                                        ContinuityModelROF *&rof_model) {
    // Create realization models by copying the water supply systems but sharing water sources.
    // Water sources should be shared across realizations, not duplicated.
    vector<WaterSource *> water_sources_realization = water_sources; // Share, don't copy
    vector<DroughtMitigationPolicy *> drought_mitigation_policies_realization =
            Utils::copyDroughtMitigationPolicyVector(
                    drought_mitigation_policies);
    
    // Extract water supply systems from utilities for realization model
    vector<WaterSupplySystems *> wss_realization;
    vector<vector<int>> water_sources_to_wss_mapping;
    
    // printf("DEBUG: Creating WSS mapping for %zu utilities\n", utilities.size());
    
    for (auto* utility : utilities) {
        // printf("DEBUG: Processing utility %d with %zu WSS\n", utility->id, utility->getWaterSupplySystems().size());
        for (const auto& wss : utility->getWaterSupplySystems()) {
            // Create copies of WSS for this realization
            wss_realization.push_back(new WaterSupplySystems(*wss));
            
            // Create mapping for this WSS - initially empty, will be populated by addWaterSource calls
            water_sources_to_wss_mapping.push_back(vector<int>());
            // printf("DEBUG: Added WSS %d to mapping, total WSS count: %zu\n", wss->getSystemId(), water_sources_to_wss_mapping.size());
        }
    }
    
    // printf("DEBUG: Final water_sources_to_wss_mapping.size() = %zu\n", water_sources_to_wss_mapping.size());
    
    // Now populate the water sources to WSS mapping based on utility-level mapping
    for (int utility_id = 0; utility_id < utilities.size(); ++utility_id) {
        auto* utility = utilities[utility_id];
        int wss_start_index = 0;
        
        // Find the starting index for this utility's WSS in the flattened list
        for (int u = 0; u < utility_id; ++u) {
            wss_start_index += utilities[u]->getWaterSupplySystems().size();
        }
        
        // printf("DEBUG: Utility %d starts at WSS index %d\n", utility_id, wss_start_index);
        // printf("DEBUG: Utility %d has %zu water sources: ", utility_id, water_sources_to_utilities[utility_id].size());
        // Removed debug printing of water source IDs
        
        // Distribute water sources from utility-level mapping to WSS-level mapping
        for (int ws_id : water_sources_to_utilities[utility_id]) {            
            int target_wss_index = wss_start_index; // Default to first WSS
            
            if (utility->getWaterSupplySystems().size() == 2) {
                // CAESB case: 2 WSS within the utility
                if (ws_id == 1 || ws_id == 3 || ws_id == 4) {
                    target_wss_index = wss_start_index + 1; // TortoSM (second WSS)
                } else {
                    target_wss_index = wss_start_index; // Descoberto (first WSS)
                }
            }
            
            // printf("DEBUG: Assigning water source %d to WSS index %d\n", ws_id, target_wss_index);
            water_sources_to_wss_mapping[target_wss_index].push_back(ws_id);
        }
    }
    
    // printf("DEBUG: Final mapping:\n");
    // Removed debug printing of WSS water source assignments
    
    vector<MinEnvFlowControl *> min_env_flow_controls_realization =
            Utils::copyMinEnvFlowControlVector(min_env_flow_controls);

    // Store realization models in vector - using WSS instead of utilities
    realization_model = new ContinuityModelRealization(
            water_sources_realization,
            water_sources_graph,
            water_sources_to_wss_mapping,  // Now properly maps water sources to individual WSS
            wss_realization,
            drought_mitigation_policies_realization,
            min_env_flow_controls_realization,
            utilities_rdm.at(realization),
            water_sources_rdm.at(realization),
            policies_rdm.at(realization),
            (int) realization);

    // Create rof models by copying the water supply systems but sharing water sources.
    // Water sources should be shared across models, not duplicated.
    vector<WaterSource *> water_sources_rof = water_sources; // Share, don't copy
    
    // Extract water supply systems from utilities for ROF model
    vector<WaterSupplySystems *> wss_rof;
    for (auto* utility : utilities) {
        for (const auto& wss : utility->getWaterSupplySystems()) {
            // Create copies of WSS for ROF model
            wss_rof.push_back(new WaterSupplySystems(*wss));
        }
    }
    
    vector<MinEnvFlowControl *> min_env_flow_controls_rof =
            Utils::copyMinEnvFlowControlVector(min_env_flow_controls);

    // Store realization models in vector - using WSS instead of utilities
    rof_model = new ContinuityModelROF(
            water_sources_rof,
            water_sources_graph,
            water_sources_to_wss_mapping,  // Use the same WSS-level mapping as realization model
            wss_rof,
            min_env_flow_controls_rof,
            utilities_rdm.at(realization),
            water_sources_rdm.at(realization),
            total_simulation_time,
            import_export_rof_tables,
            realization);

    // Initialize rof models by connecting it to realization water sources.
    rof_model->connectRealizationWaterSources(water_sources_realization);
    
    // Connect ROF model to realization WSS
    rof_model->connectRealizationWSS(wss_realization);

    // Pass ROF tables to continuity model
    if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        rof_model->setROFTablesAndShifts(
                precomputed_rof_tables->at(realization), *table_storage_shift);
    }

    // Link storage-rof tables of policies and rof models.
    for (DroughtMitigationPolicy *dmp :
            realization_model->getDrought_mitigation_policies())
        dmp->setStorage_to_rof_table_(
                rof_model->getUt_storage_to_rof_table(),
                import_export_rof_tables);
}

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

MasterDataCollector *
Simulation::runFullSimulation(unsigned long n_threads, double *vars) {
    if (rof_tables_folder.length() == 0) {
        rof_tables_folder = "rof_tables";
    }

    // Check if number of imported tables corresponds to model.
    if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        if (precomputed_rof_tables->at(0).size() != utilities.size()) {
            throw invalid_argument(
                    "Different number of utilities in model and imported ROF tables.");
        }

        auto max_realization = *max_element(realizations_to_run.begin(),
                                            realizations_to_run.end()) + 1;
        auto n_precomputed_tables = precomputed_rof_tables->size();
        if (n_precomputed_tables != max_realization) {
            string error = "There are at least " + to_string(max_realization) + 
                          " potential realizations but " + to_string(n_precomputed_tables) + 
                          " imported ROF tables.";
            throw invalid_argument(error);
        }
    }

    set<unsigned long> s(realizations_to_run.begin(),
                         realizations_to_run.end());
    vector<unsigned long> realizations_to_run_unique;
    realizations_to_run_unique.assign(s.begin(), s.end());

    // Prepare error output.
    int had_catch = 0;
    string error_m = "Error in realizations ";
    string error_file_name = "error_reals";
    string error_file_content = "#";

    // Run realizations.
#pragma omp parallel for ordered num_threads(n_threads) shared(had_catch, realizations_to_run_unique, error_m, error_file_name, error_file_content) default(none)
    for (unsigned long r = 0; r < realizations_to_run_unique.size(); ++r) {
        unsigned long realization = realizations_to_run_unique[r];
        //printf("Realization %lu\n", r);

        // Create continuity models.
        ContinuityModelRealization *realization_model = nullptr;
        ContinuityModelROF *rof_model = nullptr;
        createContinuityModels(realization, realization_model, rof_model);

        // Initialize data collector.
        // For data collection, we still need utilities (for financial data)
        // Extract utilities from WSS back-references
        vector<Utility *> utilities_for_data_collection;
        for (const auto& wss : realization_model->getContinuity_wss()) {
            utilities_for_data_collection.push_back(wss->getOwner());
        }
        
        master_data_collector->addRealization(
                realization_model->getContinuity_water_sources(),
                realization_model->getDrought_mitigation_policies(),
                utilities_for_data_collection,
                realization);

//        try {
//        double start = omp_get_wtime();
            // printf("DEBUG: Starting simulation for realization %lu, total_simulation_time = %lu\n", 
                //    realization, total_simulation_time);
            for (int w = 0; w < (int) total_simulation_time; ++w) {
                if (w % 52 == 0) {  // Print every year
                    // printf("DEBUG: Processing week %d (year %d)\n", w, w/52);
                }
//                printf("%d\n", w);
                // DO NOT change the order of the steps. This would mess up
                // important dependencies.
                // Calculate long-term risk-of-failre if current week is first week of the year.
                if (Utils::isFirstWeekOfTheYear(w)) {
                    // printf("DEBUG: Calculating long-term ROF for week %d\n", w);
                    realization_model->setLongTermROFs(
                            rof_model->calculateLongTermROF(w), w);
                }
                // Calculate short-term risk-of-failure
                // printf("DEBUG: Calculating short-term ROF for week %d\n", w);
                realization_model->setShortTermROFs(
                        rof_model->calculateShortTermROF(w,
                                import_export_rof_tables));
                // Apply drought mitigation policies
                if (import_export_rof_tables != EXPORT_ROF_TABLES) {
                    realization_model->applyDroughtMitigationPolicies(w);
                }
                // Continuity calculations for current week
                realization_model->continuityStep(w);
                // Collect system data for output printing and objective calculations.
                if (import_export_rof_tables != EXPORT_ROF_TABLES) {
                    master_data_collector->collectData(realization);
                }
            }
            // Export ROF tables for future simulations of the same problem with the same states-of-the-world.
            if (import_export_rof_tables == EXPORT_ROF_TABLES) {
                rof_model->printROFTable(rof_tables_folder);
            }
//        printf("Realization %lu took %f seconds.\n", r, omp_get_wtime() - start);

// #pragma omp critical
//         printProgress(
//                 (double) master_data_collector->getRealizations_created() /
//                 (double) realizations_to_run_unique.size());

//        } catch (...) {
//#pragma omp atomic
//            ++had_catch;
//            error_m += to_string(realization) + " ";
//            error_file_name += "_" + to_string(realization);
//            error_file_content += to_string(realization) + ",";
//            master_data_collector->removeRealization(realization);
//        }

        // Delete ROF model first since it only references shared objects
        delete rof_model;
        // Delete realization model last since it owns the shared objects
        delete realization_model;
    }
    // Handle exception from the OpenMP region and pass it up to the
    // problem class.
    if (had_catch) {
        int world_rank;
#ifdef  PARALLEL
        int mpi_initialized;
        MPI_Initialized(&mpi_initialized);
        if (mpi_initialized)
                 MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        else
            world_rank = 0;
#else
        world_rank = 0;
#endif

        // Create error file
        error_file_name += ".csv";
        error_m += ". Error data in " + error_file_name;
        ofstream error_file;
        error_file.open(error_file_name);

        // Write error rile
        error_file << error_file_content << endl;
        for (int i = 0; i < NUM_DEC_VAR - 1; ++i) {
            error_file << vars[i] << ",";
        }
        error_file << vars[NUM_DEC_VAR - 1];

        // Finalize error reporting
        error_file.close();
        printf("%s", error_m.c_str());

//	master_data_collector->cleanCollectorsOfDeletedRealizations();
//        throw_with_nested(runtime_error(error_m.c_str()));
    }
    return master_data_collector;
}

void Simulation::setRof_tables_folder(const string &rof_tables_folder) {
    Simulation::rof_tables_folder = rof_tables_folder;
}

