//
// Created by bernardo on 1/25/17.
//

#include "Simulation.h"
#include "../Utils/Utils.h"
#include <algorithm>
#include <omp.h>
#include <set>
#include <cstdio>
#include <cassert>

#ifdef  PARALLEL
#include <mpi.h>
#endif

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


Simulation::Simulation(
        vector<WaterSource *> &water_sources, Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_wss,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &wss_rdm,
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
        water_sources_to_wss(water_sources_to_wss),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        wss_rdm(wss_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm) {
    setupSimulation(
            water_sources, water_sources_graph,
            water_sources_to_wss, utilities, drought_mitigation_policies,
            min_env_flow_controls,
            wss_rdm, water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

Simulation::Simulation(
        vector<WaterSource *> &water_sources, Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_wss,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &wss_rdm,
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
        water_sources_to_wss(water_sources_to_wss),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        wss_rdm(wss_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm),
        precomputed_rof_tables(&precomputed_rof_tables),
        table_storage_shift(&table_storage_shift) {
    setRof_tables_folder(rof_tables_folder);

    setupSimulation(
            water_sources, water_sources_graph,
            water_sources_to_wss, utilities, drought_mitigation_policies,
            min_env_flow_controls,
            wss_rdm, water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

Simulation::Simulation(
        vector<WaterSource *> &water_sources,
        Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_wss,
        vector<Utility *> &utilities,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_controls,
        vector<vector<double>> &wss_rdm,
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
        water_sources_to_wss(water_sources_to_wss),
        utilities(utilities),
        drought_mitigation_policies(drought_mitigation_policies),
        min_env_flow_controls(min_env_flow_controls),
        wss_rdm(wss_rdm),
        water_sources_rdm(water_sources_rdm),
        policies_rdm(policies_rdm) {
    setRof_tables_folder(rof_tables_folder);

    setupSimulation(
            water_sources,
            water_sources_graph,
            water_sources_to_wss,
            utilities,
            drought_mitigation_policies,
            min_env_flow_controls,
            wss_rdm,
            water_sources_rdm,
            policies_rdm,
            realizations_to_run);
}

void Simulation::setupSimulation(vector<WaterSource *> &water_sources,
                                 Graph &water_sources_graph,
                                 const vector<vector<int>> &water_sources_to_wss,
                                 vector<Utility *> &utilities,
                                 const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
                                 vector<MinEnvFlowControl *> &min_env_flow_controls,
                                 vector<vector<double>> &wss_rdm,
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
            printf("The IDs of water sources %d and %d do not follow a unit progression.\n",
                   water_sources[ws]->id, water_sources[ws - 1]->id);
            throw_with_nested(
                    invalid_argument("Improper water source ID sequencing"));
        }
    }

    for (int u = 1; u < (int) utilities.size(); ++u) {
        if (utilities[u]->id != utilities[u - 1]->id + 1) {
            printf("The IDs of utilities %d and %d do not follow a unit progression.\n",
                   utilities[u]->id, utilities[u - 1]->id);
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
            if (std::find(water_sources_to_wss[u].begin(),
                          water_sources_to_wss[u].end(),
                          ws)
                == water_sources_to_wss[u].end()) {
                printf("Water source #%d is listed in the construction order for utility %d (%s) but is not present in utility's list of water sources.\n",
                       ws, utilities[u]->id, utilities[u]->name);
                throw invalid_argument("Utility's construction order and "
                                       "owned sources mismatch.");
            }

        for (int ws : water_sources_to_wss[u])
            if (find_if(water_sources.begin(),
                        water_sources.end(),
                        [&ws](
                                const WaterSource *
                                obj) { return obj->id == ws; }) ==
                water_sources.end()) {
                printf("Water source #%d not present in comprehensive water sources vector.\n", ws);
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
    // Create realization models by copying the water sources and water supply systems.
    vector<WaterSource *> water_sources_realization =
            Utils::copyWaterSourceVector(water_sources);
    vector<DroughtMitigationPolicy *> drought_mitigation_policies_realization =
            Utils::copyDroughtMitigationPolicyVector(
                    drought_mitigation_policies);
    
    // Extract water supply systems from utilities for realization model
    vector<std::unique_ptr<WaterSupplySystems>> wss_realization;
    vector<MinEnvFlowControl *> min_env_flow_controls_realization;
    vector<vector<int>> water_sources_to_wss_mapping;

    // ALL copying and construction must be serialized
    // Using single critical section to prevent race conditions when multiple threads
    // copy from shared utility objects simultaneously
    #pragma omp critical(model_creation)
    {
        // Copy WSS from utilities
        for (auto* utility : utilities) {
            for (const auto& wss : utility->getWaterSupplySystems()) {
                // Create copies of WSS for this realization using unique_ptr
                wss_realization.push_back(std::make_unique<WaterSupplySystems>(*wss));
            }
        }
        
        // Directly use the WSS connectivity matrix
        water_sources_to_wss_mapping = water_sources_to_wss;
        
        min_env_flow_controls_realization = Utils::copyMinEnvFlowControlVector(min_env_flow_controls);

        // Store realization models in vector - using WSS instead of utilities
        realization_model = new ContinuityModelRealization(
                water_sources_realization,
                water_sources_graph,
                water_sources_to_wss_mapping,  // Now properly maps water sources to individual WSS
                std::move(wss_realization),
                drought_mitigation_policies_realization,
                min_env_flow_controls_realization,
                wss_rdm.at(realization),
                water_sources_rdm.at(realization),
                policies_rdm.at(realization),
                realization);
    }

    // ROF model construction must use SAME critical section name
    // to ensure it serializes with realization model construction above
    #pragma omp critical(model_creation)
    {    
        // Create rof models by copying the water sources and WSS.
        // ROF model needs its own copies for independent ROF calculations
        vector<WaterSource *> water_sources_rof =
                Utils::copyWaterSourceVector(water_sources);
        
        // Create fresh WSS copies for ROF model (from original utilities, not realization)
        vector<std::unique_ptr<WaterSupplySystems>> wss_rof;
        for (auto* utility : utilities) {
            for (const auto& wss : utility->getWaterSupplySystems()) {
                auto wss_copy = std::make_unique<WaterSupplySystems>(*wss);
                // Mark ROF WSS as NOT used for realization BEFORE passing to constructor
                // This prevents them from trying to access realization-specific demand data they don't need
                wss_copy->setUsedForRealization(false);
                wss_rof.push_back(std::move(wss_copy));
            }
        }
        
        vector<MinEnvFlowControl *> min_env_flow_controls_rof =
                Utils::copyMinEnvFlowControlVector(min_env_flow_controls);

        // Store realization models in vector - using WSS instead of utilities
        rof_model = new ContinuityModelROF(
                water_sources_rof,
                water_sources_graph,
                water_sources_to_wss_mapping,  // Use the same WSS-level mapping as realization model
                std::move(wss_rof),
                min_env_flow_controls_rof,
                wss_rdm.at(realization),
                water_sources_rdm.at(realization),
                total_simulation_time,
                import_export_rof_tables,
                realization);

        // Initialize rof models by connecting it to realization water sources and WSS for observation
        rof_model->connectRealizationWaterSources(water_sources_realization);
        
        // Connect ROF model to realization WSS (now owned by realization_model) INSIDE critical section
        // Access the WSS from the realization model that now owns them
        rof_model->connectRealizationWSS(realization_model->getContinuity_wss());
    }

    // Pass ROF tables to continuity model
    if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        rof_model->setROFTablesAndShifts(
                precomputed_rof_tables->at(realization), *table_storage_shift);
    }

    // Link storage-rof tables of policies and rof models.
    for (DroughtMitigationPolicy *dmp :
            realization_model->getDrought_mitigation_policies())
        dmp->setStorage_to_rof_table_(
                rof_model->getWSS_storage_to_rof_table(),
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
        // Count total number of WSS across all utilities
        size_t total_wss = 0;
        for (const auto* utility : utilities) {
            total_wss += utility->getWaterSupplySystems().size();
        }
        
        if (precomputed_rof_tables->at(0).size() != total_wss) {
            char error[512];
            sprintf(error, 
                    "Different number of WSS in model (%zu) and imported ROF tables (%zu).",
                    total_wss, precomputed_rof_tables->at(0).size());
            throw invalid_argument(error);
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
    
    // VALIDATION: Check if any requested realizations have missing/corrupted ROF tables
    if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        vector<unsigned long> valid_realizations;
        for (unsigned long realization : realizations_to_run_unique) {
            if (realization >= precomputed_rof_tables->size()) {
                printf("WARNING: Realization %lu requested but only %zu ROF tables available. Skipping.\n",
                       realization, precomputed_rof_tables->size());
                continue;
            }
            
            // Simple check: if the vector of tables for this realization is not empty, consider it valid
            // More detailed validation happens in setRofTables where empty files are caught
            if (!precomputed_rof_tables->at(realization).empty()) {
                valid_realizations.push_back(realization);
            } else {
                printf("WARNING: Realization %lu has no ROF tables loaded. Skipping this realization.\n", 
                       realization);
            }
        }
        
        if (valid_realizations.empty()) {
            throw std::runtime_error("ERROR: No valid realizations to run after filtering out corrupted ROF tables!");
        }
        
        // printf("Running %zu out of %zu requested realizations (skipped %zu with corrupted ROF tables)\n",
        //        valid_realizations.size(), realizations_to_run_unique.size(), 
        //        realizations_to_run_unique.size() - valid_realizations.size());
        
        realizations_to_run_unique = valid_realizations;
    }

    // Prepare error output.
    int had_catch = 0;
    string error_m = "Error in realizations ";
    string error_file_name = "error_reals";
    string error_file_content = "#";

    // THREAD-SAFE: Clear global bond tracking BEFORE parallel region starts
    // This ensures all realizations start with clean state and prevents race conditions
    // where one thread clears globals while another thread is reading them
    Utility::clearGlobalBondTracking();

    // Run realizations.
#pragma omp parallel for ordered num_threads(n_threads) shared(had_catch, realizations_to_run_unique, error_m, error_file_name, error_file_content) default(none)
    for (unsigned long r = 0; r < realizations_to_run_unique.size(); ++r) {
        unsigned long realization = realizations_to_run_unique[r];

        // Create continuity models.
        ContinuityModelRealization *realization_model = nullptr;
        ContinuityModelROF *rof_model = nullptr;
        try {
        createContinuityModels(realization, realization_model, rof_model);

        // Initialize data collector.
        // Pass realization WSS so data collectors point to where bonds are actually issued
        master_data_collector->addRealization(
                realization_model->getContinuity_water_sources(),
                realization_model->getDrought_mitigation_policies(),
                utilities,
                realization_model->getContinuity_wss(),  // Realization WSS where bonds are issued
                realization,
                wss_rdm.at(realization));  // Pass RDM factors for deterministic discount rate

            for (int w = 0; w < (int) total_simulation_time; ++w) {
                if (w % 52 == 0) {  // Print every year
                }
//                printf("%d\n", w);
                // DO NOT change the order of the steps. This would mess up
                // important dependencies.
                // Calculate long-term risk-of-failre if current week is first week of the year.
                if (Utils::isFirstWeekOfTheYear(w)) {
                    realization_model->setLongTermROFs(
                            rof_model->calculateLongTermROF(w), w);
                }
                // Calculate short-term risk-of-failure
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

        // Check if any infrastructure was built in this realization
        #pragma omp critical
        {
            const auto& realization_wss = realization_model->getContinuity_wss();
            bool any_infrastructure_built = false;
            
            for (const auto& wss : realization_wss) {
                const auto& water_sources = wss->getWater_sources();
                for (const auto& source : water_sources) {
                    if (source && source->isOnline()) {
                        any_infrastructure_built = true;
                    }
                }
            }
        }

//        printf("Realization %lu took %f seconds.\n", r, omp_get_wtime() - start);

// #pragma omp critical
//         printProgress(
//                 (double) master_data_collector->getRealizations_created() /
//                 (double) realizations_to_run_unique.size());

        } catch (...) {
#pragma omp atomic
            ++had_catch;
#pragma omp critical
            {
            error_m += to_string(realization) + " ";
            error_file_name += "_" + to_string(realization);
            error_file_content += to_string(realization) + ",";
            }
            master_data_collector->removeRealization(realization);
        }

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

