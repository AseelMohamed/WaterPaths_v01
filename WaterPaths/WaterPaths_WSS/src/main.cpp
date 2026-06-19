#include "SystemComponents/WaterSources/Base/WaterSource.h"
#include "Utils/QPSolver/QuadProg++.h"
#include "Utils/Solutions.h"
//#include "Problem/PaperTestProblem.h"
//#include "Problem/Triangle.h"
#include "Utils/Utils.h"
#include "Problem/Caesb.h"

#ifdef  PARALLEL
#include "../Borg/borgms.h"
#include <mpi.h>
#endif

#include <sys/stat.h>
#include <algorithm>
// #include <getopt.h> // Removed for portability
#include <fstream>
#include <omp.h>


using namespace std;
using namespace Constants;
using namespace Solutions;

Caesb *problem_ptr;
//Triangle *problem_ptr;
int failures = 0;

void eval(double *vars, double *objs, double *consts) {
    try {
        for (int i = 0; i < NUM_DEC_VAR; ++i) {
            if (isnan(vars[i])) {
                string error = "Nan in decision variable " + to_string(i);
                throw invalid_argument(error);
            }
        }
        int eval_result = problem_ptr->functionEvaluation(vars, objs, consts);
        failures += eval_result;
        // Only destroy if the evaluation succeeded and master_data_collector was set.
        // On exception paths, functionEvaluation returns 1 and master_data_collector
        // is nullptr (never assigned), so we skip the destroy call to avoid noise.
        if (eval_result == 0) {
            problem_ptr->destroyDataCollector();
        }
    } catch (...) {
        int rank = 0;
#ifdef PARALLEL
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        string sol_out_name = "eval_failure_rank_" + to_string(rank) + ".csv";
        ofstream sol_out(sol_out_name.c_str(), ios::app);
        printf("\nFailure! Decision variable values (rank %d):\n", rank);
        for (int i = 0; i < NUM_DEC_VAR; ++i) {
            printf("  vars[%d] = %f\n", i, vars[i]);
            if (sol_out.is_open()) sol_out << vars[i] << ",";
        }
        if (sol_out.is_open()) sol_out << "\n";
        printf("Objectives set to 1e5 for this evaluation.\n");
        int num_objs = Constants::getNumObjectives();
        for (int i = 0; i < num_objs; ++i) objs[i] = 1e5;
    }
}

int main(int argc, char *argv[]) {
    const int c_num_dec = NUM_DEC_VAR;
    int c_num_obj = Constants::getNumObjectives(); // Use dynamic objective count
    int c_num_constr = 0;
    double c_obj[NUM_OBJECTIVES]; // Use max size for array allocation
    double c_constr[0];

    unsigned long n_realizations = 1000;
    unsigned long n_weeks = 1565;

    string system_io = DEFAULT_DATA_DIR;
    string solution_file = "-1";
    string uncertainty_file = "-1";
    string bootstrap_file = "-1";
    string utilities_rdm_file = "-1";
    string policies_rdm_file = "-1";
    string water_sources_rdm_file = "-1";
    string inflows_evap_directory_suffix = "-1";
    string rof_tables_directory = DEFAULT_ROF_TABLES_DIR;
    unsigned long standard_solution = 0;
    int n_threads = 2;
    int standard_rdm = 0;
    int first_solution = -1;
    int last_solution = -1;
    int import_export_rof_table = 0;
    bool verbose = false;
    bool tabular = false;
    bool plotting = true;
    bool run_optimization = false;
    bool print_objs_row = false;
    unsigned long n_islands = 2;
    unsigned long nfe = 1000;
    unsigned long output_frequency = 200;
    int seed = -1;
    int rdm_no = -1;
    int n_sets = -1;
    int n_bs_samples = -1;
    int target_wss_id = 0;
    vector<vector<int>> realizations_to_run;
    vector<vector<double>> utilities_rdm;
    vector<vector<double>> water_sources_rdm;
    vector<vector<double>> policies_rdm;

    // Portable argument parsing (simple version)
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-s" && i + 1 < argc) solution_file = argv[++i];
        else if (arg == "-u" && i + 1 < argc) uncertainty_file = argv[++i];
        else if (arg == "-T" && i + 1 < argc) n_threads = atoi(argv[++i]);
        else if (arg == "-r" && i + 1 < argc) n_realizations = (unsigned long) atoi(argv[++i]);
        else if (arg == "-t" && i + 1 < argc) n_weeks = (unsigned long) atoi(argv[++i]);
        else if (arg == "-d" && i + 1 < argc) system_io = argv[++i];
        else if (arg == "-f" && i + 1 < argc) first_solution = atoi(argv[++i]);
        else if (arg == "-l" && i + 1 < argc) last_solution = atoi(argv[++i]);
        else if (arg == "-m" && i + 1 < argc) standard_solution = (unsigned long) atoi(argv[++i]);
        else if (arg == "-v" && i + 1 < argc) verbose = static_cast<bool>(atoi(argv[++i]));
        else if (arg == "-c" && i + 1 < argc) tabular = static_cast<bool>(atoi(argv[++i]));
        else if (arg == "-p" && i + 1 < argc) plotting = static_cast<bool>(atoi(argv[++i]));
        else if (arg == "-b") run_optimization = true;
        else if (arg == "-i" && i + 1 < argc) n_islands = (unsigned long) atoi(argv[++i]);
        else if (arg == "-n" && i + 1 < argc) nfe = (unsigned long) atoi(argv[++i]);
        else if (arg == "-o" && i + 1 < argc) output_frequency = (unsigned long) atoi(argv[++i]);
        else if (arg == "-e" && i + 1 < argc) seed = atoi(argv[++i]);
        else if (arg == "-E" && i + 1 < argc) {
            Constants::EXPERIMENT_MODE = atoi(argv[++i]);
            if (Constants::EXPERIMENT_MODE < 1 || Constants::EXPERIMENT_MODE > 7) {
                fprintf(stderr, "Invalid experiment mode. Must be 1, 2, 3, 4, 5, 6, or 7.\n");
                return -1;
            }
            c_num_obj = Constants::getNumObjectives(); // Update objective count
        }
        else if (arg == "-w" && i + 1 < argc) {
            target_wss_id = atoi(argv[++i]);
            Constants::TARGET_WSS_ID = target_wss_id;
        }
        else if (arg == "-V" && i + 1 < argc) {
            int sev_flag = atoi(argv[++i]);
            Constants::INCLUDE_SEVERITY = (sev_flag != 0);
            c_num_obj = Constants::getNumObjectives(); // Update objective count
        }
        else if (arg == "-y" && i + 1 < argc) bootstrap_file = argv[++i];
        else if (arg == "-R" && i + 1 < argc) rdm_no = atoi(argv[++i]);
        else if (arg == "-U" && i + 1 < argc) utilities_rdm_file = argv[++i];
        else if (arg == "-P" && i + 1 < argc) policies_rdm_file = argv[++i];
        else if (arg == "-W" && i + 1 < argc) water_sources_rdm_file = argv[++i];
        else if (arg == "-I" && i + 1 < argc) inflows_evap_directory_suffix = argv[++i];
        else if (arg == "-C" && i + 1 < argc) import_export_rof_table = atoi(argv[++i]);
        else if (arg == "-S" && i + 1 < argc) n_bs_samples = atoi(argv[++i]);
        else if (arg == "-A" && i + 1 < argc) n_sets = atoi(argv[++i]);
        else if (arg == "-O" && i + 1 < argc) rof_tables_directory = argv[++i];
        else if (arg == "-?" || arg == "--help") {
            fprintf(stdout,
                    "%s\n"
                    "\t-?: print this message\n"
                    "\t-s: solutions file (hard coded solutions)\n"
                    "\t-u: uncertain factors file (hard coded values)\n"
                    "\t-T: number of threads (auto)\n"
                    "\t-r: number of realizations (%lu)\n"
                    "\t-t: total simulated time in weeks (%lu)\n"
                    "\t-d: directory for system I/O (%s)\n"
                    "\t-f: first solution number\n"
                    "\t-l: last solution number\n"
                    "\t-m: individual solution number\n"
                    "\t-v: verbose mode (false)\n"
                    "\t-c: tabular output (really big files, false)\n"
                    "\t-p: plotting output (smaller csv files)\n"
                    "\t-b: run optimization with Borg (false)\n"
                    "\t-i: number of islands if using Borg (2)\n"
                    "\t-n: number of function evaluations if using"
                    " Borg (1000)\n"
                    "\t-o: output frequency if using Borg (200)\n"
                    "\t-e: seed number (none)\n"
                    "\t-y: file with bootstrap samples\n"
                    "\t-S: number of bootstrap samples per set for bootstrap analysis.\n"
                    "\t-A: number of sets of bootstrap samples for bootstrap analysis.\n"
                    "\t-R: RDM sample number\n"
                    "\t-U: Utility RDM file\n"
                    "\t-P: Policies RDM file\n"
                    "\t-W: Water sources RDM file\n"
                    "\t-I: Inflows and evaporation folder suffix to"
                    " be added to \"inflows\" and "
                    "\"evaporation\" (e.g., _high for "
                    "inflows_high)\n"
                    "\t-O: Directory containing the pre-computed "
                    "ROF table binaries\n"
                    "\t-C: Import/export rof tables (1: export, 0:"
                    " do nothing (standard), -1: import)\n"
                    "\t-B: Export objectives for all utilities on a single line\n"
                    "\t-E: Experiment mode (1-7):\n"
                    "\t    1: 4 objs, reliability=MIN (worst case)\n"
                    "\t    2: 4 objs, reliability=AVERAGE\n"
                    "\t    3: 5 objs, reliability=MIN, affordability=MAX [DEFAULT]\n"
                    "\t    4: 5 objs, reliability=AVERAGE, affordability=AVERAGE\n"
                    "\t    5: 7 objs, per-WSS reliability and affordability (no aggregation)\n"
                    "\t    6: 5 objs, single-WSS reliability and affordability (use -w to set target WSS)\n"
                    "\t    7: 5 objs, penalty-based (diff_rel + diff_afford summed over WSS)\n"
                    "\t       diff_rel   = 100 * sum_wss[max(0.97 - rel_wss,    0)^2]\n"
                    "\t       diff_afford= 100 * sum_wss[max(afford_wss - 0.03, 0)^2]\n"
                    "\t-w: Target WSS ID for experiment 6 (default: 0)\n"
                    "\t-V: Include severity objective (1=yes [default], 0=no)\n",
                    argv[0], n_realizations, n_weeks, system_io.c_str());
            return -1;
        }
        // Unknown option handling
        else if (arg[0] == '-') {
            fprintf(stderr, "Unknown option (%s)\n", arg.c_str());
            return -1;
        }
    }

    // Print experiment mode configuration (only on master process in serial section)
    #ifndef PARALLEL
    printf("========================================\n");
    printf("Experiment Configuration:\n");
    printf("  Experiment Mode: %d\n", Constants::EXPERIMENT_MODE);
    printf("  Number of Objectives: %d\n", Constants::getNumObjectives());
    if (Constants::includePerWSSObjectives()) {
        printf("  Mode: per-WSS (no aggregation across WSS)\n");
    } else if (Constants::includeSingleWSSObjectives()) {
        printf("  Mode: single-WSS (target WSS ID: %d)\n", Constants::TARGET_WSS_ID);
    } else if (Constants::includePenaltyObjectives()) {
        printf("  Mode: penalty-based\n");
        printf("  diff_rel    = 100 * sum_wss[max(%.2f - rel_wss,    0)^2]\n", Constants::PENALTY_RELIABILITY_THRESHOLD);
        printf("  diff_afford = 100 * sum_wss[max(afford_wss - %.2f, 0)^2]\n", Constants::PENALTY_AFFORDABILITY_THRESHOLD);
    } else {
    printf("  Reliability Aggregation: %s\n", 
           Constants::getReliabilityAggregationMethod() == Constants::AVERAGE ? "AVERAGE" : "MIN");
    if (Constants::includeAffordabilityObjective()) {
        printf("  Affordability Aggregation: %s\n", 
               Constants::getAffordabilityAggregationMethod() == Constants::AVERAGE ? "AVERAGE" : "MAX");
    } else {
        printf("  Affordability: NOT INCLUDED\n");
    }
    if (Constants::includeSeverityObjective()) {
        printf("  Severity: INCLUDED (aggregation follows reliability)\n");
    } else {
        printf("  Severity: NOT INCLUDED\n");
    }
    }
    printf("========================================\n\n");
    #endif
    
    Caesb problem(n_weeks, import_export_rof_table);
//    Triangle problem(n_weeks, import_export_rof_table);
    if (seed > -1) {
        WaterSource::setSeed(seed);
        MasterDataCollector::setSeed(seed);
    }

    // Set basic realization parameters.
    problem.setN_weeks(n_weeks);
    problem.setIODirectory(system_io);
    problem.setN_threads((unsigned long) n_threads);

    // Load bootstrap samples if necessary.
    if (strlen(bootstrap_file.c_str()) > 2) {
        auto bootstrap_samples_double = Utils::parse2DCsvFile(system_io + bootstrap_file);
        for (auto &v : bootstrap_samples_double) {
            realizations_to_run.push_back(vector<int>(v.begin(), v.end()));
            if (*std::max_element(v.begin(), v.end()) >= n_realizations)
                throw invalid_argument(
                        "Number of realizations must be higher than the ID of all realizations in the bootstrap samples file.");

        }
        problem.setPrint_output_files(false);
    }

    if (!run_optimization && strlen(system_io.c_str()) == 0) {
        throw invalid_argument(
                "You must specify an input output directory.");
    }

    // Set up input/output suffix, if necessary.
    if (strlen(inflows_evap_directory_suffix.c_str()) > 2) {
        problem.setEvap_inflows_suffix(inflows_evap_directory_suffix);
    }
    // Read RDM file, if any
    if (strlen(utilities_rdm_file.c_str()) > 2) {
        printf("reading RDM file\n");
        if (rdm_no != NON_INITIALIZED) {
            auto utilities_rdm_row = Utils::parse2DCsvFile(system_io + utilities_rdm_file)[rdm_no];
            auto policies_rdm_row = Utils::parse2DCsvFile(system_io + policies_rdm_file)[rdm_no];
            auto water_sources_rdm_row = Utils::parse2DCsvFile(system_io + water_sources_rdm_file)[rdm_no];

            utilities_rdm = std::vector<vector<double>>(n_realizations, utilities_rdm_row);
            policies_rdm = std::vector<vector<double>>(n_realizations, policies_rdm_row);
            water_sources_rdm = std::vector<vector<double>>(n_realizations, water_sources_rdm_row);

            problem.setRDMReevaluation(rdm_no, utilities_rdm,
                                       water_sources_rdm, policies_rdm);

            if (strlen(inflows_evap_directory_suffix.c_str()) > 2)
                problem.setFname_sufix("_RDM" + std::to_string(rdm_no) +
                                       "_infevap" + inflows_evap_directory_suffix);
            else
                problem.setFname_sufix("_RDM" + std::to_string(rdm_no));
        } else {
            utilities_rdm = Utils::parse2DCsvFile(system_io + utilities_rdm_file);
            water_sources_rdm = Utils::parse2DCsvFile(system_io + water_sources_rdm_file);
            policies_rdm = Utils::parse2DCsvFile(system_io + policies_rdm_file);
            if (n_realizations > utilities_rdm.size()) {
                throw length_error("If no rdm number is passed, the number of realizations needs to be smaller "
                                   "or equal to the number of rows in the rdm files.");
            }
            problem.setRDMReevaluation(rdm_no, utilities_rdm,
                                       water_sources_rdm, policies_rdm);
        }
    }
    problem_ptr = &problem;

    // Set realizations to be run -- otherwise, n_realizations realizations will be run.
    if (!realizations_to_run.empty() && (n_sets <= 0 || n_bs_samples <= 0)) {
        auto realizations_to_run_ul = vector<unsigned long>(realizations_to_run[0].begin(),
                                                            realizations_to_run[0].end());
        problem.setRealizationsToRun(realizations_to_run_ul);
        problem.setN_realizations(*max_element(realizations_to_run_ul.begin(), realizations_to_run_ul.end()) + 1);
    } else {
        problem.setN_realizations(n_realizations);
    }
    problem.setImport_export_rof_tables(import_export_rof_table,
                                        system_io + rof_tables_directory);
    problem.readInputData();

    // If Borg is not called, run in simulation mode
    if (!run_optimization) {
        vector<vector<double>> solutions;
        if (strlen(solution_file.c_str()) > 2) {
            solutions = Utils::parse2DCsvFile(system_io + solution_file);
            if (standard_solution >= solutions.size())
                throw invalid_argument("Number of solutions in file <= solution ID.\n");
        } else {
            throw invalid_argument("You must specify a solutions file.\n");
        }

        vector<int> sol_range;
        // Check for basic input errors.
        if ((first_solution == -1 && last_solution != -1) ||
            (first_solution != -1 && last_solution == -1))
            throw invalid_argument("If you set a first or last solution, you "
                                   "must set the other as well.");

        // Run model
        if (first_solution == -1) {
            printf("\n\n\nRunning solution %lu%s%s\n", 
                   standard_solution,
                   (rdm_no != NON_INITIALIZED ? " RDM " : ""),
                   (rdm_no != NON_INITIALIZED ? to_string(rdm_no).c_str() : ""));
            problem.setSol_number(standard_solution);
            problem_ptr->functionEvaluation(solutions[standard_solution].data(), c_obj, c_constr);
            if (problem_ptr->getMaster_data_collector() != nullptr) {
                problem_ptr->getMaster_data_collector()->setOutputDirectory(system_io);
                problem_ptr->getMaster_data_collector()->printWeeklyReliabilityByWSSCsv(
                    "annualWSS_s" + std::to_string(standard_solution) +
                    (rdm_no == NON_INITIALIZED ? "" : "_RDM" + std::to_string(rdm_no)));
                problem_ptr->getMaster_data_collector()->printAnnualReliabilityBySourceCsv(
                    "annualSource_s" + std::to_string(standard_solution) +
                    (rdm_no == NON_INITIALIZED ? "" : "_RDM" + std::to_string(rdm_no)));
            }

            // Export pathways and objectives, otherwise, if required, run bootstrap sub-sampling.
            if (n_sets > 0 && n_bs_samples > 0) {
                problem_ptr->runBootstrapRealizationThinning(
                        (int) standard_solution, n_sets, n_bs_samples, n_threads, realizations_to_run);
            } else if (import_export_rof_table != EXPORT_ROF_TABLES) {
                if (plotting)
                    problem.printTimeSeriesAndPathways(plotting);
                auto objectives = problem_ptr->calculateAndPrintObjectives(!print_objs_row);
                //            trianglePtr->getMaster_data_collector()->printNETCDFUtilities("netcdf_output");
            }

            problem_ptr->destroyDataCollector();
        } else {
//            double time_0 = omp_get_wtime();
            ofstream objs_file;
            string file_name = system_io + "output" + BAR + "Objectives" +
                    (rdm_no == NON_INITIALIZED ? "" : "_RDM" + to_string(rdm_no)) +
                               "_sols" + to_string(first_solution) + "_to_" + to_string(last_solution) + ".csv";
            objs_file.open(file_name);
            printf("Objectives file will be printed at %s.\n", file_name.c_str());
            for (int s = first_solution; s < last_solution; ++s) {
                printf("\n\n\nRunning solution %d%s%s\n", 
                       s,
                       (rdm_no != NON_INITIALIZED ? " RDM " : ""),
                       (rdm_no != NON_INITIALIZED ? to_string(rdm_no).c_str() : ""));
                problem.setSol_number((unsigned long) s);
                problem_ptr->functionEvaluation(solutions[s].data(), c_obj, c_constr);
                if (problem_ptr->getMaster_data_collector() != nullptr) {
                    problem_ptr->getMaster_data_collector()->setOutputDirectory(system_io);
                    problem_ptr->getMaster_data_collector()->printWeeklyReliabilityByWSSCsv(
                            "annualWSS_s" + std::to_string(s) +
                            (rdm_no == NON_INITIALIZED ? "" : "_RDM" + std::to_string(rdm_no)));
                    problem_ptr->getMaster_data_collector()->printAnnualReliabilityBySourceCsv(
                            "annualSource_s" + std::to_string(s) +
                            (rdm_no == NON_INITIALIZED ? "" : "_RDM" + std::to_string(rdm_no)));
                }
                vector<double> objectives = problem_ptr->calculateAndPrintObjectives(false);
                // problem.printTimeSeriesAndPathways(plotting);
                problem.printTimeSeriesAndPathways(plotting);
                problem_ptr->destroyDataCollector();
                string line;
                for (double &o : objectives) {
                    line += to_string(o) + ",";
                }
                line.pop_back();
                objs_file << line << endl;
            }
            objs_file.close();
//            printf("Time to simulate %d solutions: %f s", last_solution - first_solution, omp_get_wtime() - time_0);
        }

        return 0;
    } else {
#ifdef  PARALLEL

        printf("Running Borg with:\n"
            "n_dec_vars: %d\n"
            "n_objectives: %d\n"
            "nfe: %lu\n"
            "output freq.: %lu\n"
            "n_weeks: %lu\n"
            "n_realizations: %lu\n\n",
            c_num_dec, c_num_obj, nfe, output_frequency, n_weeks, n_realizations);
         
        // for debugging borg, creating file to print each ranks DVs which isdone in Eval function   
        
        BORG_Algorithm_ms_startup(&argc, &argv);
        
        // Print experiment configuration from master rank only (after MPI is initialized)
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            printf("========================================\n");
            printf("Experiment Configuration:\n");
            printf("  Experiment Mode: %d\n", Constants::EXPERIMENT_MODE);
            printf("  Number of Objectives: %d\n", Constants::getNumObjectives());
            if (Constants::includePerWSSObjectives()) {
                printf("  Mode: per-WSS (no aggregation across WSS)\n");
            } else if (Constants::includeSingleWSSObjectives()) {
                printf("  Mode: single-WSS (target WSS ID: %d)\n", Constants::TARGET_WSS_ID);
            } else if (Constants::includePenaltyObjectives()) {
                printf("  Mode: penalty-based\n");
                printf("  diff_rel    = 100 * sum_wss[max(%.2f - rel_wss,    0)^2]\n", Constants::PENALTY_RELIABILITY_THRESHOLD);
                printf("  diff_afford = 100 * sum_wss[max(afford_wss - %.2f, 0)^2]\n", Constants::PENALTY_AFFORDABILITY_THRESHOLD);
            } else {
            printf("  Reliability Aggregation: %s\n", 
                   Constants::getReliabilityAggregationMethod() == Constants::AVERAGE ? "AVERAGE" : "MIN");
            if (Constants::includeAffordabilityObjective()) {
                printf("  Affordability Aggregation: %s\n", 
                       Constants::getAffordabilityAggregationMethod() == Constants::AVERAGE ? "AVERAGE" : "MAX");
            } else {
                printf("  Affordability: NOT INCLUDED\n");
            }
            if (Constants::includeSeverityObjective()) {
                printf("  Severity: INCLUDED (aggregation follows reliability)\n");
            } else {
                printf("  Severity: NOT INCLUDED\n");
            }
            }
            printf("========================================\n");
        }
        // BORG_Algorithm_ms_islands((int) n_islands);
        // BORG_Algorithm_ms_initialization(INITIALIZATION_LATIN_GLOBAL);
        BORG_Algorithm_ms_max_evaluations((int) nfe);
        BORG_Algorithm_output_frequency((int) output_frequency);

        // Define the problem.
        BORG_Problem problem = BORG_Problem_create(c_num_dec, c_num_obj,
                                                   c_num_constr,
                                                   eval);
        // Set all the parameter bounds and epsilons
        printf("setting up problem\n");
        problem_ptr->setProblemDefinition(problem);

        if (seed > -1) {
            srand(seed);
            WaterSource::setSeed(seed);
	    BORG_Random_seed(seed);
        }
        string output_directory = system_io + DEFAULT_OUTPUT_DIR;
        string output_file_name = output_directory + "NC_output_MS_E" + to_string(Constants::EXPERIMENT_MODE) + "_S" + to_string(seed) + "_N" + to_string(nfe) + ".set";
        string runtime_file = output_directory + "NC_runtime_MS_E" + to_string(Constants::EXPERIMENT_MODE) + "_S" + to_string(seed) + "_N" + to_string(nfe) + ".runtime";
        
        printf("Reference set will be in %s.\n", output_file_name.c_str());
        printf("Runtime files will be in %s.\n", runtime_file.c_str());

        BORG_Algorithm_output_runtime(const_cast<char*>(runtime_file.c_str()));

        // rank already declared above, just reuse it
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        //string rank_out_file = "diagnostic_output/DVs_rank_" + to_string(rank) + ".csv";
        //sol_out.open(rank_out_file.c_str());

        //int rank; // different seed on each processor
        //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        //BORG_Random_seed(37*seed*(rank+1));
        
        BORG_Archive result = BORG_Algorithm_ms_run(problem); // this actually runs the optimization
        //BORG_Archive result = BORG_Algorithm_run(problem, nfe);
        // If this is the master node, print out the final archive

        if (result != nullptr) {
            FILE* outputFile = fopen(output_file_name.c_str(), "w");
            printf("master node print, should see only once\n");
            if (!outputFile) {
                BORG_Debug("Unable to open final output file\n");
            }
            BORG_Archive_print(result, outputFile);
            BORG_Archive_destroy(result);
            //sol_out.close();
            fclose(outputFile);
        }

        printf("Number of failed function evaluations: %d.\n", failures);

        BORG_Algorithm_ms_shutdown();
        BORG_Problem_destroy(problem);
#else
        throw invalid_argument("This version of WaterPaths was not compiled with Borg.");
#endif

        return 0;
    }
}
