//;;
// Created by bernardoct on 8/26/17.
//

#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <numeric>
#include <random>
#include <algorithm>
#include <cmath>
#include <map>
#include "../DroughtMitigationInstruments/TransfersBilateral.h"

#ifdef NETCDF
#include <netcdf> 
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
// Return this in event of a problem.
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
static const int NC_ERR = 2;
#define DEFLATE_LEVEL_9 9
#endif

#include "MasterDataCollector.h"
#include "../Utils/ObjectivesCalculator.h"
#include "../Utils/Utils.h"
#include "../DroughtMitigationInstruments/Transfers.h"
#include "TransfersDataCollector.h"
#include "../SystemComponents/WaterSources/Quarry.h"
#include "../SystemComponents/WaterSources/WaterReuse.h"
#include "../SystemComponents/WaterSources/AllocatedReservoir.h"
#include "ReservoirDataCollector.h"
#include "IntakeDataCollector.h"
#include "QuaryDataCollector.h"
#include "WaterReuseDataCollector.h"
#include "AllocatedReservoirDataCollector.h"
#include "EmptyDataCollector.h"
#include "TransfersBilateralDataCollector.h"
#include "EmergencyTransferParanoaDataCollector.h"
#include "../DroughtMitigationInstruments/EmergencyTransferParanoa.h"

using namespace Constants;

int MasterDataCollector::seed = NON_INITIALIZED;

MasterDataCollector::MasterDataCollector(vector<unsigned long> &realizations_to_run)
        : n_realizations(*max_element(realizations_to_run.begin(), realizations_to_run.end()) + 1),
        realizations_ran(realizations_to_run) {}

MasterDataCollector::~MasterDataCollector() {
    for (vector<DataCollector *> dcs : water_source_collectors)
        for (DataCollector *dc : dcs)
            delete dc;

    for (vector<DataCollector *> dcs : drought_mitigation_policy_collectors)
        for (DataCollector *dc : dcs)
            delete dc;

    for (vector<UtilitiesDataCollector *> dcs : utility_collectors)
        for (UtilitiesDataCollector *dc : dcs)
            delete dc;
    
    for (vector<WSSDataCollector *> dcs : wss_collectors)
        for (WSSDataCollector *dc : dcs)
            delete dc;
}

/**
* @todo NetCDF file being saved with hardly any compression. Someone who knows more
* about NetCDF may be able to improve compression, which would be great given the 
* amount of output that can be generated in one run.
*/
int MasterDataCollector::printNETCDFUtilities(string base_file_name) {
#ifdef NETCDF
    unsigned long n_weeks = utility_collectors[0][realizations_ran[0]]->getCombined_storage().size();
    unsigned long n_realizations = utility_collectors[0].size();
    unsigned long n_utilities = utility_collectors.size();
    unsigned long n_vars = 6;

    try {
	string file_name = io_directory + base_file_name + ".nc";
	printf("Printing NetCDF output in %s.\n", file_name.c_str());
	int ncid, dimid, csid, strofid, retval;
	vector<int> group_ids(n_realizations);
	vector<int> utilities_group_ids(n_utilities * n_realizations);
	vector<int> vars_ids(n_utilities * n_realizations * n_vars);

        /* Create the file. The NC_NETCDF4 flag tells netCDF to
         * create a netCDF-4/HDF5 file.
	 * */
        if ((retval = nc_create(file_name.c_str(), NC_NETCDF4|NC_CLOBBER, &ncid)))
		               ERR(retval);

        nc_def_dim(ncid, "weeks", n_weeks, &dimid);
        for (int r = 0; r < (int) n_realizations; ++r) {
	    if ((retval = nc_def_grp(ncid, (string("realization_") + to_string(r)).c_str(), &group_ids[r])))
		             ERR (retval);
	    for (int u = 0; u < (int) n_utilities; ++u) {
	        if ((retval = nc_def_grp(group_ids[r], utility_collectors[u][r]->name, &utilities_group_ids[n_utilities * r + u])))
		                 ERR (retval);	    

                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "storage", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u])))
                    ERR(retval);
                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "st_rof", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u + 1])))
                    ERR(retval);
                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "lt_rof", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u + 2])))
                    ERR(retval);
                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "res_demand", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u + 3])))
                    ERR(retval);
                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "cf_size", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u+ 4])))
                    ERR(retval);
                if ((retval = nc_def_var(utilities_group_ids[n_utilities * r + u], "infr_pmts", 
						NC_FLOAT, 1, &dimid, &vars_ids[n_utilities * r * n_vars + u + 5])))
                    ERR(retval);

                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);
                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u + 1], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);
                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u + 2], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);
                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u + 3], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);
                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u + 4], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);
                if ((retval = nc_def_var_deflate(utilities_group_ids[n_utilities * r + u], 
						vars_ids[n_utilities * r * n_vars + u + 5], 1, 1, DEFLATE_LEVEL_9)))
                    ERR(retval);

                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u],
						utility_collectors[u][r]->getCombined_storage().data())))
	            ERR(retval);
                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u + 1],
						utility_collectors[u][r]->getSt_rof().data())))
	            ERR(retval);
                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u + 2],
						utility_collectors[u][r]->getLt_rof().data())))
	            ERR(retval);
                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u + 3],
						utility_collectors[u][r]->getRestricted_demand().data())))
	            ERR(retval);
                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u + 4],
						utility_collectors[u][r]->getContingency_fund_size().data())))
	            ERR(retval);
                if ((retval = nc_put_var_double(utilities_group_ids[n_utilities * r + u], vars_ids[n_utilities * r * n_vars + u + 5],
						utility_collectors[u][r]->getDebt_service_payments().data())))
	            ERR(retval);
	    }
	}
	return 0;
    } catch(NcException& e){
        e.what();
        return NC_ERR;
    }
#else
    printf("This version of WaterPaths was not compiled with NetCDF. NetCDF result files will not be printed.\n");
    return 1;
#endif
}

void MasterDataCollector::printPoliciesOutputCompact(
        int week_i, int week_f, string file_name) {
    // Print ACTUAL WSS restriction multipliers (what was applied to each WSS)
    // NOT the policy decisions (which are misleading since policies apply to all WSS)
    if (!wss_collectors.empty()) {
#pragma omp parallel for
        for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
            auto r = realizations_ran[rr];
            std::ofstream out_stream;
            out_stream.open(output_directory + file_name + "_r"
                            + std::to_string(r) + ".csv");

            // Build header: WSS restriction multipliers + any other policy data (transfers, etc.)
            string line = "";
            // First, add WSS restriction multipliers
            for (const auto& wss_vec : wss_collectors) {
                if (wss_vec[r] != nullptr) {
                    line += std::to_string(wss_vec[r]->id) + "rest_m,";
                }
            }
            // Then add non-restriction policies (transfers, insurance, etc.)
            for (vector<DataCollector *> p : drought_mitigation_policy_collectors) {
                if (p[r]->type != RESTRICTIONS) {
                    line += p[r]->printCompactStringHeader();
                }
            }
            if (!line.empty() && line.back() == ',') line.pop_back();
            out_stream << line << endl;

            // Print data for each week
            for (int w = week_i; w < week_f; ++w) {
                line = "";
                // First, print WSS restriction multipliers (actual applied values)
                for (const auto& wss_vec : wss_collectors) {
                    if (wss_vec[r] != nullptr) {
                        const auto& demand_mults = wss_vec[r]->getDemand_multiplier();
                        if (w < (int)demand_mults.size()) {
                            line += std::to_string(demand_mults[w]) + ",";
                        } else {
                            line += "1.0,";  // Default if data not available
                        }
                    }
                }
                // Then print other policy data (transfers, etc.)
                for (vector<DataCollector *> p : drought_mitigation_policy_collectors) {
                    if (p[r]->type != RESTRICTIONS) {
                        line += p[r]->printCompactString(w);
                    }
                }
                if (!line.empty() && line.back() == ',') line.pop_back();
                out_stream << line << endl;
            }

            out_stream.close();
        }
    } else if (!drought_mitigation_policy_collectors.empty()) {
        // Fallback to old behavior if no WSS collectors (shouldn't happen in WSS model)
#pragma omp parallel for
        for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
            auto r = realizations_ran[rr];
            std::ofstream out_stream;
            out_stream.open(output_directory + file_name + "_r"
                            + std::to_string(r) + ".csv");

            string line;
            for (vector<DataCollector *> p : drought_mitigation_policy_collectors)
                line += p[r]->printCompactStringHeader();
            line.pop_back();
            out_stream << line << endl;

            for (int w = week_i; w < week_f; ++w) {
                line = "";
                for (vector<DataCollector *> p : drought_mitigation_policy_collectors)
                    line += p[r]->printCompactString(w);
                line.pop_back();
                out_stream << line << endl;
            }

            out_stream.close();
        }
    }
}


void MasterDataCollector::printPoliciesOutputTabular(
        int week_i, int week_f, string file_name) {
    if (!drought_mitigation_policy_collectors.empty()) {
#pragma omp parallel for
        for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
            auto r = realizations_ran[rr];
            std::ofstream out_stream;
            out_stream.open(output_directory + file_name + "_r"
                            + std::to_string(r) + ".tab");

            out_stream << "    ";
            for (vector<DataCollector *> p : drought_mitigation_policy_collectors)
                out_stream << p[r]->printTabularStringHeaderLine1();
            out_stream << endl;

            out_stream << "Week";
            for (vector<DataCollector *> p : drought_mitigation_policy_collectors)
                out_stream << p[r]->printTabularStringHeaderLine2();
            out_stream << endl;

            for (int w = week_i; w < week_f; ++w) {
                out_stream << setw(4) << w;
                for (vector<DataCollector *> p : drought_mitigation_policy_collectors)
                    out_stream << p[r]->printTabularString(w);
                out_stream << endl;
            }

            out_stream.close();
        }
    }
}

void MasterDataCollector::printUtilitiesOutputCompact(
        int week_i, int week_f, string file_name) {
#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        std::ofstream out_stream;
        out_stream.open(output_directory + file_name + "_r"
                        + std::to_string(r) + ".csv");

        string line;
        for (vector<UtilitiesDataCollector *> &p : utility_collectors)
            line += p[r]->printCompactStringHeader();
        line.pop_back();
        out_stream << line << endl;

        for (int w = week_i; w < week_f; ++w) {
            line = "";
            for (vector<UtilitiesDataCollector *> &p : utility_collectors)
                line += p[r]->printCompactString(w);
            line.pop_back();
            out_stream << line << endl;
        }

        out_stream.close();
    }
}

void MasterDataCollector::printWCCComponentsOutput(
        int week_i, int week_f, string file_name) {
#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        std::ofstream out_stream;
        out_stream.open(output_directory + file_name + "_r"
                        + std::to_string(r) + ".csv");

        out_stream << "utility_id,utility_name,realization,year,row_type,"
                      "drought_cost_sum,gross_rev_sum,discount_rate,annual_cost"
                   << endl;

        for (vector<UtilitiesDataCollector *> &p : utility_collectors) {
            auto *u = p[r];
            if (u == nullptr) {
                continue;
            }

            const auto &drought = u->getDrought_mitigation_cost();
            const auto &gross = u->getGross_revenues();
            int max_week = std::min(week_f, (int)std::min(drought.size(), gross.size()));

            double discount_rate = u->getInfraDiscountRate();
            double drought_sum = 0.0;
            double gross_sum = 1e-6;
            int year = 0;
            double max_annual_cost = 0.0;

            for (int w = week_i; w < max_week; ++w) {
                drought_sum += drought[w];
                gross_sum += gross[w];

                if (Utils::isFirstWeekOfTheYear(w + 1)) {
                    double annual_cost = std::max(drought_sum, 0.0) /
                            (gross_sum * (1. + std::pow(1. + discount_rate, year)));

                    out_stream << u->id << ","
                               << u->name << ","
                               << r << ","
                               << year << ","
                               << "year" << ","
                               << drought_sum << ","
                               << gross_sum << ","
                               << discount_rate << ","
                               << annual_cost
                               << endl;

                    if (annual_cost > max_annual_cost) {
                        max_annual_cost = annual_cost;
                    }

                    drought_sum = 0.0;
                    gross_sum = 1e-6;
                    year++;
                }
            }

            out_stream << u->id << ","
                       << u->name << ","
                       << r << ","
                       << -1 << ","
                       << "max" << ","
                       << "" << ","
                       << "" << ","
                       << discount_rate << ","
                       << max_annual_cost
                       << endl;
        }

        out_stream.close();
    }
}


void MasterDataCollector::printUtilitesOutputTabular(
        int week_i, int week_f, string file_name) {
#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        std::ofstream out_stream;
        out_stream.open(output_directory + file_name + "_r"
                        + std::to_string(r) + ".tab");

        stringstream names;
        names << "    ";
        for (vector<UtilitiesDataCollector *> &p : utility_collectors)
            names << setw(p[realizations_ran[0]]->table_width) << p[r]->name;

        out_stream << names.str();
        out_stream << endl;

        out_stream << "    ";
        for (vector<UtilitiesDataCollector *> &p : utility_collectors)
            out_stream << p[r]->printTabularStringHeaderLine1();
        out_stream << endl;

        out_stream << "Week";
        for (vector<UtilitiesDataCollector *> &p : utility_collectors)
            out_stream << p[r]->printTabularStringHeaderLine2();
        out_stream << endl;

        for (int w = week_i; w < week_f; ++w) {
            out_stream << setw(4) << w;
            for (vector<UtilitiesDataCollector *> &p : utility_collectors)
                out_stream << p[r]->printTabularString(w);
            out_stream << endl;
        }

        out_stream.close();
    }
}

void MasterDataCollector::printWaterSourcesOutputCompact(
        int week_i, int week_f, string file_name) {
#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        try {
            std::ofstream out_stream;
            out_stream.open(output_directory + file_name + "_r"
                            + std::to_string(r) + ".csv");

            string line;
            for (vector<DataCollector *> p : water_source_collectors)
                line += p[r]->printCompactStringHeader();
            line.pop_back();
            out_stream << line << endl;

            for (int w = week_i; w < week_f; ++w) {
                line = "";
                for (vector<DataCollector *> p : water_source_collectors)
                    line += p[r]->printCompactString(w);
                line.pop_back();
                out_stream << line << endl;
            }

            out_stream.close();
        } catch (...) {
            printf("Warning: water sources data for realization %lu not saved due to error.\n", r);
        }
    }
}

void MasterDataCollector::printWaterSourcesOutputTabular(
        int week_i, int week_f, string file_name) {
#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        std::ofstream out_stream;
        out_stream.open(output_directory + file_name + "_r"
                        + std::to_string(r) + ".tab");

        stringstream names;
        names << "    ";
        for (vector<DataCollector *> p : water_source_collectors)
            names << setw(p[realizations_ran[0]]->table_width) << p[r]->name;

        out_stream << names.str();
        out_stream << endl;

        out_stream << "    ";
        for (vector<DataCollector *> p : water_source_collectors)
            out_stream << p[r]->printTabularStringHeaderLine1();
        out_stream << endl;

        out_stream << "Week";
        for (vector<DataCollector *> p : water_source_collectors)
            out_stream << p[r]->printTabularStringHeaderLine2();
        out_stream << endl;

        for (int w = week_i; w < week_f; ++w) {
            out_stream << setw(4) << w;
            for (vector<DataCollector *> p : water_source_collectors)
                out_stream << p[r]->printTabularString(w);
            out_stream << endl;
        }

        out_stream.close();
    }
}

void MasterDataCollector::printWSSOutputCompact(
        int week_i, int week_f, string file_name) {
    // Only print if WSS collectors exist
    if (wss_collectors.empty()) {
        printf("No WSS collectors available - skipping WSS output files.\\n");
        return;
    }

#pragma omp parallel for
    for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
        auto r = realizations_ran[rr];
        try {
            std::ofstream out_stream;
            out_stream.open(output_directory + file_name + "_r"
                            + std::to_string(r) + ".csv");

            // Print header
            string line;
            for (vector<WSSDataCollector *> &wss_vec : wss_collectors)
                line += wss_vec[r]->printCompactStringHeader();
            line.pop_back();  // Remove trailing comma
            out_stream << line << endl;

            // Print data for each week
            for (int w = week_i; w < week_f; ++w) {
                line = "";
                for (vector<WSSDataCollector *> &wss_vec : wss_collectors)
                    line += wss_vec[r]->printCompactString(w);
                line.pop_back();  // Remove trailing comma
                out_stream << line << endl;
            }

            out_stream.close();
        } catch (...) {
            char error_msg[256];
            sprintf(error_msg,
                    "Error printing WSS output for realization %lu\\n", r);
            throw runtime_error(error_msg);
        }
    }
}

void MasterDataCollector::printUtilityObjectivesToRowOutStream(vector<UtilitiesDataCollector *> &u,
        std::ofstream &outStream, vector<double> &objectives) {
    try {
    // Create vector with restriction policies pertaining only to the
    // utility whose objectives are being calculated.
    vector<RestrictionsDataCollector *> utility_restrictions(
            *max_element(realizations_ran.begin(), realizations_ran.end()) + 1,
            nullptr
    );
    isolateRestrictionDataCollectors(u, utility_restrictions);

    // Reliability - Use WSS-level calculation (utility = minimum reliability of its WSS)
    double reliability = 1.0;
    if (!wss_collectors.empty()) {
        // Filter WSS collectors to only include those belonging to this utility
        vector<vector<WSSDataCollector *>> utility_wss_collectors;
        isolateWSSDataCollectors(u, utility_wss_collectors);
        
        if (!utility_wss_collectors.empty()) {
            // Use configurable aggregation method based on experiment mode
            reliability = ObjectivesCalculator::calculateReliabilityObjective_WSS_Configurable(
                utility_wss_collectors, 
                realizations_ran,
                Constants::getReliabilityAggregationMethod());
        } else {
            // No WSS found for this utility, fallback to utility-level calculation
            reliability = ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran);
        }
    } else {
        // Fallback to utility-level if WSS collectors not available
        reliability = ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran);
    }
    
    /// Restriction Frequency — use WSS demand multipliers (covers all WSS, consistent with policies CSV)
    double restriction_freq = NONE;
    if (!wss_collectors.empty()) {
        vector<vector<WSSDataCollector *>> utility_wss_collectors_rf;
        isolateWSSDataCollectors(u, utility_wss_collectors_rf);
        if (!utility_wss_collectors_rf.empty())
            restriction_freq = ObjectivesCalculator::
                calculateRestrictionFrequencyObjective_WSS(utility_wss_collectors_rf, realizations_ran);
    }
    if (restriction_freq == NONE)
        restriction_freq = ObjectivesCalculator::
            calculateRestrictionFrequencyObjective(utility_restrictions, realizations_ran);
    
    /// Infrastructure NPC and Worse Case Costs - Use WSS-level calculations
    double inf_npc = 0.0;
    double worse_cost = 0.0;
    if (!wss_collectors.empty()) {
        // Use WSS-level calculations (pass utility data for safe discount rate access)
        inf_npc = ObjectivesCalculator::
        calculateNetPresentCostInfrastructureObjective_WSS(wss_collectors, u, realizations_ran);
        worse_cost = ObjectivesCalculator::
        calculateWorseCaseCostsObjective(u, realizations_ran);
    } else {
        // Fallback to utility-level if WSS collectors not available
        inf_npc = ObjectivesCalculator::
        calculateNetPresentCostInfrastructureObjective(u, realizations_ran);
        worse_cost = ObjectivesCalculator::
        calculateWorseCaseCostsObjective(u, realizations_ran);
    }
    
    /// Affordability Index - Only calculate if included in experiment
    double affordability_index = 0.0;
    if (Constants::includeAffordabilityObjective()) {
        if (!wss_collectors.empty()) {
            // Filter WSS collectors to only include those belonging to this utility
            vector<vector<WSSDataCollector *>> utility_wss_collectors_afford;
            isolateWSSDataCollectors(u, utility_wss_collectors_afford);
            
            // Use configurable aggregation method based on experiment mode
            affordability_index = ObjectivesCalculator::
            calculateAffordabilityIndexObjective_WSS_Configurable(
                utility_wss_collectors_afford, 
                realizations_ran,
                Constants::getAffordabilityAggregationMethod());
        } else {
            // No WSS data available - use a default high value
            affordability_index = 1.0;
        }
    }

    /// Failure Severity - Use same aggregation method as reliability
    double failure_severity = 0.0;
    if (Constants::includeSeverityObjective()) {
        if (!wss_collectors.empty()) {
            vector<vector<WSSDataCollector *>> utility_wss_collectors_sev;
            isolateWSSDataCollectors(u, utility_wss_collectors_sev);
            failure_severity = ObjectivesCalculator::
            calculateFailureSeverityObjective_WSS_Configurable(
                utility_wss_collectors_sev,
                realizations_ran,
                Constants::getReliabilityAggregationMethod());
        }
    }

    if (Constants::includePerWSSObjectives()) {
        // Mode 5: keep each WSS's reliability and affordability as separate objectives
        vector<double> per_wss_rel;
        vector<double> per_wss_afford;
        if (!wss_collectors.empty()) {
            vector<vector<WSSDataCollector *>> utility_wss_mode5;
            isolateWSSDataCollectors(u, utility_wss_mode5);
            if (!utility_wss_mode5.empty()) {
                per_wss_rel    = ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                    utility_wss_mode5, realizations_ran);
                per_wss_afford = ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                    utility_wss_mode5, realizations_ran);
            }
        }
        outStream << setw(COLUMN_WIDTH) << u[realizations_ran[0]]->name;
        for (double r : per_wss_rel)
            outStream << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << r;
        outStream << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << restriction_freq
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << inf_npc
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << worse_cost;
        for (double a : per_wss_afford)
            outStream << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << a;
        outStream << endl;
        for (double r : per_wss_rel)    objectives.push_back(r);
        objectives.push_back(restriction_freq);
        objectives.push_back(inf_npc);
        objectives.push_back(worse_cost);
        for (double a : per_wss_afford) objectives.push_back(a);
    } else if (Constants::includeSingleWSSObjectives()) {
        // Mode 6: single target WSS reliability and affordability
        vector<double> per_wss_rel;
        vector<double> per_wss_afford;
        if (!wss_collectors.empty()) {
            vector<vector<WSSDataCollector *>> utility_wss_mode6;
            isolateWSSDataCollectors(u, utility_wss_mode6);
            if (!utility_wss_mode6.empty()) {
                per_wss_rel    = ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                    utility_wss_mode6, realizations_ran);
                per_wss_afford = ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                    utility_wss_mode6, realizations_ran);
            }
        }
        int target = Constants::TARGET_WSS_ID;
        double target_rel    = (target < (int)per_wss_rel.size())    ? per_wss_rel[target]    : reliability;
        double target_afford = (target < (int)per_wss_afford.size()) ? per_wss_afford[target] : affordability_index;
        outStream << setw(COLUMN_WIDTH) << u[realizations_ran[0]]->name
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << target_rel
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << restriction_freq
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << inf_npc
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << worse_cost
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << target_afford
                  << endl;
        objectives.push_back(target_rel);
        objectives.push_back(restriction_freq);
        objectives.push_back(inf_npc);
        objectives.push_back(worse_cost);
        objectives.push_back(target_afford);
    } else if (Constants::includePenaltyObjectives()) {
        // Mode 7: penalty-based objectives
        // diff_rel    = sum_wss[ max(REL_THRESHOLD   - rel_wss,    0)^2 ]
        // diff_afford = 100 * sum_wss[ max(afford_wss - AFFORD_THRESHOLD, 0)^2 ]
        vector<double> per_wss_rel;
        vector<double> per_wss_afford;
        if (!wss_collectors.empty()) {
            vector<vector<WSSDataCollector *>> utility_wss_mode7;
            isolateWSSDataCollectors(u, utility_wss_mode7);
            if (!utility_wss_mode7.empty()) {
                per_wss_rel    = ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                    utility_wss_mode7, realizations_ran);
                per_wss_afford = ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                    utility_wss_mode7, realizations_ran);
            }
        }
        double diff_rel    = 0.0;
        double diff_afford = 0.0;
        for (double r : per_wss_rel) {
            double d = std::max(Constants::PENALTY_RELIABILITY_THRESHOLD - r, 0.0);
            diff_rel += d * d;
        }
        diff_rel *= 100.0;
        for (double a : per_wss_afford) {
            double d = std::max(a - Constants::PENALTY_AFFORDABILITY_THRESHOLD, 0.0);
            diff_afford += d * d;
        }
        diff_afford *= 100.0;
        outStream << setw(COLUMN_WIDTH) << u[realizations_ran[0]]->name
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << diff_rel
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << restriction_freq
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << inf_npc
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << worse_cost
                  << setw(COLUMN_WIDTH * 2) << setprecision(COLUMN_PRECISION) << diff_afford
                  << endl;
        objectives.push_back(diff_rel);
        objectives.push_back(restriction_freq);
        objectives.push_back(inf_npc);
        objectives.push_back(worse_cost);
        objectives.push_back(diff_afford);
    } else {
        // Modes 1-4: aggregated objectives
        outStream << setw(COLUMN_WIDTH) << u[realizations_ran[0]]->name
                  /// Reliability
                  << setw(COLUMN_WIDTH * 2)
                  << setprecision(COLUMN_PRECISION)
                  << reliability
                  /// Restriction Frequency
                  << setw(COLUMN_WIDTH * 2)
                  << setprecision(COLUMN_PRECISION)
                  << restriction_freq
                  /// Infrastructure NPC
                  << setw(COLUMN_WIDTH * 2)
                  << setprecision(COLUMN_PRECISION)
                  << inf_npc
                  /// Worse Case Costs
                  << setw(COLUMN_WIDTH * 2)
                  << setprecision(COLUMN_PRECISION)
                  << worse_cost;
        if (Constants::includeAffordabilityObjective()) {
            outStream << setw(COLUMN_WIDTH * 2)
                      << setprecision(COLUMN_PRECISION)
                      << affordability_index;
        }
        if (Constants::includeSeverityObjective()) {
            outStream << setw(COLUMN_WIDTH * 2)
                      << setprecision(COLUMN_PRECISION)
                      << failure_severity;
        }
        outStream << endl;
        objectives.push_back(reliability);
        objectives.push_back(restriction_freq);
        objectives.push_back(inf_npc);
        objectives.push_back(worse_cost);
        if (Constants::includeAffordabilityObjective()) {
            objectives.push_back(affordability_index);
        }
        if (Constants::includeSeverityObjective()) {
            objectives.push_back(failure_severity);
        }
    }
    } catch (const std::exception& e) {
        printf("ERROR in printUtilityObjectivesToRowOutStream: %s\n", e.what());
        throw;
    }
}

vector<double> MasterDataCollector::calculatePrintObjectives(string file_name, bool print) {
    vector<double> objectives = vector<double>();

    if (print) {
        printf("Calculating and printing Objectives\n");
        string obj_file_path = output_directory + file_name + ".out";
//        cout << obj_file_path << endl;

        std::ofstream outStream;
        outStream.open(obj_file_path);

        // Dynamic header based on experiment mode
        if (Constants::includePerWSSObjectives()) {
            // Mode 5: per-WSS headers
            outStream << setw(COLUMN_WIDTH) << "      "
                      << setw(COLUMN_WIDTH * 2) << "Reliability WSS0"
                      << setw(COLUMN_WIDTH * 2) << "Reliability WSS1"
                      << setw(COLUMN_WIDTH * 2) << "Restriction Freq."
                      << setw(COLUMN_WIDTH * 2) << "Infrastructure NPC"
                      << setw(COLUMN_WIDTH * 2) << "Worse Case Costs"
                      << setw(COLUMN_WIDTH * 2) << "Affordability WSS0"
                      << setw(COLUMN_WIDTH * 2) << "Affordability WSS1"
                      << endl;
        } else if (Constants::includeSingleWSSObjectives()) {
            // Mode 6: single target WSS headers
            string wss_label = "WSS" + to_string(Constants::TARGET_WSS_ID);
            outStream << setw(COLUMN_WIDTH) << "      "
                      << setw(COLUMN_WIDTH * 2) << ("Reliability " + wss_label)
                      << setw(COLUMN_WIDTH * 2) << "Restriction Freq."
                      << setw(COLUMN_WIDTH * 2) << "Infrastructure NPC"
                      << setw(COLUMN_WIDTH * 2) << "Worse Case Costs"
                      << setw(COLUMN_WIDTH * 2) << ("Affordability " + wss_label)
                      << endl;
        } else if (Constants::includePenaltyObjectives()) {
            // Mode 7: penalty-based headers
            outStream << setw(COLUMN_WIDTH) << "      "
                      << setw(COLUMN_WIDTH * 2) << "Diff Reliability"
                      << setw(COLUMN_WIDTH * 2) << "Restriction Freq."
                      << setw(COLUMN_WIDTH * 2) << "Infrastructure NPC"
                      << setw(COLUMN_WIDTH * 2) << "Worse Case Costs"
                      << setw(COLUMN_WIDTH * 2) << "Diff Affordability"
                      << endl;
        } else {
            outStream << setw(COLUMN_WIDTH) << "      " << setw((COLUMN_WIDTH * 2))
                      << "Reliability"
                      << setw(COLUMN_WIDTH * 2) << "Restriction Freq."
                      //              << setw(COLUMN_WIDTH * 2) << "Jordan Lake Alloc."
                      << setw(COLUMN_WIDTH * 2) << "Infrastructure NPC"
                      << setw(COLUMN_WIDTH * 2) << "Worse Case Costs";
            // Only add affordability header if included in experiment
            if (Constants::includeAffordabilityObjective()) {
                outStream << setw(COLUMN_WIDTH * 2) << "Affordability Index";
            }
            if (Constants::includeSeverityObjective()) {
                outStream << setw(COLUMN_WIDTH * 2) << "Failure Severity";
            }
            outStream << endl;
        }

        for (auto &u : utility_collectors) {
            printUtilityObjectivesToRowOutStream(u, outStream, objectives);
        }

        outStream.close();

//        for (int i = 0; i < (int) objectives.size(); ++i) {
//            double o = objectives.at(i);
//            if (o > 10e10 || o < -0.1) {
//                char error[512];
//                sprintf(error, "Objective %d has absurd value of %f. Aborting.\n", i, o);
//                throw_with_nested(runtime_error(error));
//            }
//        }
    } else {
//        cout << "Calculating Objectives" << endl;    
        for (auto &u : utility_collectors) {
            // Create vector with restriction policies pertaining only to the
            // utility whose objectives are being calculated.
            vector<RestrictionsDataCollector *> utility_restrictions(
                    *max_element(realizations_ran.begin(), realizations_ran.end()) + 1,
                    nullptr
                    );
            isolateRestrictionDataCollectors(u, utility_restrictions);

            // Reliability - Use WSS-level calculation (utility = minimum reliability of its WSS)
            if (Constants::includePerWSSObjectives()) {
                // Mode 5: push one reliability objective per WSS
                if (!wss_collectors.empty()) {
                    vector<vector<WSSDataCollector *>> utility_wss_mode5;
                    isolateWSSDataCollectors(u, utility_wss_mode5);
                    if (!utility_wss_mode5.empty()) {
                        vector<double> per_wss_rel =
                            ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                                utility_wss_mode5, realizations_ran);
                        for (double r : per_wss_rel) objectives.push_back(r);
                    } else {
                        objectives.push_back(ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran));
                        objectives.push_back(ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran));
                    }
                } else {
                    double rel = ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran);
                    objectives.push_back(rel);
                    objectives.push_back(rel);
                }
            } else if (Constants::includeSingleWSSObjectives()) {
                // Mode 6: push reliability for target WSS only
                if (!wss_collectors.empty()) {
                    vector<vector<WSSDataCollector *>> utility_wss_mode6;
                    isolateWSSDataCollectors(u, utility_wss_mode6);
                    if (!utility_wss_mode6.empty()) {
                        vector<double> per_wss_rel =
                            ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                                utility_wss_mode6, realizations_ran);
                        int target = Constants::TARGET_WSS_ID;
                        objectives.push_back(target < (int)per_wss_rel.size() ? per_wss_rel[target] : 1.0);
                    } else {
                        objectives.push_back(ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran));
                    }
                } else {
                    objectives.push_back(ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran));
                }
            } else if (Constants::includePenaltyObjectives()) {
                // Mode 7: push diff_rel = sum_wss[ max(threshold - rel_wss, 0)^2 ]
                double diff_rel = 0.0;
                if (!wss_collectors.empty()) {
                    vector<vector<WSSDataCollector *>> utility_wss_mode7;
                    isolateWSSDataCollectors(u, utility_wss_mode7);
                    if (!utility_wss_mode7.empty()) {
                        vector<double> per_wss_rel =
                            ObjectivesCalculator::calculateReliabilityObjective_WSS_PerWSS(
                                utility_wss_mode7, realizations_ran);
                        for (double r : per_wss_rel) {
                            double d = std::max(Constants::PENALTY_RELIABILITY_THRESHOLD - r, 0.0);
                            diff_rel += d * d;
                        }
                    }
                }
                objectives.push_back(diff_rel);
            } else if (!wss_collectors.empty()) {
                // Filter WSS collectors to only include those belonging to this utility
                vector<vector<WSSDataCollector *>> utility_wss_collectors;
                isolateWSSDataCollectors(u, utility_wss_collectors);
            

                
                if (!utility_wss_collectors.empty()) {
                    // Use configurable aggregation method based on experiment mode
                    double reliability = ObjectivesCalculator::calculateReliabilityObjective_WSS_Configurable(
                        utility_wss_collectors, 
                        realizations_ran,
                        Constants::getReliabilityAggregationMethod());
                    objectives.push_back(reliability);
                } else {
                    // No WSS found for this utility, fallback to utility-level calculation
                    
                    double reliability = ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran);
                    objectives.push_back(reliability);
                }
            } else {
                #ifdef PARALLEL
                printf("WARNING: wss_collectors is EMPTY, using utility-level reliability calculation\n");
                #endif
                double reliability = ObjectivesCalculator::calculateReliabilityObjective(u, realizations_ran);
                objectives.push_back(reliability);
            }
            
            // Restriction Frequency — use WSS demand multipliers (covers ALL WSS including TortoSM,
            // consistent with print path). Fallback to utility-level if WSS collectors unavailable.
            {
                double restriction_freq = NONE;
                if (!wss_collectors.empty()) {
                    vector<vector<WSSDataCollector *>> utility_wss_collectors_rf;
                    isolateWSSDataCollectors(u, utility_wss_collectors_rf);
                    if (!utility_wss_collectors_rf.empty())
                        restriction_freq = ObjectivesCalculator::
                            calculateRestrictionFrequencyObjective_WSS(utility_wss_collectors_rf, realizations_ran);
                }
                if (restriction_freq == NONE)
                    restriction_freq = ObjectivesCalculator::
                        calculateRestrictionFrequencyObjective(utility_restrictions, realizations_ran);
                objectives.push_back(restriction_freq);
            }
            
            // Use WSS-level calculations for infrastructure NPC and worse case costs
            if (!wss_collectors.empty()) {
                // Filter WSS collectors to only include those belonging to this utility
                vector<vector<WSSDataCollector *>> utility_wss_collectors_npc;
                isolateWSSDataCollectors(u, utility_wss_collectors_npc);
                
                objectives.push_back
                        (ObjectivesCalculator::calculateNetPresentCostInfrastructureObjective_WSS(utility_wss_collectors_npc, u, realizations_ran));
            } else {
                objectives.push_back
                        (ObjectivesCalculator::calculateNetPresentCostInfrastructureObjective(u, realizations_ran));
            }
            
            objectives.push_back
                    (ObjectivesCalculator::calculateWorseCaseCostsObjective(u, realizations_ran));
            
            // Calculate affordability index - only if included in experiment
            if (Constants::includeAffordabilityObjective()) {
                if (Constants::includePerWSSObjectives()) {
                    // Mode 5: push one affordability objective per WSS
                    if (!wss_collectors.empty()) {
                        vector<vector<WSSDataCollector *>> utility_wss_collectors_affordability;
                        isolateWSSDataCollectors(u, utility_wss_collectors_affordability);
                        if (!utility_wss_collectors_affordability.empty()) {
                            vector<double> per_wss_afford =
                                ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                                    utility_wss_collectors_affordability, realizations_ran);
                            for (double a : per_wss_afford) objectives.push_back(a);
                        } else {
                            objectives.push_back(1.0);
                            objectives.push_back(1.0);
                        }
                    } else {
                        objectives.push_back(1.0);
                        objectives.push_back(1.0);
                    }
                } else if (Constants::includeSingleWSSObjectives()) {
                    // Mode 6: push affordability for target WSS only
                    if (!wss_collectors.empty()) {
                        vector<vector<WSSDataCollector *>> utility_wss_collectors_affordability;
                        isolateWSSDataCollectors(u, utility_wss_collectors_affordability);
                        if (!utility_wss_collectors_affordability.empty()) {
                            vector<double> per_wss_afford =
                                ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                                    utility_wss_collectors_affordability, realizations_ran);
                            int target = Constants::TARGET_WSS_ID;
                            objectives.push_back(target < (int)per_wss_afford.size() ? per_wss_afford[target] : 1.0);
                        } else {
                            objectives.push_back(1.0);
                        }
                    } else {
                        objectives.push_back(1.0);
                    }
                } else if (Constants::includePenaltyObjectives()) {
                    // Mode 7: push diff_afford = sum_wss[ max(afford_wss - threshold, 0)^2 ]
                    double diff_afford = 0.0;
                    if (!wss_collectors.empty()) {
                        vector<vector<WSSDataCollector *>> utility_wss_collectors_affordability;
                        isolateWSSDataCollectors(u, utility_wss_collectors_affordability);
                        if (!utility_wss_collectors_affordability.empty()) {
                            vector<double> per_wss_afford =
                                ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_PerWSS(
                                    utility_wss_collectors_affordability, realizations_ran);
                            for (double a : per_wss_afford) {
                                double d = std::max(a - Constants::PENALTY_AFFORDABILITY_THRESHOLD, 0.0);
                                diff_afford += d * d;
                            }
                        }
                    }
                    objectives.push_back(diff_afford);
                } else if (!wss_collectors.empty()) {
                    // Modes 3/4: aggregated affordability
                    vector<vector<WSSDataCollector *>> utility_wss_collectors_affordability;
                    isolateWSSDataCollectors(u, utility_wss_collectors_affordability);
                    objectives.push_back
                            (ObjectivesCalculator::calculateAffordabilityIndexObjective_WSS_Configurable(
                                utility_wss_collectors_affordability,
                                realizations_ran,
                                Constants::getAffordabilityAggregationMethod()));
                } else {
                    objectives.push_back(1.0);
                }
            }

            // Calculate failure severity objective - only if included (not supported in mode 5)
            if (Constants::includeSeverityObjective() && !Constants::includePerWSSObjectives()) {
                if (!wss_collectors.empty()) {
                    vector<vector<WSSDataCollector *>> utility_wss_collectors_severity;
                    isolateWSSDataCollectors(u, utility_wss_collectors_severity);
                    
                    objectives.push_back
                            (ObjectivesCalculator::calculateFailureSeverityObjective_WSS_Configurable(
                                utility_wss_collectors_severity,
                                realizations_ran,
                                Constants::getReliabilityAggregationMethod()));
                } else {
                    objectives.push_back(0.0);
                }
            }
        }
    }
    return objectives;
}

void MasterDataCollector::isolateRestrictionDataCollectors(vector<UtilitiesDataCollector *> &u,
                                                           vector<RestrictionsDataCollector *> &utility_restrictions) const {
    for (auto &p : drought_mitigation_policy_collectors)
        if (p.at(realizations_ran.at(0))->type == RESTRICTIONS && p[realizations_ran[0]]->id == u.at(realizations_ran[0])->id)
            for (auto i : realizations_ran) {
                utility_restrictions.at(i) =
                        dynamic_cast<RestrictionsDataCollector *>(p.at(i));
            }
}

void MasterDataCollector::isolateWSSDataCollectors(vector<UtilitiesDataCollector *> &u,
                                                   vector<vector<WSSDataCollector *>> &utility_wss_collectors) const {
    // Filter WSS collectors to include only those belonging to the specified utility
    utility_wss_collectors.clear();
    int utility_id = 0;
    if (!u.empty() && !realizations_ran.empty() && u.at(realizations_ran.at(0)) != nullptr) {
        utility_id = u.at(realizations_ran.at(0))->id;
    }
    
    for (const auto &wss_realization_data : wss_collectors) {
        if (!wss_realization_data.empty() && 
            wss_realization_data[realizations_ran[0]] != nullptr) {
            
            // Get stored owner ID (safe even after WSS is deleted)
            int owner_id = wss_realization_data[realizations_ran[0]]->getOwnerId();
            
            // Check using stored ID (safe even if owner pointer is dangling)
            if (owner_id == utility_id) {
                utility_wss_collectors.push_back(wss_realization_data);
            }
        }
    }
}

unsigned long MasterDataCollector::getActualWeeksCollected() const {
    if (!utility_collectors.empty() && !utility_collectors[0].empty() && 
        !realizations_ran.empty()) {
        return utility_collectors[0][realizations_ran[0]]->getCombined_storage().size();
    }
    return 0;
}

void MasterDataCollector::printWeeklyReliabilityByWSS() const {
    if (wss_collectors.empty() || realizations_ran.empty()) {
        printf("Weekly WSS reliability: no WSS data available.\n");
        return;
    }

    printf("\nWeekly reliability by WSS (threshold ratio = %.3f)\n", STORAGE_CAPACITY_RATIO_FAIL);

    for (const auto &wss_realization_data : wss_collectors) {
        int wss_id = NON_INITIALIZED;
        unsigned long n_weeks = 0;

        for (const auto &r : realizations_ran) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                wss_id = wss_realization_data[r]->id;
                n_weeks = std::max(n_weeks, (unsigned long) wss_realization_data[r]->getCombined_storage().size());
            }
        }

        if (n_weeks == 0) {
            continue;
        }

        printf("WSS %d:\n", wss_id);

        for (unsigned long w = 0; w < n_weeks; ++w) {
            int valid_realizations = 0;
            int failed_realizations = 0;

            for (const auto &r : realizations_ran) {
                if (r >= wss_realization_data.size() || wss_realization_data[r] == nullptr) {
                    continue;
                }

                const auto &failure_flag = wss_realization_data[r]->getWeekly_failure_flag();

                if (w >= failure_flag.size()) {
                    continue;
                }

                ++valid_realizations;
                if (failure_flag[w] == 1) {
                    ++failed_realizations;
                }
            }

            if (valid_realizations == 0) {
                printf("  Week %lu: N/A (no storage-based realizations)\n", w);
            } else {
                double weekly_reliability = 1.0 - static_cast<double>(failed_realizations) / valid_realizations;
                printf("  Week %lu: %.6f (%d/%d non-fail)\n",
                       w,
                       weekly_reliability,
                       valid_realizations - failed_realizations,
                       valid_realizations);
            }
        }
    }
}

void MasterDataCollector::printWeeklyReliabilityByWSSCsv(const string &file_name) const {
    if (wss_collectors.empty() || realizations_ran.empty()) {
        printf("Weekly WSS reliability CSV: no WSS data available.\n");
        return;
    }

    if (output_directory.empty()) {
        printf("Weekly WSS reliability CSV: output directory is not set.\n");
        return;
    }

    string csv_path = output_directory + file_name + ".csv";
    ofstream out_stream;
    out_stream.open(csv_path);

    if (!out_stream.is_open()) {
        char error_msg[512];
        sprintf(error_msg, "Could not open weekly reliability CSV for writing: %s", csv_path.c_str());
        throw runtime_error(error_msg);
    }

    const vector<WSSDataCollector *> *wss0 = nullptr;
    const vector<WSSDataCollector *> *wss1 = nullptr;

    for (const auto &wss_realization_data : wss_collectors) {
        int wss_id = NON_INITIALIZED;
        for (const auto &r : realizations_ran) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                wss_id = wss_realization_data[r]->id;
                break;
            }
        }

        if (wss_id == 0) {
            wss0 = &wss_realization_data;
        } else if (wss_id == 1) {
            wss1 = &wss_realization_data;
        }
    }

    auto get_n_weeks = [&](const vector<WSSDataCollector *> *wss_data) {
        unsigned long n_weeks = 0;
        if (wss_data == nullptr) return n_weeks;
        for (const auto &r : realizations_ran) {
            if (r < wss_data->size() && wss_data->at(r) != nullptr) {
                n_weeks = std::max(n_weeks, (unsigned long) wss_data->at(r)->getCombined_storage().size());
            }
        }
        return n_weeks;
    };

    auto annual_reliability = [&](const vector<WSSDataCollector *> *wss_data, unsigned long y, bool &has_data) {
        has_data = false;
        if (wss_data == nullptr || realizations_ran.empty()) return 0.0;

        int failed_realizations = 0;
        const int n_realizations = (int) realizations_ran.size();
        int year_start = (int) round(y * WEEKS_IN_YEAR);
        int year_end = (int) round((y + 1) * WEEKS_IN_YEAR);

        for (const auto &r : realizations_ran) {
            if (r >= wss_data->size() || wss_data->at(r) == nullptr) {
                continue;
            }

            const auto &failure_flag = wss_data->at(r)->getWeekly_failure_flag();
            if (failure_flag.empty()) {
                continue;
            }

            bool realization_failed = false;
            for (int w = year_start; w < year_end; ++w) {
                if (w >= (int) failure_flag.size()) {
                    continue;
                }

                has_data = true;
                if (failure_flag[w] == 1) {
                    realization_failed = true;
                    break;
                }
            }

            if (realization_failed) {
                failed_realizations++;
            }
        }

        if (!has_data || n_realizations == 0) {
            return 0.0;
        }

        return 1.0 - static_cast<double>(failed_realizations) / n_realizations;
    };

    auto annual_affordability = [&](const vector<WSSDataCollector *> *wss_data, unsigned long y, bool &has_data) {
        has_data = false;
        if (wss_data == nullptr) return 0.0;

        vector<double> realization_max_affordability;
        int year_start = (int) round(y * WEEKS_IN_YEAR);
        int year_end = (int) round((y + 1) * WEEKS_IN_YEAR);

        for (const auto &r : realizations_ran) {
            if (r >= wss_data->size() || wss_data->at(r) == nullptr) {
                continue;
            }

            const auto &residential_prices = wss_data->at(r)->getResidential_current_price();
            const auto &restricted_demand = wss_data->at(r)->getRestricted_demand();
            const auto &demand_offset = wss_data->at(r)->getDemand_offset();
            double average_monthly_income = wss_data->at(r)->getAverage_monthly_income();
            double initial_households = wss_data->at(r)->getInitial_households();
            double weekly_average_income = (average_monthly_income > 0.0)
                                                 ? (average_monthly_income / WEEKS_IN_MONTH)
                                                 : 0.0;

            if (weekly_average_income <= 0.0 || initial_households <= 0.0) {
                continue;
            }

            int max_week = std::min(year_end, (int) std::min(residential_prices.size(), restricted_demand.size()));
            if (year_start >= max_week) {
                continue;
            }

            double max_affordability = 0.0;
            bool has_value = false;
            for (int w = year_start; w < max_week; ++w) {
                double weekly_price = residential_prices[w];
                // Add demand_offset back: transfers are paid by the utility, not residents
                double weekly_demand = restricted_demand[w] + (w < (int)demand_offset.size() ? demand_offset[w] : 0.0);
                if (weekly_price > 0.0 && weekly_demand > 0.0) {
                    double weekly_cost = (weekly_price * weekly_demand) / initial_households;
                    double affordability = weekly_cost / weekly_average_income;
                    if (!has_value || affordability > max_affordability) {
                        max_affordability = affordability;
                        has_value = true;
                    }
                }
            }

            if (has_value) {
                realization_max_affordability.push_back(max_affordability);
                has_data = true;
            }
        }

        if (realization_max_affordability.empty()) {
            return 0.0;
        }

        sort(realization_max_affordability.begin(), realization_max_affordability.end());
        size_t index_95th = (size_t) (0.95 * (realization_max_affordability.size() - 1));
        return realization_max_affordability[index_95th];
    };

    auto annual_severity = [&](const vector<WSSDataCollector *> *wss_data, unsigned long y, bool &has_data) {
        has_data = false;
        if (wss_data == nullptr || realizations_ran.empty()) return 0.0;

        vector<double> realization_failure_weeks; // failure weeks per realization for this year
        int year_start = (int) round(y * WEEKS_IN_YEAR);
        int year_end = (int) round((y + 1) * WEEKS_IN_YEAR);

        for (const auto &r : realizations_ran) {
            if (r >= wss_data->size() || wss_data->at(r) == nullptr) {
                continue;
            }

            const auto &failure_flag = wss_data->at(r)->getWeekly_failure_flag();
            if (failure_flag.empty()) {
                continue;
            }

            double realization_count = 0.0;
            for (int w = year_start; w < year_end; ++w) {
                if (w >= (int) failure_flag.size()) {
                    continue;
                }

                has_data = true;
                if (failure_flag[w] == 1) {
                    realization_count += 1.0;
                }
            }
            realization_failure_weeks.push_back(realization_count);
        }

        if (!has_data || realization_failure_weeks.empty()) {
            return 0.0;
        }

        // Sort and take 95th percentile (WCC-like approach)
        sort(realization_failure_weeks.begin(), realization_failure_weeks.end());
        unsigned long percentile_index = (unsigned long) floor(
                WORSE_CASE_COST_PERCENTILE * realization_failure_weeks.size());
        if (percentile_index >= realization_failure_weeks.size()) {
            percentile_index = realization_failure_weeks.size() - 1;
        }
        return realization_failure_weeks.at(percentile_index);
    };

    unsigned long n_weeks = std::max(get_n_weeks(wss0), get_n_weeks(wss1));
    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);
    out_stream << "0reliability,1reliability,0afford,1afford,0severity,1severity" << endl;

    for (unsigned long y = 0; y < n_years; ++y) {
        bool has_0 = false;
        bool has_1 = false;
        bool has_afford_0 = false;
        bool has_afford_1 = false;
        bool has_sev_0 = false;
        bool has_sev_1 = false;
        double rel_0 = annual_reliability(wss0, y, has_0);
        double rel_1 = annual_reliability(wss1, y, has_1);
        double afford_0 = annual_affordability(wss0, y, has_afford_0);
        double afford_1 = annual_affordability(wss1, y, has_afford_1);
        double sev_0 = annual_severity(wss0, y, has_sev_0);
        double sev_1 = annual_severity(wss1, y, has_sev_1);

        if (has_0) {
            out_stream << rel_0;
        }
        out_stream << ",";
        if (has_1) {
            out_stream << rel_1;
        }
        out_stream << ",";
        if (has_afford_0) {
            out_stream << afford_0;
        }
        out_stream << ",";
        if (has_afford_1) {
            out_stream << afford_1;
        }
        out_stream << ",";
        if (has_sev_0) {
            out_stream << sev_0;
        }
        out_stream << ",";
        if (has_sev_1) {
            out_stream << sev_1;
        }
        out_stream << endl;
    }

    out_stream.close();
}

void MasterDataCollector::printAnnualReliabilityBySourceCsv(const string &file_name) const {
    if (wss_collectors.empty() || realizations_ran.empty()) {
        printf("Annual source reliability CSV: no WSS data available.\n");
        return;
    }

    if (output_directory.empty()) {
        printf("Annual source reliability CSV: output directory is not set.\n");
        return;
    }

    string csv_path = output_directory + file_name + ".csv";
    ofstream out_stream;
    out_stream.open(csv_path);

    if (!out_stream.is_open()) {
        char error_msg[512];
        sprintf(error_msg, "Could not open annual source reliability CSV for writing: %s", csv_path.c_str());
        throw runtime_error(error_msg);
    }

    // Collect all source IDs and names across all WSS and realizations
    // Use a map to preserve source_id ordering
    map<int, string> all_source_info;  // source_id -> "wssId_sourceName"
    unsigned long n_weeks = 0;

    for (const auto &wss_realization_data : wss_collectors) {
        for (const auto &r : realizations_ran) {
            if (r < wss_realization_data.size() && wss_realization_data[r] != nullptr) {
                int wss_id = wss_realization_data[r]->id;
                const auto &src_names = wss_realization_data[r]->getSource_names();
                const auto &src_flags = wss_realization_data[r]->getPer_source_failure_flag();
                for (const auto &entry : src_names) {
                    // Only add if not already present (first WSS to claim it wins)
                    if (all_source_info.find(entry.first) == all_source_info.end()) {
                        all_source_info[entry.first] = to_string(wss_id) + "_" + entry.second;
                    }
                }
                for (const auto &entry : src_flags) {
                    n_weeks = std::max(n_weeks, (unsigned long)entry.second.size());
                }
            }
        }
    }

    if (all_source_info.empty() || n_weeks == 0) {
        out_stream << "No source data collected." << endl;
        out_stream.close();
        return;
    }

    unsigned long n_years = (unsigned long) round(n_weeks / WEEKS_IN_YEAR);

    // Write header: source columns in order of source_id
    bool first = true;
    for (const auto &entry : all_source_info) {
        if (!first) out_stream << ",";
        out_stream << entry.second << "(id" << entry.first << ")";
        first = false;
    }
    out_stream << endl;

    // For each year, calculate per-source reliability across realizations
    for (unsigned long y = 0; y < n_years; ++y) {
        int year_start = (int) round(y * WEEKS_IN_YEAR);
        int year_end = (int) round((y + 1) * WEEKS_IN_YEAR);

        // For each source, count how many realizations had a failure in this year
        map<int, double> source_reliability;
        for (const auto &entry : all_source_info) {
            int src_id = entry.first;
            int failed_realizations = 0;
            int counted_realizations = 0;

            for (const auto &wss_realization_data : wss_collectors) {
                for (const auto &r : realizations_ran) {
                    if (r >= wss_realization_data.size() || wss_realization_data[r] == nullptr) {
                        continue;
                    }

                    const auto &src_flags = wss_realization_data[r]->getPer_source_failure_flag();
                    auto it = src_flags.find(src_id);
                    if (it == src_flags.end()) continue;

                    const auto &flags = it->second;
                    counted_realizations++;
                    bool realization_failed = false;

                    for (int w = year_start; w < year_end; ++w) {
                        if (w >= (int) flags.size()) continue;
                        if (flags[w] == 1) {
                            realization_failed = true;
                            break;
                        }
                    }

                    if (realization_failed) {
                        failed_realizations++;
                    }
                }
            }

            if (counted_realizations > 0) {
                source_reliability[src_id] = 1.0 - (double) failed_realizations / counted_realizations;
            } else {
                source_reliability[src_id] = -1.0;  // No data marker
            }
        }

        // Write row
        first = true;
        for (const auto &entry : all_source_info) {
            if (!first) out_stream << ",";
            double rel = source_reliability[entry.first];
            if (rel >= 0.0) {
                out_stream << rel;
            }
            first = false;
        }
        out_stream << endl;
    }

    out_stream.close();
}

void MasterDataCollector::performBootstrapAnalysis(
		int sol_id, int n_sets, int n_samples, int n_threads, vector<vector<int>> bootstrap_samples) {
    printf("Running bootstrap samples.\n");
    vector<vector<int>> bootstrap_sample_sets((unsigned long) n_sets, vector<int>((unsigned long) n_samples));

    // Create or use specified bootstrap samples
    readOrCreateBSSamples(sol_id, n_sets, n_samples, bootstrap_samples, bootstrap_sample_sets);

    vector<vector<double>> objectives((unsigned long) n_sets);
    for (unsigned long &r : crashed_realizations) {
        for (vector<int> &bs : bootstrap_sample_sets) {
            bs.erase(remove(bs.begin(), bs.end(), r), bs.end());
        }
    }

//#pragma omp parallel for num_threads(n_threads) shared(objectives)
    for (int set = 0; set < n_sets; ++set) {
        // Calculate objectives for the set of bootstrapped realizations.
        vector<unsigned long> bootstrap_sample_set = vector<unsigned long>(
                bootstrap_sample_sets[set].begin(),
                bootstrap_sample_sets[set].end());

        for (unsigned long &r : crashed_realizations) {
            for (unsigned long &bs : bootstrap_sample_set) {
                if (bs >= r) {
                    --bs;
                }
            }
        }

        for (auto &u : utility_collectors) {
            // Create vector with restriction policies pertaining only to the
            // utility whose objectives are being calculated.
            vector<RestrictionsDataCollector *> utility_restrictions(
                    *max_element(realizations_ran.begin(), realizations_ran.end()) + 1
            );
            isolateRestrictionDataCollectors(u, utility_restrictions);

            // Populate vector of objectives for each corresponding set of bootstrap samples.
            objectives[set].push_back
                    (ObjectivesCalculator::calculateReliabilityObjective(u, bootstrap_sample_set));
            objectives[set].push_back
                    (ObjectivesCalculator::calculateRestrictionFrequencyObjective(utility_restrictions, bootstrap_sample_set));
            objectives[set].push_back
                    (ObjectivesCalculator::calculateNetPresentCostInfrastructureObjective(u, bootstrap_sample_set));
                objectives[set].push_back
                    (ObjectivesCalculator::calculateWorseCaseCostsObjective(u, bootstrap_sample_set));
        }
    }

    // Print objectives of bootstrap samples
    printObjsBSSamples(sol_id, n_sets, n_samples, objectives);

    // Print objectives of all realizations
    printObjectivesOfAllRealizationsForBSAnalysis(sol_id, n_sets, n_samples);

    // Print bootstrap samples file.
    printBSSamples(sol_id, n_sets, n_samples, bootstrap_sample_sets);

}

void MasterDataCollector::printBSSamples(int sol_id, int n_sets, int n_samples,
                                         const vector<vector<int>> &bootstrap_sample_sets) const {
    ofstream outStream_realizations; // Either read samples from file or create new ones.
    outStream_realizations.open(output_directory + "bootstrap_realizations_" +
                                to_string(n_sets) + "_" + to_string(n_samples) + "_S" +
                                to_string(sol_id) + ".csv");

    string line;
    for (int set = 0; set < n_sets; ++set) {
        // Generate one set of bootstrapped realizations, if none was specified.
        line = "";
        for (int s : bootstrap_sample_sets[set]) {
            line += to_string(s) + ",";
        }
        line.pop_back();
        outStream_realizations << line << endl;
    }

    outStream_realizations.close();
}

void MasterDataCollector::printObjectivesOfAllRealizationsForBSAnalysis(int sol_id, int n_sets, int n_samples) {
    string file_name = output_directory + "objectives_all_reals_" + to_string(n_sets) +
                       "_" + to_string(n_samples) + "_S" + to_string(sol_id) + ".csv";
    vector<double> objectives_all_reals = calculatePrintObjectives("", false);

    string line;
    line = "";
    for (double &o : objectives_all_reals) {
	    line += to_string(o) + ",";
    }
    line.pop_back();

    ofstream outStream_objs_all_reals;
    outStream_objs_all_reals.open(file_name);
    outStream_objs_all_reals << line << endl;

    outStream_objs_all_reals.close();
}

void MasterDataCollector::printObjsBSSamples(int sol_id, int n_sets, int n_samples,
                                              vector<vector<double>> &objectives) {// Print objectives.
    ofstream outStream_objs;
    string objectives_file_name = output_directory + "bootstrap_objs_" + to_string(n_sets) + "_" +
                        to_string(n_samples) + "_S" + to_string(sol_id) + ".csv";
    outStream_objs.open(objectives_file_name);
    printf("Bootstrap objectives files will be printed at %s\n", objectives_file_name.c_str());

    string line;
    for (int set = 0; set < n_sets; ++set) {
        line = "";
        for (double &o : objectives[set]) {
            line += to_string(o) + ",";
        }
        line.pop_back();
        outStream_objs << line << endl;
    }
    outStream_objs.close();
}

void MasterDataCollector::readOrCreateBSSamples(int sol_id, int n_sets, int n_samples,
                                                const vector<vector<int>> &bootstrap_samples,
                                                vector<vector<int>> &bootstrap_sample_sets) const {
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng((seed == NON_INITIALIZED ? rd() : seed));    // random-number engine used (Mersenne-Twister in this case)

    int min = 0;
    int max = (int) n_realizations - 1;
    uniform_int_distribution<int> uni(min, max); // guaranteed unbiased
    string line;
    if (!bootstrap_samples.empty()) {
	    bootstrap_sample_sets = bootstrap_samples;
    } else {
        for (int set = 0; set < n_sets; ++set) {
            // Generate one set of bootstrapped realizations, if none was specified.
            for (int &s : bootstrap_sample_sets[set]) {
                s = uni(rng);
            }
        }
    }
}

void MasterDataCollector::printPathways(string file_name) {
    std::ofstream outStream;
    outStream.open(output_directory + file_name + ".out");

    outStream << "Realization\tutility\tWSS\tweek\tinfra." << endl;

    struct PathwayRow {
        unsigned long realization;
        int utility_id;
        int wss_id;
        int week;
        int infra_id;
        double long_term_rof;
    };

    auto is_shared_paranoa = [](int infra_id) {
        return infra_id == 7 || infra_id == 8 || infra_id == 9;
    };

    vector<PathwayRow> output_rows;
    map<pair<unsigned long, int>, PathwayRow> best_paranoa_by_realization;

    int total_pathways = 0;
    // Collect pathways from WSS collectors (which point to realization WSS where infrastructure is actually built)
    for (auto &wss_collector_vec : wss_collectors) {
        for (int rr = 0; rr < (int) realizations_ran.size(); ++rr) {
            auto r = realizations_ran[rr];
            if (!wss_collector_vec[r]) {
                continue;
            }

            vector<vector<int>> pathways = wss_collector_vec[r]->getPathways();
            total_pathways += pathways.size();
            for (const vector<int>& infra : pathways) {
                // Format: realization, utility_id, wss_id, week, water_source_id
                PathwayRow row = {r, infra[0], infra[1], infra[2], infra[3], 0.0};

                if (is_shared_paranoa(row.infra_id)) {
                    const auto &lt_rof = wss_collector_vec[r]->getLong_term_rof();
                    if (row.week >= 0 && row.week < (int) lt_rof.size()) {
                        row.long_term_rof = lt_rof[row.week];
                    } else {
                        row.long_term_rof = -1.0;
                    }

                    auto key = make_pair(row.realization, row.infra_id);
                    auto it = best_paranoa_by_realization.find(key);
                    if (it == best_paranoa_by_realization.end()) {
                        best_paranoa_by_realization[key] = row;
                    } else {
                        const auto &best = it->second;
                        bool replace = false;
                        if (row.week < best.week) {
                            replace = true;
                        } else if (row.week == best.week) {
                            if (row.long_term_rof > best.long_term_rof) {
                                replace = true;
                            } else if (row.long_term_rof == best.long_term_rof &&
                                       row.wss_id < best.wss_id) {
                                replace = true;
                            }
                        }

                        if (replace) {
                            it->second = row;
                        }
                    }
                } else {
                    output_rows.push_back(row);
                }
            }
        }
    }

    for (const auto &kv : best_paranoa_by_realization) {
        output_rows.push_back(kv.second);
    }

    sort(output_rows.begin(), output_rows.end(), [](const PathwayRow &a, const PathwayRow &b) {
        if (a.realization != b.realization) return a.realization < b.realization;
        return a.week < b.week;
    });

    for (const auto &row : output_rows) {
        outStream << row.realization << "\t" << row.utility_id << "\t" << row.wss_id << "\t"
                  << row.week << "\t" << row.infra_id << endl;
    }
    outStream.close();
}

void MasterDataCollector::setOutputDirectory(string io_directory) {
    // Check if io_directory is not being set for the same io_directory it is already set. Avoids unnecessary verbose.
    if (io_directory != output_directory) {
        output_directory = io_directory + DEFAULT_OUTPUT_DIR;
        Utils::createDir(output_directory);
        printf("Output will be printed to folder %s\n", output_directory.c_str());
    }
}

DataCollector* MasterDataCollector::createPolicyDataCollector(DroughtMitigationPolicy* dmp, unsigned long r) {
    if (dmp->type == RESTRICTIONS)
        return new RestrictionsDataCollector(dynamic_cast<Restrictions *> (dmp), r);
    else if (dmp->type == TRANSFERS)
        return new TransfersDataCollector(dynamic_cast<Transfers *> (dmp), r);
    else if (dmp->type == TRANSFERS_CAESB)
        return new TransfersBilateralDataCollector(dynamic_cast<TransfersBilateral *> (dmp), r);
    else if (dmp->type == INSURANCE_STORAGE_ROF)
        return new EmptyDataCollector();
    else if (dmp->type == EMERGENCY_TRANSFER_PARANOA)
        return new EmergencyTransferParanoaDataCollector(
                dynamic_cast<EmergencyTransferParanoa *>(dmp), r);
    else
        throw invalid_argument("Drought mitigation policy not recognized. "
                                 "Did you forget to add it to the "
                                 "MasterDataCollector::addRealization"
                                 " function or to create its data collector??");
}

DataCollector* MasterDataCollector::createWaterSourceDataCollector(WaterSource* ws, unsigned long r) {
    if (ws->source_type == RESERVOIR)
        return new ReservoirDataCollector(dynamic_cast<Reservoir *> (ws), r);
    else if (ws->source_type == INTAKE)
        return new IntakeDataCollector(dynamic_cast<Intake *> (ws), r);
    else if (ws->source_type == QUARRY)
        return new QuaryDataCollector(dynamic_cast<Quarry *> (ws), r);
    else if (ws->source_type == WATER_REUSE)
        return new WaterReuseDataCollector(dynamic_cast<WaterReuse *> (ws), r);
    else if (ws->source_type == ALLOCATED_RESERVOIR)
        return new AllocatedReservoirDataCollector(dynamic_cast<AllocatedReservoir *> (ws), r);
    else if (ws->source_type ==
             RESERVOIR_EXPANSION ||
             ws->source_type ==
             NEW_WATER_TREATMENT_PLANT ||
             ws->source_type ==
             SOURCE_RELOCATION)
        return new EmptyDataCollector();
    else
        throw invalid_argument("Water source not recognized. "
                                 "Did you forget to add it to the "
                                 "MasterDataCollector::addRealization"
                                 " function?");
}

void MasterDataCollector::addRealization(
        vector<WaterSource *> water_sources_realization,
        vector<DroughtMitigationPolicy *> drought_mitigation_policies_realization,
        vector<Utility *> utilities_realization,
        const vector<std::unique_ptr<WaterSupplySystems>> &wss_realization,
        unsigned long r,
        const vector<double> &rdm_factors) {
    // If collectors vectors have not yet been initialized, initialize them.
#pragma omp critical
    {
        if (water_source_collectors.empty()) {
            water_source_collectors = vector<vector<DataCollector *>>
                    (water_sources_realization.size(), vector<DataCollector *>(n_realizations));
            drought_mitigation_policy_collectors = vector<vector<DataCollector *>>
                    (drought_mitigation_policies_realization.size(), vector<DataCollector *>(n_realizations));
            utility_collectors = vector<vector<UtilitiesDataCollector *>>
                    (utilities_realization.size(), vector<UtilitiesDataCollector *>(n_realizations));
            
            // Count total number of WSS from the realization model
            // These are the WSS where bonds are actually issued during simulation
            wss_collectors = vector<vector<WSSDataCollector *>>
                    (wss_realization.size(), vector<WSSDataCollector *>(n_realizations));
        }
        realizations_created++;
    };

    // Create utilities data collectors
    // Calculate realization-specific discount rate from RDM factors (deterministic)
    // rdm_factors[3] is the discount rate multiplier
    double base_discount_rate = 0.05; // This should match utility's base rate
    if (!utilities_realization.empty() && utilities_realization[0] != nullptr) {
        base_discount_rate = utilities_realization[0]->getBaseInfraDiscountRate();
    }
    double realization_discount_rate = base_discount_rate * (rdm_factors.size() > 3 ? rdm_factors[3] : 1.0);
    
    for (int u = 0; u < (int) utilities_realization.size(); ++u) {
        utility_collectors[u][r] = new UtilitiesDataCollector(utilities_realization[u], r, realization_discount_rate);
    }

    // Create WSS data collectors pointing to realization WSS (where bonds are actually issued)
    // This ensures we collect the correct infrastructure NPC values
    for (size_t wss_index = 0; wss_index < wss_realization.size(); ++wss_index) {
        if (wss_index >= wss_collectors.size()) {
            char error[256];
            sprintf(error, "WSS collector index %zu out of bounds (size=%zu) for realization %lu",
                    wss_index, wss_collectors.size(), r);
            throw std::out_of_range(error);
        }
        if (wss_realization[wss_index] == nullptr || wss_realization[wss_index].get() == nullptr) {
            char error[256];
            sprintf(error, "Realization WSS at index %zu is null for realization %lu", wss_index, r);
            throw std::runtime_error(error);
        }
        wss_collectors[wss_index][r] = new WSSDataCollector(wss_realization[wss_index].get(), r);
    }

    // Create drought mitigation policies data collector
    for (int dmp = 0; dmp < (int) drought_mitigation_policies_realization.size(); ++dmp)
        drought_mitigation_policy_collectors[dmp][r] =
                createPolicyDataCollector(drought_mitigation_policies_realization[dmp], r);

    // Create water sources data collectors
    for (int ws = 0; ws < (int) water_sources_realization.size(); ++ws) {
        water_source_collectors[ws][r] = createWaterSourceDataCollector(water_sources_realization[ws], r);
    }
} 

void MasterDataCollector::removeRealization(unsigned long r) {
    for (int u = 0; u < (int) utility_collectors.size(); ++u) {
        delete utility_collectors[u][r];
        utility_collectors[u][r] = nullptr;
    }
    for (int wss = 0; wss < (int) wss_collectors.size(); ++wss) {
        delete wss_collectors[wss][r];
        wss_collectors[wss][r] = nullptr;
    }
    for (int dmp = 0; dmp < (int) drought_mitigation_policy_collectors.size(); ++dmp) {
	    delete drought_mitigation_policy_collectors[dmp][r];
        drought_mitigation_policy_collectors[dmp][r] = nullptr;
    }
    for (int ws = 0; ws < (int) water_source_collectors.size(); ++ws) {
	    delete water_source_collectors[ws][r];
        water_source_collectors[ws][r] = nullptr;
    }

    realizations_ran.erase(std::remove(realizations_ran.begin(), realizations_ran.end(), r), realizations_ran.end());
    crashed_realizations.push_back(r);
}


void MasterDataCollector::cleanCollectorsOfDeletedRealizations() {

    for (auto &v : utility_collectors) {
        v.erase(remove_if(v.begin(), v.end(), [](const void *x) { return x == nullptr; }), v.end());
    }
    for (auto &v : wss_collectors) {
        v.erase(remove_if(v.begin(), v.end(), [](const void *x) { return x == nullptr; }), v.end());
    }
    for (auto &v : drought_mitigation_policy_collectors) {
        v.erase(remove_if(v.begin(), v.end(), [](const void *x) { return x == nullptr; }), v.end());
    }
    for (auto &v : water_source_collectors) {
        v.erase(remove_if(v.begin(), v.end(), [](const void *x) { return x == nullptr; }), v.end());
    }
}


void MasterDataCollector::collectData(unsigned long r) {
    for (vector<UtilitiesDataCollector *> &uc : utility_collectors)
        uc[r]->collect_data();
    for (vector<WSSDataCollector *> &wss : wss_collectors)
        wss[r]->collect_data();
    for (vector<DataCollector *> dmp : drought_mitigation_policy_collectors)
        dmp[r]->collect_data();
    for (vector<DataCollector *> ws : water_source_collectors)
        ws[r]->collect_data();
}

void MasterDataCollector::setSeed(int seed) {
    MasterDataCollector::seed = seed;
}

void MasterDataCollector::unsetSeed() {
    MasterDataCollector::seed = NON_INITIALIZED;
}

int MasterDataCollector::getRealizations_created() const {
    return realizations_created;
}