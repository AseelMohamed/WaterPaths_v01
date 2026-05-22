//
// Created by bernardoct on 5/1/17.
//

#include <fstream>
#include <sstream>
#include <stdexcept>
#include "InsuranceStorageToROF.h"
#include "../Utils/Utils.h"

// Helper to convert raw pointer vector to unique_ptr vector
static vector<std::unique_ptr<WaterSupplySystems>> convertToUniquePtr(const vector<WaterSupplySystems*>& wss_raw, bool clear_water_sources) {
    vector<std::unique_ptr<WaterSupplySystems>> result;
    for (auto* wss : wss_raw) {
        auto wss_copy = std::make_unique<WaterSupplySystems>(*wss);
        if (clear_water_sources) {
            wss_copy->clearWaterSources();
        }
        result.push_back(std::move(wss_copy));
    }
    return result;
}

InsuranceStorageToROF::InsuranceStorageToROF(const int id, vector<WaterSource *> &water_sources,
                                             const Graph &water_sources_graph,
                                             const vector<vector<int>> &water_sources_to_wss,
                                             vector<WaterSupplySystems *> &wss,
                                             vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
                                             vector<MinEnvFlowControl *> min_env_flow_controls,
                                             vector<vector<double>>& wss_rdm,
                                             vector<vector<double>>& water_sources_rdm,
                                             vector<vector<double>>& policy_rdm, vector<double> &rof_triggers,
                                             const double insurance_premium, const vector<double> &fixed_payouts,
                                             unsigned long total_simulation_time,
                                             const vector<string> &rof_freq_table_files)
        : DroughtMitigationPolicy(id, INSURANCE_STORAGE_ROF),
          ContinuityModelROF(Utils::copyWaterSourceVector(water_sources), water_sources_graph,
                             water_sources_to_wss, convertToUniquePtr(wss, true),
                             Utils::copyMinEnvFlowControlVector(min_env_flow_controls), wss_rdm.at(0),
                             water_sources_rdm.at(0), total_simulation_time, false, (unsigned int) NON_INITIALIZED),
          rof_triggers(rof_triggers),
          total_simulation_time(total_simulation_time),
          insurance_premium(insurance_premium),
          insurance_price(NONE),
          fixed_payouts(fixed_payouts),
          utilities_revenue_update(vector<double>((unsigned long) n_wss, 0.)),
          utilities_revenue_last_year(vector<double>((unsigned long) n_wss, 0.)),
          drought_mitigation_policies(Utils::copyDroughtMitigationPolicyVector(drought_mitigation_policies)),
          rof_freq_table_files(rof_freq_table_files)
{

    for (WaterSupplySystems *w : wss) wss_ids.push_back(w->system_id);

    for (auto& w : continuity_wss) {
        w->clearWaterSources();
        w->resetTotal_storage_capacity();
    }

    // Remove financial instruments, if any.
    auto it = this->drought_mitigation_policies.begin();
    while(it != this->drought_mitigation_policies.end()) {
        if((*it)->type > SUPPLY_INSTRUMENT_LIM) {
            it = this->drought_mitigation_policies.erase(it);
        }
        else ++it;
    }

    insurance_price = vector<double>((unsigned long) n_wss, 0.);
    payout_multiplier = vector<double>((unsigned long) n_wss, 1.);

    // Load the 4 ROF frequency lookup tables (one per infrastructure state).
    if (rof_freq_table_files.size() != 4)
        throw invalid_argument("InsuranceStorageToROF: exactly 4 ROF frequency table files must be provided.");
    for (const auto &filepath : rof_freq_table_files)
        rof_freq_tables.push_back(loadROFFreqTable(filepath));

    // Derive time-period dimensions from the loaded tables.
    // Each row has: [rof_threshold, fail_freq_period_0, ..., fail_freq_period_(N-1)]
    // years_per_period = 40 years / N periods.
    n_time_periods  = (int) rof_freq_tables[0][0].size() - 1;
    years_per_period = (n_time_periods > 0) ? (40 / n_time_periods) : 1;
}

InsuranceStorageToROF::InsuranceStorageToROF(
        InsuranceStorageToROF &insurance) :
        DroughtMitigationPolicy(insurance.id, insurance.type),
        ContinuityModelROF(Utils::copyWaterSourceVector(insurance.continuity_water_sources),
                           insurance.water_sources_graph,
                           insurance.water_sources_to_wss,
                           Utils::copyWSSVectorUnique(insurance.continuity_wss, true),
                           Utils::copyMinEnvFlowControlVector(insurance.min_env_flow_controls),
                           insurance.wss_rdm,
                           insurance.water_sources_rdm,
                           insurance.total_simulation_time,
                           false, insurance.realization_id),
        rof_triggers(insurance.rof_triggers),
        total_simulation_time(insurance.total_simulation_time),
        insurance_premium(insurance.insurance_premium),
        insurance_price(insurance.insurance_price),
        payout_multiplier(insurance.payout_multiplier),
        fixed_payouts(insurance.fixed_payouts),
        utilities_revenue_update(vector<double>((unsigned long) n_wss, 0.)),
        utilities_revenue_last_year(vector<double>((unsigned long) n_wss, 0.)),
        drought_mitigation_policies(Utils::copyDroughtMitigationPolicyVector(insurance.drought_mitigation_policies)),
        rof_freq_tables(insurance.rof_freq_tables),
        rof_freq_table_files(insurance.rof_freq_table_files),
        n_time_periods(insurance.n_time_periods),
        years_per_period(insurance.years_per_period),
        initial_treatment_capacity(insurance.initial_treatment_capacity),
        initial_storage_capacity(insurance.initial_storage_capacity) {
    wss_ids = insurance.wss_ids;
}

InsuranceStorageToROF::~InsuranceStorageToROF() {
    // Delete drought mitigation policies.
    for (auto dmp : drought_mitigation_policies) {
        delete dmp;
    }
}


void InsuranceStorageToROF::applyPolicy(int week) {
    // Both WSSes belong to ONE utility. Track that utility's gross revenue once
    // per week (iterating over each WSS would double-count the same owner).
    double weekly_gross_rev =
            DroughtMitigationPolicy::realization_wss[wss_ids[0]]->getOwner()->getGrossRevenue();
    utilities_revenue_update[0] += weekly_gross_rev;

    // If first week of the year, price insurance for coming year and snapshot revenue.
    if (Utils::isFirstWeekOfTheYear(week + 1)) {
        utilities_revenue_last_year = utilities_revenue_update;
        std::fill(utilities_revenue_update.begin(), utilities_revenue_update.end(), NONE);

        priceInsurance(week);
    }

    // Do not make payouts during the first year, when no insurance was purchased.
    if (week > WEEKS_IN_YEAR) {
        vector<double> wss_rof(DroughtMitigationPolicy::realization_wss.size(), NONE);

        // Get WSS ROFs
        for (unsigned long u = 0; u < DroughtMitigationPolicy::realization_wss.size(); ++u) {
            wss_rof[u] = DroughtMitigationPolicy::realization_wss[u]->getRisk_of_failure();
        }

        // If ANY WSS triggers, issue a single lump-sum payout to the utility once.
        // The payout is not per-WSS; it is one fixed amount regardless of how many WSSes triggered.
        bool any_triggered = false;
        for (size_t i = 0; i < wss_ids.size(); ++i) {
            if (wss_rof[wss_ids[i]] > rof_triggers[wss_ids[i]]) {
                any_triggered = true;
                break;
            }
        }
        double total_payout = any_triggered
            ? fixed_payouts[wss_ids[0]] * utilities_revenue_last_year[0] * payout_multiplier[0]
            : NONE;

        // Call addInsurancePayout once for the single utility owner.
        DroughtMitigationPolicy::realization_wss[wss_ids[0]]->getOwner()
                ->addInsurancePayout(total_payout);
    }
}

void InsuranceStorageToROF::addSystemComponents(vector<WaterSupplySystems *> wss,
                                                vector<WaterSource *> water_sources,
                                                vector<MinEnvFlowControl *> min_env_flow_controls) {
    DroughtMitigationPolicy::realization_wss = vector<WaterSupplySystems *>(wss.size());
    for (int i : wss_ids) {
        DroughtMitigationPolicy::realization_wss[i] = wss[i];
    }

    // Capture baseline capacities for infrastructure-state detection in priceInsurance.
    // Stored in the same order as wss_ids (index 0 = Descoberto, 1 = TortoSM).
    initial_treatment_capacity.clear();
    initial_storage_capacity.clear();
    for (int id : wss_ids) {
        initial_treatment_capacity.push_back(wss[id]->getTotal_treatment_capacity());
        initial_storage_capacity.push_back(wss[id]->getTotal_storage_capacity());
    }

    for (WaterSupplySystems *w : wss) {
        wss_base_storage_capacity.push_back(w->getTotal_storage_capacity());
    }
    current_storage_table_shift = vector<double>(wss.size());

    connectRealizationWaterSources(water_sources);
}

/**
 * Prices the insurance using pre-computed static ROF-frequency lookup tables.
 *
 * Four CSV tables (one per infrastructure state) encode the combined expected
 * trigger-weeks per year from BOTH WSSes, for a range of ROF trigger values and
 * across multiple time periods (5-year windows over a 40-year horizon).
 *
 * This reflects that:
 *  - Both WSSes belong to ONE utility (single revenue basis, single payout).
 *  - Infrastructure build-out changes trigger statistics (4 states).
 *  - Demand and operating conditions evolve over the simulation (time periods).
 *
 * Price formula:
 *   insurance_price = combined_fail_freq * fixed_payout_rate
 *                     * utility_annual_revenue * insurance_premium
 */
void InsuranceStorageToROF::priceInsurance(int week) {

    // Reset prices.
    for (int u : wss_ids) insurance_price[u] = 0;

    // Bring any newly built infrastructure online in the continuity model.
    updateOnlineInfrastructure(week);

    // Select the lookup table that matches the current infrastructure state.
    int state = getInfraState();

    // Select the time-period column: maps simulation year → 5-year window index.
    // Clamped to the last period for years beyond the 40-year pre-computed horizon.
    int current_year = (int)(week / WEEKS_IN_YEAR);
    int time_period_idx = (years_per_period > 0)
        ? min(current_year / years_per_period, n_time_periods - 1)
        : 0;

    // Combined fail_freq for both WSSes.
    // Each WSS may have a different absolute trigger (restriction_trigger + shared_offset).
    // Use the minimum of both triggers as the representative lookup value.
    double avg_trigger = rof_triggers[wss_ids[0]];
    for (int u : wss_ids) avg_trigger = min(avg_trigger, rof_triggers[u]);
    double combined_fail_freq = lookupFailFreq(state, time_period_idx, avg_trigger);

    // Single utility revenue (tracked once per week; stored at index 0).
    double utility_revenue = utilities_revenue_last_year[0];

    // Annual insurance price = expected combined payout exposure * premium loading.
    // fixed_payouts[wss_ids[0]] is the payout rate DV (same for all WSSes).
    double price = combined_fail_freq * fixed_payouts[wss_ids[0]]
                   * utility_revenue * insurance_premium;

    // Store price and multiplier at index 0 (single-utility convention).
    insurance_price[wss_ids[0]] = price;
    payout_multiplier[0] = (price == 0) ? 0. : 1.;

    // Charge the single utility once.
    DroughtMitigationPolicy::realization_wss[wss_ids[0]]->getOwner()->
            purchaseInsurance(price);
}

// ---------------------------------------------------------------------------
// Static helper: load a CSV with columns:
//   rof_threshold, fail_freq_period_0, fail_freq_period_1, ...
// The header row is skipped. Each data row must have at least 2 columns.
// Returns a vector of rows (each row is a vector<double>).
// ---------------------------------------------------------------------------
vector<vector<double>> InsuranceStorageToROF::loadROFFreqTable(const string &filepath) {
    vector<vector<double>> table;
    ifstream file(filepath);
    if (!file.is_open())
        throw runtime_error(
            "InsuranceStorageToROF: cannot open ROF frequency table: " + filepath);

    string line;
    getline(file, line); // skip header row
    while (getline(file, line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string token;
        vector<double> row;
        while (getline(ss, token, ','))
            row.push_back(stod(token));
        if (row.size() < 2)
            throw runtime_error(
                "InsuranceStorageToROF: each row must have at least 2 columns "
                "(rof_threshold, fail_freq_period_0 [, ...]): " + filepath);
        table.push_back(row);
    }
    if (table.empty())
        throw runtime_error(
            "InsuranceStorageToROF: ROF frequency table contains no data rows: " + filepath);
    return table;
}

// ---------------------------------------------------------------------------
// Determine infrastructure state (0-3) by comparing current WSS capacities
// against the baseline captured when addSystemComponents was called.
//   State 0: no new infra for either WSS
//   State 1: new infra for Descoberto (wss_ids[0]) only
//   State 2: new infra for TortoSM  (wss_ids[1]) only
//   State 3: new infra for both
// ---------------------------------------------------------------------------
int InsuranceStorageToROF::getInfraState() const {
    const double eps = 1e-6;
    auto hasNewInfra = [&](size_t idx) -> bool {
        int u = wss_ids[idx];
        return (DroughtMitigationPolicy::realization_wss[u]->getTotal_treatment_capacity()
                    > initial_treatment_capacity[idx] + eps)
            || (DroughtMitigationPolicy::realization_wss[u]->getTotal_storage_capacity()
                    > initial_storage_capacity[idx] + eps);
    };
    bool d_infra = hasNewInfra(0); // Descoberto
    bool t_infra = hasNewInfra(1); // TortoSM
    if (!d_infra && !t_infra) return 0;
    if ( d_infra && !t_infra) return 1;
    if (!d_infra &&  t_infra) return 2;
    return 3;
}

// ---------------------------------------------------------------------------
// Look up the combined expected annual trigger frequency for the utility in
// infrastructure state `state`, simulation time period `time_period_idx`, and
// ROF trigger value `rof_trigger`.
//
// The table rows are in ASCENDING threshold order, but fail_freq is
// DECREASING (lower threshold crossed more often than higher threshold).
// We scan forward and keep updating the result for every row whose threshold
// is <= rof_trigger; this yields the value at the closest threshold that does
// not exceed the trigger.  If the trigger falls below all table thresholds,
// we return the highest-frequency row (table.front()).
// ---------------------------------------------------------------------------
double InsuranceStorageToROF::lookupFailFreq(int state, int time_period_idx,
                                             double rof_trigger) const {
    const auto &table = rof_freq_tables[state];
    // Column 0 = rof_threshold; columns 1..N = fail_freq for each time period.
    int col = time_period_idx + 1;
    // Default: highest frequency (trigger below all table thresholds).
    double result = table.front()[col];
    for (const auto &row : table) {
        if (row[0] <= rof_trigger + 1e-9)
            result = row[col];  // last update before threshold exceeds trigger
        else
            break;              // thresholds are ascending; no need to continue
    }
    return result;
}

void InsuranceStorageToROF::setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                                           vector<double> &water_sources_rdm, vector<double> &policy_rdm) {
    ContinuityModel::setRealization(realization_id, wss_rdm, water_sources_rdm);

    // Create raw pointer vector for drought mitigation policies (they don't take ownership)
    vector<WaterSupplySystems*> wss_raw;
    for (const auto& wss_ptr : continuity_wss) {
        wss_raw.push_back(wss_ptr.get());
    }
    
    // Pass corresponding wss to drought mitigation instruments.
    for (DroughtMitigationPolicy *dmp : this->drought_mitigation_policies) {
        dmp->addSystemComponents(wss_raw, continuity_water_sources, min_env_flow_controls);
        dmp->setRealization(realization_id, wss_rdm, water_sources_rdm,
                            policy_rdm);
    }
}

/**
 * Runs one the full rof calculations for realization #realization_id for a
 * given week.
 * @param week for which rof is to be calculated.
 * @todo This is mostly a copy and paste from ContinuityModelROF.cpp, so
 * although I'm now on a fix this should be done in a better way at some point.
 */
vector<double> InsuranceStorageToROF::calculateShortTermROFTable(int week) {
    // vector where risks of failure will be stored.
    vector<double> risk_of_failure((unsigned long) n_wss, 0.0);
    double m;
    for (size_t u = 0; u < continuity_wss.size(); ++u) {
        // Get current stored volume for wss u.
        double wss_storage =
                continuity_wss[u]->getTotal_stored_volume();
        // Ratio of current and status-quo wss storage capacities
        //        double m = current_and_base_storage_capacity_ratio[u];
        m = continuity_wss[u]->getTotal_storage_capacity() /
            wss_base_storage_capacity[u];
        // Calculate base table tier that contains the desired ROF by
        // shifting the table around based on new infrastructure -- the
        // shift is made by the part (m - 1) * STORAGE_CAPACITY_RATIO_FAIL *
        // wss_base_storage_capacity[u] - current_storage_table_shift[u]
        double storage_convert = wss_storage +
                                 STORAGE_CAPACITY_RATIO_FAIL * wss_base_storage_capacity[u] *
                                 (1. - m) + current_storage_table_shift[u];
        int tier = min((int) (storage_convert * NO_OF_INSURANCE_STORAGE_TIERS /
                          wss_base_storage_capacity[u]), NO_OF_INSURANCE_STORAGE_TIERS - 1);
        // Mean ROF between the two tiers of the ROF table where
        // current storage is located.
//        risk_of_failure[u] = getRofFromRealizationTable(u, week, tier);
        risk_of_failure[u] = (getRofFromRealizationTable(u, week, tier) +
                getRofFromRealizationTable(u, week, min(NO_OF_INSURANCE_STORAGE_TIERS - 1, tier + 1))) / 2;
    }

    return risk_of_failure;
}

void InsuranceStorageToROF::updateOnlineInfrastructure(int week) {
    ContinuityModelROF::updateOnlineInfrastructure(week);

    // Use capacity multiplier if using pre-computed tables (which are created with such multiplier)
    if (week < WEEKS_IN_YEAR + 1 && use_imported_tables) {
        for (double &u : wss_base_storage_capacity) {
            u *= BASE_STORAGE_CAPACITY_MULTIPLIER;
        }
    } else if (!use_imported_tables) {
        for (int u = 0; u < n_wss; ++u) {
            wss_base_storage_capacity[u] = continuity_wss[u]->getTotal_storage_capacity();
        }
    }
}

