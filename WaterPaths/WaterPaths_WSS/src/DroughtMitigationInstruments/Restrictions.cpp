//
// Created by bernardo on 2/3/17.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include "Restrictions.h"
#include "../Utils/Utils.h"

/**
 * Restriction policy.
 * 
 * IMPORTANT BEHAVIOR:
 * - Restrictions are TRIGGERED by a specific WSS's Risk of Failure (ROF)
 * - But once triggered, restrictions are APPLIED to ALL WSS in the owner Utility
 * - This ensures utility-wide operational consistency during drought response
 * - Financial impacts (price surcharges) are applied at the Utility level
 * 
 * @param id WSS ID that this restriction policy monitors
 * @param stage_multipliers Demand multipliers for each restriction stage (e.g., 0.8 = 20% reduction)
 * @param stage_triggers ROF thresholds that trigger each restriction stage
 * @todo set lower ROF threshold for utilities to lift restrictions.
 * @todo implement drought surcharges.
 */
Restrictions::Restrictions(const int id, const vector<double> &stage_multipliers,
                           const vector<double> &stage_triggers)
        : DroughtMitigationPolicy(id, RESTRICTIONS),
          stage_multipliers(stage_multipliers),
          stage_triggers(stage_triggers) {
    wss_ids = vector<int>(1, id);
}

Restrictions::Restrictions(
        const int id, const vector<double> &stage_multipliers,
        const vector<double> &stage_triggers,
        const vector<vector<double>> *typesMonthlyDemandFraction,
        const vector<vector<double>> *typesMonthlyWaterPrice,
        const vector<vector<double>> *priceMultipliers)
        : DroughtMitigationPolicy(id, RESTRICTIONS),
          stage_multipliers(stage_multipliers),
          stage_triggers(stage_triggers) {
    calculateWeeklyAverageWaterPrices(typesMonthlyDemandFraction,
                                      typesMonthlyWaterPrice,
                                      priceMultipliers);
    wss_ids = vector<int>(1, id);
}

Restrictions::Restrictions(const Restrictions &restrictions) :
        DroughtMitigationPolicy(restrictions),
        stage_multipliers(restrictions.stage_multipliers),
        stage_triggers(restrictions.stage_triggers),
        restricted_weekly_average_volumetric_price(
                restrictions.restricted_weekly_average_volumetric_price) {
    wss_ids = restrictions.wss_ids;
}

Restrictions::~Restrictions() {}

void Restrictions::applyPolicy(int week) {

    current_multiplier = 1.0;
    unsigned long stage = 0;
    /// Loop through restriction stage rof triggers to see which stage should be applied, if any.
    for (unsigned long i = 0; i < stage_triggers.size(); ++i) {
        if (realization_wss[0]->getRisk_of_failure() > stage_triggers[i]) {
            /// Demand multiplier to be applied, based on the rofs.
            current_multiplier = stage_multipliers[i];
            stage = i + 1;
        } else
            break;
    }
    
    // SAFETY: Verify we have a valid triggering WSS
    if (realization_wss.empty() || realization_wss[0] == nullptr) {
        cerr << "ERROR: Restrictions policy has no valid WSS assigned. Cannot apply policy." << endl;
        return;
    }
    
    // Get the owner utility to apply restrictions utility-wide
    Utility* owner_utility = realization_wss[0]->getOwner();
    if (owner_utility == nullptr) {
        cerr << "ERROR: WSS has no owner utility. Cannot apply restrictions." << endl;
        return;
    }
    
    // Apply operational restrictions to ALL WSS in the owner utility
    // This ensures system-wide consistency when restrictions are triggered
    const vector<WaterSupplySystems*>& all_wss = owner_utility->getWSSReferences();
    if (!all_wss.empty()) {
        // Use WSS references if available (realization-specific copies)
        for (size_t i = 0; i < all_wss.size(); ++i) {
            WaterSupplySystems* wss = all_wss[i];
            if (wss != nullptr) {
                wss->setDemand_multiplier(current_multiplier);
            }
        }
    } else {
        // Fallback to owned WSS if references not set
        const auto& owned_wss = owner_utility->getWaterSupplySystems();
        for (size_t i = 0; i < owned_wss.size(); ++i) {
            const auto& wss = owned_wss[i];
            if (wss != nullptr) {
                wss->setDemand_multiplier(current_multiplier);
            }
        }
    }
    
    // Apply financial restriction (price surcharge) at utility level
    if (!restricted_weekly_average_volumetric_price.empty() && stage > 0) {
        int week_of_year = Utils::weekOfTheYear(week);
        owner_utility->setRestricted_price(
                restricted_weekly_average_volumetric_price[stage - 1][week_of_year]);
    }
}

double Restrictions::getCurrent_multiplier() const {
    return current_multiplier;
}

void Restrictions::addSystemComponents(vector<WaterSupplySystems *> systems_wss,
                                       vector<WaterSource *> water_sources,
                                       vector<MinEnvFlowControl *> min_env_flow_controls) {
    /// Get WSS whose IDs correspond to restriction policy ID.
    for (WaterSupplySystems *w : systems_wss) {
        if (w->system_id == id) {
            if (!realization_wss.empty())
                throw std::logic_error("This restriction policy already has a WSS assigned to it.");
            /// Link WSS to policy.
            this->realization_wss.push_back(w);
        }
    }
    if (systems_wss.empty())
        throw std::invalid_argument("Restriction policy ID must match WSS's ID.");
}

/**
 * Calculates average water price from consumer types and respective prices.
 * @param typesMonthlyDemandFraction
 * @param typesMonthlyWaterPrice
 */
void Restrictions::calculateWeeklyAverageWaterPrices(
        const vector<vector<double>> *typesMonthlyDemandFraction,
        const vector<vector<double>> *typesMonthlyWaterPrice, const
        vector<vector<double>> *priceMultipliers) {

    if (priceMultipliers) {
        // Validate price multipliers to ensure they're >= 1.0 (prices shouldn't decrease during restrictions)
        for (unsigned long s = 0; s < priceMultipliers->size(); ++s) {
            for (unsigned long t = 0; t < priceMultipliers->at(s).size(); ++t) {
                if ((*priceMultipliers)[s][t] < 1.0) {
                    cerr << "WARNING: Price multiplier < 1.0 detected in restriction stage " << s 
                         << ", tier " << t << ": " << (*priceMultipliers)[s][t] << endl;
                    cerr << "This may cause restricted prices to be lower than unrestricted prices." << endl;
                }
            }
        }

        unsigned long n_tiers = typesMonthlyWaterPrice->at(0).size();
        restricted_weekly_average_volumetric_price =
                std::vector<std::vector<double>>();

        for (unsigned long s = 0; s < priceMultipliers->size(); ++s) { // stages loop
            restricted_weekly_average_volumetric_price.emplace_back(
                    (int) WEEKS_IN_YEAR + 1, 0.);
            double monthly_average_price[NUMBER_OF_MONTHS] = {};

            for (int m = 0; m < NUMBER_OF_MONTHS; ++m) { // monthly loop
                for (unsigned long t = 0; t < n_tiers; ++t) { // consumer type loop
                    monthly_average_price[m] +=
                            (*typesMonthlyDemandFraction)[m][t] *
                            (*typesMonthlyWaterPrice)[m][t] *
                            (*priceMultipliers)[s][t];
                }
            }

            for (int w = 0; w < (int) (WEEKS_IN_YEAR + 1); ++w) {
                // Fix buffer overflow: ensure month index doesn't exceed array bounds
                int month_index = min((int)(w / WEEKS_IN_MONTH), NUMBER_OF_MONTHS - 1);
                restricted_weekly_average_volumetric_price[s][w] =
                        monthly_average_price[month_index] / WEEKS_IN_MONTH;
            }
        }
    }
}

void Restrictions::setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                                  vector<double> &water_sources_rdm, vector<double> &policy_rdm) {
    for (double& sm : stage_multipliers) {
        sm = min(1., (1 - (1 - sm) * policy_rdm.at((unsigned long) id)));
    }
}
