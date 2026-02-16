//
// Created by bernardo on 1/26/17.
//

// #include <iostream>
#include <algorithm>
#include <set>
#include "ContinuityModelRealization.h"

ContinuityModelRealization::ContinuityModelRealization(
        vector<WaterSource *> &water_sources,
        const Graph &water_sources_graph,
        const vector<vector<int>> &water_sources_to_wss,
        vector<std::unique_ptr<WaterSupplySystems>> &&wss,
        const vector<DroughtMitigationPolicy *> &drought_mitigation_policies,
        vector<MinEnvFlowControl *> &min_env_flow_control,
        vector<double>& wss_rdm,
        vector<double>& water_sources_rdm,
        vector<double>& policy_rdm,
        const unsigned int realization_id)
        : ContinuityModel(water_sources, std::move(wss), min_env_flow_control, water_sources_graph,
                          water_sources_to_wss, wss_rdm, water_sources_rdm,
                          realization_id),
          drought_mitigation_policies(drought_mitigation_policies) {

    // Create raw pointer vector for drought mitigation instruments (they don't take ownership)
    vector<WaterSupplySystems*> wss_raw;
    for (const auto& wss_ptr : continuity_wss) {
        wss_raw.push_back(wss_ptr.get());
    }
    
    // Pass corresponding wss to drought mitigation instruments.
    for (DroughtMitigationPolicy *dmp : this->drought_mitigation_policies) {
        dmp->addSystemComponents(wss_raw, water_sources, min_env_flow_control);
        dmp->setRealization(realization_id, wss_rdm, water_sources_rdm,
                            policy_rdm);
    }
}

ContinuityModelRealization::~ContinuityModelRealization() {
    // Delete drought mitigation policies.
    for (auto dmp : drought_mitigation_policies) {
        delete dmp;
    }
}

void ContinuityModelRealization::setShortTermROFs(const vector<double> &risks_of_failure) {
    for (unsigned long i = 0; i < continuity_wss.size(); ++i) {
        continuity_wss.at(i)->setRisk_of_failure(risks_of_failure.at(i));
    }
}

/**
 * Set long-term risk of failure for each WSS and trigger infrastructure decisions.
 * 
 * CORRECTED ARCHITECTURE:
 * - Each WSS evaluates its individual ROF and makes infrastructure decisions
 * - WSS-level decisions are based on specific local conditions, not aggregated utility-level data
 * - Utilities are notified to implement (finance, construct) the infrastructure decisions
 * - This reflects how real water systems operate: individual systems identify needs, 
 *   utilities implement and finance the solutions
 * 
 * @param risks_of_failure Vector of ROF values, one per WSS
 * @param week Current simulation week
 */
void ContinuityModelRealization::setLongTermROFs(const vector<double> &risks_of_failure, const int week) {
    vector<int> new_infra_triggered;
    int nit; // new infrastruction triggered - id.

    // Loop over utilities to see if any of them will build new infrastructure.
    for (unsigned long u = 0; u < continuity_wss.size(); ++u) {
        // Runs utility's infrastructure construction handler and get the id
        // of new source built, if any.
        nit = continuity_wss[u]->
                infrastructureConstructionHandler(risks_of_failure[u], week);
        // If new source was built, check add it to the list of sources
        // built by all utilities.
        if (nit != NON_INITIALIZED) {
            new_infra_triggered.push_back(nit);
            // Bond issuance is now handled internally by WSS infrastructureConstructionHandler
            // No separate bond issuance call needed here
        }
    }

    // Look for and remove duplicates, in the unlikely case two wss
    // build the same source at the same time. This will prevent the source
    // from being erased from a utility which will later try to build it.
    sort(new_infra_triggered.begin(),
         new_infra_triggered.end());
    new_infra_triggered.erase(unique(new_infra_triggered.begin(),
                                     new_infra_triggered.end()),
                              new_infra_triggered.end());

    // If infrastructure was built, force wss to build their share of
    // that infrastructure option (which will only happen it the listed
    // option is in the list of sources to be built for other wss.
    if (!new_infra_triggered.empty())
        for (auto& w : continuity_wss) {
            w->getOwner()->forceInfrastructureConstruction(week, new_infra_triggered);
        }
}

void ContinuityModelRealization::applyDroughtMitigationPolicies(int week) {
    // Reset drought mitigation variables at the start of the week
    // so current prices remain available for data collection after continuityStep.
    std::set<Utility*> utilities_reset;
    for (auto& wss : continuity_wss) {
        if (wss != nullptr) {
            Utility* owner = wss->getOwner();
            if (owner != nullptr && utilities_reset.insert(owner).second) {
                owner->resetDroughtMitigationVariables();
            }
        }
    }

    // Reset all WSS demand multipliers to 1.0 before applying policies
    // This ensures restrictions don't persist from previous weeks unless re-triggered
    // match Original model behavior and test impact on transfer patterns
    for (auto& wss : continuity_wss) {
        if (wss != nullptr) {
            wss->setDemand_multiplier(1.0);
        }
    }
    
    // Apply all drought mitigation policies
    for (DroughtMitigationPolicy* dmp : drought_mitigation_policies) {
        dmp->applyPolicy(week);
    }
}

const vector<DroughtMitigationPolicy *> ContinuityModelRealization::getDrought_mitigation_policies() const {
    return drought_mitigation_policies;
}
