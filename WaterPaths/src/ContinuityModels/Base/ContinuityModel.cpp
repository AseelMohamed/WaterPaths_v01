//
// Created by bernardo on 1/12/17.
//

// #include <iostream>
#include <cmath>
// #include <cstring>
#include <algorithm>
#include "ContinuityModel.h"
// #include "../../SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h"

ContinuityModel::ContinuityModel(vector<WaterSource *> &water_sources, vector<WaterSupplySystems *> &wss,
                                 vector<MinEnvFlowControl *> &min_env_flow_controls,
                                 const Graph &water_sources_graph,
                                 const vector<vector<int>> &water_sources_to_wss,
                                 vector<double> &wss_rdm,
                                 vector<double> &water_sources_rdm,
                                 unsigned long realization_id) :
        continuity_water_sources(water_sources),
        continuity_wss(wss),
        min_env_flow_controls(min_env_flow_controls),
        water_sources_graph(water_sources_graph),
        water_sources_to_wss(water_sources_to_wss),
        sources_topological_order(water_sources_graph.getTopological_order()), /// Get topological order so that mass balance is ran from up to downstream.
        wss_rdm(wss_rdm),
        water_sources_rdm(water_sources_rdm),
        n_wss((int) wss.size()),
        n_sources((int) water_sources.size()),
        realization_id(realization_id)
        {

    //FIXME: THERE IS A STUPID MISTAKE HERE IN THE SORT FUNCTION THAT IS PREVENTING IT FROM WORKING UNDER WINDOWS AND LINUX.
    std::sort(continuity_water_sources.begin(), continuity_water_sources.end(), WaterSource::compare);
    std::sort(continuity_wss.begin(), continuity_wss.end(), 
              [](WaterSupplySystems *a, WaterSupplySystems *b) { return a->system_id < b->system_id; });

    // CRITICAL: Reconnect infrastructure managers after WSS objects are copied
    // When WSS objects are copied, their internal vectors get new addresses but 
    // InfrastructureManager still points to the old addresses
    for (auto* wss_ptr : continuity_wss) {
        wss_ptr->reconnectInfrastructureManager();
    }

    // Link water sources to WSS within each utility by passing pointers of the former to each WSS.
    // printf("DEBUG: water_sources_to_wss.size() = %zu\n", water_sources_to_wss.size());
    for (unsigned long u = 0; u < wss.size(); ++u) {
        // printf("DEBUG: Processing utility %lu, wss[%lu] = %p\n", u, u, continuity_wss[u]);
        // printf("DEBUG: water_sources_to_wss[%lu].size() = %zu\n", u, water_sources_to_wss[u].size());
        for (unsigned long ws = 0; ws < water_sources_to_wss[u].size(); ++ws) {
            auto ws_id = water_sources_to_wss[u][ws];
            // printf("DEBUG: Adding water source %d to utility %lu (wss = %p)\n", ws_id, u, continuity_wss[u]);
            if (ws_id >= continuity_water_sources.size()) {
                string error = "Water source " + to_string(ws_id) + " was not added to list of water sources passed to the continuity model.";
                throw invalid_argument(error);
            }
            WaterSource *water_source =
                    continuity_water_sources.at((unsigned int) ws_id);
            this->continuity_wss[u]->addWaterSource(water_source);
        }
    }

    // Create table showing which utilities draw water from each water source.
    // printf("DEBUG: Creating wss_to_water_sources and water_sources_online_to_wss tables\n");
    // printf("DEBUG: water_sources.size() = %zu, wss.size() = %zu\n", water_sources.size(), wss.size());
    
    wss_to_water_sources.assign(water_sources.size(), vector<int>(0));
    water_sources_online_to_wss.assign(water_sources.size(), vector<int>(0));
    
    for (unsigned long u = 0; u < wss.size(); ++u) {
        // printf("DEBUG: Processing utility %lu for table creation\n", u);
        // printf("DEBUG: water_sources_to_wss[%lu].size() = %zu\n", u, water_sources_to_wss[u].size());
        for (const int &ws : water_sources_to_wss[u]) {
            // printf("DEBUG: Processing water source %d for utility %lu\n", ws, u);
            if (ws >= 0 && ws < static_cast<int>(wss_to_water_sources.size())) {
                wss_to_water_sources[ws].push_back(u);
                if (ws < static_cast<int>(water_sources.size()) && water_sources[ws]->isOnline()) {
                    // printf("DEBUG: Water source %d is online, adding to online list\n", ws);
                    water_sources_online_to_wss[u].push_back(ws);
                }
            } else {
                // printf("DEBUG: ERROR - Water source %d out of bounds (max = %zu)\n", ws, wss_to_water_sources.size());
            }
        }
    }

    // The variables below are to make the storage-ROF table calculation
    // faster by limiting the storage curve shifting to online water sources.
    for (auto water_source : water_sources) {
        bool online = false;

        for (unsigned long u = 0; u < wss.size(); ++u) {
            if (water_source->isOnline())
                online = true;
        }

        if (online) {
            water_sources_capacities.push_back(
                    water_source->getSupplyCapacity());
        } else {
            water_sources_capacities.push_back((double) NONE);
        }
    }

    // Populate vector with wss capacities and check if all wss
    // have storage capacity.
    for (WaterSupplySystems *u : continuity_wss) {
        wss_capacities.push_back(u->getTotal_storage_capacity());
        if (wss_capacities.back() == 0) {
            string error = "Water Supply System " + to_string(u->system_id) + " has no storage capacity (0 MGD), "
            "which would lead to an ROF value of 0/0. WaterPaths currently requires wss to have non-zero storage capacity";
            throw invalid_argument(error);
        }
    }

    // Populate vector indicating the downstream source from each source.
    for (vector<int> ds : water_sources_graph.getDownSources()) {
        if (ds.empty()) {
            downstream_sources.push_back(NON_INITIALIZED);
        } else {
            downstream_sources.push_back(ds[0]);
        }
    }

    // Add reference to water sources and wss so that controls can
    // access their info.
    for (MinEnvFlowControl *mef : this->min_env_flow_controls) {
        mef->addComponents(water_sources, wss);
    }

    // Set realization id on wss and water sources, so that they use the
    // right streamflow, evaporation and demand data.
    setRealization(realization_id, wss_rdm, water_sources_rdm);

    demands = std::vector<vector<double>>(
            continuity_water_sources.size(),
            vector<double>(continuity_wss.size(), 0.));
    
    // populate array delta_realization_weeks so that the rounding and casting don't
    // have to be done every time continuityStep is called, avoiding a bottleneck.
    // Variable delta_realization_weeks[0] is to be zero representing a continuity
    // step for a non-ROF simulation: the realization itself.
    delta_realization_weeks[0] = 0;
    for (int r = 0; r < NUMBER_REALIZATIONS_ROF; ++r) {
        delta_realization_weeks[r + 1] = (int) std::round((r + 1) * WEEKS_IN_YEAR);
    }
}

ContinuityModel::~ContinuityModel() {
    /// Delete water sources (with null check for shared objects)
    for (auto ws : continuity_water_sources) {
        if (ws != nullptr) {
            delete ws;
        }
    }

    /// Delete wss (with null check for shared objects)
    for (auto u : continuity_wss) {
        if (u != nullptr) {
            delete u;
        }
    }

    /// Delete min env flow controls (with null check for shared objects)
    for (auto mef : min_env_flow_controls){
        if (mef != nullptr) {
            delete mef;
        }
    }
}

/**
 * Calculates continuity for one week time step for streamflows of id_rof years
 * before current week.
 * @param week current week.
 * @param rof_realization rof realization id (between 0 and 49 inclusive).
 */
void ContinuityModel::continuityStep(
        int week, int rof_realization, bool apply_demand_buffer) {
    // For ROF calculations, we want to use PAST year data.
    // For rof_realization = 0, use current week (no offset)
    // For rof_realization = 1, use week - 52 (1 year ago)
    // For rof_realization = 2, use week - 104 (2 years ago), etc.
    int shifted_week;
    if (rof_realization == NON_INITIALIZED) {
        shifted_week = week; // Normal simulation, no offset
    } else {
        // ROF simulation: use past data (rof_realization years ago)
        shifted_week = week - (rof_realization * WEEKS_IN_YEAR_ROUND);
        if (shifted_week < 0) {
            // If we don't have enough historical data, skip this step
            return;
        }
    }

    double* upstream_spillage = new double[n_sources];
    fill_n(upstream_spillage, n_sources, 0.);
    double* wastewater_discharges = new double[n_sources];
    fill_n(wastewater_discharges, n_sources, 0.);

    // if ROF calculations use previous year unless this is the first year (to
    // simplify data input, since the first year of the simulation is probably
    // going to be as dire are latter years.
    int week_demand_rof_shift = 
	    (rof_realization != NON_INITIALIZED && 
	     week > WEEKS_IN_YEAR_ROUND ? WEEKS_IN_YEAR_ROUND : NONE);
    int week_demand = week - week_demand_rof_shift;

    /**
     * Get wastewater discharges based on previous week's demand.
     *
     * Split weekly demands among each reservoir for each utility. For each
     * water source: (1) sums the demands of each drawing utility to come up
     * with the total unrestricted_demand for that week for that water
     * source, and (2) sums the flow contributions of upstream reservoirs.
     */
    for (WaterSupplySystems *u : continuity_wss) {
        u->calculateWastewater_releases(week_demand, wastewater_discharges);
        u->splitDemands(week_demand, demands, apply_demand_buffer);
    }

    /**
     * Set minimum environmental flows for water sources based on their
     * individual controls.
     */
    for (MinEnvFlowControl *c : min_env_flow_controls) {
        continuity_water_sources[c->water_source_id]->
                setMin_environmental_outflow(c->getRelease(week));
    }

    /**
     * For all water sources, performs mass balance to update the available
     * volume. The week here is shifted back according to the rof year
     * realization (0 to 49, inclusive) so that the right flows are gotten
     * from source catchments for each rof year realization. If this is not an
     * rof calculation but an actual simulation instead, rof_realization will
     * be equal to -1 (see header file) so that there is no week shift.
     */
    auto& upstream_sources_ids = water_sources_graph.getUpstream_sources();
    for (int i : sources_topological_order) {
        // Sum spillage from all sources upstream source i.
        for (int ws : upstream_sources_ids[i]) {
            upstream_spillage[i] += continuity_water_sources.at(
                            static_cast<unsigned long>(ws))->getTotal_outflow();
        }

        // Mass balance. For ROF calculations, we want to use PAST year data.
        // For rof_realization = 0, use current week (no offset)
        // For rof_realization = 1, use week - 52 (1 year ago)
        // For rof_realization = 2, use week - 104 (2 years ago), etc.
        int week_with_rof_offset;
        if (rof_realization == NON_INITIALIZED) {
            week_with_rof_offset = week; // Normal simulation, no offset
        } else {
            // ROF simulation: use past data (rof_realization years ago)
            week_with_rof_offset = week - (rof_realization * WEEKS_IN_YEAR_ROUND);
            if (week_with_rof_offset < 0) {
                // If we don't have enough historical data, skip this step
                continue;
            }
        }
        continuity_water_sources[i]->continuityWaterSource(
                week_with_rof_offset,
                upstream_spillage[i], wastewater_discharges[i], demands[i]);
        demands[i] = vector<double>(n_wss, 0.);
    }

    // updates combined storage for Water Supply Systems.
    for (WaterSupplySystems *u : continuity_wss) {
        u->updateTotalAvailableVolume();
    }

    delete[] upstream_spillage;
    delete[] wastewater_discharges;
}

void ContinuityModel::setRealization(unsigned long realization_id, vector<double> &wss_rdm,
                                     vector<double> &water_sources_rdm) {
    if (realization_id != (unsigned) NON_INITIALIZED) {
        for (WaterSupplySystems *u : continuity_wss)
            u->setRealization(realization_id, wss_rdm);
        for (WaterSource *ws : continuity_water_sources)
            ws->setRealization(realization_id, water_sources_rdm);
        for (MinEnvFlowControl *mef : min_env_flow_controls)
            mef->setRealization(realization_id, water_sources_rdm);
    }
}

/*
 * Get list of next online downstream sources for each source.
 * @return vector with downstream water sources.
 */
vector<int> ContinuityModel::getOnlineDownstreamSources() {
    vector<int> online_downstream_sources(n_sources, NON_INITIALIZED);
    for (int ws : sources_topological_order) {
        int downstream_source_online = ws;
        do {
            downstream_source_online = downstream_sources[downstream_source_online];
        } while (downstream_source_online != NON_INITIALIZED &&
                !continuity_water_sources[downstream_source_online]->isOnline() &&
                continuity_water_sources[ws]->getSupplyCapacity() == 0);
        online_downstream_sources[ws] = downstream_source_online;
    }

    return online_downstream_sources;
}

const vector<WaterSource *> &ContinuityModel::getContinuity_water_sources() const {
    return continuity_water_sources;
}

const vector<WaterSupplySystems *> &ContinuityModel::getContinuity_wss() const {
    return continuity_wss;
}
