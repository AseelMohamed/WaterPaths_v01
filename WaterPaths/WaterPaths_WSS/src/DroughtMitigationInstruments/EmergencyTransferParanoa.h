//
// Paranoa emergency-only supply to TortoSM.
// When TortoSM's ROF exceeds rof_trigger, water is pumped from Paranoa to
// TortoSM at up to the effective pipe capacity (hm³/week). The effective
// capacity starts at base_pipe_capacity and grows by pipe_capacity_per_expansion
// for each triggered infrastructure in pipe_expansion_ids that is online.
// The total volume is also bounded by Paranoa's allocated supply volume for WSS 1.
//

#ifndef WATERPATHS_EMERGENCYTRANSFERPARANOA_H
#define WATERPATHS_EMERGENCYTRANSFERPARANOA_H

#include <vector>
#include "Base/DroughtMitigationPolicy.h"

class EmergencyTransferParanoa : public DroughtMitigationPolicy {
    const int receiver_wss_id;                  // TortoSM system_id = 1
    const int paranoa_source_id;                // Paranoa water source id = 3
    const double base_pipe_capacity;            // Initial pipeline capacity (hm³/week)
    const double pipe_capacity_per_expansion;   // Added per built expansion (hm³/week)
    const vector<int> pipe_expansion_ids;       // Water source IDs that expand the pipeline
    double rof_trigger;                         // ROF threshold to activate emergency supply
    const double transfer_cost_multiplier;      // Cost per unit volume = receiver's water price × this factor

    WaterSupplySystems* receiver_wss = nullptr;
    WaterSource* paranoa_source = nullptr;
    vector<WaterSource*> pipe_expansion_sources; // populated in addSystemComponents
    double transferred_volume = 0.0;
    bool ever_triggered_ = false;  // set true the first time the transfer actually fires

public:
    EmergencyTransferParanoa(int id,
                             int receiver_wss_id,
                             int paranoa_source_id,
                             double base_pipe_capacity,
                             double rof_trigger,
                             double transfer_cost_multiplier = 1.1,
                             vector<int> pipe_expansion_ids = {},
                             double pipe_capacity_per_expansion = 0.0);

    EmergencyTransferParanoa(const EmergencyTransferParanoa& other);

    void applyPolicy(int week) override;

    void addSystemComponents(vector<WaterSupplySystems*> wss,
                             vector<WaterSource*> water_sources,
                             vector<MinEnvFlowControl*> min_env_flow_controls) override;

    void setRealization(unsigned long realization_id,
                        vector<double>& wss_rdm,
                        vector<double>& water_sources_rdm,
                        vector<double>& policy_rdm) override;

    double getTransferredVolume() const;

    // Returns true if the emergency transfer has fired at least once this realization.
    // Used by WaterSupplySystems to gate Paranoa infrastructure construction.
    bool hasEverTriggered() const;
    const bool* getEverTriggeredPtr() const;
};

#endif //WATERPATHS_EMERGENCYTRANSFERPARANOA_H
