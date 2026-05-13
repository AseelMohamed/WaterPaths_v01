//
// Paranoa emergency-only supply to TortoSM.
// When TortoSM's ROF exceeds rof_trigger, water is pumped from Paranoa to
// TortoSM at up to pipe_capacity hm³/week. The total volume is bounded by
// Paranoa's allocated supply volume for WSS 1 (0.8% of total capacity).
//

#ifndef WATERPATHS_EMERGENCYTRANSFERPARANOA_H
#define WATERPATHS_EMERGENCYTRANSFERPARANOA_H

#include "Base/DroughtMitigationPolicy.h"

class EmergencyTransferParanoa : public DroughtMitigationPolicy {
    const int receiver_wss_id;            // TortoSM system_id = 1
    const int paranoa_source_id;          // Paranoa water source id = 3
    const double pipe_capacity;           // Max weekly transfer volume (hm³/week)
    double rof_trigger;                   // ROF threshold to activate emergency supply
    const double transfer_cost_multiplier; // Cost per unit volume = receiver's water price × this factor

    WaterSupplySystems* receiver_wss = nullptr;
    WaterSource* paranoa_source = nullptr;
    double transferred_volume = 0.0;

public:
    EmergencyTransferParanoa(int id,
                             int receiver_wss_id,
                             int paranoa_source_id,
                             double pipe_capacity,
                             double rof_trigger,
                             double transfer_cost_multiplier = 1.1);

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
};

#endif //WATERPATHS_EMERGENCYTRANSFERPARANOA_H
