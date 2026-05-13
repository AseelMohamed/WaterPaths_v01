//
// Option B implementation: Paranoa emergency supply to TortoSM.
//
// applyPolicy() is called before splitDemands() each week.
// When triggered:
//   1. paranoa->removeWater(1, vol) immediately decrements Paranoa's
//      available_allocated_volumes[1] and available_volume. AllocatedReservoir
//      continuity then processes the already-decremented volume correctly.
//   2. receiver_wss->setDemand_offset(vol, 0.0) reduces the demand that
//      TortoSM must meet from its normal sources (TortoSM + Bananal/Torto).
//

#include "EmergencyTransferParanoa.h"
#include "../Utils/Utils.h"

EmergencyTransferParanoa::EmergencyTransferParanoa(int id,
                                                   int receiver_wss_id,
                                                   int paranoa_source_id,
                                                   double pipe_capacity,
                                                   double rof_trigger,
                                                   double transfer_cost_multiplier)
        : DroughtMitigationPolicy(id, EMERGENCY_TRANSFER_PARANOA),
          receiver_wss_id(receiver_wss_id),
          paranoa_source_id(paranoa_source_id),
          pipe_capacity(pipe_capacity),
          rof_trigger(rof_trigger),
          transfer_cost_multiplier(transfer_cost_multiplier) {
    this->wss_ids = {receiver_wss_id};
}

EmergencyTransferParanoa::EmergencyTransferParanoa(const EmergencyTransferParanoa& other)
        : DroughtMitigationPolicy(other.id, EMERGENCY_TRANSFER_PARANOA),
          receiver_wss_id(other.receiver_wss_id),
          paranoa_source_id(other.paranoa_source_id),
          pipe_capacity(other.pipe_capacity),
          rof_trigger(other.rof_trigger),
          transfer_cost_multiplier(other.transfer_cost_multiplier),
          transferred_volume(0.0) {
    this->wss_ids = {receiver_wss_id};
}

void EmergencyTransferParanoa::applyPolicy(int week) {
    if (receiver_wss == nullptr || paranoa_source == nullptr) {
        transferred_volume = 0.0;
        return;
    }

    double receiver_rof = receiver_wss->getRisk_of_failure();

    if (receiver_rof > rof_trigger) {
        // Available supply allocation for TortoSM (WSS 1) in Paranoa
        double available = paranoa_source->getAvailableAllocatedVolume(receiver_wss_id);
        double transfer_vol = max(min(pipe_capacity, available), 0.0);

        if (transfer_vol > 0.0) {
            // Deduct from Paranoa's allocated volume for WSS 1.
            // AllocatedReservoir::removeWater decrements available_allocated_volumes[1]
            // and available_volume; the policy_added_demand mechanism ensures
            // continuity picks up the deduction correctly.
            paranoa_source->removeWater(receiver_wss_id, transfer_vol);

            // Inject supply into TortoSM via demand offset.
            // Cost = receiver's current water price × transfer_cost_multiplier,
            // matching the surcharge convention used by TransfersBilateral.
            int price_week = Utils::weekOfTheYear(week);
            double cost_rate = receiver_wss->getOwner()->waterPrice(price_week)
                               * transfer_cost_multiplier;
            receiver_wss->setDemand_offset(transfer_vol, cost_rate);
        }
        transferred_volume = transfer_vol;
    } else {
        transferred_volume = 0.0;
    }
}

void EmergencyTransferParanoa::addSystemComponents(
        vector<WaterSupplySystems*> wss,
        vector<WaterSource*> water_sources,
        vector<MinEnvFlowControl*> min_env_flow_controls) {

    for (auto* w : wss) {
        if (w != nullptr && w->system_id == receiver_wss_id) {
            receiver_wss = w;
        }
    }

    if (paranoa_source_id >= 0 && paranoa_source_id < (int)water_sources.size()) {
        paranoa_source = water_sources[paranoa_source_id];
    }
}

void EmergencyTransferParanoa::setRealization(unsigned long realization_id,
                                               vector<double>& wss_rdm,
                                               vector<double>& water_sources_rdm,
                                               vector<double>& policy_rdm) {
    // Reset only the per-realization counter. Do NOT null the pointers:
    // addSystemComponents is called before setRealization in the constructor,
    // so nullifying here would permanently break the policy for every realization.
    transferred_volume = 0.0;
}

double EmergencyTransferParanoa::getTransferredVolume() const {
    return transferred_volume;
}
