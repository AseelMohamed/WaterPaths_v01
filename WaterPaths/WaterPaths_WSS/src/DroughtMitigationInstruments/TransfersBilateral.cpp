//
// Created by Bernardo on 5/10/2020.
//

#include "../Utils/Utils.h"
#include "TransfersBilateral.h"

TransfersBilateral::TransfersBilateral(const int id,
                                       const vector<double> &pipe_transfer_capacities,
                                       double source_treatment_buffer,
                                       double surcharge_percentage_fee,
                                       const vector<double> &transfer_triggers,
                                       const vector<int> &wss_ids)
        : DroughtMitigationPolicy(id, TRANSFERS_CAESB),
          source_treatment_buffer(source_treatment_buffer),
          surcharge_percentage_fee(surcharge_percentage_fee),
          pipe_transfer_capacities(pipe_transfer_capacities),
          transfer_overhead(surcharge_percentage_fee),
          transfer_triggers(vector<double>(wss_ids.size(), NON_INITIALIZED)) {
        this->transfer_triggers[wss_ids[0]] = transfer_triggers[0];
        // Only set second trigger if we have a second system
        if (wss_ids.size() > 1 && transfer_triggers.size() > 1) {
            this->transfer_triggers[wss_ids[1]] = transfer_triggers[1];
        }
    this->wss_ids = wss_ids;

}

TransfersBilateral::TransfersBilateral(const TransfersBilateral &transfer_caesb):
        DroughtMitigationPolicy(transfer_caesb.id, TRANSFERS_CAESB),
        source_treatment_buffer(transfer_caesb.source_treatment_buffer),
        surcharge_percentage_fee(transfer_caesb.surcharge_percentage_fee),
        pipe_transfer_capacities(transfer_caesb.pipe_transfer_capacities),
        transfer_overhead(transfer_caesb.surcharge_percentage_fee),
        transfer_triggers(transfer_caesb.transfer_triggers) {
    this->wss_ids = transfer_caesb.wss_ids;
}

void TransfersBilateral::applyPolicy(int week) {

    // Add validation for realization_wss pointers before calling performTransfer
    if (realization_wss.size() < 2 || 
        realization_wss[0] == nullptr || 
        realization_wss[1] == nullptr ||
        realization_wss[0]->getOwner() == nullptr ||
        realization_wss[1]->getOwner() == nullptr) {
        transfered_volumes = {0.0, 0.0};
        return;
    }

    // Reset demand_offset to 0 for both WSS every week (will be set if transfer occurs)
    realization_wss[0]->setDemand_offset(0.0, 0.0);
    realization_wss[1]->setDemand_offset(0.0, 0.0);

    double transfer_volume = performTransfer(realization_wss[0],
                                      realization_wss[1],
                                      pipe_transfer_capacities[1], week);

    transfered_volumes = {-transfer_volume, transfer_volume};

    if (transfer_volume == 0) {
        transfer_volume = performTransfer(realization_wss[1],
                                          realization_wss[0],
                                          pipe_transfer_capacities[0], week);
        transfered_volumes = {transfer_volume, -transfer_volume};
    }
}

double TransfersBilateral::performTransfer(WaterSupplySystems *sender, WaterSupplySystems *receiver,
                                           double pumping_capacity,
                                           int week) const {  
    double transfer_volume = 0;
    double receiver_risk = receiver->getRisk_of_failure();
    
    // Add bounds checking for transfer_triggers access
    if (receiver->system_id < 0 || receiver->system_id >= transfer_triggers.size()) {
        printf("ERROR: system_id %d is out of bounds for transfer_triggers (size %zu)\n", 
               receiver->system_id, transfer_triggers.size());
        return 0.0;  // Return early to avoid crash
    }
    
    double trigger = transfer_triggers[receiver->system_id];

    double sender_risk = sender->getRisk_of_failure();

    
    if (receiver_risk > trigger && sender_risk == 0.) {
        double treatment_capacity = sender->getTotal_treatment_capacity();
        double unrestricted_demand = sender->getUnrestrictedDemand();
       
        // Calculate transfer volume
        double available_transfer_volume =
                (treatment_capacity - source_treatment_buffer) * PEAKING_FACTOR
                - unrestricted_demand;
        transfer_volume = max(min(available_transfer_volume,
                              pumping_capacity), 0.);

        // For within-utility transfers: only operational costs, no benefit to sender
        // Receiver pays operational cost (transfer_overhead factor)
        // Sender has zero cost/benefit (no pricing mechanism for internal transfers)
        int price_week = Utils::weekOfTheYear(week);
        double operational_cost_per_volume = sender->getOwner()->waterPrice(price_week) * 
                                             transfer_overhead;
        
        // Receiver: positive offset with operational cost rate
        receiver->setDemand_offset(transfer_volume, operational_cost_per_volume);
        
        // Sender: negative offset with SAME rate as base price (zero net cost/benefit)
        sender->setDemand_offset(-transfer_volume, 
                                 sender->getOwner()->waterPrice(price_week));
    }

    return transfer_volume;
}

void TransfersBilateral::addSystemComponents(vector<WaterSupplySystems *> wss,
                                             vector<WaterSource *> water_sources,
                                             vector<MinEnvFlowControl *> min_env_flow_controls) {
    // Validate that we have at least 2 WaterSupplySystems and they are not null
    if (wss.size() < 2) {
        throw std::invalid_argument("TransfersBilateral::addSystemComponents: Need at least 2 WaterSupplySystems");
    }
    if (wss[0] == nullptr || wss[1] == nullptr) {
        throw std::invalid_argument("TransfersBilateral::addSystemComponents: WaterSupplySystems pointers cannot be null");
    }
    
    realization_wss = {wss[0], wss[1]};
}

void TransfersBilateral::setRealization(unsigned long realization_id,
                                        vector<double> &wss_rdm,
                                        vector<double> &water_sources_rdm,
                                        vector<double> &policy_rdm) {}

const vector<double> &TransfersBilateral::getTransferedVolumes() const {
    return transfered_volumes;
}
