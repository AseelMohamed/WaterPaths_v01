//
// Created by bct52 on 6/28/17.
//

#include "MinEnvFlowControl.h"

MinEnvFlowControl::MinEnvFlowControl(int water_source_id,
                                             const vector<int> &water_sources_ids,
                                             const vector<int> &aux_WSS_ids, int type)
        : 
          water_sources_ids(water_sources_ids),
          WSS_ids(aux_WSS_ids),
	  water_source_id(water_source_id),
          type(type) {}

MinEnvFlowControl::MinEnvFlowControl(
        const MinEnvFlowControl &min_env_control) :
        water_sources(vector<WaterSource *>()),
        WSS(vector<WaterSupplySystems *>()),
        water_sources_ids(min_env_control.water_sources_ids),
        WSS_ids(min_env_control.WSS_ids),
        water_source_id(min_env_control.water_source_id),
        type(min_env_control.type) {}


void MinEnvFlowControl::addComponents(
        vector<WaterSource *> water_sources, vector<WaterSupplySystems *> WSS) {
    this->water_sources = vector<WaterSource *>(water_sources.size());

    for (int i : water_sources_ids) {
        this->water_sources[i] = water_sources[i];
    }

    this->WSS = vector<WaterSupplySystems *>(WSS.size());

    for (int i : WSS_ids) {
        this->WSS[i] = WSS[i];
    }
}

void MinEnvFlowControl::setRealization(unsigned long r, vector<double> &rdm_factors) {

}

MinEnvFlowControl::~MinEnvFlowControl() = default;
