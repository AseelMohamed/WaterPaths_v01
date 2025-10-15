//
// Created by AI Assistant on 10/13/25.
//

#ifndef TRIANGLEMODEL_WSSDATACOLLECTOR_H
#define TRIANGLEMODEL_WSSDATACOLLECTOR_H

#include "Base/DataCollector.h"
#include "../SystemComponents/Utility/WaterSupplySystems.h"

class WSSDataCollector : public DataCollector {
private:
    vector<double> combined_storage;
    vector<double> storage_capacity;
    vector<double> storage_to_capacity_ratio;
    vector<double> unrestricted_demand;
    vector<double> restricted_demand;
    vector<double> demand_multiplier;
    vector<double> demand_offset;
    vector<double> waste_water_discharge;
    vector<double> unfulfilled_demand;
    vector<double> net_stream_inflow;
    vector<double> total_treatment_capacity;
    vector<double> total_storage_treatment_capacity;
    vector<int> water_sources_count;
    const WaterSupplySystems *wss;

public:

    explicit WSSDataCollector(const WaterSupplySystems *wss, unsigned long realization);

    WSSDataCollector &operator=(const WSSDataCollector &wss_data_collector);

    string printTabularString(int week) override;

    string printCompactString(int week) override;

    void collect_data() override;

    string printTabularStringHeaderLine1() override;

    string printTabularStringHeaderLine2() override;

    string printCompactStringHeader() override;

    void checkForNans() const;

    // Getters
    const vector<double> &getCombined_storage() const;
    const vector<double> &getStorage_capacity() const;
    const vector<double> &getStorage_to_capacity_ratio() const;
    const vector<double> &getUnrestricted_demand() const;
    const vector<double> &getRestricted_demand() const;
    const vector<double> &getDemand_multiplier() const;
    const vector<double> &getDemand_offset() const;
    const vector<double> &getWaste_water_discharge() const;
    const vector<double> &getUnfulfilled_demand() const;
    const vector<double> &getNet_stream_inflow() const;
    const vector<double> &getTotal_treatment_capacity() const;
    const vector<double> &getTotal_storage_capacity() const;
    const vector<int> &getWater_sources_count() const;
    const WaterSupplySystems *getWss() const;
};

#endif //TRIANGLEMODEL_WSSDATACOLLECTOR_H