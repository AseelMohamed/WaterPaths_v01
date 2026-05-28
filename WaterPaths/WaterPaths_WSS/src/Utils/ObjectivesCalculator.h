//
// Created by bernardoct on 8/25/17.
//

#ifndef TRIANGLEMODEL_OBJECTIVESCALCULATOR_H
#define TRIANGLEMODEL_OBJECTIVESCALCULATOR_H


#include "../DataCollector/UtilitiesDataCollector.h"
#include "../DataCollector/RestrictionsDataCollector.h"
#include "../DataCollector/WSSDataCollector.h"

class ObjectivesCalculator {

public:
    static double calculateReliabilityObjective(
            const vector<UtilitiesDataCollector *>& utility_collector,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateRestrictionFrequencyObjective(
            const vector<RestrictionsDataCollector *>& restriction_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateRestrictionFrequencyObjective_WSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateNetPresentCostInfrastructureObjective(
            const vector<UtilitiesDataCollector *>& utility_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculatePeakFinancialCostsObjective(
            const vector<UtilitiesDataCollector *>& utility_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateWorseCaseCostsObjective(
            const vector<UtilitiesDataCollector *>& utility_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    // WSS-level objective calculation functions
    static double calculateNetPresentCostInfrastructureObjective_WSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            const vector<UtilitiesDataCollector *>& utility_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));
    
    static double calculateReliabilityObjective_WSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateAffordabilityIndexObjective_WSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    // New functions with configurable aggregation methods
    static double calculateReliabilityObjective_WSS_Configurable(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations,
            int aggregation_method); // 0=MIN, 1=AVERAGE

    static double calculateAffordabilityIndexObjective_WSS_Configurable(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations,
            int aggregation_method); // 0=MAX, 1=AVERAGE

    static double calculateFailureSeverityObjective_WSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static double calculateFailureSeverityObjective_WSS_Configurable(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations,
            int aggregation_method); // 0=MAX (worst case), 1=AVERAGE

    // Per-WSS functions for experiment mode 5 (no aggregation across WSS)
    static vector<double> calculateReliabilityObjective_WSS_PerWSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

    static vector<double> calculateAffordabilityIndexObjective_WSS_PerWSS(
            const vector<vector<WSSDataCollector *>>& wss_data,
            vector<unsigned long> realizations = vector<unsigned long>(0));

};


#endif //TRIANGLEMODEL_OBJECTIVESCALCULATOR_H
