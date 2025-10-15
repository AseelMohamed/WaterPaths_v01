//
// Created by bernardo on 1/13/17.
//

#ifndef TRIANGLEMODEL_UTILITY_H
#define TRIANGLEMODEL_UTILITY_H


#include <map>
#include <memory>
#include <vector>
#include "../WaterSources/Reservoir.h"
#include "../../Utils/Constants.h"
#include "../../Controls/WwtpDischargeRule.h"
#include "InfrastructureManager.h"
//#include "../Utils/Matrix3D.h"

class WaterSupplySystems; // forward declaration

class Utility {
private:
    vector<double> weekly_average_volumetric_price;
    vector<int> priority_draw_water_source;
    vector<int> non_priority_draw_water_source;
    vector<double> weekly_peaking_factor;
    // Financial and strategic variables (kept in Utility)
    double gross_revenue = 0;
    
    // Operational variables moved to WaterSupplySystems:
    // short_term_risk_of_failure, long_term_risk_of_failure, 
    // total_storage_capacity, total_available_volume, total_stored_volume,
    // total_treatment_capacity, total_storage_treatment_capacity,
    // waste_water_discharge, unfulfilled_demand, net_stream_inflow
    double price_rdm_multiplier = 1.;
    // removed available_treated_flow_rate - now owned by WaterSupplySystems
    bool used_for_realization = true;
    // removed operational variables - now in WaterSupplySystems
    vector<WaterSource *> water_sources;
    WwtpDischargeRule wwtp_discharge_rule;
    vector<vector<double>> &demands_all_realizations;
    vector<double> demand_series_realization;
    vector<double> utility_owned_wtp_capacities; /// vector with water treatment capacity shared across one or more sources.
    vector<int> water_source_to_wtp;
    InfrastructureManager infrastructure_construction_manager;
    std::vector<std::unique_ptr<WaterSupplySystems>> water_supply_systems;

    double *P_x{}, *A_x{}, *q{}, *l{}, *u{};

    /// Drought mitigation
    double fund_contribution = 0;
    double demand_multiplier = 1;
    double demand_offset = 0;
    double restricted_price = NON_INITIALIZED;
    double offset_rate_per_volume = 0;
    double contingency_fund = 0;
    double drought_mitigation_cost = 0;
    double insurance_payout = 0;
    double insurance_purchase = 0;
    double restricted_demand = 0;
    double unrestricted_demand = 0;
    int n_sources = 0;
    double infra_discount_rate;
    double bond_term_multiplier;
    double bond_interest_rate_multiplier;
    double max_capacity = 0;

    /// Infrastructure cost
    double current_debt_payment = 0;
    double infra_net_present_cost = 0;
    vector<Bond *> issued_bonds;

public:
    const int id;
    const int number_of_week_demands;
    const char *name;
    const double percent_contingency_fund_contribution;
    const double demand_buffer;

    Utility(
            const char *name, int id,
            vector<vector<double>> &demands_all_realizations,
            int number_of_week_demands,
            double percent_contingency_fund_contribution,
            const vector<vector<double>> &typesMonthlyDemandFraction,
            const vector<vector<double>> &typesMonthlyWaterPrice,
            WwtpDischargeRule wwtp_discharge_rule,
            double demand_buffer,
            vector<vector<int>> water_source_to_wtp,
            vector<double> utility_owned_wtp_capacities);

    Utility(const char *name, int id,
            vector<vector<double>> &demands_all_realizations,
            int number_of_week_demands,
            double percent_contingency_fund_contribution,
            const vector<vector<double>> &typesMonthlyDemandFraction,
            const vector<vector<double>> &typesMonthlyWaterPrice,
            WwtpDischargeRule wwtp_discharge_rule,
            double demand_buffer,
            vector<vector<int>> water_source_to_wtp,
            vector<double> utility_owned_wtp_capacities,
            const vector<int> &rof_infra_construction_order,
            const vector<int> &demand_infra_construction_order,
            const vector<double> &infra_construction_triggers,
            double infra_discount_rate,
            const vector<vector<int>> &infra_if_built_remove, double
            bond_term, double bond_interest_rate);

    Utility(const char *name, int id,
            vector<vector<double>> &demands_all_realizations,
            int number_of_week_demands,
            double percent_contingency_fund_contribution,
            const vector<vector<double>> &typesMonthlyDemandFraction,
            const vector<vector<double>> &typesMonthlyWaterPrice,
            WwtpDischargeRule wwtp_discharge_rule,
            double demand_buffer,
            vector<vector<int>> water_source_to_wtp,
            vector<double> utility_owned_wtp_capacities,
            const vector<int> &rof_infra_construction_order,
            const vector<int> &demand_infra_construction_order,
            const vector<double> &infra_construction_triggers,
            double infra_discount_rate, double bond_term,
            double bond_interest_rate);

    Utility(Utility &utility);

    ~Utility();

    Utility &operator=(const Utility &utility);

    bool operator<(const Utility *other) const;

    bool operator>(const Utility *other) const;

    static bool compById(Utility *a, Utility *b);

    void setRisk_of_failure(double risk_of_failure);

    void updateTotalAvailableVolume(); // Delegates to WaterSupplySystems

    void calculateWastewater_releases(int week, double *discharges); // Delegates to WaterSupplySystems

    void addWaterSource(WaterSource *water_source); // Delegates to WaterSupplySystems
    
    // Add new water supply system to utility
    void addWaterSupplySystem(const std::string& name, int system_id, int utility_id,
                              const std::vector<std::vector<double>>& demands_all_realizations,
                              int number_of_week_demands,
                              const std::vector<std::vector<double>>& typesMonthlyDemandFraction,
                              const WwtpDischargeRule& wwtp_discharge_rule,
                              double demand_buffer,
                              const std::vector<std::vector<int>>& water_source_to_wtp,
                              const std::vector<double>& utility_owned_wtp_capacities);

    void splitDemands(
            int week, vector<vector<double>> &demands, bool
    apply_demand_buffer = false); // Delegates to WaterSupplySystems

    void checkErrorsAddWaterSourceOnline(WaterSource *water_source); // Delegates to WaterSupplySystems

    void resetDroughtMitigationVariables();

    void issueBond(int new_infra_triggered, int week);

    void calculateWeeklyAverageWaterPrices(
            const vector<vector<double>> &typesMonthlyDemandFraction,
            const vector<vector<double>> &typesMonthlyWaterPrice);

    double waterPrice(int week);

    void
    forceInfrastructureConstruction(int week, const vector<int>& new_infra_triggered);

    int infrastructureConstructionHandler(double long_term_rof, int week);

    static void priceCalculationErrorChecking(
            const vector<vector<double>> &typesMonthlyDemandFraction,
            const vector<vector<double>> &typesMonthlyWaterPrice);

    double getTotal_storage_capacity() const;

    double getRisk_of_failure() const;

    double getStorageToCapacityRatio() const;

    double getGrossRevenue() const;

    void setDemand_multiplier(double demand_multiplier);

    void setDemand_offset(double demand_offset, double offset_rate_per_volume);

    double getTotal_treatment_capacity() const;

    void updateContingencyFundAndDebtService(
            double unrestricted_demand, double demand_multiplier,
            double demand_offset, double unfulfilled_demand, int week);

    void setWaterSourceOnline(unsigned int source_id, int week); // Delegates to WaterSupplySystems
    
    void connectWaterSources(std::vector<WaterSource*>& water_sources,
                            std::vector<int>& priority_sources,
                            std::vector<int>& non_priority_sources);

    double updateCurrent_debt_payment(int week);

    double getContingency_fund() const;

    double getUnrestrictedDemand() const;

    double getRestrictedDemand() const;

    void setRestricted_price(double restricted_price);

    double getDemand_multiplier() const;

    double getUnrestrictedDemand(int week) const;

    double getInfrastructure_net_present_cost() const;

    double getCurrent_debt_payment() const;

    double getCurrent_contingency_fund_contribution() const;

    double getDrought_mitigation_cost() const;

    void addInsurancePayout(double payout_value);

    void clearWaterSources(); // Delegates to WaterSupplySystems

    void purchaseInsurance(double insurance_price);

    double getInsurance_payout() const;

    double getInsurance_purchase() const;

    const vector<int> &getRof_infrastructure_construction_order() const;

    void setRealization(unsigned long r, vector<double> &rdm_factors);

    vector<int> getInfrastructure_built() const;

    void setNoFinaicalCalculations();

    double getLong_term_risk_of_failure() const;

    const vector<int> &getDemand_infra_construction_order() const;

    static vector<double> calculateWeeklyPeakingFactor(vector<double> *demands); // Delegates to WaterSupplySystems

    const vector<WaterSource *> &getWater_sources() const;

    double getWaste_water_discharge() const;

    double getTotal_available_volume() const;

    void setTotal_storage_capacity();

    void resetTotal_storage_capacity();

    double getUnfulfilled_demand() const;

    double getNet_stream_inflow() const;

    double getTotal_stored_volume() const;

    const InfrastructureManager &getInfrastructure_construction_manager() const;

    double getDemand_offset() const;

    double getInfraDiscountRate() const;

    void updateTreatmentAndNumberOfStorageSources(); // Delegates to WaterSupplySystems

    static bool
    idealDemandSplitUnconstrained(double *split_demands,
                                  const double *available_treated_flow_rate,
                                  double total_demand, const double *storage,
                                  double total_storage, int n_storage_sources); // Delegates to WaterSupplySystems

    static bool
    idealDemandSplitConstrained(double *split_demands, bool *over_allocated,
                                bool *has_spare_capacity,
                                const double *available_treated_flow_rate,
                                double total_demand, const double *storage,
                                double total_storage, int n_storage_sources); // Delegates to WaterSupplySystems

    void splitDemandsQP(int week, vector<vector<double>> &demands,
                        bool apply_demand_buffer);

    bool hasTreatmentConnected(int ws); // Delegates to WaterSupplySystems

    void
    unrollWaterSourceToWtpVector(
            const vector<vector<int>> &water_source_to_wtp,
            const vector<double>& utility_owned_wtp_capacities); // Delegates to WaterSupplySystems
};


#endif //TRIANGLEMODEL_UTILITY_H
