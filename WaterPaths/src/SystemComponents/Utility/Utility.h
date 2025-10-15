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

class WaterSupplySystems; // forward declaration

class Utility {
private:
    vector<double> weekly_average_volumetric_price;
    // Financial and strategic variables (kept in Utility)
    double gross_revenue = 0;
    
    double price_rdm_multiplier = 1.;
    bool used_for_realization = true;
    WwtpDischargeRule wwtp_discharge_rule;
    vector<WaterSource *> water_sources;
    vector<vector<double>> &demands_all_realizations;
    vector<double> demand_series_realization;
    vector<double> utility_owned_wtp_capacities; /// vector with water treatment capacity shared across one or more sources.
    vector<int> water_source_to_wtp;
    InfrastructureManager infrastructure_construction_manager;
    std::vector<std::unique_ptr<WaterSupplySystems>> water_supply_systems;
    std::vector<double> utility_demand_series_realization;

    /// Drought mitigation
    double fund_contribution = 0;
    double contingency_fund = 0;
    double drought_mitigation_cost = 0;
    double insurance_payout = 0;
    double insurance_purchase = 0;
    int n_sources = 0;
    double infra_discount_rate;
    double bond_term_multiplier;
    double bond_interest_rate_multiplier;
    double max_capacity = 0;
    double restricted_price = NON_INITIALIZED;
    double utility_restricted_demand = 0;
    double utility_unrestricted_demand = 0;
    double utility_total_storage_capacity = 0;
    double utility_total_available_volume = 0;
    double utility_total_stored_volume = 0;
    double offset_rate_per_volume = 0;
    double demand_offset = 0;
    double demand_multiplier = 1;

    /// Infrastructure cost
    double current_debt_payment = 0;
    double infra_net_present_cost = 0;
    vector<Bond *> issued_bonds;

    // Helper method to find system containing a specific water source
    WaterSupplySystems& systemForSource(int source_id);

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

    static bool compById(Utility *a, Utility *b); //Checked

    const vector<unique_ptr<WaterSupplySystems>> &getWaterSupplySystems() const;

    WaterSupplySystems& systemById(int system_id);
    const WaterSupplySystems& systemById(int system_id) const;

    void addWaterSupplySystem(const std::string& name, int system_id, int utility_id,
                              std::vector<std::vector<double>>& demands_all_realizations,
                              int number_of_week_demands,
                              const WwtpDischargeRule& wwtp_discharge_rule,
                              double demand_buffer,
                              const std::vector<std::vector<int>>& water_source_to_wtp,
                              const std::vector<double>& utility_owned_wtp_capacities);

    void clearWaterSupplySystems();

    void setRealization(unsigned long r, vector<double> &rdm_factors);

//     int infrastructureConstructionHandler(double long_term_rof, int week); //checked

    const vector<int> &getDemand_infra_construction_order() const; //checked

    void purchaseInsurance(double insurance_price); //checked

    double updateCurrent_debt_payment(int week); //checked

    void setNoFinancialCalculations();
    
    void priceCalculationErrorChecking(const vector<vector<double>> &typesMonthlyDemandFraction,
                                       const vector<vector<double>> &typesMonthlyWaterPrice); //checked

    void calculateWeeklyAverageWaterPrices(const vector<vector<double>> &typesMonthlyDemandFraction,
                                           const vector<vector<double>> &typesMonthlyWaterPrice); //checked

    void setDemand_offset(double demand_offset, double offset_rate_per_volume);

    void forceInfrastructureConstruction(int week, const vector<int>& new_infra_triggered); //checked

    void issueBond(int new_infra_triggered, int week); //checked

    void resetDroughtMitigationVariables(); //Checked

    double waterPrice(int week); //Checked

    double getGrossRevenue() const; //checked

    void updateContingencyFundAndDebtService(
            double unrestricted_demand, double demand_multiplier,
            double demand_offset, double unfulfilled_demand, int week); //checked

    void setRestricted_price(double restricted_price); //checked

    double getInfrastructure_net_present_cost() const; //checked

    double getCurrent_debt_payment() const; //checked

    double getCurrent_contingency_fund_contribution() const; //checked

    void addInsurancePayout(double payout_value); //checked

    double getInsurance_purchase() const; //checked

    const vector<int> &getRof_infrastructure_construction_order() const; //checked

    double getGross_revenue() const;

    double getContingency_fund() const; //checked

    double getDrought_mitigation_cost() const; //checked

    double getInsurance_payout() const; //checked

    vector<int> getInfrastructure_built() const; //checked

    double getTotal_storage_capacity() const; //check --- 

    double getUnfulfilled_demand() const; 

    double getRestrictedDemand() const; //checked

    double getUnrestrictedDemand(int week = -1) const; //checked

    double getTotal_treatment_capacity() const; 

    const InfrastructureManager &getInfrastructure_construction_manager() const;

    double getInfraDiscountRate() const;

    void
    unrollWaterSourceToWtpVector(
            const vector<vector<int>> &water_source_to_wtp,
            const vector<double>& utility_owned_wtp_capacities); // Delegates to all WaterSupplySystems
};


#endif //TRIANGLEMODEL_UTILITY_H
