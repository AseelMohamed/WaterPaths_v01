//
// Created by bernardo on 1/13/17.
//

#ifndef TRIANGLEMODEL_UTILITY_H
#define TRIANGLEMODEL_UTILITY_H


#include <map>
#include <memory>
#include <vector>
#include <set>
#include "../WaterSources/Reservoir.h"
#include "../../Utils/Constants.h"
#include "../../Controls/WwtpDischargeRule.h"
#include "InfrastructureManager.h"
#include "../Bonds/Base/Bond.h"

class WaterSupplySystems; // forward declaration

class Utility {
private:
    vector<double> weekly_average_volumetric_price;
    // NOTE: base_weekly_average_volumetric_price and price_rdm_multiplier removed to match Original model
    // Prices are no longer scaled by RDM factors

    
    // WSS-specific financial parameters (indexed by system_id)
    std::map<int, double> wss_contingency_percentages;  // Contingency fund % for each WSS
    std::map<int, vector<vector<double>>> wss_demand_fractions;  // Demand class fractions for each WSS
    std::map<int, vector<vector<double>>> wss_water_prices;  // Water prices for each WSS
    std::map<int, vector<double>> wss_weekly_average_prices;  // Calculated weekly average prices for each WSS
        std::map<int, double> wss_restricted_prices;  // Restricted price for each WSS (per week)
    
    // Financial and strategic variables (kept in Utility)
    double gross_revenue = 0;
    
    bool used_for_realization = true;
    WwtpDischargeRule wwtp_discharge_rule;
    vector<WaterSource *> water_sources;
    vector<vector<double>> &demands_all_realizations;
    vector<double> demand_series_realization;
    vector<double> utility_owned_wtp_capacities; /// vector with water treatment capacity shared across one or more sources.
    vector<int> water_source_to_wtp;
    InfrastructureManager infrastructure_construction_manager;
    std::vector<std::unique_ptr<WaterSupplySystems>> water_supply_systems;
    std::vector<WaterSupplySystems*> water_supply_systems_refs; // Non-owning references to WSS for simulation
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
    // THREAD-SAFE: Store base rates to prevent accumulation when setRealization is called multiple times
    double base_infra_discount_rate;
    double base_bond_term_multiplier;
    double base_bond_interest_rate_multiplier;
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
    
    // THREAD-SAFE: Global tracking of issued bonds across all threads/realizations
    // Key format: "utility_id:source_id:wss_id"
    static std::set<std::string> globally_issued_bonds;
    
    // THREAD-SAFE: Store NPC for each globally issued bond so other realizations can retrieve it
    // Maps bond_key -> NPC value
    static std::map<std::string, double> globally_issued_bond_npcs;
    
    // THREAD-SAFE: Track per-realization infrastructure costs to avoid shared state corruption
    std::map<unsigned long, double> realization_infra_costs;
    unsigned long current_realization_id = NON_INITIALIZED;

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
            vector<double> utility_owned_wtp_capacities,
            double infra_discount_rate,
            double bond_term = 15,
            double bond_interest_rate = 0.12); 

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
    
    void updateWSSReferences(const vector<WaterSupplySystems*>& new_wss);
    
    // Get WSS references for simulation (non-owning pointers)
    const vector<WaterSupplySystems*>& getWSSReferences() const;

    WaterSupplySystems& systemById(int system_id);
    const WaterSupplySystems& systemById(int system_id) const;

    void addWaterSupplySystem(const std::string& name, int system_id, int utility_id,
                              std::vector<std::vector<double>>& demands_all_realizations,
                              int number_of_week_demands,
                              const WwtpDischargeRule& wwtp_discharge_rule,
                              double demand_buffer,
                              const std::vector<std::vector<int>>& water_source_to_wtp,
                             const std::vector<double>& utility_owned_wtp_capacities,
                             double wss_contingency_percent,
                             const std::vector<std::vector<double>>& wss_demand_class_fractions,
                             const std::vector<std::vector<double>>& wss_water_prices_param);
    
    void clearWaterSupplySystems();
    
    // Must be called BEFORE parallel region starts, not during setRealization
    static void clearGlobalBondTracking();

//     int infrastructureConstructionHandler(double long_term_rof, int week); //checked 

    const vector<int> &getDemand_infra_construction_order() const; //checked

    void purchaseInsurance(double insurance_price); //checked

    double updateCurrent_debt_payment(int week, const std::vector<WaterSupplySystems*>& current_wss); //checked

    void setNoFinancialCalculations();
    
    void priceCalculationErrorChecking(const vector<vector<double>> &typesMonthlyDemandFraction,
                                       const vector<vector<double>> &typesMonthlyWaterPrice); //checked

    void calculateWeeklyAverageWaterPrices(const vector<vector<double>> &typesMonthlyDemandFraction,
                                           const vector<vector<double>> &typesMonthlyWaterPrice); //checked
    
    // Overloaded version for WSS-specific price calculation (output to parameter)
    void calculateWeeklyAverageWaterPrices(const vector<vector<double>> &typesMonthlyDemandFraction,
                                           const vector<vector<double>> &typesMonthlyWaterPrice,
                                           vector<double>& output_weekly_prices);

    void setDemand_offset(double demand_offset, double offset_rate_per_volume);

    void forceInfrastructureConstruction(int week, const vector<int>& new_infra_triggered); //checked

    void issueBond(int new_infra_triggered, int week); //checked

    //  Overloaded version to issue bond for infrastructure in specific WSS
    void issueBond(int new_infra_triggered, int week, WaterSupplySystems* target_wss, double discount_rate,
                   double bond_term_mult, double bond_interest_rate_mult); 

    void resetDroughtMitigationVariables(); //Checked

        double waterPrice(int week) const; //Checked
        double getCurrentWaterPrice(int week) const;

        // Per-WSS restricted price helpers
        void setRestrictedPriceForWss(int system_id, double restricted_price);
        double getRestrictedPriceForWss(int system_id) const;
        double calculateRestrictedWeeklyPriceForWss(int system_id, int stage,
                                                                                                int week_of_year,
                                                                                                const vector<vector<double>> &priceMultipliers) const;

    double getGrossRevenue() const; //checked

    void updateContingencyFundAndDebtService(
            double unrestricted_demand, double demand_multiplier,
            double demand_offset, double unfulfilled_demand, int week); //checked

    // NEW: Aggregate financial calculations across all WSS (following Original model pattern)
    void updateUtilityFinancialCalculations(int week, const std::vector<WaterSupplySystems*>& realization_wss);

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

    vector<int> getInfrastructure_built() const;
    
    vector<vector<int>> getAllAndClearInfraBuilt() const; //checked

    double getTotal_storage_capacity() const; //check --- 

    double getTotal_available_volume() const; // NEW - aggregates from all WSS
    
    double getTotal_stored_volume() const; // NEW - aggregates stored volume from all WSS

    void updateTotalAvailableVolume(); // NEW - delegates to all WSS
    
    void setRealization(unsigned long r, vector<double> &rdm_factors);

    double getUnfulfilled_demand() const; 

    double getRestrictedDemand() const; //checked

    double getUnrestrictedDemand(int week = -1) const; //checked

    double getTotal_treatment_capacity() const; 

    const InfrastructureManager &getInfrastructure_construction_manager() const;

    double getInfraDiscountRate() const;
    
    double getBaseInfraDiscountRate() const;
    
    void resetFinancialState();

    bool isUsedForRealization() const; // Getter for realization flag

    void
    unrollWaterSourceToWtpVector(
            const vector<vector<int>> &water_source_to_wtp,
            const vector<double>& utility_owned_wtp_capacities); // Delegates to all WaterSupplySystems
};


#endif //TRIANGLEMODEL_UTILITY_H
