//
// Created by bernardoct on 8/25/17.
//

#include <iomanip>
#include <sstream>
//#include <cmath>
#include "UtilitiesDataCollector.h"


UtilitiesDataCollector::UtilitiesDataCollector(const Utility *utility, unsigned long realization, double discount_rate)
    : DataCollector(utility->id, utility->name, realization, UTILITY, 13 * COLUMN_WIDTH),
          utility(utility),
          infra_discount_rate(discount_rate) {
}

string UtilitiesDataCollector::printTabularString(int week) {

    stringstream outStream;

    outStream << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << contingency_fund_size[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << insurance_payout[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << insurance_contract_cost[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << net_present_infrastructure_cost[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << debt_service_payments[week];

    return outStream.str();
}

string UtilitiesDataCollector::printCompactString(int week) {

#ifndef NDEBUG
    // Debug-mode bounds checking to catch uninitialized memory access
    auto check = [week, this](const std::vector<double> &v, const char *name) {
        if (week < 0 || (size_t)week >= v.size()) {
            std::ostringstream oss;
            oss << "Week index " << week << " out of range for " << name
                << " (size=" << v.size() << "), utility id=" << id
                << ", realization=" << realization;
            throw std::out_of_range(oss.str());
        }
    };

    check(contingency_fund_size, "contingency_fund_size");
    check(insurance_payout, "insurance_payout");
    check(insurance_contract_cost, "insurance_contract_cost");
    check(net_present_infrastructure_cost, "net_present_infrastructure_cost");
    check(debt_service_payments, "debt_service_payments");
    check(gross_revenues, "gross_revenues");
    check(contingency_fund_contribution, "contingency_fund_contribution");
    check(drought_mitigation_cost, "drought_mitigation_cost");
    check(water_price, "water_price");
#endif

    stringstream outStream;

    outStream << contingency_fund_size[week]
              << ","
              << insurance_payout[week]
              << ","
              << insurance_contract_cost[week]
              << ","
              << net_present_infrastructure_cost[week]
              << ","
              << debt_service_payments[week]
              << ","
              << gross_revenues[week]
              << ","
              << contingency_fund_contribution[week]
              << ","
              << drought_mitigation_cost[week]
              << ","
              << water_price[week]
              << ",";
    
    return outStream.str();
}

string UtilitiesDataCollector::printTabularStringHeaderLine1() {

    stringstream outStream;

    outStream << setw(COLUMN_WIDTH) << "Conting."
              << setw(COLUMN_WIDTH) << "Insurance"
              << setw(COLUMN_WIDTH) << "Insurance"
              << setw(COLUMN_WIDTH) << "Infra."
              << setw(COLUMN_WIDTH) << "Debt";

    return outStream.str();
}

string UtilitiesDataCollector::printTabularStringHeaderLine2() {

    stringstream outStream;

    outStream << setw(COLUMN_WIDTH) << "Fund"
              << setw(COLUMN_WIDTH) << "Payout"
              << setw(COLUMN_WIDTH) << "Price"
              << setw(COLUMN_WIDTH) << "NPV"
              << setw(COLUMN_WIDTH) << "Service";

    return outStream.str();
}

string UtilitiesDataCollector::printCompactStringHeader() {
    stringstream outStream;

    outStream << id << "cont_fund" << ","
              << id << "ins_pout" << ","
              << id << "ins_price" << ","
              << id << "infra_npv" << ","
              << id << "debt_serv" << ","
              << id << "gross_rev" << ","
              << id << "cont_contrib" << ","
              << id << "drought_cost" << ","
              << id << "water_price" << ",";
    
    return outStream.str();
}

void UtilitiesDataCollector::collect_data() {
    // Collect operational data for internal use (objectives, NetCDF, reliability calculations)
    // but these won't be printed in CSV output - only WSS files will show them
    double total_storage = 0.0;
    double total_capacity = 0.0;
    double total_net_inflow = 0.0;
    double total_rest_demand = 0.0;
    double total_unrest_demand = 0.0;
    double total_unfulf_demand = 0.0;
    double total_wastewater = 0.0;
    double total_treatment_cap = 0.0;
    
    // Aggregate from owned WSSs
    const auto& wss_refs = utility->getWSSReferences();
    for (const auto* wss : wss_refs) {
        if (wss != nullptr) {
            total_storage += wss->getTotal_stored_volume();
            total_capacity += wss->getTotal_storage_capacity();
            total_net_inflow += wss->getNet_stream_inflow();
            total_rest_demand += wss->getRestrictedDemand();
            total_unrest_demand += wss->getUnrestrictedDemand();
            total_unfulf_demand += wss->getUnfulfilled_demand();
            total_wastewater += wss->getWaste_water_discharge();
            total_treatment_cap += wss->getTotal_treatment_capacity();
        }
    }
    
    // Store operational data (for internal calculations, not CSV output)
    combined_storage.push_back(total_storage);
    capacity.push_back(total_capacity);
    net_stream_inflow.push_back(total_net_inflow);
    restricted_demand.push_back(total_rest_demand);
    unrestricted_demand.push_back(total_unrest_demand);
    unfulfilled_demand.push_back(total_unfulf_demand);
    waste_water_discharge.push_back(total_wastewater);
    total_treatment_capacity.push_back(total_treatment_cap);
    st_rof.push_back(0.0);  // WSS doesn't calculate ROF at utility level
    lt_rof.push_back(0.0);  // WSS doesn't calculate ROF at utility level
    
    // Collect financial data - copy atomically from shared utility to avoid race conditions
    // Each thread reads from shared utility, so serialize the reads
    double cf, npc, gr, dsp, cfc, dmc, ins_cost, ins_payout, wp;
    #pragma omp critical(utility_financial_read)
    {
        cf = utility->getContingency_fund();
        npc = utility->getInfrastructure_net_present_cost();
        gr = utility->getGrossRevenue();
        dsp = utility->getCurrent_debt_payment();
        cfc = utility->getCurrent_contingency_fund_contribution();
        dmc = utility->getDrought_mitigation_cost();
        ins_cost = utility->getInsurance_purchase();
        ins_payout = utility->getInsurance_payout();
        wp = utility->getCurrentWaterPrice((int)gross_revenues.size());
    }
    
    contingency_fund_size.push_back(cf);
    net_present_infrastructure_cost.push_back(npc);
    gross_revenues.push_back(gr);
    debt_service_payments.push_back(dsp);
    contingency_fund_contribution.push_back(cfc);
    drought_mitigation_cost.push_back(dmc);
    insurance_contract_cost.push_back(ins_cost);
    insurance_payout.push_back(ins_payout);
    water_price.push_back(wp);

//    checkForNans();

    // Note: Pathway collection moved to WSSDataCollector since those collectors
    // point to the actual realization WSS where infrastructure is built.
    // Utility data collectors point to original utilities with original WSS.
}

void UtilitiesDataCollector::checkForNans() const {
    string error = "nan collecting data for utility " + to_string(id) + " in week " + to_string(lt_rof.size
            ()) + ", realization " + to_string(realization);
    if (std::isnan(unrestricted_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(restricted_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(combined_storage.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(lt_rof.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(st_rof.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(contingency_fund_size.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(net_present_infrastructure_cost.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(gross_revenues.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(debt_service_payments.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(contingency_fund_contribution.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(drought_mitigation_cost.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(insurance_contract_cost.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(insurance_payout.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(water_price.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(capacity.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(waste_water_discharge.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(unfulfilled_demand.back()))
        throw_with_nested(runtime_error(error.c_str()));
    if (std::isnan(net_stream_inflow.back()))
        throw_with_nested(runtime_error(error.c_str()));

    error = "NPV absurdly high when collecting data for utility " + to_string(id) + " in week " + to_string(lt_rof.size
            ()) + ", realization " + to_string(realization) + "\n";
    if (net_present_infrastructure_cost.back() > 1e10) printf("%s", error.c_str());//throw_with_nested(runtime_error(error.c_str()));
    
}

const vector<double> &UtilitiesDataCollector::getCombined_storage() const {
    return combined_storage;
}

const vector<double> &UtilitiesDataCollector::getCapacity() const {
    return capacity;
}

const vector<double> &UtilitiesDataCollector::getGross_revenues() const {
    return gross_revenues;
}

const vector<double> &
UtilitiesDataCollector::getContingency_fund_contribution() const {
    return contingency_fund_contribution;
}

const vector<double> &UtilitiesDataCollector::getDebt_service_payments() const {
    return debt_service_payments;
}

const vector<double> &
UtilitiesDataCollector::getInsurance_contract_cost() const {
    return insurance_contract_cost;
}

const vector<double> &
UtilitiesDataCollector::getDrought_mitigation_cost() const {
    return drought_mitigation_cost;
}

const vector<double> &UtilitiesDataCollector::getContingency_fund_size() const {
    return contingency_fund_size;
}

const vector<vector<int>> &UtilitiesDataCollector::getPathways() const {
    return pathways;
}

const vector<double> &
UtilitiesDataCollector::getNet_present_infrastructure_cost() const {
    return net_present_infrastructure_cost;
}

const vector<double> &UtilitiesDataCollector::getSt_rof() const {
    return st_rof;
}

const vector<double> &UtilitiesDataCollector::getLt_rof() const {
    return lt_rof;
}

const vector<double> &UtilitiesDataCollector::getRestricted_demand() const {
    return restricted_demand;
}

const Utility *UtilitiesDataCollector::getUtility() const {
    return utility;
}
