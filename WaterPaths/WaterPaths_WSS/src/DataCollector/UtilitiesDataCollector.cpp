//
// Created by bernardoct on 8/25/17.
//

#include <iomanip>
#include <sstream>
//#include <cmath>
#include "UtilitiesDataCollector.h"


UtilitiesDataCollector::UtilitiesDataCollector(const Utility *utility, unsigned long realization)
        : DataCollector(utility->id, utility->name, realization, UTILITY, 15 * COLUMN_WIDTH),
          utility(utility) {
}

string UtilitiesDataCollector::printTabularString(int week) {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << combined_storage[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << capacity[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
	          << net_stream_inflow[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << st_rof[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << lt_rof[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << unfulfilled_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << restricted_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << unrestricted_demand[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << waste_water_discharge[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
              << total_treatment_capacity[week]
              << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
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

    stringstream outStream;

    outStream << combined_storage[week]
              << ","
              << capacity[week]
              << ","
	          << net_stream_inflow[week]
              << ","
              << st_rof[week]
              << ","
              << lt_rof[week]
              << ","
              << restricted_demand[week]
              << ","
              << unrestricted_demand[week]
              << ","
              << unfulfilled_demand[week]
              << ","
              << waste_water_discharge[week]
              << ","
              << total_treatment_capacity[week]
              << ","
              << contingency_fund_size[week]
              << ","
              << insurance_payout[week]
              << ","
              << insurance_contract_cost[week]
              << ","
              << net_present_infrastructure_cost[week]
              << ","
              << debt_service_payments[week]
              << ",";
    
    return outStream.str();
}

string UtilitiesDataCollector::printTabularStringHeaderLine1() {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << "Stored"
              << setw(COLUMN_WIDTH) << " "
              << setw(COLUMN_WIDTH) << "Net"
              << setw(COLUMN_WIDTH) << " "
              << setw(COLUMN_WIDTH) << " "
              << setw(COLUMN_WIDTH) << "Rest."
              << setw(COLUMN_WIDTH) << "Unrest."
              << setw(COLUMN_WIDTH) << "Unfulf."
              << setw(COLUMN_WIDTH) << "W. Water"
              << setw(COLUMN_WIDTH) << "Treatment"
              << setw(COLUMN_WIDTH) << "Conting."
              << setw(COLUMN_WIDTH) << "Insurance"
              << setw(COLUMN_WIDTH) << "Insurance"
              << setw(COLUMN_WIDTH) << "Infra."
              << setw(COLUMN_WIDTH) << "Debt";

    return outStream.str();
}

string UtilitiesDataCollector::printTabularStringHeaderLine2() {

    stringstream outStream;

    outStream << setw(2 * COLUMN_WIDTH) << "Volume"
              << setw(COLUMN_WIDTH) << "Capacity"
              << setw(COLUMN_WIDTH) << "Inflow"
              << setw(COLUMN_WIDTH) << "ST-ROF"
              << setw(COLUMN_WIDTH) << "LT-ROF"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Demand"
              << setw(COLUMN_WIDTH) << "Discharge"
              << setw(COLUMN_WIDTH) << "Capacity"
              << setw(COLUMN_WIDTH) << "Fund"
              << setw(COLUMN_WIDTH) << "Payout"
              << setw(COLUMN_WIDTH) << "Price"
              << setw(COLUMN_WIDTH) << "NPV"
              << setw(COLUMN_WIDTH) << "Service";

    return outStream.str();
}

string UtilitiesDataCollector::printCompactStringHeader() {
    stringstream outStream;

    outStream << id << "st_vol" << ","
              << id << "capacity" << ","
              << id << "net_inf" << ","
              << id << "st_rof" << ","
              << id << "lt_rof" << ","
              << id << "rest_demand" << ","
              << id << "unrest_demand" << ","
              << id << "unfulf_demand" << ","
              << id << "wastewater" << ","
              << id << "treat_capacity" << ","
              << id << "cont_fund" << ","
              << id << "ins_pout" << ","
              << id << "ins_price" << ","
              << id << "infra_npv" << ","
              << id << "debt_serv" << ",";
    
    return outStream.str();
}

void UtilitiesDataCollector::collect_data() {
    // Restore storage collection for reliability calculations
    // Use stored volume (actual water in storage) rather than available volume
    // (which includes flow-through sources that don't show storage depletion)
    combined_storage.push_back(utility->getTotal_stored_volume());
    lt_rof.push_back(0.0);  // WSS doesn't calculate LT-ROF at utility level
    st_rof.push_back(0.0);  // WSS doesn't calculate ST-ROF at utility level
    unrestricted_demand.push_back(utility->getUnrestrictedDemand());
    restricted_demand.push_back(utility->getRestrictedDemand());
    contingency_fund_size.push_back(utility->getContingency_fund());
    net_present_infrastructure_cost.push_back(utility->getInfrastructure_net_present_cost());
    
    // Track when infrastructure cost data is collected each week
    static int debug_week_counter = 0;
    double gr = utility->getGrossRevenue();
    double dsp = utility->getCurrent_debt_payment();
    double cfc = utility->getCurrent_contingency_fund_contribution();
    double dmc = utility->getDrought_mitigation_cost();
    
    debug_week_counter++;
    
    gross_revenues.push_back(gr);
    debt_service_payments.push_back(dsp);
    contingency_fund_contribution.push_back(cfc);
    drought_mitigation_cost.push_back(dmc);
    insurance_contract_cost.push_back(utility->getInsurance_purchase());
    insurance_payout.push_back(utility->getInsurance_payout());
    capacity.push_back(utility->getTotal_storage_capacity());
    waste_water_discharge.push_back(0.0);  // WSS doesn't track wastewater discharge
    unfulfilled_demand.push_back(utility->getUnfulfilled_demand());
    net_stream_inflow.push_back(0.0);  // WSS doesn't track net stream inflow
    total_treatment_capacity.push_back(utility->getTotal_treatment_capacity());

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
