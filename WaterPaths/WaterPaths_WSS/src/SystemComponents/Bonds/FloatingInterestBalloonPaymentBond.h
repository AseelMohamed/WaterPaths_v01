//
// Created by bernardo on 4/17/18.
//

#ifndef TRIANGLEMODEL_FLOATINGINTERESTBOND_H
#define TRIANGLEMODEL_FLOATINGINTERESTBOND_H


#include "Base/Bond.h"

class FloatingInterestBalloonPaymentBond : public Bond {
private:
    vector<double> interest_rate_series;
    const int begin_repayment_after_n_years;
    int n_payments_made = 0;
    int last_payment_week = NON_INITIALIZED;

public:
    FloatingInterestBalloonPaymentBond(const int id, const double cost_of_capital, const double n_payments,
                         const vector<double> interest_rate_series, vector<int> pay_on_weeks);

    FloatingInterestBalloonPaymentBond(const int id, const double cost_of_capital, double n_payments,
                         const vector<double> interest_rate_series, vector<int> pay_on_weeks,
                         const int starts_paying_after_n_years);

    double getDebtService(int week) override;

    double getNetPresentValueAtIssuance(double yearly_discount_rate, int week) const override;

    void issueBond(int week, int construction_time, double bond_term_multiplier, double bond_interest_rate_multiplier) override;
    
    void resetForRealization() override { Bond::resetForRealization(); n_payments_made = 0; last_payment_week = NON_INITIALIZED; }
};


#endif //TRIANGLEMODEL_FLOATINGINTERESTBOND_H
