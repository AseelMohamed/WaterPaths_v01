//
// Created by bernardo on 4/17/18.
//

#ifndef TRIANGLEMODEL_BALLOONPAYMENTBOND_H
#define TRIANGLEMODEL_BALLOONPAYMENTBOND_H

#include <vector>


#include "Base/Bond.h"

class BalloonPaymentBond : public Bond {
private:
    double interest_payments;
    const int begin_repayment_after_n_years;
    int n_payments_made = 0;
    int last_payment_week = NON_INITIALIZED;

public:
    BalloonPaymentBond(const int id, const double cost_of_capital, const int n_payments,
                         const double coupon_rate, vector<int> pay_on_weeks);

    BalloonPaymentBond(const int id, const double cost_of_capital, const int n_payments, const double coupon_rate,
                       vector<int> pay_on_weeks, const int begin_repayment_after_n_years);

    double getDebtService(int week) override;

    double getNetPresentValueAtIssuance(double yearly_discount_rate, int week) const override;

    void issueBond(int week, int construction_time, double bond_term_multiplier, double bond_interest_rate_multiplier) override;
    
    void resetForRealization() override { Bond::resetForRealization(); n_payments_made = 0; last_payment_week = NON_INITIALIZED; }
};


#endif //TRIANGLEMODEL_BALLOONPAYMENTBOND_H
