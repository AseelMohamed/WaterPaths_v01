//
// Criado por Bernardo em 23/11/2017 e adaptado por Andressa em 21/08/2019.
//
//DISTRITO FEDERAL MODEL
#ifndef TRIANGLEMODEL_CAESB_H
#define TRIANGLEMODEL_CAESB_H

#ifdef PARALLEL
#include "../../Borg/borgms.h"
#endif
#include "Base/Problem.h"
#include "../Simulation/Simulation.h"

using namespace std;

class Caesb : public Problem {
private:
    const int n_utilities = 2;
//obs: não inseri os streamflows do subsistema gama (provavelmente será excluído da análise)
    vector<vector<double>> streamflows_descoberto; //descoberto
    vector<vector<double>> streamflows_tortoSM; //captação no Santa Maria
    vector<vector<double>> streamflows_bananal_torto; //captação no bananal e no torto
    //vector<vector<double>> streamflows_torto; //captação no torto
    vector<vector<double>> streamflows_paranoa; //lago paranoa
    vector<vector<double>> streamflows_corumbaIV; //corumba
    vector<vector<double>> demand_caesb_descoberto;
    vector<vector<double>> demand_caesb_tortoSM;
    vector<vector<double>> evap_descoberto;
    vector<vector<double>> evap_tortoSM;
    vector<vector<double>> evap_paranoa;
    vector<vector<double>> evap_corumba;
    vector<vector<double>> demand_to_wastewater_fraction_caesb_descoberto;
    vector<vector<double>> demand_to_wastewater_fraction_caesb_tortoSM;
    vector<vector<double>> caesbDescobertoDemandClassesFractions;
    vector<vector<double>> caesbTortoSMDemandClassesFractions;
    vector<vector<double>> caesbDescobertoUserClassesWaterPrices;
    vector<vector<double>> caesbTortoSMUserClassesWaterPrices;
    vector<vector<double>> caesbPriceRestrictionMultipliers;

public:
    Caesb(unsigned long n_weeks, int import_export_rof_table);

    ~Caesb();

#ifdef PARALLEL
    void setProblemDefinition(BORG_Problem &problem);
#endif

    int functionEvaluation(double *vars, double *objs, double *consts) override;

    int simulationExceptionHander(const std::exception &e, Simulation *s,
                                  double *objs, const double *vars);

    void readInputData();

};


#endif //TRIANGLEMODEL_CAESB_H
