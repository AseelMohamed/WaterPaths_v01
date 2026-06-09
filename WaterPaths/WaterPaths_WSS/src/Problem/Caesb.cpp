// Criado por Bernardo em 23/11/2017 e adaptado por Andressa em 12/08/2019.
//Os arquivos .cpp contêm as fontes das implementações em c++. Já os arquivos .h são os chamados headers, utilizados para extrair a declaração de funções, classes e outras declarações, que podem ser compartilhados entre vários arquivos com o código fonte. São, por isso, reutilizáveis. Esses arquivos possuem códigos que o compilador precisa para compilar outras partes

#include <algorithm> //biblioteca que define uma coleção de funções a serem usadas em uma gama de elementos
#include <numeric> //biblioteca que descreve um conjunto de algoritmos para executar determinadas operações em sequências de valores numéricos
#include <iostream> //biblioteca que define os objetos padrão de fluxo de entrada / saída
#include <iterator> //iterador é qualquer objeto que, apontando para algum elemento de um intervalo de elementos (como uma matriz ou um contêiner), tem a capacidade de iterar através dos elementos desse intervalo usando um conjunto de operadores (por exemplo, o operador de incremento ++)
#include <fstream> //biblioteca que fornece classes de fluxos de arquivos
#include <omp.h> //biblioteca que permite a programação multi-processo de memória compartilhada em múltiplas plataformas, de modo a acrescentar simultaneidade aos programas escritos em C++
#include <vector>

#ifdef  PARALLEL
#include <mpi.h> // biblioteca que permite a comunicação de dados em computação paralela
#endif

#include "Caesb.h" //include é uma diretiva de compilação
#include "../Utils/Constants.h"
#include "../Controls/SeasonalMinEnvFlowControl.h"
#include "../Controls/StorageMinEnvFlowControl.h"
#include "../Controls/InflowMinEnvFlowControl.h"
#include "../Controls/FixedMinEnvFlowControl.h"
#include "../SystemComponents/WaterSources/ReservoirExpansion.h"
#include "../SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h"
#include "../Simulation/Simulation.h"
#include "../SystemComponents/Bonds/LevelDebtServiceBond.h"
#include "../SystemComponents/Bonds/BalloonPaymentBond.h"
#include "../DroughtMitigationInstruments/Transfers.h"
#include "../DroughtMitigationInstruments/TransfersBilateral.h"
#include "../DroughtMitigationInstruments/EmergencyTransferParanoa.h"
#include "../SystemComponents/WaterSources/AllocatedReservoir.h"

/**
#include "../../../../../AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/usr/include/c++/7/cmath"
#include "../../../../../AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/usr/include/c++/7/vector"
**/


#ifdef PARALLEL //#ifdef é uma diretiva que permite que uma seção de um programa seja compilado somente se a macro especificada como parâmetro tiver sido definida, não importando qual seja o seu valor.
void Caesb::setProblemDefinition(BORG_Problem &problem) //void = vazio. O tipo void permite fazer funções que não retornam nada e funções que não têm parâmetros. Algumas funções não precisam retornar nenhum valor para funcionar, apenas realizar alguma ação.
{ //BORG é o algoritmo de otimização.(problem, número de identificação da variável de decisão, limite inferior da variável, limite superior da variável)
    // The parameter bounds are the same for all formulations
    BORG_Problem_set_bounds(problem, 0, 0.001, 1.1); //Gatilho para acionar a restrição de uso da água - descoberto
    BORG_Problem_set_bounds(problem, 1, 0.001, 1.1); //Gatilho para acionar a restrição de uso da água - tortoSM
    BORG_Problem_set_bounds(problem, 2, 0.001, 1.1); //Diferença entre um estágio de maior restrição para o estágio de menor restrição - descoberto
    BORG_Problem_set_bounds(problem, 3, 0.001, 1.1); //Diferença entre um estágio de maior restrição para o estágio de menor restrição - tortoSM
    BORG_Problem_set_bounds(problem, 4, 0.001, 1.1); //Gatilho para acionar a transferência de água entre sistemas - descoberto
    BORG_Problem_set_bounds(problem, 5, 0.001, 1.1); //Gatilho para acionar a transferência de água entre sistemas - tortoSM
    BORG_Problem_set_bounds(problem, 6, 0.0, 1e-12); // Percentual da receita anual alocada para o fundo de contingência da companhia 0. Fixado em 0 na avaliação; limite não degenerado para evitar NaN no operador PM.
    BORG_Problem_set_bounds(problem, 7, 0.0, 1e-12); // Percentual da receita anual alocada para o fundo de contingência da companhia 1. Fixado em 0 na avaliação; limite não degenerado para evitar NaN no operador PM.
    BORG_Problem_set_bounds(problem, 8, 0.001, 1.1); //Gatilho para acionar construção de infraestrutura pela Caesb - descoberto
    BORG_Problem_set_bounds(problem, 9, 0.001, 1.1); //Gatilho para acionar construção de infraestrutura pela Caesb - tortoSM
    BORG_Problem_set_bounds(problem, 10, 0.0, 1.0); //Ordem de "construção" da Etapa 1 da ETA Paranoá Sul (ID 7: + 0.7 m³/s) — TortoSM
    BORG_Problem_set_bounds(problem, 11, 0.0, 1.0); //Ordem de "construção" da Etapa 2 da ETA Paranoá Sul (ID 8: + 0.7 m³/s) — TortoSM
    BORG_Problem_set_bounds(problem, 12, 0.0, 1.0); //Ordem de "construção" da Etapa 3 das ETAs Paranoá Sul e Norte (ID 9: + 0.7 m³/s) — TortoSM
    BORG_Problem_set_bounds(problem, 13, 0.0, 1.0); //Ordem de "construção" da Etapa 1 da ETA Corumbá (ID 5: + 1.4 m³/s, total 2.8 m³/s)
    BORG_Problem_set_bounds(problem, 14, 0.0, 1.0); //Ordem de "construção" da Etapa 2 da ETA Corumbá (ID 6: + 1.2 m³/s, total 4.0 m³/s)
    BORG_Problem_set_bounds(problem, 15, 0.0, 1.0); //Ordem de "construção" da expansão do Descoberto (ID 10: + 25% storage)
    BORG_Problem_set_bounds(problem, 16, 0.1, 0.5); //Buffer de infraestrutura por parte da Caesb - descoberto
    BORG_Problem_set_bounds(problem, 17, 0.1, 0.5); //Buffer de infraestrutura por parte da Caesb - tortoSM
    BORG_Problem_set_bounds(problem, 18, 0.001, 1.1); //Gatilho ROF para ativar transferência de emergência do Paranoá para TortoSM
    BORG_Problem_set_bounds(problem, 19, 0.0, 1.2); //Fração da demanda irrestrita do emissor protegida antes da transferência (0=sem proteção, 1=proteção total, 2=margem de segurança de 100%)
    BORG_Problem_set_bounds(problem, 20, 0.0, 1.0); //Ordem de "construção" da Ampliação da ETA Santa Maria (ID 12: +0.7 m³/s) — TortoSM; sequência força esta antes das ETAs do Paranoá
    BORG_Problem_set_bounds(problem, 21, 1.0, 1.4); //Multiplicador de preço para restrição de uso residencial (estágio 3, col 0) — Descoberto e TortoSM

    // Set epsilons for objectives //(problem, n° de identificação da função objetivo, valor do epsilon). O valor do epsilon indica a precisão das funções objetivo.
    if (Constants::includePerWSSObjectives()) {
        // Mode 5: 7 objectives [WSS0_rel, WSS1_rel, restr_freq, NPC, worst_cost, WSS0_afford, WSS1_afford]
        BORG_Problem_set_epsilon(problem, 0, 0.001);      // WSS0 reliability
        BORG_Problem_set_epsilon(problem, 1, 0.001);      // WSS1 reliability
        BORG_Problem_set_epsilon(problem, 2, 0.005);      // Restriction frequency
        BORG_Problem_set_epsilon(problem, 3, 10000000.);  // Infrastructure NPC
        BORG_Problem_set_epsilon(problem, 4, 0.005);      // Worst case costs
        BORG_Problem_set_epsilon(problem, 5, 0.001);      // WSS0 affordability
        BORG_Problem_set_epsilon(problem, 6, 0.001);      // WSS1 affordability
    } else {
        // Modes 1–4: [reliability, restr_freq, NPC, worst_cost (, affordability)]
        BORG_Problem_set_epsilon(problem, 0, 0.001);      // Reliability
        BORG_Problem_set_epsilon(problem, 1, 0.005);      // Restriction frequency
        BORG_Problem_set_epsilon(problem, 2, 10000000.);  // Infrastructure NPC
        BORG_Problem_set_epsilon(problem, 3, 0.005);      // Worst case costs
        if (Constants::includeAffordabilityObjective()) {
            BORG_Problem_set_epsilon(problem, 4, 0.001);  // Affordability index
        }
    }
}
#endif


/**
 * Rodar o problema do Distrito Federal
 * @param vars
 * @param n_realizations
 * @param n_weeks
 * @param sol_number
 * @param output_directory

 */
int Caesb::functionEvaluation(double *vars, double *objs, double *consts) {
    double realization_start = omp_get_wtime();

//    cout << "Building Caesb Problem." << endl;
    // ===================== SET UP DECISION VARIABLES  =====================

    Simulation *s = nullptr; //nullptr: ponteiro nulo //*operador de dereferencia: a dereferencia busca o VALOR que está contida no endereço gravado no ponteiro
//    try {
    //throw invalid_argument("Test error");

    //VARIÁVEIS DE DECISÃO

    double caesb_descoberto_restriction_trigger = vars[0]; //gatilho para acionar restrição no Descoberto
    double caesb_tortoSM_restriction_trigger = vars[1]; //gatilho para acionar restrição no TortoSM
    double delta_descoberto_restriction_trigger = vars[2];
    double delta_tortoSM_restriction_trigger = vars[3];
    double caesb_descoberto_transfer_trigger = vars[4]; //gatilho para acionar transferência de água para o Descoberto
    double caesb_tortoSM_transfer_trigger = vars[5]; //gatilho para acionar transferência de água para o Descoberto
    double caesb_descoberto_annual_payment = 0; //vars[6]; // pagamento anual ao fundo de contingência. O valor é constante (igual para todo ano).
    double caesb_tortoSM_annual_payment = 0; //vars[7]; // pagamento anual ao fundo de contingência. O valor é constante (igual para todo ano).
    
    double caesb_descoberto_inftrigger = vars[8]; //gatilho para acionar a construção de nova infraestrutura por parte da Companhia Descoberto
    double caesb_tortoSM_inftrigger = vars[9]; //gatilho para acionar a construção de nova infraestrutura por parte da Companhia Torto/SM
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        caesb_descoberto_inftrigger = 1.1;
        caesb_tortoSM_inftrigger = 1.1;
    }
    double ETA_paranoaSul_upgrade1_ranking = vars[10]; // ID 7: Etapa 1 da ETA Paranoá Sul (+0.1 m³/s)
    double ETA_paranoaSul_upgrade2_ranking = vars[11]; // ID 8: Etapa 2 da ETA Paranoá Sul (+0.1 m³/s)
    double ETA_paranoaSul_upgrade3_ranking = vars[12]; // ID 9: Etapa 3 das ETAs Paranoá Sul e Norte (+0.1 m³/s)
    double ETA_corumba_upgrade1_ranking = vars[13]; // ID 5: Etapa 1 da ETA Corumbá (+1.4 m³/s, total 2.8 m³/s)
    double ETA_corumba_upgrade2_ranking = vars[14]; // ID 6: Etapa 2 da ETA Corumbá (+1.2 m³/s, total 4.0 m³/s)
    double descoberto_expansao_ranking = vars[15]; // ID 10: Expansão do Reservatório do Descoberto (+25% storage)
    double caesb_descoberto_inf_buffer = vars[16];
    double caesb_tortoSM_inf_buffer = vars[17];
    // Paranoa is emergency-only for TortoSM. Single trigger variable replaces
    double caesb_paranoa_transfer_trigger = vars[18]; // ROF threshold to activate emergency Paranoa supply to TortoSM
    double sender_demand_protection_factor = vars[19]; // Fraction of sender's unrestricted demand protected before transferring (range 0–2)
    double ETA_santaMaria_upgrade_ranking = vars[20]; // ID 12: Ampliação da ETA Santa Maria (+0.7 m³/s) — TortoSM
    double residential_price_restriction_multiplier = vars[21]; // Multiplicador de tarifa residencial no estágio 3 (caesbPriceRestrictionMultipliers[2][0])

    //ANALISAR POSSIBILIDADE DE INCLUIR O RIO DO SAL COMO OPÇÃO DE AMPLIAÇÃO DA INFRAESTRUTURA DE OFERTA

    //IDENTIFICADOR DE CADA INFRAESTRUTURA FUTURA. Obs: as infraestruturas já existentes devem ser numeradas antes, começando do 0.

    vector<infraRank> caesb_descoberto_infra_order_raw = { // WSS Descoberto: reservatórios Descoberto e Corumbá IV
            // ID 5: ETA Corumbá Etapa 1 (+1.4 m³/s, total 2.8 m³/s)
            // ID 6: ETA Corumbá Etapa 2 (+1.2 m³/s, total 4.0 m³/s)
            // ID 10: Expansão do Descoberto (+25% storage)
            infraRank(5, ETA_corumba_upgrade1_ranking),
            infraRank(6, ETA_corumba_upgrade2_ranking),
            infraRank(10, descoberto_expansao_ranking)
    };

    // TortoSM infra: Santa Maria (ID 12) is triggered by caesb_tortoSM_inftrigger (vars[9]).
    // Paranoa ETA expansions (IDs 7-9) are emergency-only and triggered by
    // caesb_paranoa_transfer_trigger (vars[18]) — a separate, independently optimised threshold.
    vector<infraRank> caesb_tortoSM_infra_order_raw = {
            infraRank(12, ETA_santaMaria_upgrade_ranking), // Santa Maria ETA expansion (+0.7 m³/s) — TEMPORARILY DISABLED
            infraRank(7,  ETA_paranoaSul_upgrade1_ranking),
            infraRank(8,  ETA_paranoaSul_upgrade2_ranking),
            infraRank(9,  ETA_paranoaSul_upgrade3_ranking)
    };

    // GET INFRASTRUCTURE CONSTRUCTION ORDER BASED ON DECISION VARIABLES
    sort(caesb_descoberto_infra_order_raw.begin(),
         caesb_descoberto_infra_order_raw.end(),
         by_xreal());
    sort(caesb_tortoSM_infra_order_raw.begin(),
         caesb_tortoSM_infra_order_raw.end(),
         by_xreal());

    vector<int> rof_triggered_infra_order_caesb_descoberto =
            vecInfraRankToVecInt(caesb_descoberto_infra_order_raw);
    vector<int> rof_triggered_infra_order_caesb_tortoSM =
            vecInfraRankToVecInt(caesb_tortoSM_infra_order_raw);

    // Create vectors with each utility's long-term ROF values assigned to all
    // infrastructure options.
    vector<double> rofs_infra_caesb_descoberto = vector<double>
            (rof_triggered_infra_order_caesb_descoberto.size(),
             caesb_descoberto_inftrigger);
    // Santa Maria (ID 12) uses the standard TortoSM infra trigger (vars[9]).
    // Paranoa ETA expansions (IDs 7-9) use the separate emergency trigger (vars[18]).
    vector<double> rofs_infra_caesb_tortoSM;
    for (const auto& ir : caesb_tortoSM_infra_order_raw) {
        if (ir.id == 12) {
            rofs_infra_caesb_tortoSM.push_back(caesb_tortoSM_inftrigger);
        } else { // IDs 7, 8, 9 — Paranoa ETA expansions
            rofs_infra_caesb_tortoSM.push_back(caesb_paranoa_transfer_trigger);
        }
    }

    // Remove small expansions being built after big expansions that would encompass the small expansions.

    // ==================== SET UP RDM FACTORS ============================
    // RDM factors são os fatores de grande incerteza (DU factors)
    if (wss_rdm.empty()) {
        /// All matrices below have dimensions n_realizations x nr_rdm_factors
        wss_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(5, 1.));
        water_sources_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(27,
                                               1.)); //n corresponde a: 1 + 2 * N° de infraestruturas.
        policies_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(2, 1.));
    }

    // ===================== SET UP PROBLEM COMPONENTS =====================

    //cout << "BEGINNING TRIANGLE TEST" << endl << endl;
    // cout << "Using " << omp_get_num_threads() << " cores." << endl;
    //    cout << getexepath() << endl;

    /// READ STREAMFLOWS

    auto streamflow_n_weeks = (int) std::round(90 *
                                               WEEKS_IN_YEAR); //90 = 50 anos de dados históricos + 40 anos de dados sintéticos (futuros).

    /// In case a vector containing realizations numbers to be calculated is passed, set
    /// number of realizations to number of realizations in that vector.

    //    vector<double> sewageFractions = Utils::parse1DCsvFile(
    //            io_directory + "/InputFiles/sewageFractions.csv");

    //SÉRIES DE EVAPORAÇÃO DE CADA RESERVATÓRIO

    EvaporationSeries evaporation_descoberto(&evap_descoberto,
                                             streamflow_n_weeks);
    EvaporationSeries evaporation_tortoSM(
            &evap_tortoSM,
            streamflow_n_weeks);
    EvaporationSeries evaporation_paranoa(
            &evap_paranoa,
            streamflow_n_weeks);
    EvaporationSeries evaporation_corumba(&evap_corumba,
                                          streamflow_n_weeks); //evaporação obtida por meio da evap do Paranoá (multiplicada por 4.478)

    // CRIAÇÃO DOS VETORES RELACIONADOS ÀS VAZÕES AFLUENTES DE CADA RESERVATÓRIO

    // Descoberto
    Catchment descoberto_inflows(&streamflows_descoberto,
                                 streamflow_n_weeks); //as vazões dos 6 afluentes foram somadas em um único arquivo

    vector<Catchment *> bacia_descoberto;
    bacia_descoberto.push_back(&descoberto_inflows);

    // Santa Maria
    Catchment tortoSM_inflows(&streamflows_tortoSM,
                              streamflow_n_weeks); //são 3 afluentes: Santa Maria, Milho Cozido e Vargem Grande

    vector<Catchment *> bacia_tortoSM;
    bacia_tortoSM.push_back(&tortoSM_inflows);

    //Sistema Bananal e Torto (Reforçam o Santa Maria) - CAP.TOR.001 (Torto)
    Catchment bananal_torto_inflows(&streamflows_bananal_torto,
                                    streamflow_n_weeks); //esse subsistema reforça o Torto/Santa Maria -> a água captada
    //no riberião Bananal e no Torto é bombeada para a ETA Brasília
    vector<Catchment *> sistema_bananal_torto;
    sistema_bananal_torto.push_back(&bananal_torto_inflows);

//    //Captação no Torto (Torto reforça o Santa Maria) - Captação CAP.TOR.001
//    Catchment torto_inflows(&streamflows_torto,
//                            streamflow_n_weeks);
//    vector<Catchment *> sistema_torto;
//    sistema_torto.push_back(&torto_inflows);

    // Lago Paranoá (Obs: o Subsistema Lago Norte capta água no Lago Paranoá, precisamente no braço do Torto)
    Catchment paranoa_inflows(&streamflows_paranoa,
                              streamflow_n_weeks); //são 5 afluentes: Ribeirão do Torto, Bananal, Riacho Fundo, Ribeirão Gama, Cabeça de Veado

    vector<Catchment *> bacia_paranoa;
    bacia_paranoa.push_back(&paranoa_inflows);

    // Corumbá IV
    Catchment corumbaIV_inflows(&streamflows_corumbaIV,
                                streamflow_n_weeks); //são 5 afluentes: Areias, Engenho das Lajes, Alagado, Fazenda Beira e Campo Limpo

    vector<Catchment *> bacia_corumba;
    bacia_corumba.push_back(&corumbaIV_inflows);

    // CURVAS VOLUME X ÁREA DOS RESERVATÓRIOS

    //Curva do Descoberto - baseado na Nota Técnica n° 58/2016 da ADASA (volume útil em hm³)
    vector<double> descoberto_storage = {0,
                                         4.54 * table_gen_storage_multiplier,
                                         9.92 * table_gen_storage_multiplier,
                                         16.21 * table_gen_storage_multiplier,
                                         23.28 * table_gen_storage_multiplier,
                                         31.08 * table_gen_storage_multiplier,
                                         39.81 * table_gen_storage_multiplier,
                                         49.49 * table_gen_storage_multiplier,
                                         60.31 * table_gen_storage_multiplier,
                                         72.29 * table_gen_storage_multiplier};

    vector<double> descoberto_area = {412.94, 494.21, 584.10, 670.80,
                                      740.33, 825.34, 917.11, 1023.49,
                                      1141.23, 1255.34};
                                      //dados da área (hm²) do reservatório do Descoberto (correspondente a cada volume acima)

    //Curva do Santa Maria - baseado na Nota Técnica 61/2016 da ADASA (volume útil em hm³)
    vector<double> tortoSM_storage = {0,
                                      3.821 * table_gen_storage_multiplier,
                                      7.977 * table_gen_storage_multiplier,
                                      12.446 * table_gen_storage_multiplier,
                                      17.245 * table_gen_storage_multiplier,
                                      22.405 * table_gen_storage_multiplier,
                                      27.924 * table_gen_storage_multiplier,
                                      33.809 * table_gen_storage_multiplier,
                                      40.075 * table_gen_storage_multiplier,
                                      46.743 * table_gen_storage_multiplier,
                                      53.830 * table_gen_storage_multiplier,
                                      61.308 * table_gen_storage_multiplier};

    vector<double> tortoSM_area = {365.601, 399.314, 431.316, 462.723,
                                   497.954, 534.384, 569.411, 607.685,
                                   645.696, 687.915, 728.970, 765.141};
                                   //dados da area (hm²) do reservatório de Santa Maria (correspondente a cada volume acima)

    //Curva do Paranoá - baseado na batimetria da CAESB, realizada em 2003 (volume total do lago)
    vector<double> paranoa_storage = {0, 0.0278 * table_gen_storage_multiplier,
                                      0.325 * table_gen_storage_multiplier,
                                      1.581 * table_gen_storage_multiplier,
                                      7.979 * table_gen_storage_multiplier,
                                      26.169 * table_gen_storage_multiplier,
                                      62.251 * table_gen_storage_multiplier,
                                      118.494 * table_gen_storage_multiplier,
                                      193.662 * table_gen_storage_multiplier,
                                      288.837 * table_gen_storage_multiplier,
                                      460.490 * table_gen_storage_multiplier}; //dados do volume (hm³) do reservatório do Paranoá

    vector<double> paranoa_area = {0, 2.757, 15.154, 81.656, 324.785, 763.788,
                                   1321.336, 1874.680, 2429.789, 3006.939,
                                   3881.452}; //dados da area (hm²) do reservatório do Paranoá (correspondente a cada volume acima)


    DataSeries descoberto_storage_area(&descoberto_storage,
                                       &descoberto_area); //aqui ele está juntando o vetor storage e o vetor area criado para cada reservatório acima
    DataSeries tortoSM_storage_area(&tortoSM_storage,
                                    &tortoSM_area);
    DataSeries paranoa_storage_area(&paranoa_storage,
                                    &paranoa_area);


    /// REGRAS RELACIONADAS ÀS VAZÕES REMANESCENTES (RESTRIÇÕES AMBIENTAIS)
    // Valores de vazão em hm³/semana

    // Vazão remanescente do Descoberto - baseado no artigo de Rocha e Cézar (2015)
    FixedMinEnvFlowControl descoberto_min_env_control(0,
                                                      0.6e-6 * 3600 * 24 * 7);

    // Vazão remanescente do Torto/Santa Maria - não tem, apenas verte a água
    FixedMinEnvFlowControl tortoSM_min_env_control(1, 0);

    // Vazão remanescente de Corumbá IV - baseado no EIA da UHE de Corumbá IV e na ação civil pública
    FixedMinEnvFlowControl corumba_min_env_control(2, 5.3e-6 * 3600 * 24 * 7);

    // Vazão remanescente do Paranoá - baseado na resolução n° 33/2018 da ADASA
    vector<int> paranoa_weeks = {0, 18, 44,
                                 53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    vector<double> paranoa_releases = {(1.2e-6 * 3600 * 24 * 7),
                                       (0.7e-6 * 3600 * 24 * 7),
                                       (1.2e-6 * 3600 * 24 * 7)}; // mínimo de 0,7 m³/s no período de estiagem e de 1,2 m³/s no período chuvoso

    SeasonalMinEnvFlowControl paranoa_min_env_control(3, paranoa_weeks,
                                                      paranoa_releases);

    //Vazão remanescente do Ribeirão Bananal e do Torto - combinação entre as duas water sources
    vector<int> bananal_torto_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44,
                                       48,53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    
    //A vazão remanescente da water source integrada Bananal + Torto corresponde à soma das vazões remanescentes de cada um
    vector<double> bananal_torto_releases = {(0.832e-6 * 3600 * 24 * 7),
                                             (1.006e-6 * 3600 * 24 * 7),
                                             (0.924e-6 * 3600 * 24 * 7),
                                             (0.874e-6 * 3600 * 24 * 7),
                                             (0.666e-6 * 3600 * 24 * 7),
                                             (0.498e-6 * 3600 * 24 * 7),
                                             (0.432e-6 * 3600 * 24 * 7),
                                             (0.376e-6 * 3600 * 24 * 7),
                                             (0.338e-6 * 3600 * 24 * 7),
                                             (0.358e-6 * 3600 * 24 * 7),
                                             (0.484e-6 * 3600 * 24 * 7),
                                             (0.682e-6 * 3600 * 24 * 7)};
    SeasonalMinEnvFlowControl bananal_torto_min_env_control(4,
                                                            bananal_torto_weeks,
                                                            bananal_torto_releases);

//    // Vazão remanescente do Ribeirão Bananal - baseado no ato de outorga da ADASA para a captação no Bananal
//    vector<int> bananal_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44, 48,
//                                 53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
//    vector<double> bananal_releases = {(0.422e-6 * 3600 * 24 * 7), (0.446e-6 * 3600 * 24 * 7), (0.474e-6 * 3600 * 24 * 7),
//                                       (0.464e-6 * 3600 * 24 * 7), (0.406e-6 * 3600 * 24 * 7), (0.358e-6 * 3600 * 24 * 7),
//                                       (0.322e-6 * 3600 * 24 * 7), (0.286e-6 * 3600 * 24 * 7), (0.258e-6 * 3600 * 24 * 7),
//                                       (0.258e-6 * 3600 * 24 * 7), (0.314e-6 * 3600 * 24 * 7), (0.372e-6 * 3600 * 24 * 7)};
//    SeasonalMinEnvFlowControl bananal_min_env_control(4, bananal_weeks,
//                                                      bananal_releases);
//
//    // Vazão remanescente do Ribeirão do Torto - dado fornecido pela Caesb - exigência do PGIRH/DF (os valores correspondem a 20% da vazão média mínima mensal)
//    vector<int> torto_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44, 48,
//                               53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
//    vector<double> torto_releases = {(0.410e-6 * 3600 * 24 * 7), (0.560e-6 * 3600 * 24 * 7), (0.450e-6 * 3600 * 24 * 7),
//                                     (0.410e-6 * 3600 * 24 * 7), (0.260e-6 * 3600 * 24 * 7), (0.140e-6 * 3600 * 24 * 7),
//                                     (0.110e-6 * 3600 * 24 * 7), (0.90e-6 * 3600 * 24 * 7), (0.80e-6 * 3600 * 24 * 7),
//                                     (0.100e-6 * 3600 * 24 * 7), (0.170e-6 * 3600 * 24 * 7), (0.310e-6 * 3600 * 24 * 7)};
//    SeasonalMinEnvFlowControl torto_min_env_control(5, torto_weeks,
//                                                      torto_releases);

    //FixedMinEnvFlowControl captacoes_gama();

    vector<MinEnvFlowControl *> min_env_flow_controls; //criação do vetor min_env_flow_controls
    min_env_flow_controls.push_back(
            &descoberto_min_env_control); //todos esses elementos (tortoSM, paranoa..) estão sendo jogados para dentro do vetor min_env_flow_controls
    min_env_flow_controls.push_back(&tortoSM_min_env_control);
    min_env_flow_controls.push_back(&paranoa_min_env_control);
    min_env_flow_controls.push_back(&corumba_min_env_control);
    min_env_flow_controls.push_back(&bananal_torto_min_env_control);
    //min_env_flow_controls.push_back(&torto_min_env_control);
    //min_env_flow_controls.push_back(&captacoes_gama_min_env_control);


    // CRIAÇÃO DOS RESERVATÓRIOS E VETORES CORRESPONDENTES
    vector<double> construction_time_interval = {3.0,
                                                 5.0}; //o período de construção das novas infraestruturas varia entre 3 e 5 anos

//     // ========================================================
//     // =========TESTING ONLY=======: force all infrastructure online from week 0
//     // Remove these lines to restore normal behaviour
//     construction_time_interval = {0.0, 0.0};
//     caesb_descoberto_inftrigger = -1.0;
//     caesb_tortoSM_inftrigger    = -1.0;
//     // Also override the already-created rofs vectors (they were built before the trigger was overridden above)
//     fill(rofs_infra_caesb_descoberto.begin(), rofs_infra_caesb_descoberto.end(), -1.0);
//     fill(rofs_infra_caesb_tortoSM.begin(),    rofs_infra_caesb_tortoSM.end(),    -1.0);
//     // ========================================================
//     // ========================================================

    vector<double> city_infrastructure_rof_triggers = {
            caesb_descoberto_inftrigger,
            caesb_tortoSM_inftrigger};
    //todas essas são variáveis de decisão definidas lá no começo do código, gatilho para acionar a construção de nova infraestrutura

    // EXISTING SOURCES //descrição das fontes de água já existentes

    Reservoir descoberto("Descoberto",//nome do reservatório
                         0,//número de identificação
                         bacia_descoberto,//vetor criado lá em cima - contém as vazões referentes aos afluentes do Descoberto
                         72.29 * table_gen_storage_multiplier, //capacidade de armazenamento do reservatório (hm³)
                         6.0e-6 * 3600 * 24 * 7, //capacidade máxima de tratamento da ETA Descoberto (hm³/semana)
                         evaporation_descoberto,                //obs: outorga da represa de Sta Maria é de 1.478 l/s (PDSB, 2017)
                         &descoberto_storage_area);

    Reservoir tortoSM("Torto / Santa Maria", 1,
                      bacia_tortoSM,
                      61.308 * table_gen_storage_multiplier, //hm³
                      1.1e-6 * 3600 * 24 * 7, //capacidade máxima de tratamento da ETA Brasília (hm³/semana)
                      evaporation_tortoSM,
                      &tortoSM_storage_area);

    // Reservatório de Corumbá IV
    double cIV_supply_caesb_capacity = 12.849 *
                                       table_gen_storage_multiplier; //valor obtido baseado na proporção das vazões retiradas
    // (conferir arquivo excel)
    double cIV_supply_saneago_capacity = 12.849 * table_gen_storage_multiplier;
    double cIV_energy_capacity = 745.702 * table_gen_storage_multiplier;
    double cIV_wq_capacity = 2936.6 * table_gen_storage_multiplier;
    double cIV_storage_capacity = cIV_wq_capacity + cIV_energy_capacity +
                                  cIV_supply_saneago_capacity +
                                  cIV_supply_caesb_capacity; //o armazenamento de água total é igual a soma da parte destinada à abastecimento,
    // destinada à energia e da parte destinada à preservação ambiental

    //Curva de Corumbá IV - baseado nos dados do portal da ANA (volume útil)
    vector<double> corumba_storage = {0,
                                      2936.6 * table_gen_storage_multiplier, //hm³
                                      cIV_storage_capacity *
                                      table_gen_storage_multiplier};
    vector<double> corumba_area = {0, 13712.0, //hm²
                                   17330.0};
    DataSeries corumba_storage_area(&corumba_storage, &corumba_area);

    // Corumba IV allocation: storage is shared as:
    // CAESB (Descoberto WSS) = 12.849 hm³ (0.35% of 3708 hm³)
    // Saneago = 12.849 hm³ (0.35%)
    // Energy = 745.702 hm³ (20.1%)
    // Water Quality = 2936.6 hm³ (79.2%)
    vector<int> cIV_allocations_ids = {0,
                                       WATER_QUALITY_ALLOCATION};
    vector<double> cIV_allocation_fractions = {
            cIV_supply_caesb_capacity / cIV_storage_capacity,  // 0.00346 = 12.849/3708
            (cIV_wq_capacity + cIV_supply_saneago_capacity +
             cIV_energy_capacity) / cIV_storage_capacity}; 
    vector<double> cIV_treatment_allocation_fractions = {
            1.0};  // Descoberto WSS treats 100% of water withdrawn from Corumba

    // Initial treatment capacity for Corumba IV (online from start)
    // Using the first stage ETA capacity (1.4 m³/s = 1.4e-6 * 3600 * 24 * 7 hm³/week)
    double cIV_initial_treatment_capacity = 1.4e-6 * 3600 * 24 * 7;

    AllocatedReservoir corumba("Corumba IV",
                               2,
                               bacia_corumba,
                               cIV_storage_capacity,  // Total: 3708 hm³
                               cIV_initial_treatment_capacity,  // Initial: 1.4 m³/s
                               evaporation_corumba,
                               15000., //&corumba_storage_area,
                               &cIV_allocations_ids,
                               &cIV_allocation_fractions,
                               &cIV_treatment_allocation_fractions);

    // Lago Paranoá parameters: emergency source for TortoSM only.
    // TortoSM may use 10% of Paranoá's 8% supply = 0.8% of total capacity as
    // an emergency reserve. Transfer is activated by EmergencyTransferParanoa
    // policy when TortoSM ROF > caesb_paranoa_transfer_trigger.
    double lp_total_capacity = 460.490 * table_gen_storage_multiplier;
    double lp_supply_wss1 = 0.008 * lp_total_capacity;  // 0.8% emergency reserve for TortoSM
    double lp_wq_capacity = lp_total_capacity - lp_supply_wss1; // remaining 99.2% water quality
    vector<int> lp_allocations_ids = {1, WATER_QUALITY_ALLOCATION};  // TortoSM only; no Descoberto
    vector<double> lp_allocation_fractions = {
            lp_supply_wss1 / lp_total_capacity,   // WSS 1 (TortoSM): 0.8%
            lp_wq_capacity / lp_total_capacity};   // water quality pool: 99.2%
    vector<double> lp_treatment_allocation_fractions = {1.0};  // 100% of treatment capacity to WSS 1

    AllocatedReservoir paranoa("Lago Paranoa", 3,
                               bacia_paranoa,
                               lp_total_capacity,
                               0.1e-6 * 3600 * 24 * 7, // emergency pipe capacity to TortoSM (hm³/week)
                               evaporation_paranoa,
                               &paranoa_storage_area,
                               &lp_allocations_ids,
                               &lp_allocation_fractions,
                               &lp_treatment_allocation_fractions);

    Intake ribeirao_bananal_torto(
            "Captacao no Ribeirao Bananal e Ribeiro do Torto",
            4,                                                          //Obs: a série de vazão utilizada referente ao Bananal
            sistema_bananal_torto,                                          //foi retirada de uma estação fluviométrica localizada
            1.7e-6 * 3600 * 24 * 7); // hm³/semana      //a justante do ponto de captação. Não há problema,
    // pois a captação começou apenas ao final de 2017,
    // então a série é basicamente composta pela vazão natural do ribeirão.
//    Intake ribeirao_torto("Captacao no Ribeirao do Torto",
//                            5,
//                            sistema_torto,
//                            0.95e-6 * 3600 * 24 * 7); //950 m³/s corresponde à captação média no Torto.


    BalloonPaymentBond dummy_bond(11, 0., 1, 1., vector<int>(1, 0));
    Reservoir dummy_endpoint("Dummy Node", 11, vector<Catchment *>(), 1., 0,
                             evaporation_corumba, 1,
                             construction_time_interval, 0, dummy_bond);

    // OPÇÕES DE INFRAESTRUTURAS FUTURAS - Descrição das opções futuras de infraestrutura para ampliar a oferta de água

    // The capacities listed here for expansions are what additional capacity is gained relative to existing capacity,
    // NOT the total capacity after the expansion is complete. For infrastructure with a high and low option, this means
    // the capacity for both is relative to current conditions - if Lake Michie is expanded small it will gain 2.5BG,
    // and a high expansion will provide 7.7BG more water than current. if small expansion happens, followed by a large
    // expansion, the maximum capacity through expansion is 7.7BG, NOT 2.5BG + 7.7BG.

    //OBS: Assumiu-se que todas as infraestruturas já tinham licença para serem construídas (permitting period = 0)
    /* Tipos de empréstimos definidos no WaterPaths:
        - Level Debt Service
        - Ballon Payment
        - Variable-Interest Bonds */

    //Captação em Corumbá IV - Initial capacity: 1.4 m³/s (online from start)
    //                        - Etapa 1 (ID 5): + 1.4 m³/s (total 2.8 m³/s)
    //                        - Etapa 2 (ID 6): + 1.2 m³/s (total 4.0 m³/s)
    // Note: Original Etapa numbering shifted - former Etapa 2/3 are now Etapa 1/2
    
    vector<double> capacity_ETA_corumba_upgrade_1 = {1.4e-6 * 3600 * 24 * 7,
                                                     0}; // Etapa 1: ampliação da capacidade de produção (+ 1.4 m³/s) - TOTAL: 2.8 m³/s
    vector<double> capacity_ETA_corumba_upgrade_2 = {1.2e-6 * 3600 * 24 * 7,
                                                     0}; // Etapa 2: ampliação da capacidade de produção (+ 1.2 m³/s) - TOTAL: 4.0 m³/s

    // Empréstimos para Expansão da ETA Corumbá (Sistema Corumbá)
    
    vector<Bond *> debendure_expansao_ETA_corumba_1 = {
            new LevelDebtServiceBond(7, 222066142.8, 15, 0.12,
                                     vector<int>(1, 0)),
            new BalloonPaymentBond(12, 0, 15, 0.12, vector<int>(1, 0))};
    vector<Bond *> debendure_expansao_ETA_corumba_2 = {
            new LevelDebtServiceBond(8, 251383400, 15, 0.12,
                                     vector<int>(1, 0)),
            new BalloonPaymentBond(12, 0, 15, 0.12, vector<int>(1, 0))};


    // Corumbá expansions maintain sequential order: Etapa 1 (ID 5) must be built before Etapa 2 (ID 6)
    // Note: Corumbá starts with 1.4 m³/s treatment capacity online

    BalloonPaymentBond no_bond(5, 0., 1, 1, vector<int>(1, 0));
    SequentialJointTreatmentExpansion ETA_corumba_etapa1(
            "Etapa 1 de Corumba IV", 5, 2, 0, {5, 6},
            capacity_ETA_corumba_upgrade_1,
            debendure_expansao_ETA_corumba_1, construction_time_interval,
            0 * WEEKS_IN_YEAR); // Previsão: 2030-2033. ID 5

    SequentialJointTreatmentExpansion ETA_corumba_etapa2(
            "Etapa 2 de Corumba IV", 6, 2, 1, {5, 6},
            capacity_ETA_corumba_upgrade_2,
            debendure_expansao_ETA_corumba_2, construction_time_interval,
            0 * WEEKS_IN_YEAR);  //5 * WEEKS_IN_YEAR); // Previsão: depois de 2037. ID 6


    //Sistema Paranoá - Construção da ETA Paranoá Sul (0.7 m³/s), sua primeira ampliação (upgrade 2, com + 0.7 m³/s), segunda ampliação (upgrade 3, com + 0.35 m³/s) e
    // ampliação da ETA Lago Norte (upgrade 3, com + 0.35 m³/s))

    // Paranoa WTP expansions (IDs 7-9) each add +100 l/s to treatment capacity
    // AND +100 l/s to the transfer pipeline capacity (via EmergencyTransferParanoa).
    // Initial treatment capacity is 100 l/s (set on the AllocatedReservoir), and the
    // pipeline starts at 100 l/s. Each triggered expansion raises both by 100 l/s.
    //
    double paranoa_expansion_capacity = 0.1e-6 * 3600 * 24 * 7; // 100 l/s in hm³/week
    vector<double> capacities_ETA_paranoaSul_upgrade_1 = {0.0, paranoa_expansion_capacity};
    vector<double> capacities_ETA_paranoaSul_upgrade_2 = {0.0, paranoa_expansion_capacity};
    vector<double> capacities_ETAs_paranoa_upgrade_3   = {0.0, paranoa_expansion_capacity};

    // Empréstimos para a implantação e ampliação da ETA Paranoá Sul e Norte (Norte: apenas upgrade 3)

    vector<Bond *> debendure_expansao_ETA_paranoa_1 = {
            new BalloonPaymentBond(11, 0, 15, 0.12, vector<int>(1, 0)),
            new LevelDebtServiceBond(8, 60332016, 15, 0.12,
                                     vector<int>(1, 0))}; //bond term = 25 anos
    vector<Bond *> debendure_expansao_ETA_paranoa_2 = {
            new BalloonPaymentBond(11, 0, 15, 0.12, vector<int>(1, 0)),
            new LevelDebtServiceBond(9, 60332016, 15, 0.12,
                                     vector<int>(1, 0))};
    vector<Bond *> debendure_expansao_ETA_paranoa_3 = {
            new BalloonPaymentBond(11, 0, 15, 0.12, vector<int>(1, 0)),
            new LevelDebtServiceBond(10, 60332016, 15, 0.12,
                                     vector<int>(1, 0))};


    SequentialJointTreatmentExpansion etapa1_ETA_paranoaSul(
            "Etapa 1 da ETA Paranoa Sul ", 7, 3, 0, {7, 8, 9},
            capacities_ETA_paranoaSul_upgrade_1,
            debendure_expansao_ETA_paranoa_1, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2020. ID 7

    SequentialJointTreatmentExpansion etapa2_ETA_paranoaSul(
            "Etapa 2 da ETA Paranoa Sul", 8, 3, 1, {7, 8, 9},
            capacities_ETA_paranoaSul_upgrade_2,
            debendure_expansao_ETA_paranoa_2, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2022. ID 8

    SequentialJointTreatmentExpansion etapa3_ETAs_paranoa(
            "Etapa 3 da ETA Paranoa Sul e Norte", 9, 3, 2, {7, 8, 9},
            capacities_ETAs_paranoa_upgrade_3,
            debendure_expansao_ETA_paranoa_3, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2034 a 2037. ID changed to 9

    // Expansão do Reservatório do Descoberto

    LevelDebtServiceBond descoberto_exp_bond(10, 7541502, 15, 0.12,
                                             vector<int>(1,
                                                         0)); //Elevação do nível d'água da barragem do descoberto (aumento da capacidade de armazenamento em 25%)
    ReservoirExpansion descoberto_expansion(
            "Expansao da capacidade de armazenamento do Descoberto", 10, 0,
            25.216, //25.216 = aumento em hm³ da capacidade de armazenamento do Descoberto (25%)
            construction_time_interval,
            0 * WEEKS_IN_YEAR,
            descoberto_exp_bond); //previsão: 2022. ID 10

    // Ampliação da ETA Santa Maria (WSS 1 — TortoSM): +0.7 m³/s treatment capacity. ID 12.
    // Cost is a placeholder value — adjust based on project-specific estimates.
    vector<double> capacity_ETA_santaMaria = {0.0, 0.7e-6 * 3600 * 24 * 7}; // WSS0: 0, WSS1: +0.7 m³/s
    vector<Bond *> debendure_expansao_ETA_santaMaria = {
            new BalloonPaymentBond(11, 0, 15, 0.12, vector<int>(1, 0)),            // WSS0 — no cost
            new LevelDebtServiceBond(13, 60332016, 15, 0.12, vector<int>(1, 0))};  // WSS1 — same as Etapa 1 da ETA Paranoá Sul
    SequentialJointTreatmentExpansion ETA_santaMaria_expansion(
            "Ampliacao ETA Santa Maria", 12, 1, 0, {12},
            capacity_ETA_santaMaria,
            debendure_expansao_ETA_santaMaria, construction_time_interval,
            0 * WEEKS_IN_YEAR); // ID 12 — standalone expansion for TortoSM

    vector<WaterSource *> water_sources; //water_sources é um vetor comum, que comportará todas as
    // opções descritas acima de ampliação da infraestrutura de abastecimento
    water_sources.push_back(&descoberto);
    water_sources.push_back(&tortoSM);
    water_sources.push_back(&corumba);
    water_sources.push_back(&paranoa);
    water_sources.push_back(&ribeirao_bananal_torto);
//    water_sources.push_back(&ribeirao_bananal);
//    water_sources.push_back(&ribeirao_torto);

    // Corumbá now starts with 1.4 m³/s treatment capacity online
    water_sources.push_back(&ETA_corumba_etapa1);  // ID 5 at index 5
    water_sources.push_back(&ETA_corumba_etapa2);  // ID 6 at index 6
    water_sources.push_back(&etapa1_ETA_paranoaSul);
    water_sources.push_back(&etapa2_ETA_paranoaSul);
    water_sources.push_back(&etapa3_ETAs_paranoa);
    water_sources.push_back(&descoberto_expansion);
    water_sources.push_back(&ETA_santaMaria_expansion);  // ID 12 at index 12

    water_sources.push_back(&dummy_endpoint);

    //water_sources.push_back(&dummy_endpoint);

    /*
     * System connection diagram (water
     * flows from top to bottom)
     * Potential projects and expansions
     * of existing sources in parentheses
     *
     *                                  1 (Santa Maria)
     *      0(10)                       |
     *       \                          |
     *        \                         |  4(Ban)   4 (Torto) (12)
     *         \                        |   \      /
     *          \                       |    \ __/
     *           \                      |     3 (Paranoá) (7, 8, 9)
     *            \                     |     |
     *             \                    |    |
     *              \                   |   |
     *               \                  |  |
     *                \                 |_|
     *                 \               /
     *                  \             /
     *                   \           /
     *                2(5, 6)_ _    /
     *                          \  /
     *                           \/
     *                           11 (dummy)
     *
     */

    Graph g(13); //graph é uma forma de estruturar dados que consiste em dois componentes:
    // um conjunto finito de vértices denominado de nós; e um cojunto finito de pares ordenados (x,y), denominados de edges.
    g.addEdge(0,
              2); //essa conexão indica que existe uma aresta do vértice 0 ao vértice 2.
    g.addEdge(2, 11);
    g.addEdge(1, 11);
    g.addEdge(4, 3);
    //g.addEdge(11, 3);
    g.addEdge(3, 11);

    auto demand_n_weeks = (int) std::round(40 *
                                           WEEKS_IN_YEAR); //40 é o número de anos a serem simulados. A fç auto serve para
    // declarar variáveis cujo tipo vai ser inferido pelo compilador a partir da inicialização delas

    // Criação dos vetores de vazão de retorno, relacionada aos efluentes lançados em corpos hídricos utilizados também para abastecimento urbano.

    vector<int> caesb_descoberto_ws_return_id = {
            3}; //ID do reservatório onde o efluente será lançado (montante desse reservatório) - 3 (ID do Paranoá)
    WwtpDischargeRule wwtp_discharge_caesb_descoberto( //em demand_to_waste_water_fraction, deve-se colocar a % de agua tratada (que sai da ETA - demanda) que vai voltar como vazão efluente para cada corpo receptor
            demand_to_wastewater_fraction_caesb_descoberto, // n° de séries (de demanda de vazão efluente) correspondentes ao n° de reservatórios onde os efluentes são lançados. Obs: se no meu estudo de caso, apenas o Paranoá for receptor, esse arquivo csv terá uma linha e 53 colunas (vazão efluente de cada semana)
            caesb_descoberto_ws_return_id); // id de cada reservatório/corpo receptor de efluente

    vector<int> caesb_tortoSM_ws_return_id = {
            3}; //ID do reservatório onde o efluente será lançado (montante desse reservatório) - 3 (ID do Paranoá)
    WwtpDischargeRule wwtp_discharge_caesb_tortoSM(
            demand_to_wastewater_fraction_caesb_tortoSM,
            caesb_tortoSM_ws_return_id); // id de cada reservatório/corpo receptor de efluente

    //Criação das companhias de água. A descrição de cada termo está no arquivo .doc.

    // Descoberto WSS: WTP 0 (Descoberto reservoir) + WTP 1 (Corumba reservoir).
    // TortoSM WSS: WTP 0 only (TortoSM + Bananal/Torto intakes).
    // Paranoa has NO WTP slot in TortoSM — all Paranoa→TortoSM flow goes via the
    // fixed transfer pipeline (EmergencyTransferParanoa, 0.1 m³/s cap).
    vector<vector<int>> water_sources_to_wtp_caesb_1 = {{0},   // WTP 0 treats Descoberto
                                                        {2}};  // WTP 1 treats Corumba
    vector<double> wtp_capacities_caesb_1 = {6.0e-6 * 3600 * 24 * 7,    // WTP 0: Descoberto ETA (6.0 m³/s)
                                             1.4e-6 * 3600 * 24 * 7};   // WTP 1: Corumba ETA (1.4 m³/s)
    vector<vector<int>> water_sources_to_wtp_caesb_2 = {{1, 4}, {3}};  // WTP 0: TortoSM + Bananal/Torto; WTP 1: Paranoá
    vector<double> wtp_capacities_caesb_2 = {
            (1.1e-6 * 3600 * 24 * 7 + 1.7e-6 * 3600 * 24 * 7),  // WTP 0: TortoSM + Bananal ETAs (2.8 m³/s)
            0.1e-6 * 3600 * 24 * 7};  // WTP 1: ETA Lago Norte / Paranoá (0.1 m³/s)

    // Create single CAESB utility with two water supply systems
    // First system: Descoberto (system_id=0, utility_id=0)
    // Second system: TortoSM (system_id=1, utility_id=0)

    Utility caesb((char *) "CAESB", 0,
                  demand_caesb_descoberto, // Placeholder, will be managed by WSS
                  demand_n_weeks,
                  0.0,  // Placeholder - each WSS has its own contingency %
                  vector<vector<double>>(), // Empty - each WSS has its own demand fractions
                  vector<vector<double>>(), // Empty - each WSS has its own water prices
                  wwtp_discharge_caesb_descoberto, // Placeholder
                  caesb_descoberto_inf_buffer, // Placeholder
                  {}, {}, // Empty vectors for water_source_to_wtp and utility_owned_wtp_capacities
                  0.04, // Infrastructure discount rate (taxa de desconto)
                  15, 0.12); // bond_term and bond_interest_rate to match Original model

    // Add Descoberto water supply system (system_id=0)
    caesb.addWaterSupplySystem("Descoberto", 0, 0,
                               demand_caesb_descoberto, demand_n_weeks,
                               wwtp_discharge_caesb_descoberto,
                               caesb_descoberto_inf_buffer,
                               water_sources_to_wtp_caesb_1,
                               wtp_capacities_caesb_1,
                               caesb_descoberto_annual_payment,
                               caesbDescobertoDemandClassesFractions,
                               caesbDescobertoUserClassesWaterPrices);

    // Add TortoSM water supply system (system_id=1)  
    caesb.addWaterSupplySystem("TortoSM", 1, 0,
                               demand_caesb_tortoSM, demand_n_weeks,
                               wwtp_discharge_caesb_tortoSM,
                               caesb_tortoSM_inf_buffer,
                               water_sources_to_wtp_caesb_2,
                               wtp_capacities_caesb_2,
                               caesb_tortoSM_annual_payment,
                               caesbTortoSMDemandClassesFractions,
                               caesbTortoSMUserClassesWaterPrices);

    // Set infrastructure parameters for each WSS (this is where WSS gets operational control)
    // WSS_0 (Descoberto) gets Descoberto infrastructure
    auto& wss_descoberto = caesb.systemById(0);
    wss_descoberto.setInfrastructureParameters(rof_triggered_infra_order_caesb_descoberto,
                                               vector<int>(),  // Empty demand order
                                               rofs_infra_caesb_descoberto);

    // WSS_1 (TortoSM) gets TortoSM infrastructure  
    auto& wss_tortoSM = caesb.systemById(1);
    wss_tortoSM.setInfrastructureParameters(rof_triggered_infra_order_caesb_tortoSM,
                                           vector<int>(),  // Empty demand order
                                           rofs_infra_caesb_tortoSM);

    // Set average monthly income for affordability index calculation
    // Descoberto WSS (system_id = 0): average income per household = 1472 BRL/month
        wss_descoberto.setAverageMonthlyIncome(1472.0);
        // Torto/Santa Maria WSS (system_id = 1): average income per household = 4785 BRL/month
        wss_tortoSM.setAverageMonthlyIncome(4785.0);

        // Set initial households for per-household affordability scaling
        wss_descoberto.setInitialHouseholds(687805); 
        wss_tortoSM.setInitialHouseholds(207787); 

        // Set per-reservoir failure thresholds for reliability calculation
        // WSS 0 (Descoberto): existing reservoirs
        wss_descoberto.setSourceFailureThreshold(0, 0.20);  // Descoberto
        wss_descoberto.setSourceFailureThreshold(2, 0.20);  // Corumba
        // WSS 0 (Descoberto): infrastructure expansions
        wss_descoberto.setSourceFailureThreshold(5, 0.20);  // Corumba Etapa1
        wss_descoberto.setSourceFailureThreshold(6, 0.20);  // Corumba Etapa2
        // Note: Descoberto Expansion (id=10) is NOT checked independently.
        // Its capacity is added to parent Descoberto (id=0) when built, so the
        // parent's failure check already covers the expanded capacity.
        wss_descoberto.setSourceFailureThreshold(11, 0.20); // Dummy
        // WSS 1 (TortoSM): existing reservoirs
        wss_tortoSM.setSourceFailureThreshold(1, 0.20);     // SantaMaria/Torto
        wss_tortoSM.setSourceFailureThreshold(3, 0.20);     // Paranoa (now a ROF-triggered supply source for TortoSM)
        // WSS 1 (TortoSM): infrastructure expansions
        wss_tortoSM.setSourceFailureThreshold(4, 0.20);     // Bananal/Torto (Intake)
        wss_tortoSM.setSourceFailureThreshold(7, 0.20);     // Paranoa ETA Etapa 1
        wss_tortoSM.setSourceFailureThreshold(8, 0.20);     // Paranoa ETA Etapa 2
        wss_tortoSM.setSourceFailureThreshold(9, 0.20);     // Paranoa ETA Etapa 3
        wss_tortoSM.setSourceFailureThreshold(12, 0.20);    // Santa Maria ETA expansion

        vector<Utility *> utilities; //vetor que contém a única companhia CAESB
        utilities.push_back(&caesb);

    // Water-source-WSS connectivity matrix:
    // WSS 0 (Descoberto): sources {0, 2, 5, 6, 10, 11}
    // WSS 1 (TortoSM):    sources {1, 4, 3, 7, 8, 9} — Paranoa and its ETAs now ROF-triggered by TortoSM
    vector<vector<int>> reservoir_wss_connectivity_matrix = {
            {0, 2, 5, 6, 10, 11},        // Descoberto(0), Corumba(2), Corumba_E1(5), Corumba_E2(6), Descoberto_Exp(10), Dummy(11)
            {1, 4, 3, 7, 8, 9, 12}       // TortoSM(1), Bananal/Torto(4), Paranoa(3), Paranoa_E1(7), Paranoa_E2(8), Paranoa_E3(9), SantaMaria_Exp(12)
    };

//    @TODO: verificar se há necessidade de corrigir volumes de reservatórios construídos.
//    // O que table_storage_shift representa? O que são esses números (2000, 5000...) [3] [17]
    // Update table storage shift to match WSS structure (2 WSS within single utility)
    // Vector size matches water source count: 13 sources (IDs 0-12, including ETA_santaMaria ID 12)
    auto table_storage_shift = vector<vector<double>>(2,  // Two WSS
                                                      vector<double>(13, 0.));
    // No storage shift needed for current configuration
//    table_storage_shift[3][17] = 2000.; //tem a ver com a RdF
//    table_storage_shift[3][8] = 5000.;
//    table_storage_shift[1][14] = 100.;
//    table_storage_shift[1][20] = 500.;
//    table_storage_shift[1][21] = 500.;
//    table_storage_shift[1][15] = 700.;
//    table_storage_shift[1][9] = 700.;

    vector<DroughtMitigationPolicy *> drought_mitigation_policies; //criação do vetor das políticas de mitigação de seca (racionamento + transferência)

    // POLÍTICA DE RESTRIÇÃO

    vector<double> initial_restriction_triggers = {
            caesb_descoberto_restriction_trigger,
            caesb_tortoSM_restriction_trigger};
    //variáveis de decisão que representam os gatilhos para acionar o racionamento em cada companhia

    vector<double> restriction_stage_multipliers_caesb_descoberto = {0.98, 0.96,
                                                                     0.90}; //São 3 estágios de racionamento. Os fatores 0.98, 0.96, 0.90 são as restrições da demanda.
    // 0.98 significa que a demanda será restringida em 2% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_descoberto = {
            caesb_descoberto_restriction_trigger, //estágio 1 - 2% de redução (campanha de conscientização)
            caesb_descoberto_restriction_trigger +
            delta_descoberto_restriction_trigger, // estágio 2 - 4% de redução (tarifa de conting) + 2% do primeiro
            caesb_descoberto_restriction_trigger + 2. *
                                                   delta_descoberto_restriction_trigger}; //estágio 3 - 12% de redução (racionamento) + 2% do primeiro
    //Obs: o 1° estágio corresponde à campanha, o 2° à tarifa de contingência e o 3° ao racionamento.
    //Acumulou-se a campanha com os outros dois estágios.
    //Como a tarifa de contingência é uma medida bastante penosa, bem como o racionamento,
    //optou-se por não acumular ambas medidas.

    vector<double> restriction_stage_multipliers_caesb_tortoSM = {0.98, 0.96,
                                                                  0.90};

    vector<double> restriction_stage_triggers_caesb_tortoSM = {
            caesb_tortoSM_restriction_trigger, //estágio 1
            caesb_tortoSM_restriction_trigger +
            delta_tortoSM_restriction_trigger, // estágio 2
            caesb_tortoSM_restriction_trigger +
            2. * delta_tortoSM_restriction_trigger}; //estágio 3

    // Criação das políticas de restrição de uso da água (agora referenciando sistemas de abastecimento)
    // Restrictions now target specific water supply systems within the single utility

    // Override residential price restriction multiplier (row 2, col 0) with DV[21]
    caesbPriceRestrictionMultipliers[2][0] = residential_price_restriction_multiplier;

    Restrictions restrictions_descoberto(0,  // utility_id=0, targets system_id=0 (Descoberto)
                                         restriction_stage_multipliers_caesb_descoberto,
                                         restriction_stage_triggers_caesb_descoberto,
                                         &caesbDescobertoDemandClassesFractions,
                                         &caesbDescobertoUserClassesWaterPrices,
                                         &caesbPriceRestrictionMultipliers);

    Restrictions restrictions_tortoSM(1,  // utility_id=0, targets system_id=1 (TortoSM)
                                      restriction_stage_multipliers_caesb_tortoSM,
                                      restriction_stage_triggers_caesb_tortoSM,
                                      &caesbTortoSMDemandClassesFractions,
                                      &caesbTortoSMUserClassesWaterPrices,
                                      &caesbPriceRestrictionMultipliers);

    drought_mitigation_policies = {&restrictions_descoberto,
                                   &restrictions_tortoSM};

    // POLÍTICA DE TRANSFERÊNCIA
//
//    Graph transfer_graph_tortoSM_descoberto(2);
//    transfer_graph_tortoSM_descoberto.addEdge(1,
//                                              0); // Água do tortoSM para o Descoberto
//    Transfers transfer_tortoSM_descoberto(0, 1, 1, 0.1, {0},
//                                          {0.7e-6 * 3600 * 24 * 7}, //hm³/semana
//                                          {caesb_descoberto_transfer_trigger}, //TortoSM transfere até 0.7 m³/s para o Descoberto
//                                          transfer_graph_tortoSM_descoberto,
//                                          vector<double>(), vector<int>());
//
//    Graph transfer_graph_descoberto_tortoSM(2);
//    transfer_graph_descoberto_tortoSM.addEdge(0,
//                                              1); // Água do Descoberto para o TortoSM
//    Transfers transfer_descoberto_tortoSM(1, 0, 0, 0.1, {1},
//                                          {0.5e-6 * 3600 * 24 * 7},
//                                          {caesb_tortoSM_transfer_trigger}, //Descoberto transfere até 0.5 m³/s para o TortoSM
//                                          transfer_graph_descoberto_tortoSM,
//                                          vector<double>(), vector<int>());
//
//    drought_mitigation_policies.push_back(&transfer_tortoSM_descoberto);
//    drought_mitigation_policies.push_back(&transfer_descoberto_tortoSM);

    vector<double> transfer_rofs = {caesb_descoberto_transfer_trigger,  // [0] assigned to system_id 0 (Descoberto)
                                    caesb_tortoSM_transfer_trigger};  // [1] assigned to system_id 1 (TortoSM)

       //This pipline capacity for transfers was updated on May 2025 to 1.0 m3/sec                             
//     vector<double> transfer_capacities = {0.5e-6 * 3600 * 24 * 7,        // Capacity TortoSM→Descoberto
//                                           0.7e-6 * 3600 * 24 * 7};        // Capacity Descoberto→TortoSM

    vector<double> transfer_capacities = {1.0e-6 * 3600 * 24 * 7,        // Capacity TortoSM→Descoberto
                                          0};        // Capacity Descoberto→TortoSM

    vector<int> tranfers_wss_ids = {0, 1};  // system_id 0=Descoberto, system_id 1=TortoSM
    TransfersBilateral transfers(0, transfer_capacities, 0.1, 1.1,
                                 transfer_rofs, tranfers_wss_ids,
                                 sender_demand_protection_factor);
    drought_mitigation_policies.push_back(&transfers);

    // POLÍTICA DE TRANSFERÊNCIA DE EMERGÊNCIA — Paranoá → TortoSM
    // Activated when TortoSM ROF exceeds caesb_paranoa_transfer_trigger.
    // Initial pipe capacity: 0.1 m³/s (100 l/s). Each of the three Paranoa
    // expansions (IDs 7, 8, 9), when built, adds another 0.1 m³/s to the
    // effective pipeline capacity (matching the +100 l/s treatment expansion).
    // Volume also bounded by Paranoa's 0.8% supply allocation for TortoSM.
    // Must be added AFTER TransfersBilateral so that TransfersBilateral's demand
    // offset reset does not overwrite the emergency supply offset.
    // Transfer cost: receiver (TortoSM) is charged at its current water price ×
    // 1.1 (10% surcharge), consistent with the TransfersBilateral convention.
    double paranoa_base_pipe_capacity  = 0.1e-6 * 3600 * 24 * 7; // initial 100 l/s
    double paranoa_pipe_per_expansion  = 0.1e-6 * 3600 * 24 * 7; // +100 l/s per expansion
    EmergencyTransferParanoa paranoa_emergency(
            (int) drought_mitigation_policies.size(),  // unique policy id
            1,                                          // receiver: TortoSM (system_id = 1)
            3,                                          // Paranoa water source id = 3
            paranoa_base_pipe_capacity,                 // initial pipe capacity
            caesb_paranoa_transfer_trigger,             // ROF trigger
            1.1,                                        // cost multiplier (10% surcharge)
            {7, 8, 9},                                  // expansion IDs that grow the pipeline
            paranoa_pipe_per_expansion);                // capacity increment per expansion
    drought_mitigation_policies.push_back(&paranoa_emergency);

    /// Creates simulation object depending on use (or lack thereof) ROF tables
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_wss_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           wss_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run,
                           rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_wss_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           wss_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run,
                           rof_tables,
                           table_storage_shift,
                           rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else {
        s = new Simulation(water_sources,
                           g,
                           reservoir_wss_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           wss_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    }
//    double end_time = omp_get_wtime();
//	printf("Function evaluation time: %f s\n", end_time - start_time);

    double realization_end = omp_get_wtime();
    printf("Simulation took %fs\n", realization_end - realization_start);

    /// Calculate objectives and store them in Borg decision variables array.
    /// Skip objective calculation when exporting ROF tables, since no data
    /// was collected and realization models have already been deleted.
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        objectives = vector<double>(Constants::getNumObjectives(), 0.0);
    } else {
        objectives = calculateAndPrintObjectives(false);
    }

// With WSS architecture, there is only ONE utility (CAESB) with 2 WSS.
// The objectives vector layout depends on experiment mode:
// Mode 5: [WSS0_rel, WSS1_rel, restr_freq, NPC, worst_cost, WSS0_afford, WSS1_afford]
// Modes 1-4: [reliability, restr_freq, NPC, worst_cost (, affordability (, severity))]

#ifdef  PARALLEL
    if (Constants::includePerWSSObjectives()) {
        // Mode 5: copy 7 per-WSS objectives; negate both reliability objectives for minimization
        objs[0] = -objectives[0]; // WSS0 reliability (negated)
        objs[1] = -objectives[1]; // WSS1 reliability (negated)
        objs[2] = objectives[2];  // Restriction frequency
        objs[3] = objectives[3];  // Infrastructure NPC
        objs[4] = objectives[4];  // Worst case costs
        objs[5] = objectives[5];  // WSS0 affordability
        objs[6] = objectives[6];  // WSS1 affordability
    } else {
        objs[0] = -objectives[0];  // Negative reliability
        objs[1] = objectives[1];   // Restriction frequency
        objs[2] = objectives[2];   // Infrastructure NPC
        objs[3] = objectives[3];   // Worst case costs

        int obj_idx = 4;
        // Only copy affordability if it was calculated (experiments 3 & 4)
        if (Constants::includeAffordabilityObjective() && obj_idx < (int)objectives.size()) {
            objs[obj_idx] = objectives[obj_idx];   // Affordability index
            obj_idx++;
        }
        // Only copy severity if it was included
        if (Constants::includeSeverityObjective() && obj_idx < (int)objectives.size()) {
            objs[obj_idx] = objectives[obj_idx];   // Failure severity
            obj_idx++;
        }
    }
     
        if (s != nullptr) {	 // != significa "diferente de"
            delete s;
    }
    s = nullptr;
#endif
//    } catch (const std::exception& e) {
//        simulationExceptionHander(e, s, objs, vars);
//	return 1;
//    }

    delete s;

    return 0;
}

int Caesb::simulationExceptionHander(const std::exception &e,
                                     Simulation *s, // :: significa "resolução de escopo"
                                     double *objs, const double *vars) {
    int num_dec_var = 22; //número de variáveis desse estudo de caso
//        printf("Exception called during calculations. Decision variables are below:\n");
    ofstream sol;
    int world_rank;

#ifdef  PARALLEL
    // int mpi_initialized;
    // MPI_Initialized(&mpi_initialized);
    // if (mpi_initialized)
 //            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // else
        world_rank = 0;
#else
    world_rank = 0;
#endif
    string error_file = "sol_error_rank_" + to_string(world_rank) + ".csv";
    sol.open(error_file.c_str());
    for (int i = 0; i < num_dec_var; ++i) {
        sol << vars[i] << ",";
    }
    sol << flush;
    sol.close();
    printf("Error. Decision variables printed in %s\n", error_file.c_str());

#ifdef PARALLEL
    objs[0] = 0.;
    objs[1] = 1.1;
    objs[2] = 1000;
        objs[3] = 5.;
    
    // Only set affordability if included in experiment
    if (Constants::includeAffordabilityObjective()) {
                objs[4] = 1.0;  // Worst affordability
    }

    if (s != nullptr) {
        delete s;
        s = nullptr;
    }
#else
    Utils::print_exception(e);
#endif

    return 1;
}


Caesb::~Caesb() = default;

Caesb::Caesb(unsigned long n_weeks, int import_export_rof_table)
        : Problem(n_weeks) {
    if (import_export_rof_table == EXPORT_ROF_TABLES) {
        table_gen_storage_multiplier = BASE_STORAGE_CAPACITY_MULTIPLIER;
    } else {
        table_gen_storage_multiplier = 1.;
    }
}


void
Caesb::readInputData() { //A partir dessa linha serão inseridos os dados de entrada para que o modelo possa funcionar
    printf("Reading input data.\n");
    string data_dir = DEFAULT_DATA_DIR + BAR;

//#pragma omp parallel default(none) num_threads(omp_get_thread_num())
//    {
// #pragma omp single
        streamflows_descoberto = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "descoberto_inflows.csv", n_realizations);

// #pragma omp single
        streamflows_tortoSM = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "tortoSM_inflows.csv", n_realizations);
// #pragma omp single
        streamflows_bananal_torto = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "bananal_torto_inflows.csv", n_realizations);
//// #pragma omp single
//        streamflows_torto = Utils::parse2DCsvFile( //inserção dos dados de vazão
//                io_directory + DEFAULT_DATA_DIR + "inflows" +
//                evap_inflows_suffix +
//                BAR + "torto_inflows.csv", n_realizations);
        // }
// #pragma omp single
        streamflows_paranoa = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "paranoa_inflows.csv", n_realizations);

// #pragma omp single
        streamflows_corumbaIV = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "corumbaIV_inflows.csv", n_realizations);

// };
        //cout << "Reading evaporations." << endl;
// #pragma omp single
        evap_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "descoberto_evap.csv",
                n_realizations); //inserção dos dados de evaporação
// #pragma omp single
        evap_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "tortoSM_evap.csv",
                n_realizations); //inserção dos dados de evaporação
// #pragma omp single
        evap_paranoa = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "paranoa_evap.csv",
                n_realizations); //inserção dos dados de evaporação
// #pragma omp single
        evap_corumba = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                evap_inflows_suffix +
                BAR + "corumba_evap.csv",
                n_realizations); //inserção dos dados de evaporação


        //cout << "Reading demands." << endl;
// #pragma omp single
        demand_caesb_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" +
                evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_descoberto_demand.csv", n_realizations);

        demand_caesb_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" +
                evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_tortoSM_demand.csv", n_realizations);


        //cout << "Reading others." << endl;
// #pragma omp single
//        {
            demand_to_wastewater_fraction_caesb_descoberto = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "demand_to_wastewater_fraction_caesb_descoberto.csv"); //demanda de efluentes da caesb

            demand_to_wastewater_fraction_caesb_tortoSM = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "demand_to_wastewater_fraction_caesb_tortoSM.csv"); //demanda de efluentes da caesb

            caesbDescobertoDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbDescobertoDemandClassesFractions.csv"); //demanda de cada categoria de usuário da caesb

            caesbTortoSMDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbTortoSMDemandClassesFractions.csv"); //demanda de cada categoria de usuário da caesb

            caesbDescobertoUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbDescobertoUserClassesWaterPrices.csv"); //tarifa de água para cada categoria de usuário da caesb
            //OBS: inserir no valor a tarifa de esgoto também na última coluna desse arquivo (100% da de água)

            caesbTortoSMUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbTortoSMUserClassesWaterPrices.csv");

            caesbPriceRestrictionMultipliers = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbPriceRestrictionMultipliers.csv"); //% de aumento da tarifa para cada categoria durante o racionamento
//        }
//    cout << "Done reading input data." << endl;
//    }

}
