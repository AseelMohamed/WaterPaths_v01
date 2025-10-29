cmake_minimum_required(VERSION 3.6)
project(WaterPaths_Original)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP (portable way)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O0 -march=native")

set(SOURCE_FILES
        ../WaterPaths_Original/src/ContinuityModels/Base/ContinuityModel.cpp
        ../WaterPaths_Original/src/ContinuityModels/Base/ContinuityModel.h
        ../WaterPaths_Original/src/ContinuityModels/ContinuityModelRealization.cpp
        ../WaterPaths_Original/src/ContinuityModels/ContinuityModelRealization.h
        ../WaterPaths_Original/src/ContinuityModels/ContinuityModelROF.cpp
        ../WaterPaths_Original/src/ContinuityModels/ContinuityModelROF.h
        ../WaterPaths_Original/src/Controls/Base/ControlRules.h
        ../WaterPaths_Original/src/Controls/Base/MinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/Base/MinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/Custom/JordanLakeMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/Custom/JordanLakeMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/Custom/FallsLakeMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/Custom/FallsLakeMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/InflowMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/InflowMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/SeasonalMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/SeasonalMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/StorageMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/StorageMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/FixedMinEnvFlowControl.cpp
        ../WaterPaths_Original/src/Controls/FixedMinEnvFlowControl.h
        ../WaterPaths_Original/src/Controls/WwtpDischargeRule.cpp
        ../WaterPaths_Original/src/Controls/WwtpDischargeRule.h
        ../WaterPaths_Original/src/Controls/EvaporationSeries.cpp
        ../WaterPaths_Original/src/Controls/EvaporationSeries.h
        ../WaterPaths_Original/src/DataCollector/Base/DataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/Base/DataCollector.h
        ../WaterPaths_Original/src/DataCollector/AllocatedReservoirDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/AllocatedReservoirDataCollector.h
        ../WaterPaths_Original/src/DataCollector/EmptyDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/EmptyDataCollector.h
        ../WaterPaths_Original/src/DataCollector/IntakeDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/IntakeDataCollector.h
        ../WaterPaths_Original/src/DataCollector/MasterDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/MasterDataCollector.h
        ../WaterPaths_Original/src/DataCollector/QuaryDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/QuaryDataCollector.h
        ../WaterPaths_Original/src/DataCollector/ReservoirDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/ReservoirDataCollector.h
        ../WaterPaths_Original/src/DataCollector/RestrictionsDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/RestrictionsDataCollector.h
        ../WaterPaths_Original/src/DataCollector/TransfersDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/TransfersDataCollector.h
        ../WaterPaths_Original/src/DataCollector/TransfersBilateralDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/TransfersBilateralDataCollector.h
        ../WaterPaths_Original/src/DataCollector/UtilitiesDataCollector.h
        ../WaterPaths_Original/src/DataCollector/UtilitiesDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/WaterReuseDataCollector.cpp
        ../WaterPaths_Original/src/DataCollector/WaterReuseDataCollector.h
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.cpp
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.h
        ../WaterPaths_Original/src/DroughtMitigationInstruments/InsuranceStorageToROF.cpp
        ../WaterPaths_Original/src/DroughtMitigationInstruments/InsuranceStorageToROF.h
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Restrictions.cpp
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Restrictions.h
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Transfers.cpp
        ../WaterPaths_Original/src/DroughtMitigationInstruments/Transfers.h
        ../WaterPaths_Original/src/DroughtMitigationInstruments/TransfersBilateral.cpp
        ../WaterPaths_Original/src/DroughtMitigationInstruments/TransfersBilateral.h
        ../WaterPaths_Original/src/Problem/Base/Problem.cpp
        ../WaterPaths_Original/src/Problem/Base/Problem.h
        ../WaterPaths_Original/src/Problem/Caesb.cpp
        ../WaterPaths_Original/src/Problem/Caesb.h
        ../WaterPaths_Original/src/Simulation/Simulation.cpp
        ../WaterPaths_Original/src/Simulation/Simulation.h
        ../WaterPaths_Original/src/SystemComponents/Bonds/Base/Bond.cpp
        ../WaterPaths_Original/src/SystemComponents/Bonds/Base/Bond.h
        ../WaterPaths_Original/src/SystemComponents/Bonds/LevelDebtServiceBond.cpp
        ../WaterPaths_Original/src/SystemComponents/Bonds/LevelDebtServiceBond.h
        ../WaterPaths_Original/src/SystemComponents/Bonds/BalloonPaymentBond.cpp
        ../WaterPaths_Original/src/SystemComponents/Bonds/BalloonPaymentBond.h
        ../WaterPaths_Original/src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.cpp
        ../WaterPaths_Original/src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.h
        ../WaterPaths_Original/src/SystemComponents/Utility/Utility.cpp
        ../WaterPaths_Original/src/SystemComponents/Utility/Utility.h
        ../WaterPaths_Original/src/SystemComponents/Utility/InfrastructureManager.cpp
        ../WaterPaths_Original/src/SystemComponents/Utility/InfrastructureManager.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Base/WaterSource.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Base/WaterSource.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/AllocatedReservoir.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/AllocatedReservoir.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Intake.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Intake.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Quarry.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Quarry.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Relocation.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Relocation.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Reservoir.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/Reservoir.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/ReservoirExpansion.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/ReservoirExpansion.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/WaterReuse.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/WaterReuse.h
        ../WaterPaths_Original/src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.cpp
        ../WaterPaths_Original/src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h
        ../WaterPaths_Original/src/SystemComponents/Catchment.cpp
        ../WaterPaths_Original/src/SystemComponents/Catchment.h
        ../WaterPaths_Original/src/Utils/Constants.h
        ../WaterPaths_Original/src/Utils/DataSeries.cpp
        ../WaterPaths_Original/src/Utils/DataSeries.h
        ../WaterPaths_Original/src/Utils/Graph/Graph.cpp
        ../WaterPaths_Original/src/Utils/Graph/Graph.h
        ../WaterPaths_Original/src/Utils/ObjectivesCalculator.cpp
        ../WaterPaths_Original/src/Utils/ObjectivesCalculator.h
        ../WaterPaths_Original/src/Utils/QPSolver/Array.cpp
        ../WaterPaths_Original/src/Utils/QPSolver/Array.h
        ../WaterPaths_Original/src/Utils/QPSolver/QuadProg++.cpp
        ../WaterPaths_Original/src/Utils/QPSolver/QuadProg++.h
        ../WaterPaths_Original/src/Utils/Utils.cpp
        ../WaterPaths_Original/src/Utils/Utils.h
        ../WaterPaths_Original/src/main.cpp)

add_executable(WaterPaths_Original ${SOURCE_FILES})

# Add debug info & disable optimizations only for Debug builds
target_compile_options(WaterPaths_Original PRIVATE
        $<$<CONFIG:Debug>:-O0 -g -march=native>
)

# OpenMP link
target_link_libraries(WaterPaths_Original PRIVATE OpenMP::OpenMP_CXX)