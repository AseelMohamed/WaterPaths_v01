cmake_minimum_required(VERSION 3.6)
project(WaterPaths_WSS)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP (portable way)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O0 -march=native")

# # Add debugging flags for bounds checking (conditionally, since MinGW doesn't support sanitizers)
# if(NOT WIN32 OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     # Enable sanitizers only on non-Windows or with Clang
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined -fno-omit-frame-pointer -g")
#     set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address,undefined")
# else()
#     # For MinGW/Windows, just add debug symbols and bounds checking where possible
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fno-omit-frame-pointer -Wall -Wextra")
# endif()

set(SOURCE_FILES
        ../WaterPaths_WSS/src/ContinuityModels/Base/ContinuityModel.cpp
        ../WaterPaths_WSS/src/ContinuityModels/Base/ContinuityModel.h
        ../WaterPaths_WSS/src/ContinuityModels/ContinuityModelRealization.cpp
        ../WaterPaths_WSS/src/ContinuityModels/ContinuityModelRealization.h
        ../WaterPaths_WSS/src/ContinuityModels/ContinuityModelROF.cpp
        ../WaterPaths_WSS/src/ContinuityModels/ContinuityModelROF.h
        ../WaterPaths_WSS/src/Controls/Base/ControlRules.h
        ../WaterPaths_WSS/src/Controls/Base/MinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/Base/MinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/Custom/JordanLakeMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/Custom/JordanLakeMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/Custom/FallsLakeMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/Custom/FallsLakeMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/InflowMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/InflowMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/SeasonalMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/SeasonalMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/StorageMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/StorageMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/FixedMinEnvFlowControl.cpp
        ../WaterPaths_WSS/src/Controls/FixedMinEnvFlowControl.h
        ../WaterPaths_WSS/src/Controls/WwtpDischargeRule.cpp
        ../WaterPaths_WSS/src/Controls/WwtpDischargeRule.h
        ../WaterPaths_WSS/src/Controls/EvaporationSeries.cpp
        ../WaterPaths_WSS/src/Controls/EvaporationSeries.h
        ../WaterPaths_WSS/src/DataCollector/Base/DataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/Base/DataCollector.h
        ../WaterPaths_WSS/src/DataCollector/AllocatedReservoirDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/AllocatedReservoirDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/EmptyDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/EmptyDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/IntakeDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/IntakeDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/MasterDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/MasterDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/QuaryDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/QuaryDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/ReservoirDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/ReservoirDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/RestrictionsDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/RestrictionsDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/TransfersDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/TransfersDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/TransfersBilateralDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/TransfersBilateralDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/UtilitiesDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/UtilitiesDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/WaterReuseDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/WaterReuseDataCollector.h
        ../WaterPaths_WSS/src/DataCollector/WSSDataCollector.cpp
        ../WaterPaths_WSS/src/DataCollector/WSSDataCollector.h
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.cpp
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Base/DroughtMitigationPolicy.h
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/InsuranceStorageToROF.cpp
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/InsuranceStorageToROF.h
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Restrictions.cpp
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Restrictions.h
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Transfers.cpp
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/Transfers.h
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/TransfersBilateral.cpp
        ../WaterPaths_WSS/src/DroughtMitigationInstruments/TransfersBilateral.h
        ../WaterPaths_WSS/src/Problem/Base/Problem.cpp
        ../WaterPaths_WSS/src/Problem/Base/Problem.h
        ../WaterPaths_WSS/src/Problem/Caesb.cpp
        ../WaterPaths_WSS/src/Problem/Caesb.h
        ../WaterPaths_WSS/src/Simulation/Simulation.cpp
        ../WaterPaths_WSS/src/Simulation/Simulation.h
        ../WaterPaths_WSS/src/SystemComponents/Bonds/Base/Bond.cpp
        ../WaterPaths_WSS/src/SystemComponents/Bonds/Base/Bond.h
        ../WaterPaths_WSS/src/SystemComponents/Bonds/LevelDebtServiceBond.cpp
        ../WaterPaths_WSS/src/SystemComponents/Bonds/LevelDebtServiceBond.h
        ../WaterPaths_WSS/src/SystemComponents/Bonds/BalloonPaymentBond.cpp
        ../WaterPaths_WSS/src/SystemComponents/Bonds/BalloonPaymentBond.h
        ../WaterPaths_WSS/src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.cpp
        ../WaterPaths_WSS/src/SystemComponents/Bonds/FloatingInterestBalloonPaymentBond.h
        ../WaterPaths_WSS/src/SystemComponents/Utility/Utility.cpp
        ../WaterPaths_WSS/src/SystemComponents/Utility/Utility.h
        ../WaterPaths_WSS/src/SystemComponents/Utility/InfrastructureManager.cpp
        ../WaterPaths_WSS/src/SystemComponents/Utility/InfrastructureManager.h
        ../WaterPaths_WSS/src/SystemComponents/Utility/WaterSupplySystems.cpp
        ../WaterPaths_WSS/src/SystemComponents/Utility/WaterSupplySystems.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Base/WaterSource.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Base/WaterSource.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/AllocatedReservoir.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/AllocatedReservoir.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Intake.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Intake.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Quarry.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Quarry.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Relocation.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Relocation.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Reservoir.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/Reservoir.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/ReservoirExpansion.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/ReservoirExpansion.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/WaterReuse.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/WaterReuse.h
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.cpp
        ../WaterPaths_WSS/src/SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h
        ../WaterPaths_WSS/src/SystemComponents/Catchment.cpp
        ../WaterPaths_WSS/src/SystemComponents/Catchment.h
        ../WaterPaths_WSS/src/Utils/Constants.h
        ../WaterPaths_WSS/src/Utils/DataSeries.cpp
        ../WaterPaths_WSS/src/Utils/DataSeries.h
        ../WaterPaths_WSS/src/Utils/Graph/Graph.cpp
        ../WaterPaths_WSS/src/Utils/Graph/Graph.h
        ../WaterPaths_WSS/src/Utils/ObjectivesCalculator.cpp
        ../WaterPaths_WSS/src/Utils/ObjectivesCalculator.h
        ../WaterPaths_WSS/src/Utils/Constants.cpp
        ../WaterPaths_WSS/src/Utils/Constants.h
        ../WaterPaths_WSS/src/Utils/QPSolver/Array.cpp
        ../WaterPaths_WSS/src/Utils/QPSolver/Array.h
        ../WaterPaths_WSS/src/Utils/QPSolver/QuadProg++.cpp
        ../WaterPaths_WSS/src/Utils/QPSolver/QuadProg++.h
        ../WaterPaths_WSS/src/Utils/Utils.cpp
        ../WaterPaths_WSS/src/Utils/Utils.h
        ../WaterPaths_WSS/src/main.cpp)

add_executable(WaterPaths_WSS ${SOURCE_FILES})

# Add debug info & disable optimizations only for Debug builds
target_compile_options(WaterPaths_WSS PRIVATE
        $<$<CONFIG:Debug>:-O0 -g -march=native>
)

# OpenMP link
target_link_libraries(WaterPaths_WSS PRIVATE OpenMP::OpenMP_CXX)