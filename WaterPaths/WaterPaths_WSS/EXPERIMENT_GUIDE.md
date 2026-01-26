# Computational Experiments Guide

## Overview
The WaterPaths_WSS model now supports 4 different experiment configurations to test different objective function formulations. The experiments vary in:
1. Number of objectives (5 or 6)
2. How reliability is aggregated across Water Supply Systems (MIN vs AVERAGE)
3. How affordability is aggregated across Water Supply Systems (MAX vs AVERAGE)

## Experiment Modes

### Experiment 1: 5 Objectives, MIN Reliability
- **Number of Objectives**: 5 (no affordability)
- **Reliability Aggregation**: MIN (worst-case WSS)
- **Use Case**: Baseline without affordability, conservative reliability

### Experiment 2: 5 Objectives, AVERAGE Reliability
- **Number of Objectives**: 5 (no affordability)
- **Reliability Aggregation**: AVERAGE across all WSS
- **Use Case**: Baseline without affordability, average reliability performance

### Experiment 3: 6 Objectives, MIN Reliability, MAX Affordability [DEFAULT]
- **Number of Objectives**: 6 (includes affordability)
- **Reliability Aggregation**: MIN (worst-case WSS)
- **Affordability Aggregation**: MAX (worst-case WSS)
- **Use Case**: Current configuration, most conservative approach

### Experiment 4: 6 Objectives, AVERAGE Reliability, AVERAGE Affordability
- **Number of Objectives**: 6 (includes affordability)
- **Reliability Aggregation**: AVERAGE across all WSS
- **Affordability Aggregation**: AVERAGE across all WSS
- **Use Case**: Average performance across all metrics

## Usage

### Command Line Flag
Use the `-E` flag to specify experiment mode:

```bash
# Experiment 1
./WaterPaths_WSS -E 1 -s solutions.csv -d InputFiles/

# Experiment 2
./WaterPaths_WSS -E 2 -s solutions.csv -d InputFiles/

# Experiment 3 (default if -E not specified)
./WaterPaths_WSS -E 3 -s solutions.csv -d InputFiles/

# Experiment 4
./WaterPaths_WSS -E 4 -s solutions.csv -d InputFiles/
```

### With Borg Optimization
```bash
# Run optimization with experiment mode 2
mpirun -np 4 ./WaterPaths_WSS -b -E 2 -n 100000 -o 1000 -e 42
```

## Technical Details

### Modified Files
1. **Constants.h/cpp**: Added experiment configuration and helper functions
2. **ObjectivesCalculator.h/cpp**: Added configurable aggregation methods
3. **MasterDataCollector.cpp**: Uses experiment mode to select calculation method
4. **main.cpp**: Parses `-E` flag and dynamically sets objective count
5. **CMakeLists.txt**: Added Constants.cpp to build

### Aggregation Methods

#### Reliability
- **MIN (Experiments 1 & 3)**: `reliability = min(reliability_wss1, reliability_wss2, ...)`
  - Utility is only as reliable as its weakest WSS
  - Conservative approach, ensures all systems meet reliability standards
  
- **AVERAGE (Experiments 2 & 4)**: `reliability = mean(reliability_wss1, reliability_wss2, ...)`
  - Average performance across all WSS
  - Allows trade-offs between different WSS

#### Affordability  
- **MAX (Experiment 3)**: `affordability = max(afford_wss1, afford_wss2, ...)`
  - Takes worst-case affordability across all WSS
  - Ensures no WSS has unacceptable affordability burden
  
- **AVERAGE (Experiment 4)**: `affordability = mean(afford_wss1, afford_wss2, ...)`
  - Average affordability across all WSS
  - Allows some WSS to have higher costs if balanced by others

### Objective Ordering
All experiments maintain the same objective ordering:
1. Reliability (maximize → minimize 1-reliability)
2. Restriction Frequency (minimize)
3. Net Present Cost Infrastructure (minimize)
4. Peak Financial Cost (minimize)
5. Worst Case Costs (minimize)
6. Affordability Index (minimize) - *only in experiments 3 & 4*

## Stability Considerations

### Safe Defaults
- Experiment mode defaults to 3 (current behavior) if `-E` not specified
- All existing code continues to work without modifications
- Backward compatible with existing solution files

### Error Handling
- Invalid experiment mode (not 1-4) triggers error message and exit
- All calculation functions include bounds checking
- Empty WSS data handled gracefully with warnings

### Testing Recommendations
1. **Start with small test**: Run single solution with each experiment mode
2. **Compare outputs**: Verify objective values are reasonable
3. **Check consistency**: Run same solution multiple times, ensure repeatability
4. **Validate extremes**: Test with edge cases (few/many realizations, short/long timelines)

### Performance
- Negligible overhead from experiment configuration
- All aggregation methods have O(n) complexity where n = number of WSS
- Memory usage unchanged (arrays sized for maximum 6 objectives)

## Example Workflow

```bash
# 1. Test single solution with all experiments
for exp in 1 2 3 4; do
    ./WaterPaths_WSS -E $exp -m 0 -s solutions.csv -d InputFiles/
    mv output/Objectives_s0.out output/Objectives_s0_exp$exp.out
done

# 2. Compare results
# Check that experiments 1 & 2 have 5 objectives
# Check that experiments 3 & 4 have 6 objectives
# Verify reliability values differ between MIN and AVERAGE methods

# 3. Run full optimization for selected experiment
mpirun -np 16 ./WaterPaths_WSS -b -E 2 -n 1000000 -o 10000 -e 123
```

## Troubleshooting

### Compilation Issues
- Ensure Constants.cpp is added to CMakeLists.txt
- Clean build directory: `rm -rf cmake-build-*`
- Rebuild: `cmake . && make`

### Runtime Issues
- **"Invalid experiment mode"**: Use -E with values 1-4 only
- **Objective count mismatch**: Ensure all code uses `getNumObjectives()`
- **Empty WSS data**: Check that WSS data collectors are properly initialized

### Result Interpretation
- Lower reliability in AVERAGE vs MIN is expected (less conservative)
- Affordability values should be lower in AVERAGE vs MAX
- Total infrastructure costs should be similar across experiments (not aggregated)

## Future Extensions
To add more experiments:
1. Update `EXPERIMENT_MODE` range check in main.cpp
2. Add logic to `get*AggregationMethod()` functions in Constants.h
3. Document new experiment in this guide
4. Test thoroughly before production runs
