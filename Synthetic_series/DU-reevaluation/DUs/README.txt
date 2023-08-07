The DU samplings to be used in WaterPaths are segregated in three main categories (and .csv files): water sources, utilities and infrastructures.
Each of these files have a specific column for each POSSIBLE Deep Uncertainty originally set for the Sedento Valley test case of WaterPaths.
If a certain DU (column) should not be considered for the Study Case, its value is set to zero.
- Water Sources: contains a columun for Evaporation (deep) uncertainty, and 26 columns grouped in pairs. Each pair corresponds to the uncertainties associated 
to a specific infrastructure option - permitting time and construction costs;
- Utilities: the deep uncertainties associated to this category are demand growth, bond interest rate, bond term and discount rate, and water tariffs variation;
- Policies: the DUs associated to this category are the effectiveness of drought mitigation instruments (one column for Santa Maria service area and another to
Descoberto service area).
The values for each column are calculated as multiplicative factors as described below, and were executed in the 'DU calculation' file.
1) First, LHS values between 0 and 1 must be generated for each DU (column), considering the number of scenarios that will be evaluated in WaterPaths simulations. 
For example, for Water Sources DUs a study case will only consider Evaporation as a DU, and 500 scenarios will be generated using this DU, LHS application should 
generate 500 values, placed in one single column and 500 rows. In FDB Study Case, 7 DUs were considering for generating 1,000 scenarios, so LHS matrix was generated
on Python with 7 columns and 1,000 rows with LHS sampled values between 0 and 1. The Python file that generated LHS matrix is THE 'DU_generation_LHS' file.
2) The final multiplicative factors for each DU are gauged with the follow operation: for each specific DU, its upper and lower boundaries are subtracted, then
multiplied by the LHS respective value, then summed to lower boundary. The result is that all DU values, in multiplicative factors, are inside the DU range defined
according DU definitions (Table S6 of Supporting Information).
3) Then, the values are distributed in specific tabs for Water Sources, Utilities and Policies. 
4) Lastly, each of these tabs are separately saved in specific .csv files, to be used in WaterPaths. 