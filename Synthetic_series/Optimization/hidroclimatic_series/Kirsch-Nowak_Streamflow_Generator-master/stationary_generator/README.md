# Kirsch-Nowak_Streamflow_Generator
This directory contains MATLAB code for generating correlated synthetic daily streamflow time series at multiple sites assuming stationary hydrology. Monthly flows are generated using Cholesky decomposition (see [Kirsch et al. 2013](http://ascelibrary.org/doi/abs/10.1061/(ASCE)WR.1943-5452.0000287)) and then disaggregated to daily flows by proportionally scaling daily flows from a randomly selected historical month +/- 7 days as in [Nowak et al. (2010)](http://onlinelibrary.wiley.com/doi/10.1029/2009WR008530/full).

The historical data are presented as inflows for all water sources (Descoberto, Corumba, Santa Maria, Paranoa resevoirs, and Torto-Bananal stream) and evaporation for the four resevoirs. To run the stationary streamflow generator, the historical data must be in `../data` first have to be processed into a 2D array where each row is the data on a different date, and each column is the flow/evaporation rate at a different site. Leap days are also removed from the time series to create a prediodic signal. This array can be generate by running `clean_data.m`, which writes `../data/Qdaily.txt`. `clean_data.m` also writes files to `../validation/historical` which are used to make the validation figures.

Once the data are pre-processed, the synthetic generator can be run with `script_example.m`. This will write files to `../validation/synthetic` for validation.

Contents:  
`chol_corr.m`: Computes Cholesky decomposition of correlation matrix and attempts to repair non-positive definite matrices

`clean_data.m`: Pre-processes historical data into format read by generator and validation codes

`combined_generator.m`: Calls `monthly_main.m` to generate monthly streamflows and then `KNN_identification.m` and `KNN_sampling.m` to disaggregate them to daily flows

`convert_data_to_monthly.m`: Converts daily data (volume/sec) to monthly totals (volume/month)

`KNN_identification.m`: Finds the k-nearest neighbors for disaggregation

`KNN_sampling.m`: Randomly selects one of the k-nearest neighbors for disaggregation

`monthly_gen.m`: Generates monthly streamflows

`monthly_main.m`: Calls `monthly_gen.m` to generate monthly streamflows for a specified number of realizations

`script_example.m`: Shows how to generate synthetic streamflows on a sample data set (the Susquehanna River Basin)
