The reference period for this study in order to obtain 40 year long series was the period between 2018 and 2060 (with two years more). Thus, the demand between 2020 to 2060 was obtained 
multiplying population projection for this horizon (based upon IBGE estimatives) by average per capita consumption obtained for the period of July, 2013 to June,2017 
provided by Federal District Sanitation Plan (both decomposed to weekly time scale). The procedure used to calculate population projection, per capita consumption and 
1,000 demand series generations are described below.
Hypotesis considered for synthetic series generation:
- Urban population will keep growing, at a 0.045% rate;
- Per capita consumption will keep constant along all planning period (40 years or 2086 weeks);
- Percentage losses will keep constant along all planning period (32.7%);
- Weekly growth population is treated as an uncertainty factor (thus, 1,000 different population projection are generated)
Population projection until 2060
Obs.: the procedures executed in this stage are in Excel file "Pop. projections (2020 a 2060) + per capita consumption – Santa Maria and Descoberto"
1) The data of this first tab was extracted from annual population projection of IBGE (the Brazilian Institute for Geography and Statistics), from 2010 to 2060;
2) Then, based on the urban population percentage of 2010 (IBGE) and 2018 (taken from the Household Sampling Districtal Plan - PDAD), an estimation of urban 
population growth was calculated until 2060;
3) Having the estimated population growth until 2060, the urban population was obtained year by year until 2060;
4) Territorial Planning Units (TPUs, which is a set of Administrative Regions, AR) - initially, the relative contribution (%) of each TPU in urban population of FDB 
was calculated, considering the 2010-2018 period. To do so, the data from 2010 (IBGE) and 2011-13-15-18 (PDAD) was applied (TPU tabs). Next, the average percentage 
obtain for these five years (aritmetic average) was calculated, and used as a basis to calculate relative contribution of each TPU between 2018 and 2060 (tab Urban 
Pop. Evolution by TPU);
5) With the average percentage, the projection of relative contribution of each TPU in FDB urban population was calculated (tab Population projection by TPU);
6) After, the urban population projection was calculated for each TPU until 2060, by multiplying urban population projection for each year (IBGE) and respective
percentage of TPU relative contribution for that year (tab Population projection by TPU);
7) ARs - initially, the percentage of relative contribution of each RA in its respective TPU was calculated, considering the 2010-2018 period. To do so, the data
from 2010 (IBGE) and 2011-13-15-18 (PDAD, tab ARs Urban Pop.) was used. After, the average percentage of relative contribution of each RA in its respective TPU was 
calculated, in the corresponding period (2010-11-13-15-18). The average of this five-year period was used as a reference for the evolution of relative contribution
of each AR in the 2018-2060 period;
8) With the average percentage obtained, the projection of relative contribution of ARs in total urban population of its respective TPU was calculated until 2060 
(tabs TPU I ARs urban projections, TPU II ARs urban projections and so on);
9) The projection of urban population of each AR was calculated until 2060, by multiplying urban population projection of each TPU and the respective percentage of 
relative contribution of the AR in TPU population for each year (same tabs mentioned on item 8);
10) Lastly the ARs projections (2018-2060) were all aggregated in one single tab for better visualization. Besides, 6 ARs were excluded, since they don't oficially 
bellong to Santa Maria or Descoberto service areas (Brazlândia, Sobradinho I, Sobradinho II, Planaltina, Fercal and Sao Sebastiao, tab Urb proj of all 25 ARs);
11) After calculating the population projection until 2060 for all AR, the weekly population for was calculated for each year, from 2018 to 2060 (tabs 2018, 2019,...,
2060).
For weekly population of each year, the following operations were made:
	a. The annual growth rate (A) was calculated between the year involved and the previous year: [(population of present year-population of previous year)/population of previous year];
	b. The weekly growth rate (W) was calculated: W=(1+A)^(1/52)-1;
	c. The population of week 1 of year t corresponds to the multiplication of last week of year t-1 and weekly growth rate. Population of week 2 is the result of 
	multiplication of week 1 population and weekly growth rate. And so on, for each week of year t in all tabs 2018 until 2060;
12) Weekly population projection was segregated in two groups: population projection attended by Santa Maria service area and the population projection attended by 
Descoberto service area. The result is in tab Pop Proj (2020-2060) - 2 Syst;
13) After, this reference population projection for each service area was multiplied by random scalars, to obtain 1,000 different scenarios of population projection. 
The procedure to obtain weekly per capita consumption by service area (Santa Maria/Descoberto) is described below.
Per capita consumption
Obs: The procedures taken place in this stage are in file Pop. projections (2020 a 2060) + per capita consumption – Santa Maria and Descoberto
1) For weekly per capita calculation, the data from each AR per capita consumption was used, in the period from July 2013 to June,2016 (tab Weekly per capita 
(2013-2016) and Monthly per capita (2013-2016)). This per capita demand doesn't include net water losses. This data was extracted from FDB Sanitation Plan of
2017, provided by CAESB;
2) Initially, per capita consumption was converted from (L/hab.dia) to (m³/inhab.week);
3) The weekly per capita data was then distributed in the 1-52 weeks period, as follows: for january, february, april, may, july, august, october and november
there was an extension of 4 weeks, while the remaining months were assumpted to have a 5-week extension, summing up 52 weeks (tab Weekly per capita (2013-2016));
4) With this, the per capita weekly averages were calculated between July, 2013 and June, 2016, by taking the average of each week (1 to 52) of each year;
5) Finally, average weekly per capita was calculated between all AR that belongs to the same service area (Santa Maria or Descoberto). Then, the per capita unit
was turned from m³/inhab.week to hm³/inhab.week. The result is in tab Average per cap. per serv. area.
Report- 1,000 population projection series generation
1) Initially, 999 samples were generated through LHS method (random samples between 0 and 1), whose boundaries were adjusted to 0.77-1.9 range. This range was set
because it generated, by the end of population projection series generation, a decrease of until 5% and one increase up until 20% if compared to the original 
population projection. This variation could be associated to uncertainties related to populational growth rates of IBGE to FDB;
2) This random numbers were then multiplied to the population projection associated to each service area, previously calculated (File Random scalar matrix + 1000 
synthetic demand series, tab 1000 series of pop. projection) so that 1,000 different demand scenarios were generated. One of them (the first line) corresponds to the reference scenario (obtained
from file Pop. projections (2020 a 2060) + per capita consumption – Santa Maria and Descoberto, tab Pop Proj (2020-2060) - 2 Syst) while the remaining 999
correspond to demand values between the range of 77% and 190% of weekly growth rate variation;
3) To multiply the random numbers to the reference population projection, the 999 random numbers obtained were organized in a 999 row matrix, with (52*40) columns (Tab LHS(0.77-1.9)), since the 
reference population projection has 2080 columns;
4) After, the weekly per capita (calculated in file Pop. projections (2020 a 2060) + per capita consumption – Santa Maria and Descoberto) was organized in a matrix of one row and (52*40) columns
(tab Average per cap.-format adjust), so that per capita values repeat every 52 weeks (since per capita was assumed constant along all planning horizon). This was made so the 1,000 series of
population projection could be multiplied by weekly per capita values;
5) Then, the 1,000 population projection series with weekly per capita values (tab Demand (NO NET LOSSES)), so demands without net losses were obtained;
6) After, the net loss percentage (1.327 or 32.7% of net losses) was multiplied by the 1,000 previous demand series (Demand)NO NET LOSSES)), so demands with net losses were obtained, whose
results are presented in tab Demand(WITH NET LOSSES). This net loss percentage (32.7%) was obtained by calculating index averages of net losses presented in Caesb report of 2018;
7) Finally, an adjustment in the week numbers was made in the demand series (with net losses) for both service areas. Since one year has 52.14285 weeks (not exact 52 weeks), every 312 weeks 
(6 years) an one week demand value was added to the last week of the sixth year (so that every 6 years, 53 are complete). The result is in tab Demand-Adjust.in week number. The result is that
all synthetic series have 2086 errks (40 years of planning horizon);
8) The result in tab Demand-Adjust.in week number was then saved in .csv format (caesb_descoberto_demand and caesb_torto_SM_demand). 


























