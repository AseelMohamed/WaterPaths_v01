# Exploring Equity Challenges within Deeply Uncertain Water Supply Investment Pathways in the Federal District of Brazil
Bruna M. Araujo<sup>1</sup>, David F. Gold<sup>2</sup>, Lillian B. Lau<sup>2</sup>, Patrick M. Reed<sup>2</sup>, Conceição M. A. Alves<sup>1</sup>

<sup>1</sup>Department of Civil and Environmental Engineering, National University of Brasilia, Federal District, Brazil

<sup>2</sup>Department of Civil and Environmental Engineering, Cornell University, Ithaca, NY

*corresponding author: bruna.mattos.araujo@gmail.com

## Replicate our simulations and recreate our figures
This GitHub repository contains the code and data needed to replicate the computational experiments and recreate all figures for Araujo et al. (2022). Note that this experiment was run using high performance computing and cannot easily be replicated on a personal computer.

*Note: All filepaths in the code files provided should be modified to reflect current individual data and information locations. Only bootstrapped realization is provided due to memory constraints. Please contact the corresponding author at bruna.mattos.araujo@gmail.com for full set of realizations.*

## Contributing software

|  Software/Model  |  Version  |  Get the software  |  DOI  | 
| ---------------- | --------- | ---------------- |------ |
| Borg MOEA | NA | http://borgmoea.org/#contact | 10.1162/EVCO_a_00075 |
| MOEAFramework | 3.3 | http://moeaframework.org/downloads.html | NA |
| SALib | 1.4.7 | https://github.com/SALib/SALib | doi:10.18174/sesmo.18155 <br> doi:10.21105/joss.00097 |
| WaterPaths | NA | https://github.com/bernardoct/WaterPaths.git | doi.org/10.1016/j.envsoft.2020.104772 |

## Folders
- `Sythetic_series` Contains all database, codes and Excel files used for syntethic series generation
- `WaterPaths` The WaterPaths code and structure
- `post_processing` Contains all python code files used for post-processing the output of the DU Re-Evaluation
- `figure_generation` Contains all python code files used for generating figures. Note that the figures presented in the paper were generated using python code and consolidated on Adobe Illustrator to improve visualization in pdf format.

## Setup
### Download and compile WaterPaths
1. Clone this repository and unzip all files. In the command line, enter the directory where WaterPaths repository is stored.
2. Type `make gcc` into the command line to compile WaterPaths.
3. Refresh the directory. Confirm that the `triangleSimulation` file is created.
4. Running `triangleSimulation -?` will give you a list of flags to call WaterPaths with.

*To obtain the exact version of WaterPaths used to model the Federal District of Brazil, please contact the corresponding author.*

### Compiling WaterPaths with Borg MS 
To run WaterPaths in optimization mode, WaterPaths needs to be re-compiled with Borg and the appropriate flags must be passed when calling it from the command line. To do so, follow the following steps in WaterPaths folder:
1. In order to run WaterPaths with Borg, the `libborgms.a` file must be provided in the `/lib` folder. To do so, visit [this website](borgmoea.org) to request a license.
2. Once obtained, move Borg files to the folder `/borg`. Compile Borg by running `make mpi` while in the `/borg` folder.
3. Finally, move the file  `libborgms.a` to the `/lib` folder. Compile WaterPaths with `make borg`

### Generate ROF tables for realizations
To run WaterPaths in optimization mode, risk-of-failure (ROF) tables need to be built. To do so, run the `table_gen.sh` script provided in the WaterPaths folder into your terminal, previously changing the flags needed to your system specific details, location and configurations. The ROF tables are placed in `rof_tables/` folder. Keep in mind that if using a different set of inflow, evaporation, of RDM files will warrant new tables.

### Running WaterPaths in Optimization Mode
Optimization mode for FDB study case was executed running `borg_caesb.sh` file in WaterPaths folder, and 84 hours can be used as a reference of the required optimization period. The optimization results are placed in the output folder, containing in each .csv file the policy performance for all five performance objectives in each service area (Descoberto and Santa Maria) throughout all SOWs. To build SOWs during optimization, WaterPaths uses synthetic series of deep uncertainties, demand, natural inflows and evaporation, placed in  `TestFiles` folder.

### Running WaterPaths in DU-Reevaluation Mode
DU-Reevaluation mode for FDB study case was executed running `du_reeval_submission.sh` file in WaterPaths folder. The results are placed in the output_du folder, containing in each .csv file the policy performance for all five Objectives in each service area (Descoberto and Santa Maria) throughout all DU SOWs. To build SOWs during DU-Reevaluation, WaterPaths also uses synthetic series of deep uncertainties, demand, natural inflows and evaporation. 

## Synthetic Series Generation
SOWs are comprised of synthetic series of deep uncertainties, demand, natural inflows and evaporation, and need to be built for the Optimization and to DU-Reevaluation steps of DU Pathways framework. The Synthetic_series folder contains databases, code files, and explanations on how to build each type of synthetic series. One realization of deeply uncertain demand, natural inflow, and evaporation, is provided for demonstration purposes. Please contact the corresponding author at [bruna.mattos.araujo@gmail.com](bruna.mattos.araujo@gmail.com) for the full set of realizations.

Demand series were generated using official reference population growth projection UNTIL 2060 from the Brazilian Institute of Geography and Statistics (IBGE), which was disaggregated by Administrative Region using its contribution to the total population of FDB. This reference projection was multiplied by water consumption data extracted from FDB Sanitation Plan, resulting in a reference demand projection until 2060. The demand synthetic series were then obtained by multiplying this projection with various multiplicative factors derived from Latin Hypercube Sampling applications.

Natural Inflow and Evaporation series were built using the modified Fractional Gaussian Noise (mFGN) method, developed by Kirsch et al. (2013), coupled with the framework structured by Matteo Giuliani, Jon Herman and Julianne Quinn. The data used as historical observations were obtained for all affluents and main water sources for Descoberto and Santa Maria service areas in the National Water Agency of Brazil (ANA) database. 

The vectors of deep uncertainty factors were obtained by sampling multiplicative factors with LHS between 0 and 1, and adjusted according each DU reference value, lower and upper boundaries. 

## Reproducing figures
The results figures can be reproduced by executing the .py files in the `figure_generation` folder. The user should be aware that the figures will be different from those presented in the paper, since the data provided in this repository is presented for demonstration purposes. The full set of realizations and data can be provided by contacting corresponding author at [bruna.mattos.araujo@gmail.com](bruna.mattos.araujo@gmail.com).

