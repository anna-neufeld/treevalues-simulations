# treevalues-simulations
Code / data / other files to reproduce figures in Neufeld, A., Gao, L.L., and Witten, D. (2021+), treevalues: Selective Inference for Regression Trees.

# Organization

The folder "Code to run simulations for paper" contains the code to actually run the simualations from Section 5 of the paper. 
- The file "Null_testing_sims.R" runs the simulation for Figure 4. Results from running this file are stored in 
"null_res_1-18-2021-5000.csv". 
- The file "run_power_rand.R" runs the simulation for Figure 5 (relies on code in "rand_power_indices.R"). Results from running this file are stored in 
"Full_Rand_Comps_1-18-2021-maxdepth3.csv". 
- The file "non_null_sums_RUN.R" runs the simulation for Table 2 and Figure 6 (relies on code in "non_null_CI_sims.R"). Results from running this file are stored in 
"non_null_CI_main-1-18-2021-maxdepth3.csv". 

Code to reproduce each figure using these stored results is stored in a separate, labeled file in the main directory. 
