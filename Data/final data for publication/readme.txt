
This dataset includes consists of two .csv (Comma Separated Values) files, two .rds (R Data Serialization) files, and an R script file.

File 1: data_resubmission.csv
This file contains raw measurements from individual birds (N = 2,057) included for analysis in this study. Variables include:
sample number, bird ID (band number), population (municipality), sampling year, sex, latitude, longitude, and color measurements. The color measurements are abbreviated as follows:
t= throat, r = breast, b = belly, v = vent; avg.bright = brightness, hue = hue, chrom = chroma

File 2: pairwise_distance_correctpopulation.csv
This file contains the pairwise genetic distance, calculated as described in the Methods section. The variables are the two focal populations (pop1 and pop2), and mean genetic distance between those populations. 


File 3: results_10k_bootstraps_correctpopulations.RDS
This file contains the results of analyses using 10,000 bootstrapped samples. In each bootstrap iteration, we resampled the same number of individuals from each population. This was done to assess the robustness of our results to unequal sampling across populations.
List item 1 ($boot_n): number of bootstrap iterations
List item 2 ($boot_summ): The raw data of the bootstrapped samples for each iteration
List item 3 ($mean_traits): The mean trait values and the PINTc values for each population and sex based on the bootstrap replicates


File 4: PINTanalysis_results_w_Fst_scaled.rds
This file contains the results of linear models applied to each of the 10,000 bootstrapped replicates and is used to generate the mean and confidence intervals for the results shown in Table 1. The model specifications are shown in annotations in the script.

File 5: analysis_EVOLUTION.R
This is R script for reproducing all results presented in the paper. 