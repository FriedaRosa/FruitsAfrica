# FruitsAfrica


 
R-scripts:

1) 0_generateTDWGData.R:
	- requires:
		- PalmTraits1.0
		- kew_merged_Arecaceae.rds (subset WCVP data)
		- palms_in_tdwg3.csv (from Lim2020)
		- palms_tdwg_realms.csv (list of tdwg3 units and biogeographic realms)
		- BBM_sim_.rds	(created in: 1_OUwieModelsCluster.R)
		- BBMsim_Sp_mean_100trees.csv (created in 2_BBM_spAvg_100trees.R)
		- TDWG_Environment_AllData_2019Feb.csv (from Lim 2020)
		- Phylacine_Trait_data.csv (from Lim 2020)
		- mammal_curr_occ.csv (from Lim 2020)
		- mammal_presnat_occ.csv (from Lim 2020)
		- frugivoreClassification.csv (from Lim 2020)
	- produces:
		- palms_final = dataset that was used in OUwie models
		- tdwg_final = dataset that was used in ANOVA, LMs, spatial maps

2) 1_OUwieModelsCluster.R:
	- requires: 
		- palms_final.csv 
		- set of 100 phylogenetic trees 
		- MCC tree
	- produces: 
		- 101x bestfit.rds
		- 101x raw.rds
		- 101x processed.rds
		- 1 simmap reconstruction for regimes
		- 101x BBM simulated traits (for averaging)
		- 1 file with models that show convergence issues 
		- 2 merged model tables for empirical and simulated traits (for 4_OUwieFigs.R)
	- Steps:
		- match species in trees and data
		- reorder trait data by order of species in trees
		- simmap reconstruction (1000x per tree)
		- averaging over 1000 simmap simulations per tree
		- BBM trait simulations
		- OUwie models for empirical and simulated data 
		- Formating the output (create model fit table) 
		- remove models with convergence problems and evaluate problematic models
		- model comparison with aicc weights
		- save raw, processed, bestfit tables

3) 2_BBM_spAvg_100trees.R
	- requires:
		- 101x BBM simulations 
	- produces: 
		- BBMsim_Sp_mean_100trees.csv

4) 3_BBM_rootstates_100trees.R
	- requires:
		- 101x BBM simulations
	- produces: 
		- list with all simulations, all root states and their probability densities
		(FS_at_root_100sims.csv)

5) 4_OUwieFigs.R
	- requires:
		- 2 merged tables for empirical and simulated traits (from 1_OUwieModelsCluster.R)
	- produces: 
		- master_file_long_sep_emp_sim.txt (for 5_Calc_OUwie_CI.R)
	- Steps:
		- read both tables in and transform dataframe to long-format
		- create ggplot dataset with split columns for variables and estimates 
		to enable grouping for plotting
		- Plot selected models
		- Plot dot-line plots for theta and sigma

6) 5_Calc_OUwie_CI.R
	- requires: 
		- master_file_long_sep_emp_sim.txt (from 4_OUwieFigs.R)
	- produces: 
		- table with confidence intervals around thetas


7) 6_ ANOVA_LinearModels.R:
	- requires: 
		- tdwg_final.csv
	- Summary of analyses:
		- Data transformations
		- ANOVA FL (emp) ~ Africa
		- ANOVA FL (sim) ~ Africa
		- ANOVA BS change ~ Realms
		- Termplots LM Predictors
		- Data transformation and density plots
		- Violin plots, Regression plots,  Effect sizes
		- Linear models comparison: with or without interaction
		- Africa / elsewhere subset analysis
		- correlation test for predictors 

8) 7_3palms_maps.R
	- requires:
		- tdwg_final.csv
		- shapefile: TDWG_level3_Coordinates
	- Steps:
		- Body mass change ANOVA and map
		- Fruit length maps (empirical, simulated)

9) Sensitivity_Extinctions:
	- Folder with scripts for sensitivity analysis where 10% of species were randomly removed
