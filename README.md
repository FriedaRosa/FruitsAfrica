# Fruits Africa

Find below file descriptions and workflows to reproduce the study. For running the OUwie models, we relied on the HPC Cluster from UFZ and worked with SLURM scheduler. The R scripts can be run on a local machine as well, some adaptations are necessary to read in the data (in our analyses this was set up in the submit scripts that accompanied the R script).

Feel free to contact me for any questions and problems with reproducing the code via email: friederike.woelke\@gmail.com.

## R-scripts:

### 0_generateTDWGData.R:

	- requires (Data directory):

		- PalmTraits1.0.txt

		- average_fruit_width_Frieda_v2.csv (vertebrate-dispersed Angiosperm trait data)

		- BBM_summary_all.csv (BBM simulated fruit length for palms. 100 simulations on MCC tree, averaged)

		- kew_merged_Arecaceae.rds (subset WCVP data)

		- kew_merged_Angios.rds (subset WCVP data)

		- palms_in_tdwg3.csv (from Lim2020)

		- palms_tdwg_realms.csv (list of tdwg3 units and biogeographic realms)

		- TDWG_Environment_AllData_2019Feb.csv (from Lim 2020)

		- Phylacine_Trait_data.csv (from Lim 2020)

		- mammal_curr_occ.csv (from Lim 2020)

		- mammal_presnat_occ.csv (from Lim 2020)

		- frugivoreClassification.csv (from Lim 2020)

	- produces (Results directors):

		- palms_final = dataset that was used in OUwie models

		- tdwg_final = dataset that was used in ANOVA, LMs, spatial maps

### 1_OUwieModelsCluster.R: (was run on the HPC Cluster of the UFZ Leipzig)

	- requires:

		- palms_final.csv

		- set of 100 phylogenetic trees

		- MCC tree

	- produces:

		- 101x bestfit.rds

		- 101x raw.rds

		- 101x processed.rds

		- 1 simmap reconstruction for regimes

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

### 1b_OUwieModelsConvergenceAnalysis.R: 

	- requires:

		- HPC Cluster output (folder with 100 files, 1 for each phylogenetic tree. Not provided here)

		- issues.txt (produced in the first part, can be read in from the folder)

	- produces:

		- isses.txt (file with list of failed models)

		- merged_data_AVGL.txt (results from OUwie models excluding those that produced convergence issues)

		- Plot for convergence issues

### 2_OUwieFigs.R

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

### 3_ANOVA_LinearModels.R:

	- requires:

		- tdwg_final.csv

		- shapefile: TDWG_level3_Coordinates

	- Summary of analyses:

		- Data transformations

		- ANOVA FL (emp) \~ Africa

		- ANOVA FL (sim) \~ Africa

		- ANOVA FW (angios) \~ Africa

		- ANOVA BS change \~ Realms

		- Termplots LM Predictors

		- Data transformation and density plots

		- Violin plots, Regression plots, Effect sizes

		- Linear models comparison: with or without interaction

		- Africa / elsewhere subset analysis

		- correlation test for predictors

		- included world maps with plotted traits

### 4_SAR models

	- requires:

		- tdwg_final.csv

## Folder: Additional_Analyses

### 	- OUwie Americas and OUwie Asia

### 	- Sensitivity_Extinctions (10% of species were randomly removed)

### 	- Sensitivity_Extinctions2 (25, 50, 75% of fruits \> 4 cm from Africa removed)

### 	- SI_Figure_Phylogeny (code and data for SI figure)

##  
