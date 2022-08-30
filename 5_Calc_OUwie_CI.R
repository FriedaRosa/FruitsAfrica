### Script to calculate CIs for OUwie-estimated values ###


dd_master_sep <- read.table(file = "master_file_long_sep_emp_sim.txt", sep="\t", header=T)
str(dd_master_sep)
dd_master_sep$df

library(rstatix); library(dplyr); library(Rmisc); library(tidyr); library(qwraps2)

dd_master_sep %>% group_by(df) %>% get_summary_stats(type="full")

afr_emp <-dd_master_sep %>% filter(df == "observed", regime == "Africa")
afr_sim <-dd_master_sep %>% filter(df == "simulated", regime == "Africa")

CI(afr_emp$v_logTheta, ci=0.95)
CI(afr_sim$v_logTheta, ci=0.95)

library(tidyverse)

CI_tab <- dd_master_sep %>%
  group_by(df, regime, g_logTheta) %>%
  dplyr::summarise(N=n(),
            mean.ci = list(mean_ci(v_logTheta)),
            "Percent"= n_perc(v_logTheta > 0)) %>% 
  unnest_wider(mean.ci)




