### =================================================================================================== ###
## New script for simulating 'neutral' fruit size evolution in palms and subsequent analyses and figures ##
# Part 0: Simulating bounded brownian motion trait evolution 100-times on the MCC tree for palms.
# Part I: Creating a new TDWG3-data set with the simulated fruit sizes.
# Part II: New world maps with the same scale as the empirical fruit size data
# Part III: Linear models
#         a) with the newly simulated data, averaged from 100 simulations on the MCC tree
#         b) with the old data, where we simulated once for 100 phylogenetic trees
### =================================================================================================== ###
rm(list=ls())


## Part 0: Simulating bounded brownian motion trait evolution 100-times on the MCC tree for palms. ==========
setwd("C:/Users/Friederike/Desktop/BBM_100sim_MCC_cluster")

# libraries ===================================
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(evobiR))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(BBMV))

# read tree ===================================
tree <- read.nexus("TREE")

# read data ===================================
dd <- read.csv("Step1_duplicated_regimes_checked.csv", stringsAsFactors = T)

dd2 <- dd[,c("SpecName","accGenus","PalmTribe","PalmSubfamily","AverageFruitLength_cm","Americas",
             "Asia","Africa","Madagascar","Pacific","accRealm", "log_BBM_mean", "log_BBM_MCC")]
dd <- dd2 %>% distinct()

dd$SpecName <- as.character(dd$SpecName)
dd$AverageFruitLength_cm1 <- as.numeric(as.character(dd$AverageFruitLength_cm))

dd2 <- dd[!is.na(dd["AverageFruitLength_cm1"]),]
dd <- as.data.frame(dd2)

dd$log_FS <- log(dd$AverageFruitLength_cm1) #log-transformation
row.names(dd) <- dd$SpecName


# Match species in data and tree  ================================
setdiff(tree$tip.label, dd$SpecName) #502 sp.
tree2 <- drop.tip(tree, setdiff(tree$tip.label, row.names(dd)))
tree <- ladderize(tree2) # 2037 tips
tree <- tree2
data2 <- ReorderData(tree, dd, taxa.names = "row names")

## log-transformed trait
log_avgl <- data2$log_FS
names(log_avgl) <- row.names(data2)

## Sort data like phylogeny to prevent errors
log_avgl <- ReorderData(tree, log_avgl, taxa.names="names")

## BBM simulations =================================== 
# = bounded brownian motion between hard reflective bounds

# step 1. estimate z0, sigma2 and bounds from empirical data using BBM package:
## 1.1 create likelihood function from empirical data
BBM <- lnL_BBMV(tree, log_avgl, Npts=100, bounds=c(min(log_avgl),max(log_avgl)), a=0,b=0,c=0)

## 1.2. find maximum likelihood function
# mle_BBM <- find.mle_FPK(BBM, safe=F) # safe = F ~ 20min computing time; identical results as safe=T
# saveRDS(mle_BBM, "out/mle_BBM.rds")
mle_BBM <- readRDS("out/mle_BBM.rds")

# 1.3 extract estimated sigsq and transform to sigma (needed for BBM trait simulations)
BBMsigsq <- mle_BBM$par$sigsq       # 0.04855282
BBMsigma <- sqrt(BBMsigsq)          # 0.220347

# 1.4 extract estimated trait at the root (z0), find the one with the maximum density and save as 'root' for trait simulations
root_prob <- mle_BBM$root
root_state <- root_prob %>% slice(which.max(density)) # root state with highest probability:0.9823201, density: 0.01837733
root <- root_state$root             # 0.9823201

# 1.5 extract bounds
bounds <- mle_BBM$par_fixed$bounds  # -1.203973  3.401197

# step 2. simulate traits for species at the tips of the tree based on model fitted parameter estimates
# Loop for 100 simulations ==========================

# runs <- c(1:100)
# simulations_list <- list()
# 
# for (i in 1:length(runs)){
#   print(i)
#   BBM_sim <- Sim_FPK(tree, x0 = root, V = rep(0, 100), bounds=bounds, sigma = BBMsigma)
#   simulations_list[[i]] <- BBM_sim
# }
# 
# saveRDS(simulations_list, "out/BBM_simulations_MCC_100sim_list.rds")

simulations_list <- readRDS("out/BBM_simulations_MCC_100sim_list.rds")

# reduce list to dataframe
BBM_df <- reduce(simulations_list, bind_rows)

# average over simulations for species
BBM_MCC_mean_100sims <- apply(BBM_df, 2, mean)
range(BBM_MCC_mean_100sims)               # 0.7790414 1.5768057

# save to file
# write.csv(BBM_MCC_mean_100sims, "out/BBM_MCC_mean_100sim_log.csv")

### =================================================================================================== ###

## Part I: Creating a new TDWG3-data set with the simulated fruit sizes ==========

## Directories ========================
rm(list = ls())
setwd("~/GitHub/")
main.dir <- "FruitsAfrica"
data.dir <- file.path(main.dir, "Data")
res.dir <- file.path(main.dir, "Results/")
src.dir <- file.path(main.dir, "src/")
fig.dir <- file.path(main.dir, "Figures/")
options(stringsAsFactors =FALSE)

## Packages ========================
library(plyr); library(dplyr); library(reshape2)
library(sp); library(rgdal)
library(car); library(MuMIn)
library(readxl)
library(picante)

## Data import ========================
# 1.1 Palm traits ========================
palm_trait0 <- read.csv(file.path(data.dir, "PalmTraits_1.0.txt"), header=T, sep="\t")
palm_trait0$SpecName <- gsub(palm_trait0$SpecName, pattern= " ", replacement = "_")
palm_trait <- palm_trait0 %>% dplyr::select(SpecName, accGenus, PalmTribe, PalmSubfamily, AverageFruitLength_cm) %>% filter_at(vars(AverageFruitLength_cm), all_vars(!is.na(.)))

# 1.2 Distribution data on botanical country (TDWG3) scale ================

# Skip this part and use this .rds file to save computing time:
dd_merge <- readRDS(file.path(data.dir, "kew_merged_Arecaceae.rds"))

# =================================================================================================================================

# # Names data =====
# names <- read.table(file.path(data.dir, "wcvp_names.txt"), sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# names <- names %>% filter(taxon_rank == "Species" & family == "Arecaceae")
# names <- names %>% dplyr::select(plant_name_id, taxon_status, family, genus, species, taxon_name, accepted_plant_name_id)
# names$SpecName <- paste(names$genus, names$species, sep="_")
# 
# # short comparison:
# names_inclSyn <- names %>% select(accepted_plant_name_id, taxon_status, plant_name_id, genus, species, taxon_name, SpecName) %>% filter(!is.na(accepted_plant_name_id))
# names_accepted <- names %>% select(accepted_plant_name_id, taxon_status, plant_name_id, genus, species, taxon_name, SpecName) %>% filter(taxon_status == "Accepted")
# 
# # Distribution data ======
# dist <- read.table(file.path(data.dir, "wcvp_distribution.txt"), sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# dist <- dist %>% dplyr::select(plant_name_id, continent, region, area_code_l3, area)
# 
# # Merge ====
# dd_merge <- merge(names, dist, by="plant_name_id", all.x=TRUE)
# dd_merge <- unique(dd_merge)
# 
# #x <-dd_merge %>% filter(!is.na(area_code_l3) & !is.na(SpecName))
# 
# # Counts =====
# length(unique(dd_merge$family)) # Arecaceae
# length(unique(dd_merge$genus))  # 453 palm genera
# length(unique(dd_merge$area_code_l3)) #230 tdwg3
# length(unique(dd_merge$SpecName))# 7459 species
# length(unique(dd_merge$accepted_plant_name_id)) #2832
# length(unique(dd_merge$plant_name_id)) #7550
# 
# # Save to .RDS ====
# saveRDS(dd_merge, file.path(data.dir, "kew_merged_Arecaceae.rds"))

# =================================================================================================================================

dd_merge <- dd_merge %>% select(plant_name_id, taxon_status, family, genus, species, taxon_name, accepted_plant_name_id, SpecName, area_code_l3, area, region, continent)

regions <- dd_merge %>% 
  dplyr::select(area_code_l3, area, region, continent) %>% 
  mutate_all(na_if,"") %>% 
  filter_at(vars(area_code_l3), any_vars(!is.na(.))) %>% 
  distinct(.)


# reduced to those with distribution data (!= NA)
dd_merge <- dd_merge %>%
  filter_at(vars(area_code_l3), all_vars(!is.na(.)))

length(setdiff(palm_trait$SpecName, dd_merge$SpecName)) 
length(setdiff(dd_merge$SpecName, palm_trait$SpecName))


# 1.3 Palm occurrences ========================

palms_occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
length(unique(palms_occ$SpecName)) #2557 species

palms_realms <- read.csv(file.path(data.dir, "palms_tdwg_realms.csv")) #1885 species
palms_realms <- palms_realms %>% select(SpecName, accAfrica, area_code_l3, Americas, Asia, Africa, Madagascar, Pacific, accRealm)

palms_traits_v2 <- merge(palm_trait, palms_occ, by="SpecName", all.x=T)
palms_traits_v2 <- merge(palms_traits_v2, palms_realms, by.x=c("SpecName","Area_code_L3"), by.y=c("SpecName","area_code_l3"), all.x=T)
palms_traits_v2 <- merge(palms_traits_v2, regions, by.x=c("Area_code_L3"), by.y=c("area_code_l3"), all.x=T)
palms_traits_v2 <- palms_traits_v2 %>% filter_at(vars(AverageFruitLength_cm, Area_code_L3), all_vars(!is.na(.))) # 2051 species

# 2.1 Brownian-motion simulated palm traits ========================

BBM_summary <- read.csv(file.path(data.dir, "BBM_summary_all.csv"))
str(BBM_summary)

## a) New data set ============
newBBM_MCC_mean <- data.frame(log_BBM_MCC_avg = BBM_summary$new_BBM_MCC_mean_100sim_log, SpecName = BBM_summary$X)

## b) Old: Species-Mean from 100 trees (1 simulation each) ============
oldBBM_MCC_trees_mean <- data.frame(log_BBM_MCC_100trees_avg = BBM_summary$old_BBM_100trees_mean_1sim_log, SpecName = BBM_summary$X)

## c) Old: One simulation for fruit size for each tree (100 columns, one for each tree) ============
old_BBM_100trees_sep <- data.frame(BBM_summary[,c(1,5:104)])

## For a) and b) =============
palm_traits_v3 <- merge(palms_traits_v2, newBBM_MCC_mean, by =c("SpecName"))
palm_traits_v3 <- merge(palm_traits_v3, oldBBM_MCC_trees_mean, by = c("SpecName"))

length(unique(palm_traits_v3$SpecName)) # 2037 species with fruit data

# some corrections for NA values based on other columns:
palm_traits_v3$accRealm[palm_traits_v3$continent == "AFRICA"] <- "Africa"
palm_traits_v3$accRealm[palm_traits_v3$continent == "ASIA-TROPICAL" | palm_traits_v3$continent == "ASIA-TEMPERATE" | palm_traits_v3$continent =="EUROPE"] <- "Asia"
palm_traits_v3$accRealm[palm_traits_v3$continent == "SOUTHERN AMERICA" | palm_traits_v3$continent == "NORTHERN AMERICA"] <- "Americas"
palm_traits_v3$accRealm[palm_traits_v3$continent == "AUSTRALASIA" | palm_traits_v3$continent == "PACIFIC"] <- "Pacific"
palm_traits_v3$accRealm[palm_traits_v3$Area_code_l3 == "MDG"] <- "Madagascar"

palm_traits_v3$accAfrica[palm_traits_v3$continent == "AFRICA"] <- "1"
palm_traits_v3$accAfrica[palm_traits_v3$continent != "AFRICA"] <- "0"
palm_traits_v3$accAfrica[palm_traits_v3$Area_code_L3 == "MDG"] <- "0"

palm_traits_v3$Americas[palm_traits_v3$accRealm == "Americas"] <- ("1")
palm_traits_v3$Americas[palm_traits_v3$accRealm != "Americas"] <- ("0")
palm_traits_v3$Pacific[palm_traits_v3$accRealm == "Pacific"] <- ("1")
palm_traits_v3$Pacific[palm_traits_v3$accRealm != "Pacific"] <- ("0")
palm_traits_v3$Asia[palm_traits_v3$accRealm == "Asia"] <- ("1")
palm_traits_v3$Asia[palm_traits_v3$accRealm != "Asia"] <- ("0")
palm_traits_v3$Africa[palm_traits_v3$accRealm == "Africa"] <- ("1")
palm_traits_v3$Africa[palm_traits_v3$accRealm != "Africa"] <- ("0")
palm_traits_v3$Madagascar[palm_traits_v3$Area_code_L3 == "MDG"] <- ("1")
palm_traits_v3$Madagascar[palm_traits_v3$Area_code_L3 != "MDG"] <- ("0")

# save final dataframe as .rds and .csv
#saveRDS(palm_traits_v3, file=paste0(res.dir, "/palms_final.rds"))
#write.csv(palm_traits_v3, file=paste0(res.dir, "/palms_final.csv"))

palms <- palm_traits_v3

# 3. Environmental data ========================

tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))

env.vars <- c("LEVEL_3_CO", "LAT", "LONG",  "CONTINENT", "CH_Mean", "CH_PercCover", "bio1_mean", "bio4_mean", "bio15_mean", "bio12_mean")
tdwg_env2<- tdwg_env[,c(env.vars)]
rownames(tdwg_env2) <- tdwg_env2$LEVEL_3_CO
palms_env <- unique(merge(palm_traits_v3, tdwg_env2, by.x=c("Area_code_L3", "continent"), by.y=c("LEVEL_3_CO", "CONTINENT")))

# 4. Mammal / Frugivore data ========================	

phylacine_trait <- read.csv(file.path(data.dir, "Phylacine_Trait_data.csv")) ## Phylacine trait = for body masses of mammals
mammal_curr_comb_occ <- read.csv(file.path(data.dir, "mammal_curr_occ.csv")) ## mammal curr = list of botanic countries and extant mammals (from IUCN range maps)
mammal_presnat_comb_occ <- read.csv(file.path(data.dir, "mammal_presnat_occ.csv")) ## mammal presnat = list of botanic countries and extant+extinct mammals
frugivoreClass <- read.csv(file.path(data.dir, "frugivoreClassification.csv")) ## frugivoreClass = mammal list with indication of proportion frugivory (family, order)

## step 1. make mammal country lists (1. current, 2. extinct)

mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_comb_occ$LEVEL_3_CO),   # 191 botanic countries
                                         unique(mammal_curr_comb_occ$LEVEL_3_CO)))     # 190 botanic countries

## 2. make list of countries where mammals and species overlap

mammal_plant_intersect <- intersect(unique(palms$Area_code_L3), mammal_countrylist)  #191 botanic countries
sort(mammal_plant_intersect, decreasing = FALSE)


## 3. aggregate maximum + median angiosperm and (empirical, simulated) palm traits on tdwg3; merge both datasets together
tdwg_palms_summary <- ddply(.data = subset(palms, Area_code_L3 %in% mammal_plant_intersect),
                            .variables = .(Area_code_L3),
                            .fun = dplyr::summarise,
                            sp_richness_palms = length(SpecName),
                            accRealm = accRealm,
                            accAfrica = accAfrica,
                            area = area,
                            
                            # No gap filling
                            meanFL_palms = mean(AverageFruitLength_cm, na.rm = T),
                            medianFL_palms = median(AverageFruitLength_cm, na.rm = T),
                            max95FL_palms = exp(quantile(log(AverageFruitLength_cm), probs = 0.95, na.rm = T)),
                            
                            # species averages from 100 simulations (NEW, a)
                            log_meanFL_palms_BBMmean = mean(log_BBM_MCC_100trees_avg, na.rm = T),
                            log_medianFL_palms_BBMmean = median(log_BBM_MCC_100trees_avg, na.rm = T),
                            log_max95FL_palms_BBMmean = quantile(log_BBM_MCC_100trees_avg, 0.95),
                            
                            # MCC results from 1 simulation (OLD, b)
                            log_meanFL_palms_BBM_MCC = mean(log_BBM_MCC_avg, na.rm = T),
                            log_medianFL_palms_BBM_MCC = median(log_BBM_MCC_avg, na.rm = T),
                            log_max95FL_palms_BBM_MCC = quantile(log_BBM_MCC_avg, 0.95))
tdwg_palms_summary <- strings2factors(tdwg_palms_summary)
str(tdwg_palms_summary)
tdwg_palms_summary <- unique(tdwg_palms_summary)
#saveRDS(tdwg_palms_summary, file.path(res.dir, "tdwg_palms_summary.rds"))

# reduced to tdwg3 with >3 palms ==================

palms3 <- tdwg_palms_summary %>% filter(sp_richness_palms > 3)
palms3tdwg3 <- unique(palms3$Area_code_L3) #131 botanic countries

# add to environmental data ===============0
tdwg_final <- read.csv(file.path(data.dir, "tdwg_final.csv"))
names(tdwg_final)

tdwg_final$log_max95FL_palms_BBMmean <- NULL
tdwg_final$log_medianFL_palms_BBM_MCC <- NULL 
tdwg_final$log_max95FL_palms_BBMmean <- NULL 

tdwg_final2 <- merge(tdwg_final, palms3[,c("Area_code_L3", "log_max95FL_palms_BBM_MCC")], by.x=c("LEVEL_3_CO"), by.y=c("Area_code_L3")) #new BBM simulations (100 simulations on MCC)
tdwg_final2 <- merge(tdwg_final2, palms3[,c("Area_code_L3", "log_max95FL_palms_BBMmean")], by.x=c("LEVEL_3_CO"), by.y=c("Area_code_L3")) #old BBM simulations (100 trees)


# write.csv(tdwg_final2, file=paste0(data.dir, "/tdwg_final_newBBM.csv"))

### For c) ============== (merging data and running linear models)

## Loop to run over 100 trees-simulation:

runs <- c(1:ncol(old_BBM_100trees_sep))
output <- list()
coeff_output <- list()
list_significant <- list()
for (i in 2:length(runs)){
  print(i)
  BBM_trees_sep <- data.frame(BBM_sim = old_BBM_100trees_sep[,i], SpecName = old_BBM_100trees_sep$X)
  palm_traits_v3 <- merge(palms_traits_v2, BBM_trees_sep, by =c("SpecName"))
  
  
  # some corrections for NA values based on other columns:
  palm_traits_v3$accRealm[palm_traits_v3$continent == "AFRICA"] <- "Africa"
  palm_traits_v3$accRealm[palm_traits_v3$continent == "ASIA-TROPICAL" | palm_traits_v3$continent == "ASIA-TEMPERATE" | palm_traits_v3$continent =="EUROPE"] <- "Asia"
  palm_traits_v3$accRealm[palm_traits_v3$continent == "SOUTHERN AMERICA" | palm_traits_v3$continent == "NORTHERN AMERICA"] <- "Americas"
  palm_traits_v3$accRealm[palm_traits_v3$continent == "AUSTRALASIA" | palm_traits_v3$continent == "PACIFIC"] <- "Pacific"
  palm_traits_v3$accRealm[palm_traits_v3$Area_code_l3 == "MDG"] <- "Madagascar"
  
  palm_traits_v3$accAfrica[palm_traits_v3$continent == "AFRICA"] <- "1"
  palm_traits_v3$accAfrica[palm_traits_v3$continent != "AFRICA"] <- "0"
  palm_traits_v3$accAfrica[palm_traits_v3$Area_code_L3 == "MDG"] <- "0"
  
  palm_traits_v3$Americas[palm_traits_v3$accRealm == "Americas"] <- ("1")
  palm_traits_v3$Americas[palm_traits_v3$accRealm != "Americas"] <- ("0")
  palm_traits_v3$Pacific[palm_traits_v3$accRealm == "Pacific"] <- ("1")
  palm_traits_v3$Pacific[palm_traits_v3$accRealm != "Pacific"] <- ("0")
  palm_traits_v3$Asia[palm_traits_v3$accRealm == "Asia"] <- ("1")
  palm_traits_v3$Asia[palm_traits_v3$accRealm != "Asia"] <- ("0")
  palm_traits_v3$Africa[palm_traits_v3$accRealm == "Africa"] <- ("1")
  palm_traits_v3$Africa[palm_traits_v3$accRealm != "Africa"] <- ("0")
  palm_traits_v3$Madagascar[palm_traits_v3$Area_code_L3 == "MDG"] <- ("1")
  palm_traits_v3$Madagascar[palm_traits_v3$Area_code_L3 != "MDG"] <- ("0")
  
  # save final dataframe as .rds and .csv
  #saveRDS(palm_traits_v3, file=paste0(res.dir, "/palms_final.rds"))
  #write.csv(palm_traits_v3, file=paste0(res.dir, "/palms_final.csv"))
  
  palms <- palm_traits_v3
  # 3. Environmental data ========================
  
  tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))
  env.vars <- c("LEVEL_3_CO", "LAT", "LONG",  "CONTINENT", "CH_Mean", "CH_PercCover", "bio1_mean", "bio4_mean", "bio15_mean", "bio12_mean")
  tdwg_env2<- tdwg_env[,c(env.vars)]
  rownames(tdwg_env2) <- tdwg_env2$LEVEL_3_CO
  
  
  palms_env <- unique(merge(palm_traits_v3, tdwg_env2, by.x=c("Area_code_L3", "continent"), by.y=c("LEVEL_3_CO", "CONTINENT")))
  
  
  # 4. Mammal / Frugivore data ========================	
  
  phylacine_trait <- read.csv(file.path(data.dir, "Phylacine_Trait_data.csv")) ## Phylacine trait = for body masses of mammals
  mammal_curr_comb_occ <- read.csv(file.path(data.dir, "mammal_curr_occ.csv")) ## mammal curr = list of botanic countries and extant mammals (from IUCN range maps)
  mammal_presnat_comb_occ <- read.csv(file.path(data.dir, "mammal_presnat_occ.csv")) ## mammal presnat = list of botanic countries and extant+extinct mammals
  frugivoreClass <- read.csv(file.path(data.dir, "frugivoreClassification.csv")) ## frugivoreClass = mammal list with indication of proportion frugivory (family, order)
  
  ## step 1. make mammal country lists (1. current, 2. extinct)
  mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_comb_occ$LEVEL_3_CO),   # 191 botanic countries
                                           unique(mammal_curr_comb_occ$LEVEL_3_CO)))     # 190 botanic countries
  ## 2. make list of countries where mammals and species overlap
  mammal_plant_intersect <- intersect(unique(palms$Area_code_L3), mammal_countrylist)  #191 botanic countries
  sort(mammal_plant_intersect, decreasing = FALSE)
  
  
  ## 3. aggregate maximum + median angiosperm and (empirical, simulated) palm traits on tdwg3; merge both datasets together
  tdwg_palms_summary <- ddply(.data = subset(palms, Area_code_L3 %in% mammal_plant_intersect),
                              .variables = .(Area_code_L3),
                              .fun = dplyr::summarise,
                              sp_richness_palms = length(SpecName),
                              accRealm = accRealm,
                              accAfrica = accAfrica,
                              area = area,
                              
                              # No gap filling
                              meanFL_palms = mean(AverageFruitLength_cm, na.rm = T),
                              medianFL_palms = median(AverageFruitLength_cm, na.rm = T),
                              max95FL_palms = exp(quantile(log(AverageFruitLength_cm), probs = 0.95, na.rm = T)),
                              
                              # species averages from 100 simulations 
                              log_meanFL_palms_BBMmean = mean(BBM_sim, na.rm = T),
                              log_medianFL_palms_BBMmean = median(BBM_sim, na.rm = T),
                              log_max95FL_palms_BBMmean = quantile(BBM_sim, 0.95))
  tdwg_palms_summary <- strings2factors(tdwg_palms_summary)
  str(tdwg_palms_summary)
  tdwg_palms_summary <- unique(tdwg_palms_summary)
  #saveRDS(tdwg_palms_summary, file.path(res.dir, "tdwg_palms_summary.rds"))
  
  # reduced to tdwg3 with >3 palms ==================
  palms3 <- tdwg_palms_summary%>% filter(sp_richness_palms > 3)
  palms3tdwg3 <- unique(palms3$Area_code_L3) #131 botanic countries
  
  
  # add to environmental data ===============0
  
  tdwg_final <- read.csv(file.path(data.dir, "tdwg_final.csv"))
  names(tdwg_final)
  
  tdwg_final$log_max95FL_palms_BBMmean <- NULL
  tdwg_final$log_medianFL_palms_BBM_MCC <- NULL 
  tdwg_final$log_max95FL_palms_BBMmean <- NULL 
  
  tdwg_final2 <- merge(tdwg_final, palms3[,c("Area_code_L3", "log_meanFL_palms_BBMmean", "log_medianFL_palms_BBMmean", "log_max95FL_palms_BBMmean" )], by.x=c("LEVEL_3_CO"), by.y=c("Area_code_L3")) #new BBM simulations (100 simulations on MCC)
  ## Data handling  ========================
  tdwg_final2$accRealm[tdwg_final2$LEVEL_3_CO == "MDG"] <- "Madagascar"
  tdwg_final2$LEVEL_3_CO <- as.factor(tdwg_final2$LEVEL_3_CO)
  #tdwg_final2$accAfrica <- factor(tdwg_final2$accAfrica, levels = c("1", "0"))
  
  # Set negative values (i.e., positive body mass changes) to zero before log-transformation:
  tdwg_final_0 <- tdwg_final2
  tdwg_final_0$mean_abs_BSchange[tdwg_final_0$mean_abs_BSchange < 0] <- 0
  tdwg_final_0$mean_perc_BSchange[tdwg_final_0$mean_perc_BSchange < 0] <- 0
  
  tdwg_final_0$mean_abs_BSchange <- tdwg_final_0$mean_abs_BSchange/1000
  
  # sqrt transform
  tdwg_final_0$sqrt_change <- sqrt(tdwg_final_0$mean_abs_BSchange)
  
  # Scale:
  normalized<-function(y) {
    x<-y[!is.na(y)]
    x<-(x - min(x)) / (max(x) - min(x))
    y[!is.na(y)]<-x
    return(y)
  }
  
  
  lm_full1 <- lm(normalized(exp(log_max95FL_palms_BBMmean)) ~ factor(accAfrica)*
                   normalized(sqrt(CH_Mean)) +
                   normalized((sqrt_change)) +
                   normalized(sqrt(curr_max95BodySize)) +
                   normalized(I(MAT^2)) +
                   normalized(sqrt(TempSeas)) +
                   normalized(sqrt(AnnualPrec)) +
                   normalized(sqrt(PrecSeas)),
                 data = tdwg_final_0, na.action = na.exclude)
  summary(lm_full1)
  
  output[[i]] <- summary(lm_full1)
  data.frame(log_BBM_MCC_avg = BBM_summary$new_BBM_MCC_mean_100sim_log, SpecName = BBM_summary$X)
  sign_vars <- data.frame(variable = rownames(summary(lm_full1)$coefficients), 
                          P_sign = summary(lm_full1)$coefficients[,4] < 0.05,
                          P = round(summary(lm_full1)$coefficients[,4],4),
                          coeff = round(summary(lm_full1)$coefficients[,1],3),
                          R2 = summary(lm_full1)[8],
                          adj.R2 = summary(lm_full1)[9])
  
  coeff_output[[i]] <- sign_vars
  significant <- sign_vars %>% filter(P_sign == T)
  significant$tree <- i
  
  list_significant[[i]] <- significant
}

lms_all_100trees <- do.call(rbind, list_significant)

#write.csv(lms_all_100trees, file.path(res.dir, "lms_all_100trees.csv"))

## Count how many times each variable was significant =====

lms_all_100trees %>% group_by(variable) %>% dplyr::summarise(coeff_avg =mean(coeff),
                                                      coeff_sd = sd(coeff),
                                                      p_avg = mean(P),
                                                      p_sd = sd(P),
                                                      n = n())

# variable                                       coeff_avg coeff_sd   p_avg   p_sd     n
# <chr>                                              <dbl>    <dbl>   <dbl>  <dbl> <int>
#   1 (Intercept)                                    0.450    0.532 0.00913 0.0118    51
# 2 factor(accAfrica)1                              -0.281    0.356 0.0185  0.0170    21
# 3 factor(accAfrica)1:normalized(sqrt(CH_Mean))     0.221    0.687 0.0173  0.0150    16
# 4 normalized((sqrt_change))                        0.184    0.342 0.0102  0.0136    50
# 5 normalized(I(MAT^2))                             0.303    0.274 0.0130  0.0103    17
# 6 normalized(sqrt(AnnualPrec))                     0.405    0.214 0.0141  0.0147    31
# 7 normalized(sqrt(CH_Mean))                       -0.113    0.511 0.0146  0.0159     8
# 8 normalized(sqrt(curr_max95BodySize))             0.110    0.414 0.0104  0.0141    38
# 9 normalized(sqrt(PrecSeas))                       0.641    0.219 0.00774 0.0123    39
# 10 normalized(sqrt(TempSeas))                      -0.201   0.444 0.00762 0.0116    49


## Plotting some histograms of the coefficient estimates ==============

par(mfrow=c(3,3))
lms_all_100trees %>% filter(variable == "factor(accAfrica)1" ) %>% with(hist(coeff, main = "accAfrica"))
lms_all_100trees %>% filter(variable == "factor(accAfrica)1:normalized(sqrt(CH_Mean))" ) %>% with(hist(coeff, main = "factor(accAfrica)1:normalized(sqrt(CH_Mean))"))
lms_all_100trees %>% filter(variable == "normalized((sqrt_change))" ) %>% with(hist(coeff, main = "normalized((sqrt_change))"))
lms_all_100trees %>% filter(variable == "normalized(I(MAT^2))" ) %>% with(hist(coeff, main = "normalized(I(MAT^2))"))
lms_all_100trees %>% filter(variable == "normalized(sqrt(AnnualPrec))" ) %>% with(hist(coeff, main = "normalized(sqrt(AnnualPrec))"))
lms_all_100trees %>% filter(variable == "normalized(sqrt(CH_Mean))" ) %>% with(hist(coeff, main = "normalized(sqrt(CH_Mean))"))
lms_all_100trees %>% filter(variable == "normalized(sqrt(curr_max95BodySize))" ) %>% with(hist(coeff, main = "normalized(sqrt(curr_max95BodySize))"))
lms_all_100trees %>% filter(variable == "normalized(sqrt(PrecSeas))" ) %>% with(hist(coeff, main = "normalized(sqrt(PrecSeas))"))
lms_all_100trees %>% filter(variable == "normalized(sqrt(TempSeas))" ) %>% with(hist(coeff, main = "normalized(sqrt(TempSeas))"))

### =================================================================================================== ###

## Part II: Linear models ================================

## Data import ========================
tdwg_final3 <- read.csv(file.path(data.dir, "tdwg_final_newBBM.csv"), header=T, sep=",")

## Data handling  ========================
tdwg_final3$accRealm[tdwg_final3$LEVEL_3_CO == "MDG"] <- "Madagascar"
tdwg_final3$LEVEL_3_CO <- as.factor(tdwg_final3$LEVEL_3_CO)
#tdwg_final3$accAfrica <- factor(tdwg_final3$accAfrica, levels = c("1", "0"))

# Set negative values (i.e., positive body mass changes) to zero before log-transformation:
tdwg_final_0 <- tdwg_final3
tdwg_final_0$mean_abs_BSchange[tdwg_final_0$mean_abs_BSchange < 0] <- 0
tdwg_final_0$mean_perc_BSchange[tdwg_final_0$mean_perc_BSchange < 0] <- 0

tdwg_final_0$mean_abs_BSchange <- tdwg_final_0$mean_abs_BSchange/1000

# sqrt transform
tdwg_final_0$sqrt_change <- sqrt(tdwg_final_0$mean_abs_BSchange)


# Scale:
normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
}


lm_full1 <- lm(normalized(log_max95FL_palms_BBM_MCC) ~ factor(accAfrica)*
                 normalized(sqrt(CH_Mean)) +
                 normalized((sqrt_change)) +
                 normalized(sqrt(curr_max95BodySize)) +
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) +
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)),
               data = tdwg_final_0, na.action = na.exclude)
summary(lm_full1)

# Call:
#   lm(formula = normalized(log_max95FL_palms_BBM_MCC) ~ factor(accAfrica) * 
#        normalized(sqrt(CH_Mean)) + normalized((sqrt_change)) + normalized(sqrt(curr_max95BodySize)) + 
#        normalized(I(MAT^2)) + normalized(sqrt(TempSeas)) + normalized(sqrt(AnnualPrec)) + 
#        normalized(sqrt(PrecSeas)), data = tdwg_final_0, na.action = na.exclude)
# 
# Residuals:
#   Min       1Q         Median       3Q      Max 
# -0.38429  -0.07167    0.00944   0.08577   0.29665 
# 
# Coefficients:
#                                               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   0.55441    0.13215   4.195  5.6e-05 ***
# factor(accAfrica)1                           -0.15578    0.10046  -1.551 0.123927    
# normalized(sqrt(CH_Mean))                    -0.01012    0.12122  -0.083 0.933629    
# normalized((sqrt_change))                     0.00469    0.06541   0.072 0.942967    
# normalized(sqrt(curr_max95BodySize))         -0.27573    0.07973  -3.458 0.000779 ***
# normalized(I(MAT^2))                          0.18995    0.09697   1.959 0.052711 .  
# normalized(sqrt(TempSeas))                    0.13815    0.08594   1.608 0.110832    
# normalized(sqrt(AnnualPrec))                 -0.14610    0.09695  -1.507 0.134753    
# normalized(sqrt(PrecSeas))                    0.01912    0.12000   0.159 0.873721    
# factor(accAfrica)1:normalized(sqrt(CH_Mean))  0.17362    0.16764   1.036 0.302652    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1287 on 108 degrees of freedom
# (12 Beobachtungen als fehlend gelöscht)
# Multiple R-squared:  0.3985,	Adjusted R-squared:  0.3484 
# F-statistic:  7.95 on 9 and 108 DF,  p-value: 6.065e-09

## Simulated fruit sizes: perfect layout!

my_colors <- setNames(c("blue", "#bdbdbd"), c("1", "0"))

pLMsim <- plot_model(lm_full1, 
                     colors = c("#440154FF","#21908CFF"), 
                     type= "est", 
                     show.intercept = F, 
                     show.values = T, 
                     value.offset = .3,
                     show.p = T, 
                     vline.color = "grey", 
                     order.terms = c(1,2,3,6,7,8,9, 4,5),
                     title = "", 
                     axis.labels=c( "Mean annual temperature", "Current body mass", "Africa * Canopy Height", "Precipitation Seasonality", "Annual Precipitation", "Temperature Seasonality", 
                                   "Body mass decrease", "Canopy Height", 
                                   "Distribution in Africa")) + 
  theme_classic() + 
  theme(text=element_text(size=13, color = "black"),
        axis.text=element_text(size=13), 
        axis.title=element_text(size=13))


pLMsim <- pLMsim +  theme(axis.text=element_text(size=11.5, color = "black"), #change font size of axis text
                          axis.title=element_text(size=13), #change font size of axis titles
                          plot.title=element_text(size=20),
                          text = element_text(size = 13), #change font size of plot title
                          legend.text=element_text(size=20), #change font size of legend text
                          legend.title=element_text(size=20),
                          plot.margin = margin(t = 12,  # Top margin
                                               r = 12,  # Right margin
                                               b = 12,  # Bottom margin
                                               l = 12)) # Left margin) #change font size of legend title 
pLMsim


### =================================================================================================== ###

## Part III: Worl Map =============
cleanCuts <- function(x){
  # Cleans up vector of factor levels generated using `cut` function for nicer plotting
  x <- gsub(",", " - ",x)
  x <- gsub("\\[|\\(|\\]|\\)", "", x)
  num <- str_split_fixed(x, pattern = " - ", n=2)
  new <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(abs(as.numeric(num[i,1])) < 0.001){
      #start <- changeSciNot(scales::scientific(as.numeric(num[i,1]), 0))
      start <- 0 #scales::scientific(as.numeric(num[i,1]), 1)
    } else {
      start <- round(as.numeric(num[i,1]), 3 )
    }  
    if(abs(as.numeric(num[i,2])) < 0.001){
      #end <- changeSciNot(scales::scientific(as.numeric(num[i,2]), 0))
      end <- 0 #scales::scientific(as.numeric(num[i,2]), 1)
    } else {
      end <- round(as.numeric(num[i,2]), 3 )
    }
    new[i] <- paste(start,end, sep = " - ")# em-dash
    #new[i] <- paste(start,end, sep = " \u2013 ")# em-dash
  }
  return(new)
}
### World map shape file ============

library(rgeos); library(stringr); library(rgdal)

shp <- readOGR(dsn = "FruitsAfrica/Data/shp", layer = "TDWG_level3_Coordinates")
shp.df <- gSimplify(shp, tol=0.05, topologyPreserve = TRUE)
shp.df <- fortify(shp.df)
county_ids <- as.data.frame(cbind(id = rownames(shp@data), LEVEL_3_CO = shp@data$LEVEL_3_CO, region = shp@data$LEVEL_NAME))
shp.df2 <- inner_join(shp.df, county_ids, by="id")

map_theme <- theme(panel.background = element_blank(),
                   panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.text = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   legend.position = "bottom",
                   legend.justification = "center",
                   legend.title.align = 0,
                   strip.background = element_blank(),
                   legend.key = element_blank(),
                   strip.text = element_text(size = 18))


### 

dd <- tdwg_final_0

dd <- dd[complete.cases(dd$max95FL_palms),]
dd$accAfrica <- factor(dd$accAfrica)
dd$accRealm <- as.factor(dd$accRealm)

dd$sqrtFL <- sqrt(dd$max95FL_palms)
dd$sqrt_BBM_mean <- sqrt(exp(dd$log_max95FL_palms_BBMmean))
dd$sqrt_change <- sqrt(dd$mean_abs_BSchange)

### a) NEW Simulated fruit size =========
my_colors <- setNames(c("blue", "#bdbdbd"), c("1", "0"))
dd$max95_BBM <- exp(dd$log_max95FL_palms_BBM_MCC) # new

### For scale = empirical:
# Fruit size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FL_palms")
tdwg_maxFS <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FL_palms", na.rm=T)
q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_maxFS$value, probs = q_probs, na.rm = T, type = 8))

## Mapping ==============
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95_BBM")
tdwg_maxBBM <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95_BBM", na.rm=T)

max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxBBM$maxFS_q <- cut(tdwg_maxBBM$value, q_values, include.lowest = T)
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFSsim_p_new <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white", data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT, fill = maxFS_q, size = maxFS_q), data =tdwg_maxBBM, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]", override.aes = list(size = c(1.285714,1.571429)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = c(1.285714, 1.571429))+
  scale_color_manual(values = c("#46327e","#365c8d"),
                     aesthetics = c("colour", "fill"),
                     #breaks = levels(tdwg_maxBBM$maxFS_q),
                     #labels = levels(tdwg_maxBBM$maxFS_q)
  ) +
  theme(legend.position = "right",
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) 

maxFSsim_p_new

## quick t-test ====
library(rstatix)
stat.test <- dd %>% wilcox_test(max95_BBM ~ accAfrica)
stat.test 

#   .y.       group1 group2    n1    n2   statistic        p
#  max95_BBM  0      1         88    42     2910.     0.000000128

dd %>%  wilcox_effsize(max95_BBM ~ accAfrica) # moderate, effsize = 0.463

library(EnvStats)
sim_violin <- ggpubr::ggviolin(dd, x = "accAfrica", y = "max95_BBM", fill="accAfrica", color ="#636363", show.legend=F, trim = T) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(labels = c(
    "1" = "Africa",
    "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y=expression(sqrt("max 95-percentile fruit length [cm]")), x=NULL)+
  stat_n_text(y.pos = 1.35) + 
  theme(legend.position="none")+
  theme_classic()+
  theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=11.5), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 20,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title   
sim_violin <- sim_violin + labs(subtitle = get_test_label(stat.test, detailed = TRUE, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.25))
#sim_violin <- sim_violin + labs(subtitle = expression(paste("Africa:"~beta~"= -0.04, se = 0.11, p = n.s."))) + theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015)) +coord_fixed(ratio=3.25)   #  Add p-value



## b) OLD Simulated Fruit size ===========

dd$max95_BBM <- exp(dd$log_max95FL_palms_BBMmean) #old


## Mapping ==============
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95_BBM")
tdwg_maxBBM <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95_BBM", na.rm=T)

max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxBBM$maxFS_q <- cut(tdwg_maxBBM$value, q_values, include.lowest = T)
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFSsim_p_old <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white", data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT, fill = maxFS_q, size = maxFS_q), data =tdwg_maxBBM, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]", override.aes = list(size = c(1.285714,1.571429)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = c(1.285714, 1.571429))+
  scale_color_manual(values = c("#46327e","#365c8d"),
                     aesthetics = c("colour", "fill"),
                     #breaks = levels(tdwg_maxBBM$maxFS_q),
                     #labels = levels(tdwg_maxBBM$maxFS_q)
  ) +
  theme(legend.position = "right",
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) 

maxFSsim_p_old

## quick t-test ====
stat.test <- dd %>% wilcox_test(max95_BBM ~ accAfrica)
stat.test # n.s.
#   .y.       group1 group2      n1    n2   statistic        p
# max95_BBM     0      1         88    42      2096       0.218


dd %>%  wilcox_effsize(max95_BBM ~ accAfrica) # small, effsize = 0.108



library(EnvStats)
sim_violin <- ggpubr::ggviolin(dd, x = "accAfrica", y = "max95_BBM", fill="accAfrica", color ="#636363", show.legend=F, trim = T) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(labels = c(
    "1" = "Africa",
    "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y=expression(sqrt("max 95-percentile fruit length [cm]")), x=NULL)+
  stat_n_text(y.pos = 1.35) + 
  theme(legend.position="none")+
  theme_classic()+
  theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=11.5), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 20,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title   
sim_violin <- sim_violin + labs(subtitle = get_test_label(stat.test, detailed = TRUE, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.25))
#sim_violin <- sim_violin + labs(subtitle = expression(paste("Africa:"~beta~"= -0.04, se = 0.11, p = n.s."))) + theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015)) +coord_fixed(ratio=3.25)   #  Add p-value
sim_violin
