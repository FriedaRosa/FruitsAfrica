## Directories ========================
rm(list = ls())
setwd("~/share/groups/ea/Frieda/Spatial_analyses_May22/Spatial_analyses_RServer/")
main.dir <- "ProjectAfricaMegafaunaOnly/Sensitivity_Angiosperms"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results/")
src.dir <- file.path(main.dir, "src/")
fig.dir <- file.path(main.dir, "figures/")
options(stringsAsFactors =T)

## Packages ========================
library(plyr); library(dplyr); library(reshape2)
library(sp); library(rgdal)
library(car); library(MuMIn)
library(readxl)
library(picante)

## Data import ========================

fdb <- read.csv(file.path(data.dir, "fdb_trait_final.csv"), stringsAsFactors = T)
str(fdb)
fdb <- fdb %>% filter_at(vars(Fruit_width_mm_combi, area_code_l3), any_vars(!is.na(.))) %>% filter(area_code_l3 != "")

length(unique(fdb$Genus_Species_Accepted)) # 4663 species
length(unique(fdb$genus)) #1074
length(unique(fdb$family)) # 173 accepted families
length(unique(fdb$Family)) # 201 in total 

fdb %>% filter(accAfrica == "1") %>% dplyr::summarize(n_Sp = length(unique(Genus_Species_Accepted)),
                                                      n_G = length(unique(genus)),
                                                      n_F = length(unique(family))) 

# n_Sp n_G n_F
# 1  837 425 112

fdb %>% filter(accAfrica == "0") %>% dplyr::summarize(n_Sp = length(unique(Genus_Species_Accepted)),
                                                      n_G = length(unique(genus)),
                                                      n_F = length(unique(family)))
# n_Sp n_G n_F
# 1 4232 961 171


fdb$accAfrica[fdb$continent == "AFRICA"] <- ("1")
fdb$accAfrica[fdb$continent != "AFRICA"] <- ("0")
fdb$accAfrica[fdb$area_code_l3 == "MDG"] <- ("0")
fdb$Africa[fdb$accRealm == "Africa"] <- ("1")
fdb$Africa[fdb$accRealm != "Africa"] <- ("0")

fdb$accAfrica <- as.factor(fdb$accAfrica)
fdb$accRealm <- as.factor(fdb$accRealm)
fdb$area_code_l3 <- as.factor(fdb$area_code_l3)
fdb$Africa <- as.factor(fdb$Africa)

str(fdb)
# 3. Environmental data ========================

tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))
env.vars <- c("LEVEL_3_CO", "LAT", "LONG",  "CONTINENT", "CH_Mean", "CH_PercCover", "bio1_mean", "bio4_mean", "bio15_mean", "bio12_mean")
tdwg_env2<- tdwg_env[,c(env.vars)]
rownames(tdwg_env2) <- tdwg_env2$LEVEL_3_CO


#fdb_env <- unique(merge(fdb, tdwg_env2, by.x=c("area_code_l3", "continent"), by.y=c("LEVEL_3_CO", "CONTINENT")), all.x=T)

# 4. Mammal / Frugivore data ========================	

phylacine_trait <- read.csv(file.path(data.dir, "Phylacine_Trait_data.csv")) ## Phylacine trait = for body masses of mammals
mammal_curr_comb_occ <- read.csv(file.path(data.dir, "mammal_curr_occ.csv")) ## mammal curr = list of botanic countries and extant mammals (from IUCN range maps)
mammal_presnat_comb_occ <- read.csv(file.path(data.dir, "mammal_presnat_occ.csv")) ## mammal presnat = list of botanic countries and extant+extinct mammals
frugivoreClass <- read.csv(file.path(data.dir, "frugivoreClassification.csv")) ## frugivoreClass = mammal list with indication of proportion frugivory (family, order)

## step 1. make mammal country lists (1. current, 2. extinct)
mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_comb_occ$LEVEL_3_CO),   # 191 botanic countries
                                         unique(mammal_curr_comb_occ$LEVEL_3_CO)))     # 190 botanic countries

palms <- read.csv("ProjectAfricaMegafaunaOnly/results/tdwg_final.csv") # to match countries lists to palms

## 2. make list of countries where mammals and species overlap
mammal_plant_intersect <- intersect(unique(palms$LEVEL_3_CO), mammal_countrylist)  #191 botanic countries
sort(mammal_plant_intersect, decreasing = FALSE)


## 3. Aggregate fdb data:
tdwg_angios_summary <- ddply(.data = subset(fdb, area_code_l3 %in% mammal_plant_intersect),
                             .variables = .(area_code_l3),
                             .fun = dplyr::summarise,
                             sp_richness_angios = length(Genus_Species_Accepted),
                             accRealm = accRealm,
                             accAfrica = accAfrica,
                             area = area,
                             max95FW_angios = exp(quantile(log(Fruit_width_mm_combi), probs = 0.95, na.rm = T)))
tdwg_angios_summary <- unique(tdwg_angios_summary)
str(tdwg_angios_summary)

## 4. aggregate mammal bodysizes (extant, extinct) on tdwg3; merge both together =================

# Calculate maximum and median body sizes of present natural mammal assemblages
mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
mammal_presnat_occ_trait <- merge(mammal_presnat_occ_trait, frugivoreClass, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)

# Extant species are coded as NAs but they are included under all definitions
mammal_presnat_occ_trait$Liberal[is.na(mammal_presnat_occ_trait$Liberal)] <- "Y"
mammal_presnat_occ_trait$Default[is.na(mammal_presnat_occ_trait$Default)] <- "Y"
mammal_presnat_occ_trait$Conservative[is.na(mammal_presnat_occ_trait$Conservative)] <- "Y"

tdwg_presnat_summary <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = T),
        presNat_medianBodySize = median(Mass.g[Default == "Y"], na.rm = T),
        presNat_max95BodySize = exp(quantile(log(Mass.g[Default == "Y"]), probs = 0.95, na.rm = T)),
        presNat_sdBodySize = exp(sd(log(Mass.g), na.rm = T)),
        area = LEVEL_NAME,
        
        sp_richness_PresNat = length(unique(SpecName)),
        sp_richness_presNat_meso = length(unique(SpecName[Mass.g > 10000])),
        sp_richness_presNat_mega = length(unique(SpecName[Mass.g > 44000])))


tdwg_presnat_summary <- unique(tdwg_presnat_summary)
#write.csv(subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect), file.path(res.dir, "mammal_presnat_occ_trait.csv"), row.names = F)

mammal_presnat_occ_trait_ls <- subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect)

# Current the maximum and median body sizes of current mammal assemblages
mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)

tdwg_curr_summary <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = T),
        curr_medianBodySize = median(Mass.g, na.rm = T),
        curr_max95BodySize = exp(quantile(log(Mass.g), probs = 0.95, na.rm = T)),
        curr_sdBodySize = exp(sd(log(Mass.g), na.rm = T)),
        area = LEVEL_NAME,
        sp_richness_curr = length(unique(SpecName)),
        sp_richness_curr_meso = length(unique(SpecName[Mass.g > 10000])),
        sp_richness_curr_mega = length(unique(SpecName[Mass.g > 44000]))
  )
tdwg_curr_summary <- unique(tdwg_curr_summary)
# NOTE: Pteropus niger (Mascarene fruit bat) is found on both Mauritius and Reunion but is only on Mauritius for the current dataset, as a result, there are no mammals on Reunion in the current case (Pteropus) but 2 in the present-natural dataset

mammal_curr_occ_trait <- subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect)

# Calculate differences in body sizes across tdwg3 ====================
tdwg_mammals_summary <- merge(tdwg_curr_summary, tdwg_presnat_summary)
tdwg_mammals_summary2<-unique(tdwg_mammals_summary)
rownames(tdwg_mammals_summary2) <- NULL
str(tdwg_mammals_summary2)

tdwg_mammals_summary_3palms <- tdwg_mammals_summary2 %>% mutate(mean_perc_BSchange = round((1 - curr_meanBodySize / presNat_meanBodySize) * 100, 2),
                                                                  median_perc_BSchange = round((1 - curr_medianBodySize / presNat_medianBodySize) * 100, 2),
                                                                  max95_perc_BSchange = round((1 - curr_max95BodySize / presNat_max95BodySize) * 100, 2),
                                                                  mean_abs_BSchange = presNat_meanBodySize - curr_meanBodySize,
                                                                  median_abs_BSchange = presNat_medianBodySize - curr_medianBodySize,
                                                                  max95_abs_BSchange = presNat_max95BodySize - curr_max95BodySize)
rownames(tdwg_mammals_summary_3palms) <- tdwg_mammals_summary_3palms$LEVEL_3_CO
rownames(tdwg_angios_summary) <- tdwg_angios_summary$area_code_l3

# Combine all data to one dataframe
tdwg_final <- merge(tdwg_mammals_summary_3palms, tdwg_angios_summary, by = 0, all = T)
tdwg_final2 <- merge(tdwg_final, tdwg_env2, by = "LEVEL_3_CO", all.x = TRUE)

# rename bioclim variables:
colnames(tdwg_final2)[colnames(tdwg_final2)      # Rename two variable names
                      %in% c("bio1_mean", "bio4_mean", "bio12_mean", "bio15_mean")] <- c("MAT", "TempSeas", "AnnualPrec", "PrecSeas")

# reduce to needed columns:

tdwg_final3 <- tdwg_final2 %>% filter_at(vars(LEVEL_3_CO, max95FW_angios, curr_max95BodySize, mean_abs_BSchange), any_vars(!is.na(.))) %>%
  dplyr::select(LEVEL_3_CO, LAT, LONG, accAfrica, max95FW_angios, curr_max95BodySize, mean_abs_BSchange, mean_perc_BSchange, "MAT", "TempSeas", "AnnualPrec", "PrecSeas", "CH_Mean", accRealm ) %>% distinct(.)

write.csv(tdwg_final3, file=paste0(res.dir, "tdwg_final_angios.csv"))
save.image(file="ProjectAfricaMegafaunaOnly/Sensitivity_Angiosperms/src/Sensitivity_Angiosperms_generateTDWGdata_ws.RData")
