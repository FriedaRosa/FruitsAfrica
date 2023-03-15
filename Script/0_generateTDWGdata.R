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

palm_traits_v3 <- merge(palms_traits_v2, newBBM_MCC_mean, by =c("SpecName"))

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
write.csv(palm_traits_v3, file=paste0(res.dir, "/palms_final.csv"))

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

                            log_meanFL_palms_BBM_MCC = mean(log_BBM_MCC_avg, na.rm = T),
                            log_medianFL_palms_BBM_MCC = median(log_BBM_MCC_avg, na.rm = T),
                            log_max95FL_palms_BBM_MCC = quantile(log_BBM_MCC_avg, 0.95))

tdwg_palms_summary <- strings2factors(tdwg_palms_summary)
str(tdwg_palms_summary)
tdwg_palms_summary <- unique(tdwg_palms_summary)
saveRDS(tdwg_palms_summary, file.path(res.dir, "tdwg_palms_summary.rds"))

# reduced to tdwg3 with >3 palms ==================
palms3 <- tdwg_palms_summary%>% filter(sp_richness_palms > 3)
palms3tdwg3 <- unique(palms3$Area_code_L3) #131 botanic countries



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
write.csv(subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect), file.path(res.dir, "mammal_presnat_occ_trait.csv"), row.names = F)

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
write.csv(subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_plant_intersect), file.path(res.dir, "mammal_curr_occ_trait.csv"), row.names = F)


# Megafauna in Africa?
curr_mammals_africa <- mammal_curr_occ_trait %>% subset(LEVEL_3_CO %in% mammal_plant_intersect) %>% filter(CONTINENT=="AFRICA" & Mass.g > 50000)
mammal_sp <- unique(curr_mammals_africa$SpecName)

# Calculate differences in body sizes across tdwg3 ====================
tdwg_mammals_summary <- merge(tdwg_curr_summary, tdwg_presnat_summary)
tdwg_mammals_summary2<-unique(tdwg_mammals_summary)
rownames(tdwg_mammals_summary2) <- NULL
str(tdwg_mammals_summary2)


## Reduce to >3 palm species:
tdwg_curr_summary_3 <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% palms3tdwg3),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = T),
        curr_medianBodySize = median(Mass.g, na.rm = T),
        curr_max95BodySize = exp(quantile(log(Mass.g), probs = 0.95, na.rm = T)),
        curr_sdBodySize = exp(sd(log(Mass.g), na.rm = T)),
        region = unique(LEVEL_NAME),
        sp_richness_curr = length(unique(SpecName)),
        sp_richness_curr_meso = length(unique(SpecName[Mass.g > 10000])),
        sp_richness_curr_mega = length(unique(SpecName[Mass.g > 44000]))
  )

tdwg_presnat_summary_3 <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% palms3tdwg3),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = T),
        presNat_medianBodySize = median(Mass.g[Default == "Y"], na.rm = T),
        presNat_max95BodySize = exp(quantile(log(Mass.g[Default == "Y"]), probs = 0.95, na.rm = T)),
        presNat_sdBodySize = exp(sd(log(Mass.g), na.rm = T)),
        area = unique(LEVEL_NAME),
        sp_richness_PresNat = length(unique(SpecName)),
        sp_richness_presNat_meso = length(unique(SpecName[Mass.g > 10000])),
        sp_richness_presNat_mega = length(unique(SpecName[Mass.g > 44000])))

tdwg_mammals_summary_3 <- merge(tdwg_curr_summary_3, tdwg_presnat_summary_3)
tdwg_mammals_summary2_3<-unique(tdwg_mammals_summary_3)
rownames(tdwg_mammals_summary2_3) <- NULL

tdwg_mammals_summary_3palms <- tdwg_mammals_summary2_3 %>% mutate(mean_perc_BSchange = round((1 - curr_meanBodySize / presNat_meanBodySize) * 100, 2),
                                                                  median_perc_BSchange = round((1 - curr_medianBodySize / presNat_medianBodySize) * 100, 2),
                                                                  max95_perc_BSchange = round((1 - curr_max95BodySize / presNat_max95BodySize) * 100, 2),
                                                                  mean_abs_BSchange = round(presNat_meanBodySize - curr_meanBodySize,2),
                                                                  median_abs_BSchange = round(presNat_medianBodySize - curr_medianBodySize,2),
                                                                  max95_abs_BSchange = round(presNat_max95BodySize - curr_max95BodySize,2))
rownames(tdwg_mammals_summary_3palms) <- tdwg_mammals_summary_3palms$LEVEL_3_CO
rownames(palms3) <- palms3$Area_code_L3

## Angiosperms =================================

angios <- read.csv(file.path(data.dir, "average_fruit_width_Frieda_v2.csv"))
sp_angios <- angios$Genus_Species_Accepted_LCVP
length(unique(sp_angios))

# 1.2 Distribution data on botanical country (TDWG3) scale ================

dd_merge <- readRDS(file.path(data.dir, "kew_merged_Angios.rds"))

# Names data =====
# names <- read.table(file.path(data.dir, "wcvp_names.txt"), sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# names <- names %>% filter(taxon_rank == "Species" & family != "Arecaceae")
# names <- names %>% dplyr::select(plant_name_id, taxon_status, family, genus, species, taxon_name, accepted_plant_name_id)
# names$SpecName <- paste(names$genus, names$species, sep="_")
# #
# names_angio <- names %>% filter(SpecName %in% sp_angios)
# 
# # Distribution data ======
# dist <- read.table(file.path(data.dir, "wcvp_distribution.txt"), sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# dist <- dist %>% dplyr::select(plant_name_id, continent, region, area_code_l3, area)
# 
# # Merge ====
# dd_merge <- merge(names_angio, dist, by="plant_name_id", all.x=TRUE)
# dd_merge <- unique(dd_merge)
# 
# 
# # Counts =====
# length(unique(dd_merge$family)) # 168 families
# length(unique(dd_merge$genus))  # 1016 genera
# length(unique(dd_merge$area_code_l3)) #363 tdwg3
# length(unique(dd_merge$SpecName))# 4400 species
# length(unique(dd_merge$accepted_plant_name_id)) #4596
# length(unique(dd_merge$plant_name_id)) #4677
# 
# # Save to .RDS ====
# saveRDS(dd_merge, file.path(data.dir, "kew_merged_Angios.rds"))

# =================================================================================================================================

dd_merge <- dd_merge %>% select("SpecName", "genus", "family", "area_code_l3", "area", "region", "continent")

regions <- dd_merge %>% 
  dplyr::select(area_code_l3, area, region, continent) %>% 
  mutate_all(na_if,"") %>% 
  filter_at(vars(area_code_l3), any_vars(!is.na(.))) %>% 
  distinct(.)


# reduced to those with distribution data (!= NA)
dd_merge <- dd_merge %>%
  filter_at(vars(area_code_l3), all_vars(!is.na(.)))

str(dd_merge)

angios_FW_kew <- merge(angios, dd_merge, by.x="Genus_Species_Accepted_LCVP", by.y="SpecName", all.x=T)

angios_env <- unique(merge(angios_FW_kew, tdwg_env2, by.x=c("area_code_l3", "continent"), by.y=c("LEVEL_3_CO", "CONTINENT")))
angios_env <- unique(angios_env)
angios_env$FW_cm <- round(angios_env$Fruit_width_mm_combi /10, 3)

## 3. aggregate maximum + median angiosperm and (empirical, simulated) palm traits on tdwg3; merge both datasets together
tdwg_angios_summary <- ddply(.data = subset(angios_env, area_code_l3 %in% mammal_plant_intersect),
                             .variables = .(area_code_l3),
                             .fun = dplyr::summarise,
                             area = area,
                             
                             meanFW_angios = mean(FW_cm, na.rm = T),
                             medianFW_angios = median(FW_cm, na.rm = T),
                             max95FW_angios = exp(quantile(log(FW_cm), probs = 0.95, na.rm = T)))

tdwg_angios_summary <- strings2factors(tdwg_angios_summary)
str(tdwg_angios_summary)
tdwg_angios_summary <- unique(tdwg_angios_summary)



## Combine all data to one dataframe ## ======================================
tdwg_final <- merge(palms3, tdwg_angios_summary, by.x=c("Area_code_L3", "area"), by.y=c("area_code_l3", "area"), all.x = T)
rownames(tdwg_final) <- tdwg_final$Area_code_L3
tdwg_final2a <- merge(tdwg_mammals_summary_3palms, tdwg_final, by = 0, all = T)
tdwg_final2b <- merge(tdwg_final2a, tdwg_env2, by = "LEVEL_3_CO", all.x = TRUE)
tdwg_final2 <- tdwg_final2b
# rename bioclim variables:
colnames(tdwg_final2)[colnames(tdwg_final2)      # Rename two variable names
                      %in% c("bio1_mean", "bio4_mean", "bio12_mean", "bio15_mean")] <- c("MAT", "TempSeas", "AnnualPrec", "PrecSeas")

# reduce to needed columns:

tdwg_final3 <- tdwg_final2 %>% 
  filter_at(vars(LEVEL_3_CO, 
                 max95FL_palms,
                 log_max95FL_palms_BBM_MCC, 
                 curr_max95BodySize, 
                 mean_abs_BSchange), 
            any_vars(!is.na(.))) %>%
  select(LEVEL_3_CO, LAT, LONG, accAfrica, max95FL_palms, max95FW_angios,
         log_max95FL_palms_BBM_MCC, 
         curr_max95BodySize, mean_abs_BSchange, presNat_max95BodySize, mean_perc_BSchange, 
         "MAT", "TempSeas", "AnnualPrec", "PrecSeas", "CH_Mean", "CH_PercCover", 
         accRealm ) %>% 
  distinct(.)




write.csv(tdwg_final3, file.path(res.dir, "tdwg_final.csv"))
save.image(file="generateTDWGdata_ws.RData")
