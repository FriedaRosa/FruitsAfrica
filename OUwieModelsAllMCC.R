## Simulate neutral fruit sizes: different ways. Compare here: ##

rm(list=ls())
### Africa x Average Length ###

################################################################################

# full R script for Masterthesis by Frieda #

# Tasks in this script:
## 1. read in a tree and trait data for Palms
## 2. match tree and data 
## 3. make trait vector for regimes and fruit length or simulate BM evolution across the tips of the tree
## 4. reconstruction of ancestral regime states (Africa / moistdry)
## 5. make data frames for OUwie trait models (Species, regime, trait)
## 6. run OUwie trait models (1. BM, 2. BMS, 3. OU1, 4. OUM, 5. OUMV, 6. OUMA, 7. OUMVA)
## 7. shape OUwie output for analysis (calculate AICc weights, model selection, write to file)

### note : please comment out code blocks for different analysis ###

################################################################################


# libraries ===================================
suppressPackageStartupMessages(library(OUwie))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(evobiR))
suppressPackageStartupMessages(library(BBMV))

# define arguments ===================================
args <- commandArgs(trailingOnly=TRUE)
nex_file <- args[1]
output_dir <- args[2]
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# read tree ===================================
tree <- read.nexus(nex_file)

# read data ===================================
dd <- read.csv2("/home/woelke/PalmDataPhyloRegimes.csv")
str(dd)

# correct data format ===================================
dd$SpecName <- as.character(dd$SpecName)
dd$AverageFruitLength_cm1 <- as.numeric(as.character(dd$AverageFruitLength_cm))
dd2 <- dd[!is.na(dd["AverageFruitLength_cm1"]),]
dd <- as.data.frame(dd2)
dd$AverageFruitLength_cm <- log(dd$AverageFruitLength_cm1)
dd$accAfrica <- as.factor(dd$accAfrica)
dd$moist0_dry1 <- as.factor(dd$moist0_dry1)

# match data and tree ===================================
row.names(dd) <- dd$SpecName
setdiff(tree$tip.label, dd$SpecName)

# Drop these species from the tree, because we need matching data between tree and traits
tree2 <- drop.tip(tree, setdiff(tree$tip.label, dd$SpecName))
tree2 <- ladderize(tree2)

matches <- match(dd$SpecName, tree2$tip.label, nomatch = 0)
data2 <- subset(dd, matches != 0)

# Assign species names as rownames for identification
row.names(dd) <- dd$SpecName
row.names(data2) <- data2$SpecName
tree <- tree2

data2 <- ReorderData(tree, data2, taxa.names = "SpecName")
head(data2)
head(tree$tip.label)


# Make vectors with regimes ===================================
## Africa / non-Africa
africa <- numeric(length = nrow(data2))
names(africa) <- rownames(data2)
africa[data2[,"accAfrica"]=="0"] <- "1_world"
africa[data2[,"accAfrica"]=="1"] <- "2_africa"

## moist / dry
moistdry <- numeric(length = nrow(data2))
names(moistdry) <- rownames(data2)
moistdry[data2[,"moist0_dry1"]=="0"] <- "1_moist"
moistdry[data2[,"moist0_dry1"]=="1"] <- "2_dry"

# trait vectors ===================================

## log-transformed
log_avgl <- data2$AverageFruitLength_cm
names(log_avgl) <- row.names(data2)

## Sort data like phylogeny to prevent errors
moistdry <- ReorderData(tree, moistdry, taxa.names="names")
africa <- ReorderData(tree, africa, taxa.names="names")
log_avgl <- ReorderData(tree, log_avgl, taxa.names="names")



BBM_sim <- readRDS("/home/woelke/BBM_sim.rds")
BBM_sim <- ReorderData(tree, BBM_sim, taxa.names = "names")



## Africa
simmap_afr <-readRDS("/home/woelke/simmap_afr.rds")
#saveRDS(simmap_afr, file.path(output_dir, paste0("simmap_afr", task_id, ".rds")))

# summarize 100 simulations from simmap
x <- getStates(simmap_afr,"nodes")
sumr <- apply(x, 1, function(x){
  mfv <- table(x)
  return(names(mfv[mfv==max(mfv)]))
})
simmap_afr[[1]]$node.label <- sumr
simmap_AFR <- simmap_afr[[1]]




## moist/dry
simmap_md <-readRDS("/home/woelke/simmap_md.rds")
#saveRDS(simmap_md, file.path(output_dir, paste0("simmap_md", task_id, ".rds")))

# summarize 100 simulations from simmap
x <- getStates(simmap_md,"nodes")
sumr <- apply(x, 1, function(x){
  mfv <- table(x)
  return(names(mfv[mfv==max(mfv)]))
})
simmap_md[[1]]$node.label <- sumr
simmap_MD <- simmap_md[[1]]


# make OUwie frames ===================================

## moistdry x avgl
data_md_avgl <- data.frame(species=names(log_avgl), regime = moistdry, trait = as.data.frame(log_avgl))
data_md_avgl$species <- as.character(data_md_avgl$species)
data_md_avgl$regime <- as.factor(data_md_avgl$regime)
data_md_avgl$log_avgl <- as.numeric(as.character(data_md_avgl$log_avgl))
print(str(data_md_avgl))

## Africa x avgl 
data_afr_avgl <- data.frame(species=names(log_avgl), regime = africa, trait = as.data.frame(log_avgl))
data_afr_avgl$species <- as.character(data_afr_avgl$species)
data_afr_avgl$regime <- as.factor(data_afr_avgl$regime)
data_afr_avgl$log_avgl <- as.numeric(as.character(data_afr_avgl$log_avgl))
print(str(data_afr_avgl))

## Africa x BM
data_Afr_BM <- data.frame(species=names(BBM_sim), regime = africa, trait = as.data.frame(BBM_sim))
data_Afr_BM$species <- as.character(data_Afr_BM$species)
data_Afr_BM$regime <- as.factor(data_Afr_BM$regime)
data_Afr_BM$BBM_sim <- as.numeric(as.character(data_Afr_BM$BBM_sim))
str(data_Afr_BM)

## moistdry x BM
data_md_BM <- data.frame(species=names(BBM_sim), regime = moistdry, trait = as.data.frame(BBM_sim))
data_md_BM$species <- as.character(data_md_BM$species)
data_md_BM$regime <- as.factor(data_md_BM$regime)
data_md_BM$BBM_sim <- as.numeric(as.character(data_md_BM$BBM_sim))
str(data_md_BM)

# run OUwie models ===================================
## Africa x empirical
BM_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_Afr_avgl <- OUwie(simmap_AFR, data_afr_avgl, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
## save output files
ouwie_Afr_avgl <- list(BM_Afr_avgl, BMS_Afr_avgl, OU1_Afr_avgl, OUM_Afr_avgl, OUMV_Afr_avgl, OUMA_Afr_avgl, OUMVA_Afr_avgl)
saveRDS(ouwie_Afr_avgl, file.path(output_dir, paste0("ouwie_Afr_avgl_", task_id, ".rds")))

## Moist/Dry x empirical 
BM_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_md_avgl <- OUwie(simmap_MD, data_md_avgl, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
## save output files
ouwie_md_avgl <- list(BM_md_avgl, BMS_md_avgl, OU1_md_avgl, OUM_md_avgl, OUMV_md_avgl, OUMA_md_avgl, OUMVA_md_avgl)
saveRDS(ouwie_md_avgl, file.path(output_dir, paste0("ouwie_md_avgl_", task_id, ".rds")))



## Africa x BBM:
BM_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_Afr_BM <- OUwie(simmap_AFR, data_Afr_BM, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
## save output files
ouwie_Afr_BM <- list(BM_Afr_BM, BMS_Afr_BM, OU1_Afr_BM, OUM_Afr_BM, OUMV_Afr_BM, OUMA_Afr_BM, OUMVA_Afr_BM)
saveRDS(ouwie_Afr_BM, file.path(output_dir, paste0("ouwie_Afr_BM_", task_id, ".rds")))

## moistdry x BM
BM_md_BM <- OUwie(simmap_MD, data_md_BM, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_md_BM <- OUwie(simmap_MD, data_md_BM, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_md_BM <- OUwie(simmap_MD, data_md_BM, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_md_BM <- OUwie(simmap_MD, data_md_BM, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_md_BM <- OUwie(simmap_MD, data_md_BM, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_md_BM <- OUwie(simmap_MD, data_md_BM, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_md_BM <- OUwie(simmap_MD, data_md_BM, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
## save output files
ouwie_md_BM <- list(BM_md_BM, BMS_md_BM, OU1_md_BM, OUM_md_BM, OUMV_md_BM, OUMA_md_BM, OUMVA_md_BM)
saveRDS(ouwie_md_BM, file.path(output_dir, paste0("ouwie_md_BM_", task_id, ".rds")))