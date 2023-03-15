rm(list=ls())

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
tree <- read.nexus("TREES.nex")

# read data ===================================
dd <- read.csv("Step1_duplicated_regimes_checked.csv")
str(dd)

# selecting columns and removing duplicated species rows
dd2 <- dd[,c("SpecName","accGenus","PalmTribe","PalmSubfamily","AverageFruitLength_cm","Americas",
             "Asia","Africa","Madagascar","Pacific","accRealm","log_BBM_MCC","log_BBM_mean")]
dd <- dd2 %>% distinct()

# correct data format ===================================
dd$SpecName <- as.character(dd$SpecName)
dd$AverageFruitLength_cm1 <- as.numeric(as.character(dd$AverageFruitLength_cm))
dd2 <- dd[!is.na(dd["AverageFruitLength_cm1"]),]
dd <- as.data.frame(dd2)
dd$AverageFruitLength_cm <- log(dd$AverageFruitLength_cm1) #log-transformation
dd$Americas <- as.factor(dd$Americas)

# match data and tree ===================================
row.names(dd) <- dd$SpecName
setdiff(tree$tip.label, dd$SpecName)

# Drop these species from the tree, because we need matching data between tree and traits
tree2 <- list()
for (i in 1:length(tree)){
  tree2[[i]] <- drop.tip(tree[[i]], setdiff(tree[[i]]$tip.label, dd$SpecName))
  tree2[[i]] <- ladderize(tree2[[i]])
  }
for (i in 1:length(tree)){
  matches <- match(dd$SpecName, tree2[[i]]$tip.label, nomatch = 0)
  data2 <- subset(dd, matches != 0)
}

# Assign species names as rownames for identification
row.names(dd) <- dd$SpecName
row.names(data2) <- data2$SpecName
tree <- tree2

for (i in 1:length(tree)){
  data2 <- ReorderData(tree[[1]], data2, taxa.names = "row names")
}
head(data2)
head(tree[[1]]$tip.label)

# Make vectors with regimes ===================================
## Americas / non-Americas
americas <- numeric(length = nrow(data2))
names(americas) <- rownames(data2)
americas[data2[,"Americas"]=="0"] <- "1_world"
americas[data2[,"Americas"]=="1"] <- "2_americas"

## log-transformed trait
log_avgl <- data2$AverageFruitLength_cm
names(log_avgl) <- row.names(data2)

## Sort data like phylogeny to prevent errors
for (i in 1:length(tree)){
  americas <- ReorderData(tree[[i]], americas, taxa.names="names")
}
for (i in 1:length(tree)){
  log_avgl <- ReorderData(tree[[i]], log_avgl, taxa.names="names")
}

## Mk models for transition rates - compare fit
aic.discrete <- list()
weights <- list()
for (i in 1:length(tree)){
  equal <- fitDiscrete(tree[[i]], americas, model = "ER")
  sym <- fitDiscrete(tree[[i]], americas, model = "SYM")
  ard <- fitDiscrete(tree[[i]], americas, model = "ARD")
  aic.discrete[[i]] <- setNames(c(equal$opt$aic, sym$opt$aic, ard$opt$aic), c("equal", "sym", "different"))
  weights[[i]] <- aicw(aic.discrete[[i]])
}

#Finding which model has the lowest AIC -> comparing results from the 100 trees
result_mk_fit <- vector()
for (e in 1:length(weights)){
  result_mk_fit[e] <- which(weights[[e]]$fit == min(weights[[e]]$fit))
}#All of them indicated "different" (row 3) as the best model

save.image("mk_models_americas_22Sep2022.rda")




