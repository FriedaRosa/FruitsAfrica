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
tree2 <- drop.tip(tree, setdiff(tree$tip.label, dd$SpecName))
tree2 <- ladderize(tree2)

matches <- match(dd$SpecName, tree2$tip.label, nomatch = 0)
data2 <- subset(dd, matches != 0)

# Assign species names as rownames for identification
row.names(dd) <- dd$SpecName
row.names(data2) <- data2$SpecName
tree <- tree2

data2 <- ReorderData(tree, data2, taxa.names = "row names")
head(data2)
head(tree$tip.label)

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
americas <- ReorderData(tree, americas, taxa.names="names")
log_avgl <- ReorderData(tree, log_avgl, taxa.names="names")


## Simmap reconstructions of regimes ===================================
# model fitting of Mk models revealed best support for ARD model (see output from step 3)

## Americas
simmap_americ <- make.simmap(tree, americas, model="ARD", nsim=100, pi=c(1,0))

# summarize 100 simulations from simmap
x <- getStates(simmap_americ,"nodes")
sumr <- apply(x, 1, function(x){
  mfv <- table(x)
  return(names(mfv[mfv==max(mfv)]))
})
simmap_americ[[1]]$node.label <- sumr
simmap_AMERIC <- simmap_americ[[1]]


## BBM simulations =================================== 
# = bounded brownian motion between hard reflective bounds

# step 1. estimate z0, sigma2 and bounds from empirical data using BBM package:
## 1.1 create likelihood function from empirical data
#BBM <- lnL_BBMV(tree, log_avgl, Npts=100, bounds=c(min(log_avgl),max(log_avgl)), a=0,b=0,c=0)

## 1.2. find maximum likelihood function
#mle_BBM <- find.mle_FPK(BBM, safe=F) # safe = F ~ 20min computing time; identical results as safe=T

# 1.3 extract estimated sigsq and transform to sigma (needed for BBM trait simulations)
#BBMsigsq <- mle_BBM$par$sigsq
#BBMsigma <- sqrt(BBMsigsq)

# 1.4 extract estimated trait at the root (z0), find the one with the maximum density and save as 'root' for trait simulations
#root_prob<-mle_BBM$root
#root_state<-root_prob %>% slice(which.max(density)) # root state with highest probability: 1.07358869; density: 0.02105593
#root <- root_state$root

# 1.5 extract bounds
#bounds <- mle_BBM$par_fixed$bounds

# Overview / Print variables:
#print(BBMsigsq)
#print(root_state)
#print(bounds)

# step 2. simulate traits for species at the tips of the tree based on model fitted parameter estimates
#BBM_sim <- Sim_FPK(tree, x0 = root, V = rep(0, 100), bounds=bounds, sigma = BBMsigma)


# make OUwie frames [,Species] [,Regime] [,Trait] ===================================

## Americas x BBM:
#data_Americ_BM <- data.frame(species=data2$SpecName, regime = americas, trait = data2$log_BBM_MCC)
#data_Americ_BM$species <- as.character(data_Americ_BM$species)
#data_Americ_BM$regime <- as.factor(data_Americ_BM$regime)
#data_Americ_BM$trait <- as.numeric(as.character(data_Americ_BM$trait))
#str(data_Americ_BM)

## Americas x avgl
data_Americ_avgl <- data.frame(species=data2$SpecName, regime = americas, trait = data2$log_BBM_mean)
data_Americ_avgl$species <- as.character(data_Americ_avgl$species)
data_Americ_avgl$regime <- as.factor(data_Americ_avgl$regime)
data_Americ_avgl$trait <- as.numeric(as.character(data_Americ_avgl$trait))
print(str(data_Americ_avgl))

# run OUwie models ===================================

## Americas x avgl
BM_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_Americ_avgl <- OUwie(simmap_AMERIC, data_Americ_avgl, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
# save output files
ouwie_Americ_avgl <- list(BM_Americ_avgl, BMS_Americ_avgl, OU1_Americ_avgl, OUM_Americ_avgl, OUMV_Americ_avgl, OUMA_Americ_avgl, OUMVA_Americ_avgl)


## Americas x BBM:
#BM_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="BM1", simmap.tree = TRUE, diagn=TRUE)
#BMS_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="BMS", simmap.tree = TRUE, diagn=TRUE)
#OU1_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="OU1", simmap.tree = TRUE, diagn=TRUE)
#OUM_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="OUM", simmap.tree = TRUE, diagn=TRUE)
#OUMV_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
#OUMA_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
#OUMVA_Americ_BM <- OUwie(simmap_AMERIC, data_Americ_BM, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
# save output files
#ouwie_Americ_BM <- list(BM_Americ_BM, BMS_Americ_BM, OU1_Americ_BM, OUM_Americ_BM, OUMV_Americ_BM, OUMA_Americ_BM, OUMVA_Americ_BM)



# Format model output  ===================================

### choose: 
ouwie1 <- ouwie_Americ_avgl
# ouwie1 <- ouwie_Americ_BM 


## extract data for each model: 

### BM1 ===================================
model_fit <- as.data.frame(ouwie1[[1]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[1]][7])
names(alpha_sig2) <- c("estimate")
alpha_sig2$x <- rownames(alpha_sig2)
data_wide <- spread(alpha_sig2, x,  estimate)
alpha_sig2 <- data_wide
# theta + se
theta <- as.data.frame(ouwie1[[1]][8])
theta <- theta[1,]
names(theta) <- c("theta", "theta_se")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[1]][26]))
names(eigvals) <- c("eigvals")
rownames(eigvals) <- NULL
BM1 <- cbind(model_fit, eigvals, alpha_sig2, theta)

## BMS ===================================
model_fit <- as.data.frame(ouwie1[[2]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[2]][7])
names(alpha_sig2) <- c("0", "1")
alpha_sig2$x <- rownames(alpha_sig2)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = c("0", "1")
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[2]][8])
theta <- theta[1,]
names(theta) <- c("theta", "theta_se")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[2]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
BMS <- cbind(model_fit, eigvals, alpha_sig2, theta)


## OU1 ===================================
model_fit <- as.data.frame(ouwie1[[3]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[3]][7])
names(alpha_sig2) <- c("estimate")
rownames(alpha_sig2) <- c("0_alpha", "0_sigma.sq")
alpha_sig2$x <- rownames(alpha_sig2)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = estimate
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[3]][8])
names(theta) <- c("theta", "theta_se")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[3]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
OU1 <- cbind(model_fit, eigvals, alpha_sig2, theta)


### OUM ===================================
model_fit <- as.data.frame(ouwie1[[4]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[4]][7])
names(alpha_sig2) <- c("0", "1")
alpha_sig2$x <- rownames(alpha_sig2)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = c("0", "1")
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[4]][8])
names(theta) <- c("theta", "theta_se")
theta$x <- rownames(theta)
df_wide<-theta %>% pivot_wider(
  names_from = x,
  values_from = c(theta, theta_se)
)
df_wide2<-as.data.frame(df_wide)
names(df_wide2) <- c("theta_0", "theta_1", "theta_se_0", "theta_se_1")
theta <- df_wide2
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[4]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
OUM <- cbind(model_fit, eigvals, alpha_sig2, theta)


### OUMV ===================================
model_fit <- as.data.frame(ouwie1[[5]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[5]][7])
names(alpha_sig2) <- c("0", "1")
alpha_sig2$x <- rownames(alpha_sig2)
#data_wide <- spread(alpha_sig2, x,  estimate)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = c("0","1")
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[5]][8])
names(theta) <- c("theta", "theta_se")
theta$x <- rownames(theta)
df_wide<-theta %>% pivot_wider(
  names_from = x,
  values_from = c(theta, theta_se)
)
df_wide2<-as.data.frame(df_wide)
theta <- df_wide2
names(theta) <- c("theta_0", "theta_1", "theta_se_0", "theta_se_1")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[5]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
OUMV <- cbind(model_fit, eigvals, alpha_sig2, theta)


### OUMA ===================================
model_fit <- as.data.frame(ouwie1[[6]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[6]][7])
names(alpha_sig2) <- c("0", "1")
alpha_sig2$x <- rownames(alpha_sig2)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = c("0","1")
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[6]][8])
names(theta) <- c("theta", "theta_se")
theta$x <- rownames(theta)
df_wide<-theta %>% pivot_wider(
  names_from = x,
  values_from = c(theta, theta_se)
)
df_wide2<-as.data.frame(df_wide)
theta <- df_wide2
names(theta) <- c("theta_0", "theta_1", "theta_se_0", "theta_se_1")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[6]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
OUMA <- cbind(model_fit, eigvals, alpha_sig2, theta)

## OUMVA ===================================
model_fit <- as.data.frame(ouwie1[[7]][1:6])
# alpha sigma2
alpha_sig2 <- as.data.frame(ouwie1[[7]][7])
names(alpha_sig2) <- c("0", "1")
alpha_sig2$x <- rownames(alpha_sig2)
#data_wide <- spread(alpha_sig2, x,  estimate)
df_wide<-alpha_sig2 %>% pivot_wider(
  names_from = x,
  values_from = c("0","1")
)
df_wide2<-as.data.frame(df_wide)
alpha_sig2 <- df_wide2
# theta + se
theta <- as.data.frame(ouwie1[[7]][8])
names(theta) <- c("theta", "theta_se")
theta$x <- rownames(theta)
df_wide<-theta %>% pivot_wider(
  names_from = x,
  values_from = c(theta, theta_se)
)
df_wide2<-as.data.frame(df_wide)
theta <- df_wide2
names(theta) <- c("theta_0", "theta_1", "theta_se_0", "theta_se_1")
# eigvals
eigvals <- as.data.frame(unlist(ouwie1[[5]][26]))
names(eigvals) <- c("eigvals")
eigvals$x <- rownames(eigvals)
df_wide <- eigvals %>% pivot_wider(
  names_from = x,
  values_from = eigvals
)
eigvals <- as.data.frame(df_wide)
OUMVA <- cbind(model_fit, eigvals, alpha_sig2, theta)


# Merge and exclude odd fit and parameter estimates =============================

mods <- merge(OUMVA, OUMA, all=T)
mods <- merge(mods, OUMV, all=T)
mods <- merge(mods, OUM, all=T)
mods<- merge(mods, OU1, all=T)
mods<-merge(mods, BMS, all=T)
mods<- merge(mods, BM1, all=T)
mods$id <- task_id

# merging of columns with NAs (where only one parameter was estimated)
mods$theta0 <- ifelse(!is.na(mods$theta), mods$theta, mods$theta_0)
mods$theta_se0 <- ifelse(!is.na(mods$theta_se), mods$theta_se, mods$theta_se_0)
mods$sigma0 <- ifelse(is.na(mods$sigma.sq), mods$`0_sigma.sq`, mods$sigma.sq)
mods$eigval1 <- ifelse(is.na(mods$eigval1), mods$eigvals, mods$eigval1)

# re-ordering data
mods2 <- mods[,c("id", "model","loglik", "AIC", "AICc", "BIC", "param.count", "eigval1", "eigval2", "eigval3", "0_alpha", "1_alpha", "sigma0", "1_sigma.sq", "theta0", "theta_se0", "theta_1", "theta_se_1")]
mods2[is.na(mods2)] <- 0 # set NAs = 0 

# negative eigenvalues, odd estimates, negative AICc (convergence problems)
mods3 <- mods2[!is.na(mods2$eigval1) & mods2$eigval1 >= 0 & !is.na(mods2$eigval2) & mods2$eigval2 >= 0 & !is.na(mods2$eigval3) & mods2$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights
fit<-geiger::aicw(mods4$AICc)
mods5<-cbind(mods4, fit)

# select best model based on models without convergence problems
bestfit <- mods5[mods5$delta == 0, ]


# write output to file ===============================
saveRDS(bestfit, file.path(output_dir, paste0("Bestfit_", task_id, ".rds")))
saveRDS(mods2, file.path(output_dir, paste0("Raw_", task_id, ".rds")))
saveRDS(mods5, file.path(output_dir, paste0("Processed_", task_id, ".rds")))





## The following was run on a local computer ===================

setwd("~/Afr_x_avgl/")
setwd("~/Afr_x_BM/")

# Function to merge data frames in list
multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}  


# Evaluate convergence problems ===================
out <- list()
input_files <-gsub("\\.rds$","", list.files(pattern=".*Raw.*\\.rds$"))  
n_input<-(1:length(input_files))

for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  out[[i]] <- output
}

out2 <- list()
input_files <-gsub("\\.rds$","", list.files(pattern=".*Processed.*\\.rds$"))  
for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  out2[[i]] <- output
}

# Merge data frames 
merged_data <- multi_join(out, full_join) # raw
merged_data_raw <- multi_join(out2, full_join) # without convergence problems

head(merged_data)
str(merged_data)
str(merged_data_raw)
summary(merged_data)
n_input

# Convergence Problems:
unique(merged_data_raw$model)
unique(merged_data$model)

issues<-setdiff(merged_data_raw,merged_data[,1:18])
issues

write.table(issues, file = "issues.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)
sum(table(issues$model))
n_input

# write to .txt
write.table(merged_data, file = "merged_data_BM.txt", sep = "\t", row.names = F, col.names = T)
write.table(merged_data, file = "merged_data_AVGL.txt", sep = "\t", row.names = F, col.names = T)