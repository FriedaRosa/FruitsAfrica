## Sensitivity Analysis: Removal of 50, 70, 90 percentile largest fruits from Africa, see if it has an effect on OUwie outcomes

suppressPackageStartupMessages(library(OUwie))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(evobiR))
#suppressPackageStartupMessages(library(BBMV))

# define arguments ===================================
args <- commandArgs(trailingOnly=TRUE)
nex_file <- args[1]
output_dir <- args[2]
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# read tree ===================================
tree <- read.nexus(nex_file)

# read data ===================================
dd <- read.csv2("/home/woelke/palms_final.csv")
str(dd)
length(unique(dd$SpecName)) # 2037

sum(is.na(dd$AverageFruitLength_cm))
dd2 <- subset(dd, !is.na(AverageFruitLength_cm))
dd2$logFS <- log(dd2$AverageFruitLength_cm) #log-transformation


library(plyr)
dd3 <- plyr::ddply(.data = dd2, .variables = .(SpecName), .fun = summarise, 
                   logFS = logFS[1],
                   accAfrica = sum(accAfrica == 1)>0,
                   Asia = sum(Asia == 1)>0,
                   Americas = sum(Americas == 1)>0)

rownames(dd3)<- dd3$SpecName
nrow(dd3) #2037 species

# [JY NOTES] only 2 species in the dataset are found in more than 2 continent
table(rowSums(dd3[c("accAfrica", "Americas", "Asia")])) 


tree2 <- drop.tip(tree, setdiff(tree$tip.label, dd3$SpecName))
tree2 <- ladderize(tree2)
length(tree2$tip.label) # 2037

dd3 <- ReorderData(tree2, dd3)

# Number of species with megafaunal fruits per realm
nrow(dd3[dd3$logFS >= log(4),])
nrow(dd3[dd3$logFS >= log(4) & dd3$accAfrica == TRUE,])
# nrow(dd3[dd3$logFS >= log(4) & dd3$Asia == TRUE,])
# nrow(dd3[dd3$logFS >= log(4) & dd3$Americas == TRUE,])

# Create named vector of Africa distribution
africa <- ifelse(dd3$accAfrica == 0, "1_world", "2_africa")
names(africa) <- rownames(dd3)

palm_phylo <- tree2

## ANCESTRAL STATE RECONSTRUCTION ====================
simmap_afr <- make.simmap(palm_phylo, africa, model="ARD", nsim=100, pi=c(1,0))

# summarize 100 simulations from simmap
x <- getStates(simmap_afr,"nodes")
sumr <- apply(x, 1, function(x){
  mfv <- table(x)
  return(names(mfv[mfv==max(mfv)]))
})

simmap_afr[[1]]$node.label <- sumr
simmap_AFR <- simmap_afr[[1]]


#### WIP



## Remove 25%, 50%, 75% of fruits > 4 cm in Africa from data

mega_afr <- dd3 %>% filter(accAfrica == "TRUE" & logFS >= log(4)) # fruits > 4 cm in Africa

# calculate 25, 50, 75 percentiles of fruits > 4 cm and make species lists for removal

threshold25 <- quantile(mega_afr$logFS, probs=c(0.25))
extinctions25 <- mega_afr$SpecName[mega_afr$accAfrica == 1 & mega_afr$logFS >= threshold25  & mega_afr$Asia == FALSE & mega_afr$Americas== FALSE]

#t hreshold50 <- quantile(mega_afr$logFS, probs=c(0.50))
# extinctions50 <- mega_afr$SpecName[mega_afr$accAfrica == 1 & mega_afr$logFS >= threshold50  & mega_afr$Asia == FALSE & mega_afr$Americas== FALSE]

# threshold75 <- quantile(mega_afr$logFS, probs=c(0.75))
# extinctions75 <- mega_afr$SpecName[mega_afr$accAfrica == 1 & mega_afr$logFS >= threshold75  & mega_afr$Asia == FALSE & mega_afr$Americas== FALSE]

# Prune simmap trees:

simmap25 <- drop.tip.simmap(simmap_AFR, extinctions25)
# simmap50 <- drop.tip.simmap(simmap_AFR, extinctions50)
# simmap75 <- drop.tip.simmap(simmap_AFR, extinctions75)


# Create named vector of log-transformed trait values
log_avgl <- dd3$logFS
names(log_avgl) <- rownames(dd3)

data_Afr_avgl <- data.frame(species=names(log_avgl), regime = factor(africa), trait = log_avgl)

data_Afr_avgl_25 <- subset(data_Afr_avgl, !species %in% extinctions25)
# data_Afr_avgl_50 <- subset(data_Afr_avgl, !species %in% extinctions50)
# data_Afr_avgl_75 <- subset(data_Afr_avgl, !species %in% extinctions75)


BM_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="BM1", simmap.tree = TRUE, diagn=TRUE)
BMS_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="BMS", simmap.tree = TRUE, diagn=TRUE)
OU1_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="OU1", simmap.tree = TRUE, diagn=TRUE)
OUM_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="OUM", simmap.tree = TRUE, diagn=TRUE)
OUMV_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
OUMA_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
OUMVA_Afr_avgl <- OUwie(simmap25, data_Afr_avgl_25, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)

ouwie_Afr_avgl_25 <- list(BM_Afr_avgl, BMS_Afr_avgl, OU1_Afr_avgl, OUM_Afr_avgl, OUMV_Afr_avgl, OUMA_Afr_avgl, OUMVA_Afr_avgl)



# BM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="BM1", simmap.tree = TRUE, diagn=TRUE)
# BMS_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="BMS", simmap.tree = TRUE, diagn=TRUE)
# OU1_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="OU1", simmap.tree = TRUE, diagn=TRUE)
# OUM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="OUM", simmap.tree = TRUE, diagn=TRUE)
# OUMV_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
# OUMA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
# OUMVA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_50, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)

# ouwie_Afr_avgl_50 <- list(BM_Afr_avgl, BMS_Afr_avgl, OU1_Afr_avgl, OUM_Afr_avgl, OUMV_Afr_avgl, OUMA_Afr_avgl, OUMVA_Afr_avgl)


# BM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="BM1", simmap.tree = TRUE, diagn=TRUE)
# BMS_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="BMS", simmap.tree = TRUE, diagn=TRUE)
# OU1_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="OU1", simmap.tree = TRUE, diagn=TRUE)
# OUM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="OUM", simmap.tree = TRUE, diagn=TRUE)
# OUMV_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
# OUMA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
# OUMVA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl_75, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)

# ouwie_Afr_avgl_75 <- list(BM_Afr_avgl, BMS_Afr_avgl, OU1_Afr_avgl, OUM_Afr_avgl, OUMV_Afr_avgl, OUMA_Afr_avgl, OUMVA_Afr_avgl)

# Format model output  ===================================

### choose: 
ouwie1 <- ouwie_Afr_avgl_25
# ouwie1 <- ouwie_Afr_avgl_50 
# ouwie1 <- ouwie_Afr_avgl_75 

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
saveRDS(bestfit, file.path(output_dir, paste0("Bestfit_25", task_id, ".rds")))
saveRDS(mods2, file.path(output_dir, paste0("Raw_25", task_id, ".rds")))
saveRDS(mods5, file.path(output_dir, paste0("Processed_25", task_id, ".rds")))
