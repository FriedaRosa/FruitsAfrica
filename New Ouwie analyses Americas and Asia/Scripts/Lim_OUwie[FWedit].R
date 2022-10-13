suppressPackageStartupMessages(library(OUwie))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(evobiR))
#suppressPackageStartupMessages(library(BBMV))
library(devtools)
#install_github("fcboucher/BBMV")
#"fcboucher/BBMV"
#setwd("~/Dropbox/manuscripts/2022_onstein_fruit-size-evolution/analysis/")

setwd("~/GitHub/FruitsAfrica/New Ouwie analyses Americas and Asia/Data")


# define arguments ===================================
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# read tree ===================================
tree <- read.nexus("TREE")

# [JY NOTES] Let's just reproduce the example with just the MCC for now

# read data ===================================
dd <- read.csv("palms_final.csv", stringsAsFactors = TRUE)

head(subset(dd, SpecName == "Chamaedorea_ernesti-augusti"))
subset(dd, SpecName == "Nypa_fruticans")
length(unique(dd$SpecName)) # 2037

# [JY NOTES]
# 1) The dataframe is actually a merge of BC distribution by palm species
# 2) Some species already removed (e.g., Nypa)

dd$Africa == dd$accAfrica 
dd[800,]
# 3) Madagascar was not considered part of Africa

sum(is.na(dd$AverageFruitLength_cm))
dd2 <- subset(dd, !is.na(AverageFruitLength_cm))
dd2$logFS <- log(dd2$AverageFruitLength_cm) #log-transformation


# JY: Let's make the dataset species level
library(plyr)
dd3 <- plyr::ddply(.data = dd2, .variables = .(SpecName), .fun = summarise, 
                   logFS = logFS[1],
                   accAfrica = sum(accAfrica == 1)>0,
                   Asia = sum(Asia == 1)>0,
                   Americas = sum(Americas == 1)>0)
nrow(dd3) #2037 species

# [JY NOTES] only 2 species in the dataset are found in more than 2 continent
table(rowSums(dd3[c("accAfrica", "Americas", "Asia")])) 

# match data and tree ===================================
row.names(dd3) <- dd3$SpecName
setdiff(tree$tip.label, dd3$SpecName)

# [JY NOTES] Quite a lot of species missing!
# A lot of calamus especially = possibly due to taxonomic revision by Henderson
# Also the code does not work because row names need to be unique and there are duplicates

# Drop these species from the tree, because we need matching data between tree and traits
tree2 <- drop.tip(tree, setdiff(tree$tip.label, dd3$SpecName))
tree2 <- ladderize(tree2)
length(tree2$tip.label) # 2037

sum(dd3$SpecName %in% tree2$tip.label) # 2037

dd3 <- ReorderData(tree2, dd3)
head(dd3$SpecName)
head(tree2$tip.label)

# Create named vector of Africa distribution
africa <- ifelse(dd3$accAfrica == 0, "1_world", "2_africa")
names(africa) <- rownames(dd3)

## Create named vector of log-transformed trait values
log_avgl <- dd3$logFS
names(log_avgl) <- rownames(dd3)

## Sort data like phylogeny to prevent errors
africa <- ReorderData(tree2, africa, taxa.names = "names")
log_avgl <- ReorderData(tree2, log_avgl, taxa.names = "names")

# [JY NOTES] Reordering not really necessary since it's already done on line 66
# I also don't like overwriting variable names but will do it here so the rest of
# the code is reproducible without modification
palm_phylo <- tree2

## ALTERNATIVE ANALYSIS USING PGLS ====================
nrow(dd3[dd3$logFS >= log(4),])
nrow(dd3[dd3$logFS >= log(4) & dd3$accAfrica == TRUE,])
nrow(dd3[dd3$logFS >= log(4) & dd3$Asia == TRUE,])
nrow(dd3[dd3$logFS >= log(4) & dd3$Americas == TRUE,])

library(phylolm)
afr_BM <- phylolm(logFS ~ accAfrica, data = dd3, phy = palm_phylo, model = "BM")
afr_OU <- phylolm(logFS ~ accAfrica, data = dd3, phy = palm_phylo, model = "OUfixedRoot")
AIC(afr_BM); AIC(afr_OU)

asia_BM <- phylolm(logFS ~ Asia, data = dd3, phy = palm_phylo, model = "BM")
asia_OU <- phylolm(logFS ~ Asia, data = dd3, phy = palm_phylo, model = "OUfixedRoot")
AIC(asia_BM); AIC(asia_OU)

amer_BM <- phylolm(logFS ~ Americas, data = dd3, phy = palm_phylo, model = "BM")
amer_OU <- phylolm(logFS ~ Americas, data = dd3, phy = palm_phylo, model = "OUfixedRoot")
AIC(amer_BM); AIC(amer_OU)

extractModResid <- function(mod){
  EffectSize <- coefficients(mod)[2]
  LowerCI <- confint(mod)[2,][1]
  UpperCI <- confint(mod)[2,][2]
  data.frame(EffectSize, LowerCI, UpperCI)
}

all_mod_list <- list(afr_OU, asia_OU, amer_OU)
all_mod_coef <- do.call("rbind",lapply(all_mod_list, extractModResid))
all_mod_coef$Ext <- 1
all_mod_coef$Continent <- c("Africa", "Asia", "Americas")

extinctionModel <- function(x, p, continent){
  # x = data frame
  # tree
  # p = percentile to remove
  # x = dd3
  # tree = palm_phylo
  # p = 0.9
  if(continent == "Africa"){
    threshold <- quantile(x$logFS[x$accAfrica == 1 & x$logFS >= log(4)], probs = p)  
    extirpations <- x$SpecName[x$accAfrica == 1 & x$logFS >= threshold & (x$Asia == TRUE | x$Americas== TRUE) ]
    extinctions <- x$SpecName[x$accAfrica == 1 & x$logFS >= threshold & x$Asia == FALSE & x$Americas== FALSE]
    if(length(extirpations) > 0){
      x$accAfrica[x$SpecName %in% extirpations] <- FALSE
    }
  }
  if(continent == "America"){
    threshold <- quantile(x$logFS[x$Americas == 1 & x$logFS >= log(4)], probs = p)
    extirpations <- x$SpecName[x$Americas == 1 & x$logFS >= threshold & (x$Asia == TRUE | x$accAfrica == TRUE) ]
    extinctions <- x$SpecName[x$Americas == 1 & x$logFS >= threshold & x$Asia == FALSE & x$accAfrica == FALSE]
    if(length(extirpations) > 0){
      x$Americas[x$SpecName %in% extirpations] <- FALSE
    }
  }
  if(continent == "Asia"){
    threshold <- quantile(x$logFS[x$Americas == 1 & x$logFS >= log(4)], probs = p)  
    extirpations <- x$SpecName[x$Asia == 1 & x$logFS >= threshold & (x$Americas == TRUE | x$accAfrica == TRUE) ]
    extinctions <- x$SpecName[x$Asia == 1 & x$logFS >= threshold & x$Americas == FALSE & x$accAfrica == FALSE]
    if(length(extirpations) > 0){
      x$Asia[x$SpecName %in% extirpations] <- FALSE
    }
  }
  
  x <- subset(x, !SpecName %in% extinctions)
  return(x)
}

pruneTree <- function(x, tree){
  tree <- drop.tip(tree, tip = setdiff(tree$tip.label, x$SpecName))
  return(tree)
}

fitPhylolm <- function(x, tree, continent){
  if(continent == "Africa"){
    mod_form <- "logFS ~ accAfrica"
  }
  if(continent == "America"){
    mod_form <- "logFS ~ Americas"
  }
  if(continent == "Asia"){
    mod_form <- "logFS ~ Asia"
  }
  phylolm(as.formula(mod_form), data = x, phy = tree, model = "OUfixedRoot")
}

africa_ext_list <- lapply(c(0.9, 0.75, 0.5, 0.25, 0.1), FUN = extinctionModel, continent = "Africa", x = dd3)
americas_ext_list <- lapply(c(0.9, 0.75, 0.5, 0.25, 0.1), FUN = extinctionModel, continent = "America", x = dd3)
asia_ext_list <- lapply(c(0.9, 0.75, 0.5, 0.25, 0.1), FUN = extinctionModel, continent = "Asia", x = dd3)

africa_ext_tree_list <- lapply(africa_ext_list, FUN = pruneTree, tree = palm_phylo)
americas_ext_tree_list <- lapply(americas_ext_list, FUN = pruneTree, tree = palm_phylo)
asia_ext_tree_list <- lapply(asia_ext_list, FUN = pruneTree, tree = palm_phylo)

lapply(africa_ext_list, FUN = function(x){nrow(x)})
lapply(africa_ext_list, FUN = function(x){sum(x$accAfrica)})
lapply(africa_ext_tree_list, FUN = function(x){length(x$tip.label)})
lapply(americas_ext_list, FUN = function(x){nrow(x)})
lapply(americas_ext_list, FUN = function(x){sum(x$Americas)})
lapply(americas_ext_tree_list, FUN = function(x){length(x$tip.label)})
lapply(asia_ext_list, FUN = function(x){nrow(x)})
lapply(asia_ext_list, FUN = function(x){sum(x$Asia)})
lapply(asia_ext_tree_list, FUN = function(x){length(x$tip.label)})

africa_ext_mod_list <- mapply(africa_ext_list, africa_ext_tree_list, FUN = fitPhylolm, continent = "Africa", SIMPLIFY = FALSE)
america_ext_mod_list <- mapply(americas_ext_list, americas_ext_tree_list, FUN = fitPhylolm, continent = "America", SIMPLIFY = FALSE)
asia_ext_mod_list <- mapply(asia_ext_list, asia_ext_tree_list, FUN = fitPhylolm, continent = "Asia", SIMPLIFY = FALSE)

all_ext_mod_list <- c(lapply(africa_ext_mod_list, FUN = extractModResid),
                      lapply(america_ext_mod_list, FUN = extractModResid),
                      lapply(asia_ext_mod_list, FUN = extractModResid))
all_ext_mod_df <- do.call("rbind",all_ext_mod_list)
all_ext_mod_df$Ext <- rep(c(0.9, 0.75, 0.5, 0.25, 0.1), 3)
all_ext_mod_df$Continent <- rep(c("Africa", "Americas", "Asia"), each = 5)

res <- rbind(all_mod_coef, all_ext_mod_df)


library(ggplot2)
group.colors <- c(Africa = "#ecd70d",Americas = "#21609a", Asia = "#843c54")


palm_ext_mod_plot <- ggplot(res) + 
  geom_point(aes(y = EffectSize, x = factor(Ext), color = Continent, size = 1.2), position = position_dodge(width = 0.5)) + 
  scale_color_manual("Continent", labels = c("Africa", "Americas", "Asia"), values = group.colors)+
  geom_errorbar(aes(ymax = UpperCI, ymin = LowerCI, x = factor(Ext), color = Continent),
                size =1.4,
                position = position_dodge(width = 0.5),
                width = 0) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=11.5), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  theme(plot.subtitle = element_text(vjust = -6, hjust=0.015))+
  guides(size = "none")+
  xlab("Degree of Extinction")
palm_ext_mod_plot

ggsave(palm_ext_mod_plot, filename = "~/Desktop/palm_ext_mod_plot.pdf")

## ANCESTRAL STATE RECONSTRUCTION ====================
#simmap_afr <- make.simmap(tree, africa, model="ARD", nsim=100, pi=c(1,0))
# saveRDS(simmap_afr, file = "simmap_afr.rds")
# simmap_afr <- readRDS(file = "simmap_afr.rds")

#plotSimmap(simmap_afr[[1]])
# [JY NOTES] Using an Mk model to model biogeography (probably ok of just looking at Africa
# but if using other continents separately, might make more sense to do a model with all possible
# states (Africa, Asia, Americas and potentially "others", species not in any of three bins)

# summarize 100 simulations from simmap
# x <- getStates(simmap_afr,"nodes")
# sumr <- apply(x, 1, function(x){
#   mfv <- table(x)
#   return(names(mfv[mfv==max(mfv)]))
# })
# simmap_afr[[1]]$node.label <- sumr
# simmap_AFR <- simmap_afr[[1]]



# [JY NOTES] each column = 1 simulation, each row = node states (n - 1)
# Apply code calculates the most frequent node state across simulations
# Assigns to the first tree; using the first one as THE one
# But simmap also does assignments to edges, so how does this 


## FIT TRAIT MODELS ====================
# [JY NOTES] Original code calls for "avgl" object, which doesn't exist
# data_Afr_avgl <- data.frame(species=names(log_avgl), regime = factor(africa), trait = log_avgl)
# 
# BM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="BM1", simmap.tree = TRUE, diagn=TRUE)
# BMS_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="BMS", simmap.tree = TRUE, diagn=TRUE)
# OU1_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="OU1", simmap.tree = TRUE, diagn=TRUE)
# OUM_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="OUM", simmap.tree = TRUE, diagn=TRUE)
# OUMV_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="OUMV", simmap.tree = TRUE, diagn=TRUE)
# OUMA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="OUMA", simmap.tree = TRUE, diagn=TRUE)
# OUMVA_Afr_avgl <- OUwie(simmap_AFR, data_Afr_avgl, model="OUMVA", simmap.tree = TRUE, diagn=TRUE)
# 
# ouwie_Afr_avgl <- list(BM_Afr_avgl, BMS_Afr_avgl, OU1_Afr_avgl, OUM_Afr_avgl, OUMV_Afr_avgl, OUMA_Afr_avgl, OUMVA_Afr_avgl)
# saveRDS(ouwie_Afr_avgl, file = "ouwie_Afr_avgl.rds")

#ouwie_Afr_avgl <- readRDS(ouwie_Afr_avgl, file = "ouwie_Afr_avgl.rds")

## SIMULATED EXTINCTION ====================






