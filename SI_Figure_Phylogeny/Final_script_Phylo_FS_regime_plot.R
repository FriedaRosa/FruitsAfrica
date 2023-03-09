setwd("../Desktop/Phylo_FS_regime_plot/")

# I. Read in libraries  ===================================
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(evobiR))


# II. Read in Simmap reconstructions for geographical area ===================================

# a) Africa ####
simmap_afr <- readRDS("out/simmap_africa_500.rds")

# simmap_afr <- make.simmap(tree, africa, model="ARD", nsim=500, pi=c(1,0))
# Q =         1_world     2_africa
# 1_world  -0.001475705  0.001475705
# 2_africa  0.002463797 -0.002463797
# (estimated using likelihood);
# and (mean) root node prior probabilities
# pi = 1_world 2_africa 
#         1        0

# b) Americas ####
simmap_amer <- readRDS("out/simmap_americas_500.rds")

# simmap_americas <- make.simmap(tree, americas, model="ARD", nsim=500, pi=c(1,0))
# Q =           1_world    2_americas
# 1_world    -0.0005803574  0.0005803574
# 2_americas  0.0020967723 -0.0020967723
# (estimated using likelihood);
# and (mean) root node prior probabilities
# pi = 1_world 2_americas 
#       1          0 

# c) Asia ####
simmap_asia <- readRDS("out/simmap_asia_500.rds")

# simmap_asia <- make.simmap(tree, asia, model="ARD", nsim=500, pi=c(1,0))
# Q =       1_world      2_asia
# 1_world -0.00187808  0.00187808
# 2_asia   0.00720775 -0.00720775
# (estimated using likelihood);
# and (mean) root node prior probabilities
# pi = 1_world  2_asia 
#       1       0 

# III. Left plots: Posterior Densities of geographical area from 500 simulations

# a) Africa ####

cols_simmap<-setNames(c("black","#FDE725FF"), c("1_world", "2_africa"))
#density_afr<-densityMap(simmap_afr)
density_afr2 <- setMap(density_afr,colors=c("grey", "black", "#FDE725FF"))
plot(density_afr2, lwd=2, ftype="off", xlim=c(0,1.25*max(nodeHeights(simmap_afr[[1]]))), ftype="off", outline=F, legend=F,  mar=c(4.1,1.1,1.1,0) );
add.color.bar(50,cols=density_asi2$cols,title="",prompt=F, x=0, y=30)
# geo.legend(colors=c("grey", "black", "#FDE725FF"), alpha=0.2, cex=0.9);
# axisPhylo()
title("Africa")

save.image("out/simmap_afr_plot.RData")

# b) Americas ####
cols_simmap<-setNames(c("black","#FDE725FF"), c("1_world", "2_americas"))
#density_amer<-densityMap(simmap_amer)
density_amer2 <- setMap(density_amer,colors=c("grey", "black", "#FDE725FF"))
plot(density_amer2, lwd=2, ftype="off", xlim=c(0,1.25*max(nodeHeights(simmap_amer[[1]]))), ftype="off", outline=F,  legend=F,  mar=c(4.1,1.1,1.1,0)  );
add.color.bar(50,cols=density_asi2$cols,title="",prompt=F, x=0, y=30)
# geo.legend(colors=c("grey", "black", "#FDE725FF"), alpha=0.2, cex=0.9);
# axisPhylo()
title("Americas")

save.image("out/simmap_afr_amer_plot.RData")

# c) Asia ####
cols_simmap<-setNames(c("black","#FDE725FF"), c("1_world", "2_asia"))
#density_asi<-densityMap(simmap_asia)
density_asi2 <- setMap(density_asi,colors=c("grey", "black", "#FDE725FF"))
plot(density_asi2, lwd=2, ftype="off", xlim=c(0,1.25*max(nodeHeights(simmap_asia[[1]]))), ftype="off", outline=F,legend=F,  mar=c(4.1,1.1,1.1,0));
add.color.bar(50,cols=density_asi2$cols,title="",prompt=F, x=0, y=30)
title("Asia")

# geo.legend(colors=c("grey", "black", "#FDE725FF"), alpha=0.2, cex=0.9);
# axisPhylo()

save.image("out/simmap_afr_amer_asi_plot.RData")

# IV. Right plots: Reconstruction of fruit size with tip-bars for species at the tips

# Read in Fruit Size data
dd <- read.csv("data/Step1_duplicated_regimes_checked.csv", stringsAsFactors = T)

dd2 <- dd[,c("SpecName","accGenus","PalmTribe","PalmSubfamily","AverageFruitLength_cm","Americas",
             "Asia","Africa","Madagascar","Pacific","accRealm")]
dd <- dd2 %>% distinct()
dd$SpecName <- as.character(dd$SpecName)
dd$AverageFruitLength_cm1 <- as.numeric(as.character(dd$AverageFruitLength_cm))
dd2 <- dd[!is.na(dd["AverageFruitLength_cm1"]),]
dd <- as.data.frame(dd2)
dd$log_FS <- log(dd$AverageFruitLength_cm1) #log-transformation
row.names(dd) <- dd$SpecName
# Read in tree ================================
tree <- read.nexus("data/TREE")

setdiff(tree$tip.label, dd$SpecName) #502 sp.
tree2 <- drop.tip(tree, setdiff(tree$tip.label, row.names(dd)))
tree <- ladderize(tree2) # 2037 tips
tree <- tree2
data2 <- ReorderData(tree, dd, taxa.names = "row names")

# Trait vector ==================
log_avgl <- as.numeric(data2$log_FS)
names(log_avgl) <- row.names(data2)

# Plotting ==============
layout(matrix(c(1,2),1,2),c(0.3,0.6))
par(bg="white",
    fg="black",
    col.lab="white",
    col.axis="white",
    xpd=TRUE)

library(viridis); library(scales)
colors5 <- viridis_pal()(6)
cols <- viridis_pal()(2)
par(mar=c(4.1,0,1.1,1.1))
barplot(log_avgl[sim_pc$tree$tip.label],horiz=TRUE,width=2,space=0,names="")


#sim_p <- contMap(tree, log_avgl, ftype="off", res=100, type="phylogram")
#sim_pc <- setMap(sim_p, colors=c("grey","black", "#FDE725FF"))
#sim_pc2 <- setMap(sim_p, colors=colors5)

plot.contMap(sim_pc, ftype="off", res=200, type="phylogram", direction = "leftwards", outline=F, mar=c(4.1,1.1,1.1,0), legend=F)
add.color.bar(50,cols=sim_pc$cols,lims=sim_pc$lims,title="", prompt=F, x=50, y=0)

plot.contMap(sim_pc2, ftype="off", res=200, type="phylogram", direction = "leftwards", outline=F, mar=c(4.1,1.1,1.1,0), legend=F)
add.color.bar(50,cols=sim_pc2$cols,lims=sim_pc2$lims,title="", prompt=F, x=50, y=0)




