library(RColorBrewer)

tree <- read.nexus("../MCC/TREE")


tree2 <- drop.tip(tree, setdiff(tree$tip.label, rownames(X2)))
tree <- ladderize(tree2)


data2 <- read.csv("../../Frieda/Overview_simulated_traits.csv")
data2 <- data2[-1]
columns <- c("SpecName", "genus", "Americas", "Asia", "Africa", "Madagascar", "Pacific", "accRealm", "accAfrica", "moist0_dry1", "log_BBM_mean")

sim <- data2 %>% dplyr::select(any_of(columns))
sim <- unique(sim)


#rownames(sim) <- sim$SpecName

dups_names <- sim[which(duplicated(sim$SpecName)),]$SpecName

dups2 <- sim %>% filter(SpecName %in% c(dups_names)) # 41 species with duplicated rows
non_dups2 <- sim %>% filter(!SpecName %in% c(dups_names)) # 1826 species with unique rows
#write.csv(dups2, "../duplicated_afr.csv")

sim_nodups <- non_dups2 %>% filter(!SpecName %in% c(non_dups2)) %>% select(-c(accRealm))

unique <- read.csv("../duplicated_realms_unique.csv", sep=";") # read corrected data back in (merged rows in excel) --> 20 unique entries left
unique <- unique[-1]

sim2 <- rbind(sim_nodups, unique)
unique(sim2)
dups2 <- sim2[which(duplicated(sim2$SpecName)),]$SpecName

str(sim2)
#sim2[2:9] <- factor(sim2[2:9])
sim2[2:9] <- lapply(sim2[2:9], factor)

data2 <-sim2
rownames(data2) <- data2$SpecName

X <- data.frame(data2)
tree2 <- drop.tip(tree, setdiff(tree$tip.label, rownames(X)))
tree <- ladderize(tree2)


X3 <- X[c(1,10)]
matches <- match(X3$SpecName, tree$tip.label, nomatch = 0)
X3 <- subset(X3, matches != 0)

BBM <- setNames(X3$log_BBM_mean, rownames(X3))
ace <- readRDS("../BBM_ACE.rds")


min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

BBM <- min_max_norm(BBM)
ace <- min_max_norm(ace)
hist(BBM)
hist(ace)
ace2 <- fastAnc(tree,BBM,CI=F)
ace2 <- min_max_norm(ace2)
hist(ace2)

pdf("BM10FS_Reconst_thick20.pdf", height=20, width=11)
par(bg="white",
    fg="black",
    col.lab="white",
    col.axis="white",
    xpd=TRUE);

sim_p <- contMap(tree, BBM, method="user", anc.states=ace2, ftype="off", res=100, type="fan")
sim_p <- contMap(tree, BBM, ftype="off", res=100, type="fan")
sim_pc <- setMap(sim_p, col=colors5)
plot.contMap(sim_pc, ftype="off", res=200, type="fan")

###
X2 <- X[c(3:9)]

md <- as.factor(setNames(X2[,7],rownames(X2)))
amer <- as.factor(setNames(X2[,1], rownames(X2)))
asi <<- as.factor(setNames(X2[,2], rownames(X2)))
afr <- as.factor(setNames(X2[,3], rownames(X2)))
mad <- as.factor(setNames(X2[,4], rownames(X2)))
pac <- as.factor(setNames(X2[,5], rownames(X2)))
accAfr <-as.factor(setNames(X2[,6], rownames(X2)))

##

colors5 <- viridis_pal()(6)

col.hab<-setNames(c("darkgreen","grey"),levels(md))
col.accAfr <- setNames(c("white", colors5[1]), levels(accAfr))
col.afr<-setNames(c("white",colors5[2]),levels(afr))
col.amer<-setNames(c("white",colors5[3]),levels(amer))
col.asi<-setNames(c("white",colors5[4]),levels(asi))
col.pac<-setNames(c("white",colors5[5]),levels(pac))
col.mad<-setNames(c("white",colors5[6]),levels(mad))


## 1. load dataset for clade annotations

tax <- read.csv2("~/share/groups/ea/Frieda/MCC/tips_nodes_clades.csv")
tax <- tax[1:5]

matches <- match(tax$x, tree$tip.label, nomatch = 0)
tax2 <- subset(tax, matches != 0)
tax <- tax2

str(tax)
tax$x <- as.character(tax$x)
tax$PalmSubfamily <- as.character(tax$PalmSubfamily)
tax$accGenus <- as.character(tax$accGenus)
tax$PalmTribe <- as.character(tax$PalmTribe)


length(unique(tax$PalmSubfamily))
length(unique(tax$accGenus))
length(unique(tax$PalmTribe)) #29

# Subfamilies
Arecoidea <- subset(tax, PalmSubfamily == "Arecoideae")
Calamoideae <-subset(tax, PalmSubfamily == "Calamoideae")
Ceroxyloideae <- subset(tax, PalmSubfamily == "Ceroxyloideae")
Coryphoideae <-subset(tax, PalmSubfamily == "Coryphoideae")

# function to find mrca node for clades
fmrca=geiger:::.mrca  

#Palm tribes for annotation:
tribes <- c(unique(tax$PalmTribe))
tribes_overlap <-c("Phytelepheae","Cyclosphatheae", "Roystoneeae", "Reinhardtieae", "Sclerospermeae", "Podococceae", "Oranieae", "Manicarieae", "Leopoldinieae", "Euterpeae")
tribes<- tribes [-29]
tribes2 <- tribes[!tribes %in% tribes_overlap]  
list_tribes <- list()

for (i in tribes2) {
  print (i)
  x <- subset(tax, PalmTribe == i)
  x$label <- x$PalmSubfamily
  list_tribes[[i]] <- x
}

# merge overlapping labels:
list_tribes_OL <- list()
for (i in tribes_overlap) {
  print (i)
  x <- subset(tax, PalmTribe == i)
  x$label <- x$PalmTribe
  list_tribes_OL[[i]] <- x
}

list_tribes_OL[[1]]$label <- " "
list_tribes_OL[[2]]$label <- "Phytelepheae, Cyclosphatheae"
list_tribes_OL[[3]]$label <- " "
list_tribes_OL[[4]]$label <- "Roystoneeae, Reinhardtieae"
list_tribes_OL[[5]]$label <- " "
list_tribes_OL[[6]]$label <- "POS"
list_tribes_OL[[7]]$label <- " "
list_tribes_OL[[8]]$label <- " "
list_tribes_OL[[9]]$label <- "Manicarieae, Euterpeae, Leopoldinieae"
list_tribes_OL[[10]]$label <- " "





tribes3 <- tribes[!tribes %in% c("Cyclosphatheae", "Reinhardtieae", "Sclerospermeae", "Oranieae", "Manicarieae", "Euterpeae")]



#plotTree(tree,type="fan",lwd=1,ftype="off")
pdf(file="FanPlot_Realms.pdf",  width = 23.386, height = 12.402, pointsize = 2)
  plot.contMap(sim_pc, xlim=c(-100,1.25*max(nodeHeights(tree))), ylim=c(-140, 1.25*max(nodeHeights(tree))) ,ftype="off", res=40, type="fan", oma=c(3,2,1,1))
  par(lend=3)
  
  for(i in 1:Ntip(tree)){
    cc<-if(obj$xx[i]>0) 5 else -5
    th<-atan(obj$yy[i]/obj$xx[i])
    segments(obj$xx[i],obj$yy[i],
             obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),lwd=2,
             col=col.hab[md[tree$tip.label[i]]])
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+7*cc*cos(th),
             obj$yy[i]+7*cc*sin(th),lwd=2,
             col=col.afr[afr[tree$tip.label[i]]])  
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+6*cc*cos(th),
             obj$yy[i]+6*cc*sin(th),lwd=2,
             col=col.accAfr[accAfr[tree$tip.label[i]]])
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+5*cc*cos(th),
             obj$yy[i]+5*cc*sin(th),lwd=2,
             col=col.amer[amer[tree$tip.label[i]]])
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+4*cc*cos(th),
             obj$yy[i]+4*cc*sin(th),lwd=2,
             col=col.asi[asi[tree$tip.label[i]]])
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+3*cc*cos(th),
             obj$yy[i]+3*cc*sin(th),lwd=2,
             col=col.pac[pac[tree$tip.label[i]]])
    segments(obj$xx[i]+cc*cos(th),
             obj$yy[i]+cc*sin(th),
             obj$xx[i]+2*cc*cos(th),
             obj$yy[i]+2*cc*sin(th),lwd=2,
             col=col.mad[mad[tree$tip.label[i]]])
  }
  
  legend("topleft",leg.txt, 
         c("moist", "dry", "madagascar", "pacific",  "eurasia", "americas", "africa",  "accAfrica"),
         #c(levels(md),levels(mad), levels(pac), levels(asi), levels(amer), levels(afr), levels(accAfr)),
         pch=15,col=c(col.hab,col.mad[2], col.pac[2], col.asi[2], col.amer[2], col.afr[2], col.accAfr[2]),
         pt.cex=1.5,cex=0.8,bty="n")
  arc.cladelabels(tree=tree, "Arecoidea", node=fmrca(Arecoidea$x, tree), orientation="horizontal", cex=1.3, lwd=3, colors="black", ln.offset=1.65, lab.offset = 1.5)
  arc.cladelabels(tree=tree, "Coryphoideae", node=fmrca(Coryphoideae$x, tree), orientation="horizontal", cex=1.3, lwd=3, col="black", ln.offset=1.65, lab.offset = 1.5)
  arc.cladelabels(tree=tree, "Ceroxyloideae", node=fmrca(Ceroxyloideae$x, tree), orientation="horizontal", cex=1.3, lwd=3, col="darkgrey", ln.offset=1.65, lab.offset = 1.5)
  arc.cladelabels(tree=tree, "Calamoideae", node=fmrca(Calamoideae$x, tree), orientation="horizontal", cex=1.3, lwd=3, col="darkgrey", ln.offset=1.65, lab.offset=1.5)
  #dev.off()
  
  dev.off()
  

pdf(file="FanPlot_Realms.pdf", width=10, height=8)



plot.contMap(sim_pc, ftype="off", res=40, type="fan", oma=c(3,2,1,1))

#plot.contMap(sim_pc, ftype="off", res=40, type="fan", mar=c(5.1,2.1,2.1,0.1), oma=c(0,1,1,1))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
n<-Ntip(tree)



par(lend=3)

for(i in 1:Ntip(tree)){
  cc<-if(obj$xx[i]>0) 5 else -5
  th<-atan(obj$yy[i]/obj$xx[i])
  segments(obj$xx[i],obj$yy[i],
           obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),lwd=4,
           col=col.hab[md[tree$tip.label[i]]])
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+7*cc*cos(th),
           obj$yy[i]+7*cc*sin(th),lwd=4,
           col=col.afr[afr[tree$tip.label[i]]])  
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+6*cc*cos(th),
           obj$yy[i]+6*cc*sin(th),lwd=4,
           col=col.accAfr[accAfr[tree$tip.label[i]]])
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+5*cc*cos(th),
           obj$yy[i]+5*cc*sin(th),lwd=4,
           col=col.amer[amer[tree$tip.label[i]]])
  
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+4*cc*cos(th),
           obj$yy[i]+4*cc*sin(th),lwd=4,
           col=col.asi[asi[tree$tip.label[i]]])
  
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+3*cc*cos(th),
           obj$yy[i]+3*cc*sin(th),lwd=4,
           col=col.pac[pac[tree$tip.label[i]]])
  
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+2*cc*cos(th),
           obj$yy[i]+2*cc*sin(th),lwd=4,
           col=col.mad[mad[tree$tip.label[i]]])
}

legend("topleft",leg.txt, 
       c("moist", "dry", "madagascar", "pacific",  "eurasia", "americas", "africa",  "accAfrica"),
       #c(levels(md),levels(mad), levels(pac), levels(asi), levels(amer), levels(afr), levels(accAfr)),
       pch=15,col=c(col.hab,col.mad[2], col.pac[2], col.asi[2], col.amer[2], col.afr[2], col.accAfr[2]),
       pt.cex=1.5,cex=0.8,bty="n")

dev.off()






















####
moist <- data2 %>% filter(!is.na(AverageFruitLength_cm), moist0_dry1=="0") %>% group_by(accRealm) %>% dplyr::summarize(n=n())
moist <- as.data.frame(moist)
moist$moist <- moist$n

sum2 <- cbind(sum, moist$moist)

dry <- data2 %>% filter(!is.na(AverageFruitLength_cm), moist0_dry1!="0") %>% group_by(accRealm) %>% dplyr::summarize(n=n())
dry <- as.data.frame(dry)
dry$dry <- dry$n

sum3 <- cbind(sum2, dry$dry)


africa <- africa
moistdry <- moistdry
log_avgl
avgl 
BBM_sim

moistdry <- ReorderData(tree, moistdry, taxa.names="names")
africa <- ReorderData(tree, africa, taxa.names="names")
accAfrica <- ReorderData(tree, accAfrica, taxa.names="names")
accRealm <- ReorderData(tree, accRealm, taxa.names="names")
avgl <- ReorderData(tree, avgl, taxa.names="names")
log_avgl <- ReorderData(tree, log_avgl, taxa.names="names")
BBM_sim <- ReorderData(tree, BBM_sim, taxa.names="names")

X <- data.frame(accRealm = accRealm, accAfrica = accAfrica, moistdry = moistdry, sim_log = as.data.frame(BBM_sim), obs_log = as.data.frame(log_avgl), obs = as.data.frame(avgl))
X2 <- X[1:3]

####






## works
BBM_obj <- contMap(tree, BBM_sim, plot=F, res=200)
BBM_obj<-setMap(BBM_obj,c("#440154FF","#21908CFF","#FDE725FF"))
plot(BBM_obj, leg.txt="Simulated fruit length (log, cm)", type="fan", ftype="off")


obs_obj <- contMap(tree, log_avgl, plot=F, res=200)
obs_obj<-setMap(obs_obj,c("#440154FF","#21908CFF","#FDE725FF"))

cols5 <- setNames(c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"), c("AFRICA", "AMERICAS", "PACIFIC", "MADAGASCAR", "EURASIA"))
cols2 <- setNames(c("#440154FF", "#FDE725FF"), c("1_moist", "2_dry"))

md<-as.factor(setNames(X2[,1],rownames(X2)))
md <- as.factor(setNames(X2[,3], rownames(X2)))
afr <- as.factor(setNames(X2[,2], rownames(X2)))


plotTree(tree,type="fan",lwd=1,ftype="off")


obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
n<-Ntip(tree)
col.hab<-setNames(c("darkgreen","grey"),levels(md))
col.afr<-setNames(cols5,levels(afr))
par(lend=3)
for(i in 1:Ntip(tree)){
  cc<-if(obj$xx[i]>0) 5 else -5
  th<-atan(obj$yy[i]/obj$xx[i])
  segments(obj$xx[i],obj$yy[i],
           obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           lwd=4,
           col=col.hab[md[tree$tip.label[i]]])
  segments(obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           obj$xx[i]+2*cc*cos(th),
           obj$yy[i]+2*cc*sin(th),lwd=4,
           col=col.afr[afr[tree$tip.label[i]]])
}
legend("topleft",c(levels(md),levels(afr)),
       pch=15,col=c(col.hab,col.afr),
       pt.cex=1.5,cex=0.8,bty="n")















###

X<-X[tree$tip.label,]
plotTree(tree,plot=FALSE)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.3),
         ftype="off")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
h<-max(obj$xx)
fsize<-0.6
for(i in 1:Ntip(tree)){
  lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
  text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=4)
}
s<-max(fsize*strwidth(tree$tip.label))
start.x<-1.05*h+s
palettes<-c("OrRd","PuOr","RdYlBu","Spectral")
cols<-list()
for(i in 1:ncol(X)){
  text(start.x,max(obj$yy)+1,paste("trait",colnames(X)[i]),pos=4,srt=60,
       cex=0.8,offset=0)
  cols[[i]]<-setNames(sample(brewer.pal(max(3,length(levels(X[[i]]))),
                                        palettes[i]),length(levels(X[[i]]))),levels(X[[i]]))
  for(j in 1:nrow(X)){
    xy<-c(start.x,obj$yy[j])
    y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
    asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
      par()$pin[2]/par()$pin[1]
    x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
    polygon(x,y,col=cols[[i]][as.character(X[[i]][j])])
  }
  start.x<-start.x+2*asp
}
start.y<-max(obj$yy)
for(i in 1:ncol(X)){
  text(start.x,start.y,paste("trait",colnames(X)[i]),pos=4,cex=0.9,
       offset=0)
  add.simmap.legend(colors=cols[[i]],shape="square",prompt=FALSE,
                    x=start.x,y=start.y-2*strheight("W")*0.9,fsize=0.9)
  start.y<-start.y-1.5*0.9*strheight("W")*(length(cols[[i]])-1)-6
}
