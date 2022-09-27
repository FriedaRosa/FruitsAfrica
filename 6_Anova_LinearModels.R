## Directories ========================
rm(list = ls())
setwd("~/share/groups/ea/Frieda/Spatial_analyses_May22/Spatial_analyses_RServer/")
main.dir <- "ProjectAfricaMegafaunaOnly"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src/")
fig.dir <- file.path(main.dir, "figures/")
options(stringsAsFactors =FALSE)

## Packages ========================
library(plyr); library(dplyr); library(reshape2); library(EnvStats);
library(ggplot2); library(rstatix); library(EnvStats); library(ggpubr)
library(performance); library(sjPlot)

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

## Data import ========================
tdwg_final3 <- read.csv(file.path(res.dir, "tdwg_final.csv"), header=T, sep=",")
str(tdwg_final3)
tdwg_final3$accAfrica <- factor(tdwg_final3$accAfrica, levels = c("0", "1"))


### World map shape file ============

library(rgeos); library(stringr)

shp <- readOGR(dsn = "./shp/", layer = "TDWG_level3_Coordinates")
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

# log1p transform
x <- log(tdwg_final_0$mean_abs_BSchange+1)
plotNormalDensity(x, main = "log1p-transformed")





## 1. ANOVA / T-tests  ========================
## 1.1 T-Test Africa / elsewhere
dd <- tdwg_final_0

dd <- dd[complete.cases(dd$max95FL_palms),]
dd$accAfrica <- factor(dd$accAfrica)
dd$accRealm <- as.factor(dd$accRealm)

dd$sqrtFL <- sqrt(dd$max95FL_palms)
dd$sqrt_BBM_mean <- sqrt(exp(dd$log_max95FL_palms_BBMmean))
dd$sqrt_change <- sqrt(dd$mean_abs_BSchange)

### 1.1a Empirical ==============

stat.test <- dd %>% wilcox_test(sqrtFL ~ accAfrica)
stat.test # p < 0.0001
dd %>%  wilcox_effsize(sqrtFL ~ accAfrica) # large

my_colors <- setNames(c("#bdbdbd", "blue"), c("0", "1"))
dd$accAfrica <- factor(dd$accAfrica, levels = c("1", "0"))

FL_violin <- ggpubr::ggviolin(dd, x = "accAfrica", y = "sqrtFL", fill="accAfrica", color ="#636363", show.legend=F, trim = T) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  #scale_y_sqrt()+
  scale_x_discrete(labels = c(
      "1" = "Africa",
      "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y=expression(sqrt("max 95-percentile fruit length [cm]")), x=NULL)+
  stat_n_text(y.pos = 0.5) + 
  theme(legend.position="none")+
  theme_classic()+
  ylim(0, 5)+ theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=11.5), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 12,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title

pdf(file=paste0(fig.dir, "Anova_emp.pdf"), width = 6, height = 4)
#FL_violin <- FL_violin+ labs(subtitle = get_test_label(stat.test, detailed = T, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -7, hjust=1)) +coord_fixed(ratio=0.4) #  Add p-value
FL_violin <- FL_violin + labs(subtitle = expression(paste("Africa:"~beta~"= 0.24, se = 0.11, p = 0.036"))) + theme(legend.position="none", plot.subtitle = element_text(vjust = -7, hjust=0.1)) 
FL_violin
dev.off()

# Fruit size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FL_palms")
tdwg_maxFS <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FL_palms", na.rm=T)
q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_maxFS$value, probs = q_probs, na.rm = T, type = 8))
max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxFS$maxFS_q <- cut(tdwg_maxFS$value, q_values, include.lowest = T)
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")


maxFS_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white",
               data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT,
                 fill = maxFS_q, size = maxFS_q),
             data =tdwg_maxFS, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]",
                                          override.aes = list(size = seq(1,3, length.out = 8)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = seq(1,3, length.out = 8))+
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "right",
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) 
map_vio_FL <- FL_violin +maxFS_p
map_vio_FL
ggsave(file.path(fig.dir, "3_palms_maxFS_p.pdf"), map_vio_FL, width = 14, height = 10)


### 1.1b Simulated fruit size =========

stat.test <- dd %>% wilcox_test(sqrt_BBM_mean ~ accAfrica)
stat.test # n.s.
dd %>%  wilcox_effsize(sqrt_BBM_mean ~ accAfrica) # small

sim_violin <- ggpubr::ggviolin(dd, x = "accAfrica", y = "sqrt_BBM_mean", fill="accAfrica", color ="#636363", show.legend=F, trim = T) +
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
#  Add p-value
pdf(file=paste0(fig.dir, "Anova_sim.pdf"), width = 6, height = 4)
#sim_violin <- sim_violin + labs(subtitle = get_test_label(stat.test, detailed = TRUE, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.25)) +coord_fixed(ratio=3.25)   #  Add p-value
sim_violin <- sim_violin + labs(subtitle = expression(paste("Africa:"~beta~"= -0.04, se = 0.11, p = n.s."))) + theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015)) +coord_fixed(ratio=3.25)   #  Add p-value
sim_violin
dev.off()


# Simulated fruit size mapped  ========================
dd$max95_BBM <- exp(dd$log_max95FL_palms_BBMmean)

target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95_BBM")
tdwg_maxBBM <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95_BBM", na.rm=T)
q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_maxBBM$value, probs = q_probs, na.rm = T, type = 8))
max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxBBM$maxFS_q <- cut(tdwg_maxBBM$value, q_values, include.lowest = T)
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxBBM$maxFS_q) <- gsub(levels(tdwg_maxBBM$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFSsim_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white", data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT, fill = maxFS_q, size = maxFS_q), data =tdwg_maxBBM, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]", override.aes = list(size = seq(1,3, length.out = 8)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = seq(1,3, length.out = 8))+
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "right",
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) 
map_vio_FL_sim <- sim_violin +maxFSsim_p
map_vio_FL_sim
ggsave(file.path(fig.dir, "3_palms_maxFS_BBM_p.pdf"), map_vio_FL_sim, width = 14, height = 10)


### 1.1d Mean Body size change =========

## T-Test Africa / elsewhere

dd <- dd[complete.cases(dd$sqrt_change),]

stat.test <- dd %>% wilcox_test(sqrt_change ~ accAfrica, detailed=T)
stat.test # n.s.
dd %>%  wilcox_effsize(sqrt_change ~ accAfrica) # small

Bschange_violin <- ggpubr::ggviolin(dd, x = "accAfrica", y = "sqrt_change", fill="accAfrica", color ="#636363", show.legend=F, trim =T) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(labels = c(
    "1" = "Africa",
    "0" = "Elsewhere"))+
  labs(y=expression(sqrt("Mean decrease in frugivore body size [kg]")), x=NULL)+
  stat_n_text() + 
  theme(legend.position="none")+
  theme_classic()+
  theme(
    axis.text=element_text(size=13, color = "black"), 
    axis.title=element_text(size=11.5), 
    plot.title=element_text(size=20),
    text = element_text(size = 13), 
    legend.text=element_text(size=20), 
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  
                         r = 12,  
                         b = 12,  
                         l = 12))  
#  Add p-value
pdf(file=paste0(fig.dir, "Anova_absBSchange_binary.pdf"), width = 6, height = 4)
Bschange_violin <- Bschange_violin + labs(subtitle = get_test_label(stat.test, detailed = TRUE, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.15))  #  Add p-value
Bschange_violin + coord_fixed(ratio=0.05)
dev.off()


# Mean body size decreases mapped =================
maxBSchange_zCuts <- quantile(dd$mean_abs_BSchange, probs = seq(0,1,length.out = 6), na.rm = TRUE)
dd$maxBSchange_categ <-  cut(dd$mean_abs_BSchange, breaks = maxBSchange_zCuts, include.lowest = T, ordered_result = T)
maxBSchange_zmidpt <- vector()
for(i in 1:(length(maxBSchange_zCuts)-1)){ maxBSchange_zmidpt[i] <- (maxBSchange_zCuts[i] + maxBSchange_zCuts[i+1])/2 }
levels(dd$maxBSchange_categ) <- cleanCuts(levels(dd$maxBSchange_categ))


abs_decr_BS_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = shp.df2) +
  geom_point(aes(y = LAT, x= LONG, size = maxBSchange_categ, fill = maxBSchange_categ), data = dd, pch = 21, colour = "grey80") + 
  guides(size = "none", fill = guide_legend(title = "Mean decrease \nin mammalian frugivore \nbody size [kg]", override.aes = list(size = seq(1,3, length.out = 5))))+
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  theme(legend.position = 'right',
        legend.key = element_blank(),
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])

map_BSchange_abs<-Bschange_violin+ coord_fixed(ratio=0.05)+ abs_decr_BS_map 
ggsave(file.path(fig.dir, "3_palms_meanBSchange_map.pdf"), map_BSchange_abs, width = 14, height = 10, device = "pdf")



### Percentage body size change ===============

dd <- dd[complete.cases(dd$mean_perc_BSchange),]

stat.test <- dd %>% wilcox_test(mean_perc_BSchange ~ accAfrica, detailed=T)
stat.test # p < 0.0001
dd %>%  wilcox_effsize(mean_perc_BSchange ~ accAfrica) # moderate

percBSchange <- ggpubr::ggviolin(dd, x = "accAfrica", y = "mean_perc_BSchange", fill="accAfrica", color ="#636363", show.legend=F, trim = T) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(labels = c(
    "1" = "Africa",
    "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y="Decrease in mammalian frugivore body size [%]", x=NULL)+
  stat_n_text(y.pos = -5) + 
  theme(legend.position="none")+
  theme_classic()+
  ylim(-5, 105)+ 
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
#  Add p-value
pdf(file=paste0(fig.dir, "Anova_percBSchange_binary.pdf"), width = 6, height = 4)
perc_Bschange_violin <- percBSchange + labs(subtitle = get_test_label(stat.test, detailed = TRUE, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=.2)) #  Add p-value
perc_Bschange_violin
dev.off()

# Percentage body size change map  =================
maxBSchange_zCuts <- quantile(dd$mean_perc_BSchange, probs = seq(0,1,length.out = 6), na.rm = TRUE)
dd$maxBSchange_categ <-  cut(dd$mean_perc_BSchange, breaks = maxBSchange_zCuts, include.lowest = T, ordered_result = T)
maxBSchange_zmidpt <- vector()
for(i in 1:(length(maxBSchange_zCuts)-1)){ maxBSchange_zmidpt[i] <- (maxBSchange_zCuts[i] + maxBSchange_zCuts[i+1])/2 }
levels(dd$maxBSchange_categ) <- cleanCuts(levels(dd$maxBSchange_categ))

perc_decr_BS_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = shp.df2) +
  geom_point(aes(y = LAT, x= LONG, size = maxBSchange_categ, fill = maxBSchange_categ), data = dd, pch = 21, colour = "grey80") + 
  guides(size = "none",
         fill = guide_legend(title = "Mean relative decrease\nin mammalian frugivore \nbody size [%]",
                             override.aes = list(size = seq(1,3, length.out = 5))))+
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  theme(legend.position = 'right',
        legend.key = element_blank(),
        text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])

perc_decr_BS_map 
map_BSchange_perc <- perc_Bschange_violin + coord_fixed(ratio=0.015) + perc_decr_BS_map 
ggsave(file.path(fig.dir, "3_palms_mean_percBSchange_map.pdf"), map_BSchange_perc, width = 14, height = 10, device = "pdf")








## Some Test plots (not for manuscript) ========================
my_colors <- setNames(c("blue", "#bdbdbd"), c("1", "0"))

# a) max FL ~ sqrt Change

x1 <- ggplot(tdwg_final_0, aes(x=(sqrt_change), y=sqrt(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  labs(colour = NULL, 
       x = "(Sqrt) Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# b) max FL ~ raw change
x2 <- ggplot(tdwg_final_0, aes(x=(mean_abs_BSchange), y=(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +  
  labs(colour = NULL, 
       x = "Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# c) max FL ~ sqrt max BS
x3 <- ggplot(tdwg_final_0, aes(x=(sqrt(curr_max95BodySize)), y=sqrt(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa'))  +
  labs(colour = NULL, 
       x = "Maximum 95-percentile current frugivore body mass (sqrt, g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

grid<-x1+ x2+x3 
grid
library(rcompanion)

##  transform data and check plots  ========================
pdf(file = paste0(fig.dir, "Raw_data_skewness.pdf"))
par(mfrow=c(4,2))
plotNormalDensity(tdwg_final3$max95FL_palms, main = "raw \n max95 FL")
plotNormalDensity(tdwg_final3$curr_max95BodySize, main ="max95 BS")
plotNormalDensity(tdwg_final3$MAT, main ="MAT")
plotNormalDensity(tdwg_final3$TempSeas, main ="TempSeas")
plotNormalDensity(tdwg_final3$AnnualPrec, main ="AnnualPrec") 
plotNormalDensity(tdwg_final3$PrecSeas, main ="PrecSeas")
plotNormalDensity(tdwg_final3$CH_PercCover, main ="PercCover")
plotNormalDensity(tdwg_final3$CH_Mean, main ="CH mean")
dev.off()

pdf(file = paste0(fig.dir, "log_data_skewness.pdf"))
par(mfrow=c(4,2))
plotNormalDensity(log(tdwg_final3$max95FL_palms), main = "log \n max95 FL")
plotNormalDensity(log(tdwg_final3$curr_max95BodySize), main ="max95 BS")
plotNormalDensity(log(tdwg_final3$MAT), main ="MAT")
plotNormalDensity(log(tdwg_final3$TempSeas), main ="TempSeas")
plotNormalDensity(log(tdwg_final3$AnnualPrec), main ="AnnualPrec")
plotNormalDensity(log(tdwg_final3$PrecSeas), main ="PrecSeas")
plotNormalDensity(log(tdwg_final3$CH_PercCover), main ="PercCover")
plotNormalDensity(log(tdwg_final3$CH_Mean), main ="CH mean")
dev.off()

pdf(file = paste0(fig.dir, "Sqrt_data_skewness.pdf"))
par(mfrow=c(4,2))
plotNormalDensity(sqrt(tdwg_final3$max95FL_palms), main = "sqrt \n max95 FL")
plotNormalDensity(sqrt(tdwg_final3$curr_max95BodySize), main ="max95 BS")
plotNormalDensity(((tdwg_final3$MAT)^2), main ="MAT")
plotNormalDensity(sqrt(tdwg_final3$TempSeas), main ="TempSeas")
plotNormalDensity(sqrt(tdwg_final3$AnnualPrec), main ="AnnualPrec")
plotNormalDensity(sqrt(tdwg_final3$PrecSeas), main ="PrecSeas")
plotNormalDensity(sqrt(tdwg_final3$CH_PercCover), main ="PercCover")
plotNormalDensity(sqrt(tdwg_final3$CH_Mean), main ="CH mean")
dev.off()

plotNormalDensity(sqrt(exp(tdwg_final3$log_max95FL_palms_BBMmean))) #?
plotNormalDensity(sqrt(exp(tdwg_final3$log_medianFL_palms_BBM_MCC))) #?

### Observed fruit sizes: perfect Layout ======================================================================

p0ax <- ggplot(tdwg_final_0, aes(factor(accAfrica), sqrt(max95FL_palms))) +
  geom_violin(aes(fill=accAfrica, col=accAfrica)) +
  scale_fill_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere')) +
  scale_color_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere'))+
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
  scale_x_discrete(limits = c("1", "0"), labels = c("1" = "Africa", "0" = "Elsewhere")) +
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=NULL, color = "Distribution", fill = "Distribution")

p0ax

FL_violin

p1a <- ggplot(tdwg_final_0, aes(sqrt_change, sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.08, se = 0.08, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=expression(sqrt("Mean decrease in  body mass [g]")), color = NULL)

p2a <- ggplot(tdwg_final_0, aes(sqrt(curr_max95BodySize), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
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
  theme(legend.position="none")+
  labs(subtitle = expression(paste(beta~"= 0.28, se = 0.1, p = 0.01"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Current frugivore body mass [g]")), color = NULL)

## do not use in main text: 
p3a <- ggplot(tdwg_final_0, aes((MAT)^2, sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.02, se = 0.12, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression("Mean annual temperature"^2), color = NULL)

p4a <-ggplot(tdwg_final_0, aes(sqrt(TempSeas), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
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
  labs(subtitle = expression(paste(beta~"= -0.32, se = 0.1, p = 0.002"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Temperature seasonality")), color = NULL)

p5a <- ggplot(tdwg_final_0, aes(sqrt(AnnualPrec), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.23, se = 0.11, p = 0.04"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Annual precipitation")), color = NULL)

p6a <- ggplot(tdwg_final_0, aes(sqrt(PrecSeas), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
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
  labs(subtitle = expression(paste(beta~"= 0.17, se = 0.14, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Precipitation seasonality")), color = NULL)

p7a <- ggplot(tdwg_final_0, aes(sqrt(CH_Mean), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, aes(col = accAfrica))+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  #labs(subtitle = expression(paste(beta~"= -0.15, se = 0.12, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Canopy height")), color = NULL)


pdf(file=paste0(fig.dir, "allPlotsLMsFormulasObserved.pdf"), height = 12, width = 12.75 )
legend <- get_legend(p0ax + theme(legend.box.margin = margin(0, -10, 0, 0)))
x <- cowplot::plot_grid(FL_violin + coord_fixed(ratio=0.2), p1a, p2a, p3a, p4a, p5a, p6a, p7a, legend, hjust = -1, align = "hv", ncol= 3,axis = c("lb"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

### Perfect layout! =================================================================================================

## Simulated fruit sizes: perfect layout!
p0bx <- ggplot(tdwg_final_0, aes(factor(accAfrica), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_violin(aes(fill=accAfrica, col=accAfrica)) +
  scale_fill_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere')) +
  scale_color_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere'))+
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
  scale_x_discrete(limits = c("1", "0"), labels = c("1" = "Africa", "0" = "Elsewhere")) +
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x="Distribution", color = "Distribution", fill = "Distribution")

p0b <- ggplot(tdwg_final_0, aes(factor(accAfrica), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_violin(aes(fill=accAfrica, col=accAfrica), show.legend = F) +
  scale_fill_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere')) +
  scale_color_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere'))+
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
  theme(legend.position="none")+
  scale_x_discrete(limits = c("1", "0"), labels = c("1" = "Africa", "0" = "Elsewhere")) +
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x="Distribution", color = "Distribution", fill = "Distribution")

p1b <- ggplot(tdwg_final_0, aes(sqrt_change, sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.30, se = 0.08, p = 0.0004"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=expression(sqrt("Mean decrease in frugivore body mass [g]")), color = "Distribution")

p2b <- ggplot(tdwg_final_0, aes(sqrt(curr_max95BodySize), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
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
  theme(legend.position="none")+
  labs(subtitle = expression(paste(beta~"= 0.11, se = 0.1, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Current frugivore body mass [g]")), color = "Distribution")


## do not use in main text: 
p3b <- ggplot(tdwg_final_0, aes((MAT)^2, sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.04, se = 0.12, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression("Mean annual temperature"^2), color = "Distribution")


p4b <-ggplot(tdwg_final_0, aes(sqrt(TempSeas), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
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
  labs(subtitle = expression(paste(beta~"= -0.25, se = 0.1, p = 0.01"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Temperature seasonality")), color = "Distribution")



p5b <- ggplot(tdwg_final_0, aes(sqrt(AnnualPrec), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= 0.29, se = 0.11, p = 0.009"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Annual precipitation")), color = "Distribution")

p6b <- ggplot(tdwg_final_0, aes(sqrt(PrecSeas), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
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
  labs(subtitle = expression(paste(beta~"= -0.02, se = 0.14, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Precipitation seasonality")), color = "Distribution")

p7b <- ggplot(tdwg_final_0, aes(sqrt(CH_Mean), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  #stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
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
  labs(subtitle = expression(paste(beta~"= -0.05, se = 0.12, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Canopy height")), color = "Distribution")

pdf(file=paste0(fig.dir, "allPlotsLMsFormulasSimulated.pdf"), height = 12, width = 12.75 )

legend <- get_legend(p0bx + theme(legend.box.margin = margin(0, -10, 0, 0)))
x <- cowplot::plot_grid(sim_violin, p1b, p2b, p3b, p4b, p5b, p6b, p7b, legend, hjust = -1, align = c("hv"), ncol= 3,axis = c("lb"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

### Perfect layout! ======

# Save relevant subplots for maintext:
pdf(file=paste0(fig.dir, "finalPlotLMsFormulas.pdf"), width = 16, height = 5)
x<- cowplot::plot_grid(FL_violin, p2a, p1a, legend, nrow=1, labels=c('a)','b)', 'c)'), align = "hv",axis = c("lbt"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

pdf(file=paste0(fig.dir, "finalPlotLMsFormulasSimulated.pdf"), width = 16, height = 5)
x <- cowplot::plot_grid(sim_violin, p2b, p1b, legend, nrow=1, labels=c('a)','b)','c)'), align="hv",axis = c("lbt"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

### Perfect layout! ======




## Preparations for linear models ========================
full_df <- tdwg_final_0
Afr <- tdwg_final_0 %>% filter(accAfrica == "1")
Glob <- tdwg_final_0 %>% filter(accAfrica == "0")

# Scale:
normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
}

cols <- c("#636363", "blue")


# 1) Full model: (w/o interactions)

lm_full1 <- lm(normalized(sqrt(max95FL_palms)) ~ factor(accAfrica)+ 
                 normalized((sqrt_change)) + 
                 normalized(sqrt(curr_max95BodySize)) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
                      data = full_df, na.action = na.exclude)
summary(lm_full1)
plot(normalized(sqrt(max95FL_palms)) ~ 
       normalized((sqrt_change)), data = full_df)

plot_model(lm_full1, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Maximum fruit length \n(Palms) w/o interactions", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease", "Distribution in Africa"))

# 2) Full model: (with interactions)

lm_full2 <- lm(normalized(sqrt(max95FL_palms)) ~ 
                 normalized((sqrt_change))*factor(accAfrica) + 
                 normalized(sqrt(curr_max95BodySize))*factor(accAfrica) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
              data = full_df, na.action = na.exclude)
summary(lm_full2)

# models have different DFs, therefore cannot be compared using R-squared
# here I use the performance package:

compare_performance(lm_full1, lm_full1, lm_full2, rank = TRUE, metrics = "common")
# model 2 has a better fit (with interactions) 


pdf(file=paste0(fig.dir, "Effecs_lms.pdf"), width = 6, height = 4)
pLM <- plot_model(lm_full2, 
           colors = c("#440154FF","#21908CFF"), 
           type= "est", 
           show.intercept = F, 
           show.values = T, 
           value.offset = .3,
           show.p = T, 
           vline.color = "grey", 
           order.terms = c(3, 1,2, 10, 9, 4, 5, 6, 7, 8),
           title = "",
           axis.labels=c("Canopy Height", "Precipitation Seasonality", "Annual Precipitation", "Temperature Seasonality", "Mean annual temperature", 
                         "Body mass decrease * Africa", "Current body mass * Africa", 
                         "Distribution in Africa", "Body mass decrease", "Current body mass")) + theme_classic() +   theme(text=element_text(size=13, color = "black"), #change font size of all text
                                                                                                                           axis.text=element_text(size=13), #change font size of axis text
                                                                                                                           axis.title=element_text(size=13)) 
pLM <- pLM +  theme(axis.text=element_text(size=11.5, color = "black"), #change font size of axis text
             axis.title=element_text(size=13), #change font size of axis titles
             plot.title=element_text(size=20),
             text = element_text(size = 13), #change font size of plot title
             legend.text=element_text(size=20), #change font size of legend text
             legend.title=element_text(size=20),
             plot.margin = margin(t = 12,  # Top margin
                                  r = 12,  # Right margin
                                  b = 12,  # Bottom margin
                                  l = 12)) # Left margin) #change font size of legend title 
pLM
dev.off()

library(cowplot)
pdf(file=paste0(fig.dir, "finalPlotFormulas2.pdf"), width = 14, height = 10)
cowplot::plot_grid(pLM, FL_violin, p2a + theme(plot.margin = margin(t = 12,  # Top margin
                                                                    r = 12,  # Right margin
                                                                    b = 12,  # Bottom margin
                                                                    l = 100)), p1a + theme(plot.margin = margin(t = 12,  # Top margin
                                                                                                                r = 12,  # Right margin
                                                                                                                b = 12,  # Bottom margin
                                                                                                                l = 80)), nrow = 2, ncol=2, labels=c('a)','b)','c)','d)'), align = "hv", title = "Maximum 95-percentile fruit length", axis="l")
dev.off()


# 2b) Simulated full model with interactions:

lm_full1b <- lm(normalized(sqrt(exp(log_max95FL_palms_BBMmean))) ~ factor(accAfrica) +
                 normalized((sqrt_change)) + 
                 normalized(sqrt(curr_max95BodySize)) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
               data = full_df, na.action = na.exclude)
summary(lm_full1b)
plot(normalized(sqrt(exp(log_max95FL_palms_BBMmean))) ~ 
       normalized((sqrt_change)), data = full_df)

plot_model(lm_full1b, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Maximum fruit length \n(Palms) w/o interactions", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease", "Distribution in Africa"))


lm_full2b <- lm(normalized(sqrt(exp(log_max95FL_palms_BBMmean))) ~ 
                 normalized((sqrt_change))*factor(accAfrica) + 
                 normalized(sqrt(curr_max95BodySize))*factor(accAfrica) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
               data = full_df, na.action = na.exclude)
summary(lm_full2b)
compare_performance(lm_full1b, lm_full2b, rank = TRUE, metrics = "common")

cols <- c("#636363", "blue")
pdf(file=paste0(fig.dir, "Simulated_Effecs_lms.pdf"), width = 6, height = 4)
pLMsim <- plot_model(lm_full2b, 
                  colors = c("#440154FF","#21908CFF"), 
                  type= "est", 
                  show.intercept = F, 
                  show.values = T, 
                  value.offset = .3,
                  show.p = T, 
                  vline.color = "grey", 
                  order.terms = c(3, 1,2, 10, 9, 4, 5, 6, 7, 8),
                  title = "",
                  axis.labels=c("Canopy Height", "Precipitation Seasonality", "Annual Precipitation", "Temperature Seasonality", "Mean annual temperature", 
                                "Body mass decrease * Africa", "Current body mass * Africa", 
                                "Distribution in Africa", "Body mass decrease", "Current body mass")) + theme_classic() +   theme(text=element_text(size=13, color = "black"), #change font size of all text
                                                                                                                                  axis.text=element_text(size=13), #change font size of axis text
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
dev.off()

library(cowplot)
pdf(file=paste0(fig.dir, "sim_finalPlotFormulas2.pdf"), width = 14, height = 10)
cowplot::plot_grid(pLMsim, sim_violin, p2b + theme(plot.margin = margin(t = 12,  # Top margin
                                                                    r = 12,  # Right margin
                                                                    b = 12,  # Bottom margin
                                                                    l = 100)), p1b + theme(plot.margin = margin(t = 12,  # Top margin
                                                                                                                r = 12,  # Right margin
                                                                                                                b = 12,  # Bottom margin
                                                                                                                l = 80)), nrow = 2, ncol=2, labels=c('a)','b)','c)','d)'), align = "hv", title = "Maximum 95-percentile fruit length", axis="l")
dev.off()



## Not for maintext ==========================


# 3) Africa without interactions

lm_afr1 <- lm(normalized(sqrt(max95FL_palms)) ~ 
                 normalized((sqrt_change)) + 
                 normalized(sqrt(curr_max95BodySize)) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(CH_Mean), 
               data = Afr)
summary(lm_afr1)
p3 <- plot_model(lm_afr1, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Africa: Maximum fruit length \n(Palms)", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease"))

# 4) Elswhere without interactions

lm_wrld1 <- lm(normalized(sqrt(max95FL_palms)) ~ 
                normalized((sqrt_change)) + 
                normalized(sqrt(curr_max95BodySize)) + 
                normalized(I(MAT^2)) +
                normalized(sqrt(TempSeas)) + 
                normalized(sqrt(AnnualPrec)) +
                normalized(sqrt(PrecSeas)) +
                normalized(CH_Mean), 
              data = Glob)
summary(lm_wrld1)
p4 <- plot_model(lm_wrld1, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Elswhere: Maximum fruit length \n(Palms)", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease"))

p1+p2+p3+p4


### Correlation graph ============================


c<-tdwg_final_0 %>% select_if(is.numeric) %>% na.omit()
c$X <- NULL
c$log_medianFL_palms_BBM_MCC <- NULL
c$PercentagCo

#saveRDS(c, "corr.rds")
### problems with correlation package; the following was run on the local computer:

# library(correlation)
# c <- readRDS("../Downloads/corr.rds")
# c$CH_PercCover <- NULL
# #c$presNat_max95BodySize <- NULL
# 
# cc <- correlation::correlation(c, corr.method = "spearman", na.deletion = "pairwise", sort.corr =T, show.p = T, decimals = 2)
# print(cc)
# summary(cc)
# write.csv(cc, "LinuxBackup/FruitSizeEvo_Project/OnlyAfrica/corr_tab.csv")
# 
# ####
# untransformed data has no linear relationship between parameters
# that's why we use Spearman rank correlation



cc_mat <- corr.test(c)$p # save p-values as matrix
cc <- cor(c, use="pairwise.complete.obs") # calculate correlations
symnum(cc)



library(psych)
corPlot(cc, cex = 1.2, pval=T, p.mat = cc_mat, insign = "blank") # Canopy Cover ~ Annual Precipitation = 0.78


library(corrplot); library(RColorBrewer)
corrplot(cc, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"), insig = "blank", method="circle", p.mat = cc_mat)

corrplot(cc, order = 'hclust', addrect = 2, insig = "pch", method = "color", type="")

