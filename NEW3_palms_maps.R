## New script to cosntruct maps based on at least 3 palms ======== august 2022


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
library(plyr); library(dplyr); library(reshape2)
library(sp); library(rgdal)
library(car); library(MuMIn)
library(readxl)
library(picante)


## Data import ========================
tdwg_final3 <- read.csv(file.path(res.dir, "tdwg_final.csv"), header=T, sep=",")
str(tdwg_final3)

## Data handling  ========================
tdwg_final3$LEVEL_3_CO <- as.factor(tdwg_final3$LEVEL_3_CO)
tdwg_final3$accRealm <- as.factor(tdwg_final3$accRealm)

# Set negative values (i.e., positive body mass changes) to zero before log-transformation:
tdwg_final_0 <- tdwg_final3
tdwg_final_0$mean_abs_BSchange[tdwg_final_0$mean_abs_BSchange < 0] <- 0



tdwg_mammals_summary_3palms <- tdwg_final_0

######### ============================
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

tdwg_final_3palms <- tdwg_final3
tdwg_final_3palms <- tdwg_final_3palms %>% filter_at(vars( mean_abs_BSchange), any_vars(!is.na(.)))


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

# Absolute change:

# ANOVA:

library(viridis); library(rstatix); library(dplyr); library(ggpubr); library(ggplot2); library(tidyr); library(EnvStats)
cols <- "light grey"

tdwg_final_3palms$accRealm[tdwg_final_3palms$LEVEL_3_CO == "MDG"] <- "Madagascar"
tdwg_final_3palms$accRealm <- as.factor(tdwg_final_3palms$accRealm)

tdwg_final_3palms %>%
  group_by(accRealm) %>%
  get_summary_stats(mean_abs_BSchange, type="common")

### 

# ANOVA:

tdwg_final <- tdwg_final3
boxplot(sqrt(mean_abs_BSchange) ~ accRealm, data = tdwg_final)
res.kruskal <- tdwg_final %>% kruskal_test(mean_abs_BSchange ~ accRealm)
(res.kruskal)

# Dunn's test (post-hoc)

pwc <- tdwg_final %>% 
  dunn_test(mean_abs_BSchange ~ accRealm, p.adjust.method = "bonferroni", detailed=T) 
pwc #

tdwg_final %>% kruskal_effsize((mean_abs_BSchange) ~ accRealm)

y <- seq(from=1100, to=1500, length.out=10)
pwc <- pwc %>% mutate(y.position = y)

cols <- "light grey"
figIa <- ggboxplot(tdwg_final, x = "accRealm", y = "mean_abs_BSchange", fill=cols) +
  stat_pvalue_manual(pwc, hide.ns = T, tip.length = 0.01, y.position = 782, step.increase = 0.07) +  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(
    "Realm",
    labels = c(
      "Africa" = "Africa",
      "Americas" = "Americas",
      "Pacific" = "Pacific realm",
      "Asia" = "Eurasia",
      "Madagascar" = "Madagascar")) +
  labs(
    #title = "a) Observed Fruit width (Angiosperms)", 
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc),
    y="Absolute decrease \nin frugivore body size [kg]", x="Realm")+
  geom_abline(intercept=0, alpha=0.3, col="red")+
  stat_n_text()+
  theme_classic()+
  ylim(-200, 1000)+ 
  theme(
    axis.text=element_text(size=12), #change font size of axis text
    axis.title=element_text(size=18), #change font size of axis titles
    plot.title=element_text(size=20), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    text = element_text(size = 19),
    plot.margin = margin(t = 20,  # Top margin
                         r = 160,  # Right margin
                         b = 20,  # Bottom margin
                         l = 10)) # Left margin) #change font size of legend title   
figIa

abs_decr_BS_p <- figIa
ggsave(figIa, filename = paste(fig.dir, "BSdecrease_realm.pdf"), height = 6, width = 8.5, bg = "transparent")



# MAP:
maxBSchange_zCuts <- quantile(tdwg_final_3palms$mean_abs_BSchange/1000, probs = seq(0,1,length.out = 6), na.rm = TRUE)
tdwg_final_3palms$maxBSchange_categ <-  cut(tdwg_final_3palms$mean_abs_BSchange/1000, breaks = maxBSchange_zCuts, include.lowest = T, ordered_result = T)
maxBSchange_zmidpt <- vector()
for(i in 1:(length(maxBSchange_zCuts)-1)){ maxBSchange_zmidpt[i] <- (maxBSchange_zCuts[i] + maxBSchange_zCuts[i+1])/2 }
levels(tdwg_final_3palms$maxBSchange_categ) <- cleanCuts(levels(tdwg_final_3palms$maxBSchange_categ))


abs_decr_BS_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = shp.df2) +
  geom_point(aes(y = LAT, x= LONG, size = maxBSchange_categ, fill = maxBSchange_categ), data = tdwg_final_3palms, pch = 21, colour = "grey80") + 
  guides(size = "none",
         fill = guide_legend(title.position = "top",
                             title.hjust = 0.5,
                             title = "Mean decrease in mammalian frugivore body size [kg]",
                             override.aes = list(size = seq(1,3, length.out = 5))))+
  coord_fixed() +
  map_theme +
  #ggtitle("Absolute decreases in frugivore body size [kg]") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box = 'vertical',
        legend.key = element_blank(),
        text = element_text(size = 20),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 15)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +#rescale(maxFSchange_zmidpt, to = c(1,5))) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
abs_decr_BS_map 

ggsave(file.path(fig.dir, "3_palms_meanBSchange_map.pdf"), abs_decr_BS_map, width = 9, height = 12, device = "pdf")

### fruit sizes mapped ==============

# Fruit size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FL_palms")
tdwg_maxFS <- melt(tdwg_final_3palms[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FL_palms", na.rm=T)
q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_maxFS$value, probs = q_probs, na.rm = T, type = 8))
max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxFS$maxFS_q <- cut(tdwg_maxFS$value, q_values, include.lowest = T)
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")
xx <- sqrt(max_q_mdpt)


maxFS_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white",
               data = shp.df2) +
  geom_point(aes(x = LONG, y = LAT,
                 fill = maxFS_q, size = maxFS_q),
             data =tdwg_maxFS, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]",
                                         override.aes = list(size=sqrt(max_q_mdpt)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = sqrt(max_q_mdpt)) +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "right",
        text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) 
  #labs(caption="TDWG3 units reduced to those with > 3 palms")
maxFS_p
ggsave(file.path(fig.dir, "3_palms_maxFS_p.pdf"), maxFS_p, width = 9, height = 10)

## Simulated: Species averages from 100 trees
# fruit length mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "log_max95FL_palms_BBMmean")
tdwg_BBMmean <- melt(tdwg_final_3palms[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "log_max95FL_palms_BBMmean", na.rm=T)
tdwg_BBMmean$value <- round(exp(tdwg_BBMmean$value),2)

q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_BBMmean$value, probs = q_probs, na.rm = T, type = 8))
max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_BBMmean$maxFS_q <- cut(tdwg_BBMmean$value, q_values, include.lowest = T)
levels(tdwg_BBMmean$maxFS_q) <- gsub(levels(tdwg_BBMmean$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_BBMmean$maxFS_q) <- gsub(levels(tdwg_BBMmean$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")


BBMmean_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white",
               data = shp.df) +
  geom_point(aes(x = LONG, y = LAT,
                 fill = maxFS_q, size = maxFS_q),
             data =tdwg_BBMmean, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length [cm]",
                                          override.aes = list(size=xx))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = xx) +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "right",
        text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) 
  #labs(caption="TDWG3 units reduced to those with > 3 palms") 
BBMmean_p

ggsave(file.path(fig.dir, "3_palms_BBMmean_p.pdf"), BBMmean_p, width = 9, height = 10)

