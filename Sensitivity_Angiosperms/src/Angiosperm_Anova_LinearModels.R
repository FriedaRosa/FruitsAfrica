## Directories ========================
rm(list = ls())
setwd("~/share/groups/ea/Frieda/Spatial_analyses_May22/Spatial_analyses_RServer/")
main.dir <- "ProjectAfricaMegafaunaOnly/Sensitivity_Angiosperms/"
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
tdwg_final3 <- read.csv(file.path(res.dir, "tdwg_final_angios.csv"), header=T, sep=",")
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

# Set negative values (i.e., positive body mass changes) to zero before log-transformation:
tdwg_final_0 <- tdwg_final3
tdwg_final_0$mean_abs_BSchange[tdwg_final_0$mean_abs_BSchange < 0] <- 0
tdwg_final_0$mean_perc_BSchange[tdwg_final_0$mean_perc_BSchange < 0] <- 0

tdwg_final_0$mean_abs_BSchange <- tdwg_final_0$mean_abs_BSchange/1000

# sqrt transform
tdwg_final_0$sqrt_change <- sqrt(tdwg_final_0$mean_abs_BSchange)

## 1. ANOVA / T-tests  ========================
## 1.1 T-Test Africa / elsewhere
dd <- tdwg_final_0

dd <- dd[complete.cases(dd$max95FW_angios),]
dd$max95FW_angios <- dd$max95FW_angios/10

dd$accAfrica <- factor(dd$accAfrica)
dd$accRealm <- as.factor(dd$accRealm)

dd$sqrtFL <- sqrt(dd$max95FW_angios)
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
    y=expression(sqrt("max 95-percentile fruit width [cm]")), x=NULL)+
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

#pdf(file=paste0(fig.dir, "Anova_emp.pdf"), width = 6, height = 4)
#FL_violin <- FL_violin+ labs(subtitle = get_test_label(stat.test, detailed = T, type="expression"))+theme(legend.position="none", plot.subtitle = element_text(vjust = -7, hjust=1)) +coord_fixed(ratio=0.4) #  Add p-value
FL_violin <- FL_violin + labs(subtitle = expression(paste("Africa:"~beta~"= 0.22, se = 0.07, p < 0.001"))) + theme(legend.position="none", plot.subtitle = element_text(vjust = -7, hjust=0.1)) 
FL_violin
dev.off()

# Fruit size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FW_angios")
tdwg_maxFS <- melt(dd[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FW_angios", na.rm=T)
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
  guides(size="none", fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit width [cm]",
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


## Some Test plots (not for manuscript) ========================
my_colors <- setNames(c("blue", "#bdbdbd"), c("1", "0"))

# a) max FL ~ sqrt Change

x1 <- ggplot(tdwg_final_0, aes(x=(sqrt_change), y=sqrt(max95FW_angios)))+  
  theme_classic()+
  #geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  labs(colour = NULL, 
       x = "(Sqrt) Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# b) max FL ~ raw change
x2 <- ggplot(tdwg_final_0, aes(x=(mean_abs_BSchange), y=(max95FW_angios)))+  
  theme_classic()+
  #geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +  
  labs(colour = NULL, 
       x = "Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# c) max FL ~ sqrt max BS
x3 <- ggplot(tdwg_final_0, aes(x=(sqrt(curr_max95BodySize)), y=sqrt(max95FW_angios)))+  
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
#pdf(file = paste0(fig.dir, "Raw_data_skewness.pdf")) #8x6
par(mfrow=c(4,2))
plotNormalDensity(tdwg_final3$max95FW_angios, main = "raw \n max95 FL")
plotNormalDensity(tdwg_final3$curr_max95BodySize, main ="max95 BS")
plotNormalDensity(tdwg_final3$MAT, main ="MAT")
plotNormalDensity(tdwg_final3$TempSeas, main ="TempSeas")
plotNormalDensity(tdwg_final3$AnnualPrec, main ="AnnualPrec") 
plotNormalDensity(tdwg_final3$PrecSeas, main ="PrecSeas")
plotNormalDensity(tdwg_final3$CH_PercCover, main ="PercCover")
plotNormalDensity(tdwg_final3$CH_Mean, main ="CH mean")
dev.off()

#pdf(file = paste0(fig.dir, "log_data_skewness.pdf"))
par(mfrow=c(4,2))
plotNormalDensity(log(tdwg_final3$max95FW_angios), main = "log \n max95 FL")
plotNormalDensity(log(tdwg_final3$curr_max95BodySize), main ="max95 BS")
plotNormalDensity(log(tdwg_final3$MAT), main ="MAT")
plotNormalDensity(log(tdwg_final3$TempSeas), main ="TempSeas")
plotNormalDensity(log(tdwg_final3$AnnualPrec), main ="AnnualPrec")
plotNormalDensity(log(tdwg_final3$PrecSeas), main ="PrecSeas")
plotNormalDensity(log(tdwg_final3$CH_PercCover), main ="PercCover")
plotNormalDensity(log(tdwg_final3$CH_Mean), main ="CH mean")
dev.off()

#pdf(file = paste0(fig.dir, "Sqrt_data_skewness.pdf"))
par(mfrow=c(4,2))
plotNormalDensity(sqrt(tdwg_final3$max95FW_angios), main = "sqrt \n max95 FL")
plotNormalDensity(sqrt(tdwg_final3$curr_max95BodySize), main ="max95 BS")
plotNormalDensity(((tdwg_final3$MAT)^2), main ="MAT")
plotNormalDensity(sqrt(tdwg_final3$TempSeas), main ="TempSeas")
plotNormalDensity(sqrt(tdwg_final3$AnnualPrec), main ="AnnualPrec")
plotNormalDensity(sqrt(tdwg_final3$PrecSeas), main ="PrecSeas")
plotNormalDensity(sqrt(tdwg_final3$CH_PercCover), main ="PercCover")
plotNormalDensity(sqrt(tdwg_final3$CH_Mean), main ="CH mean")
dev.off()


### Observed fruit sizes: perfect Layout ======================================================================

p0ax <- ggplot(tdwg_final_0, aes(factor(accAfrica), sqrt(max95FW_angios))) +
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

p1a <- ggplot(tdwg_final_0, aes(sqrt_change, sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= 0.03, se = 0.05, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=expression(sqrt("Mean decrease in  body mass [kg]")), color = NULL)

p2a <- ggplot(tdwg_final_0, aes(sqrt(curr_max95BodySize), sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= -0.11, se = 0.06, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Current frugivore body mass [kg]")), color = NULL)

## do not use in main text: 
p3a <- ggplot(tdwg_final_0, aes((MAT)^2, sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= 0.25, se = 0.07, p < 0.001"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression("Mean annual temperature"^2), color = NULL)

p4a <-ggplot(tdwg_final_0, aes(sqrt(TempSeas), sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= -0.03, se = 0.6, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Temperature seasonality")), color = NULL)

p5a <- ggplot(tdwg_final_0, aes(sqrt(AnnualPrec), sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= 0.14, se = 0.07, p = 0.05"))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Annual precipitation")), color = NULL)

p6a <- ggplot(tdwg_final_0, aes(sqrt(PrecSeas), sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= 0.04, se = 0.09, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Precipitation seasonality")), color = NULL)

p7a <- ggplot(tdwg_final_0, aes(sqrt(CH_Mean), sqrt(max95FW_angios))) +
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
  labs(subtitle = expression(paste(beta~"= -0.13, se = 0.08, p = n.s."))) + 
  theme(legend.position="none", plot.subtitle = element_text(vjust = -6, hjust=0.015))+ 
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Canopy height")), color = NULL)


pdf(file=paste0(fig.dir, "allPlotsLMsFormulasObserved.pdf"), height = 12, width = 12.75 )
legend <- get_legend(p0ax + theme(legend.box.margin = margin(0, -10, 0, 0)))
x <- cowplot::plot_grid(FL_violin + coord_fixed(ratio=0.2), p1a, p2a, p3a, p4a, p5a, p6a, p7a, legend, hjust = -1, align = "hv", ncol= 3,axis = c("lb"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
plot(legend)
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

lm_full1 <- lm(normalized(sqrt(max95FW_angios)) ~ factor(accAfrica)+ 
                 normalized((sqrt_change)) + 
                 normalized(sqrt(curr_max95BodySize)) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
                      data = full_df, na.action = na.exclude)
summary(lm_full1)
plot(normalized(sqrt(max95FW_angios)) ~ 
       normalized((sqrt_change)), data = full_df)

plot_model(lm_full1, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Maximum fruit width \n(Angiosperms) w/o interactions", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease", "Distribution in Africa"))

# 2) Full model: (with interactions)

lm_full2 <- lm(normalized(sqrt(max95FW_angios)) ~ 
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


#pdf(file=paste0(fig.dir, "Effecs_lms.pdf"), width = 6, height = 4)
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

## Not for maintext ==========================


# 3) Africa without interactions

lm_afr1 <- lm(normalized(sqrt(max95FW_angios)) ~ 
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

lm_wrld1 <- lm(normalized(sqrt(max95FW_angios)) ~ 
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

