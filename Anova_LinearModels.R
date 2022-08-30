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

## Data import ========================
tdwg_final3 <- read.csv(file.path(res.dir, "tdwg_final.csv"), header=T, sep=",")
str(tdwg_final3)
tdwg_final3$accRealm[tdwg_final3$LEVEL_3_CO == "MDG"] <- "Madagascar"

## Data handling  ========================
tdwg_final3$LEVEL_3_CO <- as.factor(tdwg_final3$LEVEL_3_CO)
tdwg_final3$accAfrica <- as.factor(tdwg_final3$accAfrica)
plotNormalDensity(tdwg_final3$mean_abs_BSchange)

# Set negative values (i.e., positive body mass changes) to zero before log-transformation:
tdwg_final_0 <- tdwg_final3
tdwg_final_0$mean_abs_BSchange[tdwg_final_0$mean_abs_BSchange < 0] <- 0
plotNormalDensity(tdwg_final_0$mean_abs_BSchange, main = "untransformed")

# sqrt transform
tdwg_final_0$sqrt_change <- sqrt(tdwg_final_0$mean_abs_BSchange)
plotNormalDensity(tdwg_final_0$sqrt_change, main = "sqrt transformed")

x <- log(tdwg_final_0$mean_abs_BSchange+1)
plotNormalDensity(x, main = "log1p-transformed")

## 1. ANOVA / T-tests  ========================
## 1.1 T-Test Africa / elsewhere
dd <- tdwg_final_0

dd <- dd[complete.cases(dd$max95FL_palms),]
dd$accAfrica <- as.factor(dd$accAfrica)
dd$accRealm <- as.factor(dd$accRealm)

dd$sqrtFL <- sqrt(dd$max95FL_palms)
dd$sqrt_BBM_mean <- sqrt(exp(dd$log_max95FL_palms_BBMmean))
dd$sqrt_change <- sqrt(dd$mean_abs_BSchange)


### 1.1a Empirical ==============
my_colors <- setNames(c("#bdbdbd", "blue"), c("0", "1"))

dd$accAfrica <- factor(dd$accAfrica, levels = c("1", "0"))


pdf(file=paste0(fig.dir, "Anova_emp.pdf"), width = 6, height = 4)
anova <- ggpubr::ggboxplot(dd, x = "accAfrica", y = "sqrtFL", fill="accAfrica", color ="#636363", show.legend=F) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  #scale_y_sqrt()+
  scale_x_discrete(
    "Distribution",labels = c(
      "1" = "Africa",
      "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y=expression(sqrt("max 95-percentile fruit length [cm]")), x="Realm")+
  stat_n_text() + 
  theme(legend.position="none")+
  theme_classic()+
  ylim(0, 4)+ theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=13), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 12,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title 
#  Add p-value
anova + ggpubr::stat_compare_means(hjust = 0.75, show.legend=F) +theme(legend.position="none") +stat_compare_means(aes(label = paste0(", ", ..p.signif..)), hjust = -1)
dev.off()


### 1.1b Simulated =========

## T-Test Africa / elsewhere

pdf(file=paste0(fig.dir, "Anova_sim.pdf"), width = 6, height = 4)
anova <- ggpubr::ggboxplot(dd, x = "accAfrica", y = "sqrt_BBM_mean", fill="accAfrica", color ="#636363", show.legend=F) +
  scale_fill_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  stat_summary(fun="mean", color="white", shape=20) +
  scale_x_discrete(
    "Distribution", labels = c(
      "1" = "Africa",
      "0" = "Elsewhere"))+
  labs(
    #title = "c) Observed Fruit length [cm] (Palms)", 
    y=expression(sqrt("max 95-percentile fruit length [cm]")), x="Realm")+
  stat_n_text() + 
  theme(legend.position="none")+
  theme_classic()+
  ylim(1.4, 1.9)+ theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=13), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 12,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title   
#  Add p-value
anova + ggpubr::stat_compare_means(hjust = 0.75, show.legend=F) +theme(legend.position="none") +stat_compare_means(aes(label = paste0(", ", ..p.signif..)), hjust = -1.25)
dev.off()

## 1.2 ANOVA per biogeographic region

cols <- setNames(c("blue", "#bdbdbd","#bdbdbd","#bdbdbd","#bdbdbd"), c("Africa", "Americas", "Pacific", "Asia", "Madagascar"))
dd$accRealm <- factor(dd$accRealm , levels = c("Africa", "Americas", "Pacific", "Asia", "Madagascar"))
dd2 <- dd[complete.cases(dd$sqrt_change),]

res.kruskal <- dd2 %>% kruskal_test(sqrt_change ~ accRealm)
res.kruskal

# Dunn's test (post-hoc)
pwc <- dd2 %>% 
  dunn_test(sqrt_change ~ accRealm, p.adjust.method = "bonferroni", detailed=T) 
pwc #

# Effect size
dd2 %>% kruskal_effsize((sqrt_change) ~ accRealm)

pdf(file=paste0(fig.dir, "Anova_BSchange.pdf"), width = 6, height = 4)
ggboxplot(dd2, x = "accRealm", y = "sqrt_change", fill="accRealm", color ="#636363", show.legend=F) +
  stat_pvalue_manual(pwc, hide.ns = T, tip.length = 0.01, y.position = 1050, step.increase = 0.07) +  stat_summary(fun="mean", color="white", bg = "black",shape=20) +
  scale_fill_manual("Realms", values = cols)+
  scale_x_discrete(
    "Realm",
    labels = c(
      "Africa" = "Africa",
      "Americas" = "Americas",
      "Pacific" = "Pacific",
      "Asia" = "Asia",
      "Madagascar" = "Madagascar")) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = T),
    caption = get_pwc_label(pwc),
    y=expression(sqrt("Decrease in  frugivore body mass [g]")), x="Realm")+
  #geom_abline(intercept=0, alpha=0.3, col="red")+
  stat_n_text()+
  theme_classic()+
  ylim(-100, 1300)+ 
  theme(legend.position="none")+
  theme(
    axis.text=element_text(size=13, color = "black"), #change font size of axis text
    axis.title=element_text(size=13), #change font size of axis titles
    plot.title=element_text(size=20),
    text = element_text(size = 13), #change font size of plot title
    legend.text=element_text(size=20), #change font size of legend text
    legend.title=element_text(size=20),
    plot.margin = margin(t = 12,  # Top margin
                         r = 12,  # Right margin
                         b = 12,  # Bottom margin
                         l = 12)) # Left margin) #change font size of legend title 
dev.off()

## Some Test plots (not for manuscript) ========================
my_colors <- setNames(c("blue", "#bdbdbd"), c("1", "0"))

# a) max FL ~ sqrt Change
par(mfrow=c(2,2))
ggplot(tdwg_final_0, aes(x=(sqrt_change), y=sqrt(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Africa', 'Elsewhere')) +
  labs(colour = "Region", 
       x = "(Sqrt) Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# b) max FL ~ raw change
ggplot(tdwg_final_0, aes(x=(mean_abs_BSchange), y=(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Elswhere', 'Africa')) +
  labs(colour = "Region", 
       x = "Decrease in  body mass (g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

# c) max FL ~ sqrt max BS
ggplot(tdwg_final_0, aes(x=(sqrt(curr_max95BodySize)), y=sqrt(max95FL_palms)))+  
  theme_classic()+
  geom_smooth(method="auto", colour="#636363", alpha = 0.1) +
  geom_point(aes(color=accAfrica))+ 
  scale_color_manual(values=my_colors, labels=c('Africa', 'Elsewhere')) +
  labs(colour = "Region", 
       x = "Maximum 95-percentile current frugivore body mass (sqrt, g) \nper botanical country",
       y = "Maximum 95-percentile fruit length (cm)")

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
        axis.title=element_text(size=13), #change font size of axis titles
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

p0a <- ggplot(tdwg_final_0, aes(factor(accAfrica), sqrt(max95FL_palms))) +
  geom_violin(aes(fill=accAfrica, col=accAfrica), show.legend = F) +
  scale_fill_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere')) +
  scale_color_manual(values=my_colors, limits = c("1", "0"), labels=c('Africa', 'Elsewhere'))+
  theme_classic()+  
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
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

p1a <- ggplot(tdwg_final_0, aes(sqrt_change, sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=expression(sqrt("Mean decrease in  body mass [g]")), color = "Distribution")

p2a <- ggplot(tdwg_final_0, aes(sqrt(curr_max95BodySize), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
  theme_classic()+ 
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  theme(legend.position="none")+
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Current frugivore body mass [g]")), color = "Distribution")

## do not use in main text: 
p3a <- ggplot(tdwg_final_0, aes((MAT)^2, sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression("Mean annual temperature"^2), color = "Distribution")

p4a <-ggplot(tdwg_final_0, aes(sqrt(TempSeas), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Temperature seasonality")), color = "Distribution")

p5a <- ggplot(tdwg_final_0, aes(sqrt(AnnualPrec), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Annual precipitation")), color = "Distribution")

p6a <- ggplot(tdwg_final_0, aes(sqrt(PrecSeas), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 4.4) +
  stat_regline_equation(label.x = 3, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Precipitation seasonality")), color = "Distribution")

p7a <- ggplot(tdwg_final_0, aes(sqrt(CH_Mean), sqrt(max95FL_palms))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 1.2, label.y = 4.4) +
  stat_regline_equation(label.x = 1.2, label.y = 4.2)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Canopy height")), color = "Distribution")


pdf(file=paste0(fig.dir, "allPlotsLMsFormulasObserved.pdf"), height = 12, width = 12.75 )
legend <- get_legend(p0ax + theme(legend.box.margin = margin(0, -10, 0, 0)))
x <- cowplot::plot_grid(p0a, p1a, p2a, p3a, p4a, p5a, p6a, p7a, legend, hjust = -1, align = "h", ncol= 3,axis = c("lr"))
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
        axis.title=element_text(size=13), #change font size of axis titles
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
        axis.title=element_text(size=13), #change font size of axis titles
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
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x=expression(sqrt("Mean decrease in frugivore body mass [g]")), color = "Distribution")

p2b <- ggplot(tdwg_final_0, aes(sqrt(curr_max95BodySize), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors,limits = c("1", "0"), labels=c("0" ='Elsewhere', "1" = 'Africa')) +
  theme_classic()+ 
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  theme(legend.position="none")+
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Current frugivore body mass [g]")), color = "Distribution")


## do not use in main text: 
p3b <- ggplot(tdwg_final_0, aes((MAT)^2, sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression("Mean annual temperature"^2), color = "Distribution")


p4b <-ggplot(tdwg_final_0, aes(sqrt(TempSeas), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Temperature seasonality")), color = "Distribution")



p5b <- ggplot(tdwg_final_0, aes(sqrt(AnnualPrec), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Annual precipitation")), color = "Distribution")

p6b <- ggplot(tdwg_final_0, aes(sqrt(PrecSeas), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme(legend.position="none")+
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 3, label.y = 2.4) +
  stat_regline_equation(label.x = 3, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Precipitation seasonality")), color = "Distribution")

p7b <- ggplot(tdwg_final_0, aes(sqrt(CH_Mean), sqrt(exp(log_max95FL_palms_BBMmean)))) +
  geom_point(aes(color=accAfrica), show.legend = F) +
  stat_smooth(method = lm, col="#636363")+
  scale_color_manual(values=my_colors, labels=c('Elsewhere', 'Africa')) +
  theme_classic()+
  theme(axis.text=element_text(size=13, color = "black"), #change font size of axis text
        axis.title=element_text(size=13), #change font size of axis titles
        plot.title=element_text(size=20),
        text = element_text(size = 13), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20),
        plot.margin = margin(t = 12,  # Top margin
                             r = 12,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12))+ # Left margin) #change font size of legend title 
  stat_cor(label.x = 1.2, label.y = 2.4) +
  stat_regline_equation(label.x = 1.2, label.y = 2.32)+
  labs(y=expression(sqrt("Maximum 95-percentile fruit length [cm]")), x = expression(sqrt("Canopy height")), color = "Distribution")

pdf(file=paste0(fig.dir, "allPlotsLMsFormulasSimulated.pdf"), height = 12, width = 12.75 )

legend <- get_legend(p0bx + theme(legend.box.margin = margin(0, -10, 0, 0)))
x <- cowplot::plot_grid(p0b, p1b, p2b, p3b, p4b, p5b, p6b, p7b, legend, hjust = -1, align = "h", ncol= 3,axis = c("lr"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

### Perfect layout! ======

# Save relevant subplots for maintext:
pdf(file=paste0(fig.dir, "finalPlotLMsFormulas.pdf"), width = 16, height = 5)
x<- cowplot::plot_grid(p0a, p2a, p1a, legend, nrow=1, labels=c('a)','b)', 'c)'), align = "v",axis = c("lr"))
x +  theme(plot.margin = margin(12, 12, 12, 12))
dev.off()

pdf(file=paste0(fig.dir, "finalPlotLMsFormulasSimulated.pdf"), width = 16, height = 5)
x <- cowplot::plot_grid(p0b, p2b, p1b, legend, nrow=1, labels=c('a)','b)','c)'), align="v",axis = c("lr"))
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

lm_full1 <- lm(normalized(sqrt(max95FL_palms)) ~ 
                 normalized((sqrt_change)) + 
                 normalized(sqrt(curr_max95BodySize)) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
                      data = full_df)
summary(lm_full1)
plot(normalized(sqrt(max95FL_palms)) ~ 
       normalized((sqrt_change)), data = full_df)

plot_model(lm_full1, col= cols, type= "est", show.intercept = F, show.values = T, show.p = T, title = "Maximum fruit length \n(Palms) w/o interactions", axis.labels = c("Canopy height", "PrecSeas", "Annual Prec", "TempSeas", "MAT", "body mass", "body mass decrease"))

# 2) Full model: (with interactions)

lm_full2 <- lm(normalized(sqrt(max95FL_palms)) ~ 
                 normalized((sqrt_change))*factor(accAfrica) + 
                 normalized(sqrt(curr_max95BodySize))*factor(accAfrica) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
              data = full_df)
summary(lm_full2)

# models have different DFs, therefore cannot be compared using R-squared
# here I use the performance package:

plot(compare_performance(lm_full1, lm_full1, lm_full2, rank = TRUE, metrics = "common"))
# model 2 has a better fit (with interactions) 


pdf(file=paste0(fig.dir, "Effecs_lms.pdf"), width = 6, height = 4)
pLM <- plot_model(lm_full2, 
           col = cols, 
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
pLM 
dev.off()

library(cowplot)
pdf(file=paste0(fig.dir, "finalPlotFormulas2.pdf"), width = 10, height = 10)
cowplot::plot_grid(NULL, p0a, p2a, p1a, nrow = 2, ncol=2, labels=c('a)','b)','c)','d)'), align = "h", title = "Maximum 95-percentile fruit length", axis="l")
dev.off()


# 2b) Simulated full model with interactions:


lm_full2b <- lm(normalized(sqrt(exp(log_max95FL_palms_BBMmean))) ~ 
                 normalized((sqrt_change))*factor(accAfrica) + 
                 normalized(sqrt(curr_max95BodySize))*factor(accAfrica) + 
                 normalized(I(MAT^2)) +
                 normalized(sqrt(TempSeas)) + 
                 normalized(sqrt(AnnualPrec)) +
                 normalized(sqrt(PrecSeas)) +
                 normalized(sqrt(CH_Mean)), 
               data = full_df)
summary(lm_full2b)


cols <- c("#636363", "blue")
pdf(file=paste0(fig.dir, "Simulated_Effecs_lms.pdf"), width = 6, height = 4)
pLMsim <- plot_model(lm_full2b, 
                  col = cols, 
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
                                "Distribution in Africa", "Body mass decrease", "Current body mass")) + theme_classic() +   theme(text=element_text(size=13), #change font size of all text
                                                                                                                                  axis.text=element_text(size=13), #change font size of axis text
                                                                                                                                  axis.title=element_text(size=13)) 
pLMsim
dev.off()

pLM + pLMsim

pdf(file=paste0(fig.dir, "finalPlotFormulas2simulated.pdf"), width = 10, height = 10)
cowplot::plot_grid(NULL, p0b, p2b, p1b, nrow = 2, ncol=2, labels=c('a)','b)','c)','d)'), title = "Simulated maximum 95-percentile fruit length")
dev.off()

library(cowplot)
pdf(file=paste0(fig.dir, "finalPlotFormulas2.pdf"), width = 10, height = 10)
cowplot::plot_grid(NULL, p0a, p2a, p1a, nrow = 2, ncol=2, labels=c('a)','b)','c)','d)'), title = "Maximum 95-percentile fruit length")
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

