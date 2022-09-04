## Models with convergence problems ##

setwd("~/GitHub/FruitsAfrica/Data/OUwie/out")

Afr_avgl <- read.table("issues_AVGL.txt", header=T)
Afr_BM <- read.table("issues_BM.txt", header=T)


Afr_avgl$regime <- "Africa"
Afr_avgl$dat <- "empirical"

Afr_BM$regime <- "Africa"
Afr_BM$dat <- "simulated"


Afr_avgl$group <- "Afr_emp"
Afr_BM$group <- "Afr_BM"


Africa <- merge(Afr_avgl, Afr_BM, all=T)
str(Africa)

###############################################################################
library(dplyr)
library(rstatix)
library(ggplot2)

Africa %>% filter(dat == "empirical") %>%
ggplot(aes(x = factor(model)))+ 
  geom_bar(width=.7, position="dodge2")+
  scale_x_discrete(NULL, limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
  theme(text= element_text(color="black"),
        title=element_text(size=15, color = "black"),
        axis.text=element_text(size=11.5, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15, color="black"), 
        axis.title.x = element_text(size=15, color="black"),
        axis.text.y = element_text(size=15, color="black"),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+ 
  geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
  ylim(0, 105)+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of convergence failures",
       title ="a) Observed fruit sizes ")+
  theme_classic()


Africa %>% filter(dat == "simulated") %>%
  ggplot(aes(x = factor(model)))+ 
  geom_bar(width=.7, position="dodge2")+
  scale_x_discrete(NULL, limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
  theme(text= element_text(color="black"),
        title=element_text(size=15, color = "black"),
        axis.text=element_text(size=11.5, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15, color="black"), 
        axis.title.x = element_text(size=15, color="black"),
        axis.text.y = element_text(size=15, color="black"),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+ 
  geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
  ylim(0, 105)+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of convergence failures",
       title ="a) Simulated fruit sizes ")+
  theme_classic()



