## New script: plots from OUwie models ~ MCC tree =================================

# 1. Line plots (y = AICweight x = data, color = model)
## 1.1 best fitting models (1. filter for delta = 0, 2. calculate average aicw, 3. plot)
## 1.2 all models (1. group by model, 2. calculate average aicw, 3. plot)

# 2. Line plots (y = theta, x = regime, color = data)
## 2.1 best fitting models (1. filter for delta = 0, 2. group by model (exclude BM1, BMS, OU1?), 3. plot MCC results as means)
## 2.2 all models (1. group by model, 2. plot MCC results as means)

# 2. Line plots (y = sigma, x = regime, color = data)
## 2.1 best fitting models (1. filter for delta = 0, 2. filter for OUMV, OUMVA, 3. plot MCC results as means)
## 2.2 all models (1. group by model, 2. filter for OUMV, OUMVA, 3. plot MCC results as means)

# Libraries  =======================
library(tidyr); library(dplyr); library(scales); library(viridis); library(viridisLite); library(ggplot2)
setwd(dir = "~/GitHub/FruitsAfrica/Sensitivity_Extinctions/")
# =======================
rm(list = ls())

# =======================

# 100 trees results:
Afr_emp <- read.table("raw_avgl_Afr.txt", header=T)
Afr_BBM <- read.table("raw_BM_Afr.txt", header=T)

# negative eigenvalues, odd estimates, negative AICc (convergence problems)
mods3 <- Afr_emp[!is.na(Afr_emp$eigval1) & Afr_emp$eigval1 >= 0 & !is.na(Afr_emp$eigval2) & Afr_emp$eigval2 >= 0 & !is.na(Afr_emp$eigval3) & Afr_emp$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights
fit<- mods4 %>% group_by(id) %>% geiger::aicw(mods4$AICc)
mods5<-cbind(mods4, fit)

avgl_processed <- mods5

####

mods3 <- Afr_BBM[!is.na(Afr_BBM$eigval1) & Afr_BBM$eigval1 >= 0 & !is.na(Afr_BBM$eigval2) & Afr_BBM$eigval2 >= 0 & !is.na(Afr_BBM$eigval3) & Afr_BBM$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights
fit<-geiger::aicw(mods4$AICc)
mods5<-cbind(mods4, fit)

BBM_processed <- mods5

####

head(avgl_processed)
str(avgl_processed)
str(Afr_emp)

# Convergence Problems:
unique(Afr_emp$model)
unique(avgl_processed$model)

issues<-setdiff(Afr_emp,avgl_processed[,1:18])
issues

write.table(issues, file = "issues_AVGL.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)

# BMS   OU1   OUM  OUMA  OUMV OUMVA 
# 1     1     4   100     7    81

sum(table(issues$model)) # 194 in total

# write to .txt
write.table(avgl_processed, file = "processed_AVGL.txt", sep = "\t", row.names = F, col.names = T)

####

head(BBM_processed)
str(BBM_processed)
str(Afr_BBM)

# Convergence Problems:
unique(Afr_BBM$model)
unique(BBM_processed$model)

issues<-setdiff(Afr_BBM,BBM_processed[,1:18])
issues

write.table(issues, file = "issues_BM.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)

# OUM  OUMA  OUMV OUMVA 
# 13   101    14    27 

sum(table(issues$model)) # 155

# write to .txt
write.table(BBM_processed, file = "processed_BM.txt", sep = "\t", row.names = F, col.names = T)

#####

Afr_merged <- read.csv("dd_Africa_all.csv") # best fit

# transform parameters to long-format
## megafauna:
data_long_Afr <- gather(Afr_merged, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long_Afr)
head(data_long_Afr)

# round values
data_long_Afr$value <- round(as.numeric(data_long_Afr$value),3)
data_long_Afr$w <- round(data_long_Afr$w, 2)
data_long_Afr$delta <- round(data_long_Afr$delta, 2)


# merge to masterfile ============================================================
dd_master <- data_long_Afr
head(dd_master)
summary(dd_master)

str(dd_master)
dd_master$data <- as.factor(dd_master$data)
dd_master$regime <- as.factor(dd_master$regime)
dd_master$tree <- as.factor(dd_master$tree)

#write.table(dd_master, file = "master_file_long_emp_sim.txt", sep = "\t", row.names = F, col.names = T)
# ================================================================================

# dataframes grouped by variables:
## megafauna
data_long_Afr_sep <- gather(Afr_merged, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_logThetaSE, v_logThetaSE, c(theta_se0, theta_se_1), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_sigma, v_sigma, c(sigma0, X1_sigma.sq), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_alpha, v_alpha, c(X0_alpha, X1_alpha), factor_key=TRUE)
head(data_long_Afr_sep)


# merge to parameter seperated masterfile ============================================================
dd_master_sep <- data_long_Afr_sep
head(dd_master_sep)
summary(dd_master_sep)

str(dd_master_sep)
dd_master_sep$data <- as.factor(dd_master_sep$data)
dd_master_sep$regime <- as.factor(dd_master_sep$regime)
dd_master_sep$tree <- as.factor(dd_master_sep$tree)

#write.table(dd_master_sep, file = "master_file_long_sep_emp_sim.txt", sep = "\t", row.names = F, col.names = T)
# ================================================================================

dd_master_bestfit_sep <- subset(dd_master_sep, delta==0)
head(dd_master_bestfit_sep)
multi_theta <- dd_master_bestfit_sep %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"), data == "empirical")
multi_sigma <- dd_master_bestfit_sep %>% filter(model %in% c("OUMV", "OUMVA", "BMS"), data == "empirical")

afr_multi_theta <- multi_theta %>% filter(regime == "Megafauna")
afr_multi_theta_MCC <- afr_multi_theta %>% filter(tree == "MCC")

afr_multi_sigma <- multi_sigma %>% filter(regime == "Megafauna")
afr_multi_sigma_MCC <- afr_multi_sigma %>% filter(tree == "MCC")

# =================================


# ================================================================================
MCC <- MCC[-c(19:21)]
MCC_emp <- MCC %>% filter(data == "empirical")
MCC_sim <- MCC %>% filter(data == "simulated")

fit_emp<-geiger::aicw(MCC_emp$AICc)
MCC_emp<-cbind(MCC_emp, fit_emp)

fit_sim<-geiger::aicw(MCC_sim$AICc)
MCC_sim<-cbind(MCC_sim, fit_sim)

MCC <- rbind(MCC_emp, MCC_sim)

all_merged_wide <- all_merged_wide %>% filter(tree == "set_trees" | delta==0)

all_best <- merge(all_merged_wide, MCC, all=T)

multi_theta <- all_best %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
multi_sigma <- all_best %>% filter(model %in% c("OUMV", "OUMVA"))


## Fig 1: proportion of selected models ===============
dd_Afr <- Afr_merged  
group.colors <- c( OUM = "#cdcdcd", OUMV = "#cdcdcd", OUMVA = "#cdcdcd")

empty_models <-c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")

aicw_emp <- dd_Afr %>%
  filter(delta==0 & data == "empirical" & tree == "set_trees") %>% group_by(model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw_emp$mean <- round(aicw_emp$mean,2)
aicw_emp$y <- c(3, 82, 18, 7)


aicw_sim <- dd_Afr %>%
  filter(delta==0 & data == "simulated" & tree == "set_trees") %>% group_by(model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw_sim$mean <- round(aicw_sim$mean,2)
aicw_sim
aicw_sim$y <- c(65, 25, 9, 12)



# =============================================
library("ggrepel")       
p2_1<-dd_Afr %>%
  filter(delta==0 & data == "empirical") %>%
  ggplot(aes( x=factor(model)))+
  geom_bar(width=.7, position="dodge2")+
  scale_x_discrete(limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
  theme(text= element_text(color="black"),
        title=element_text(size=15, color = "black"),
        axis.text=element_text(size=11.5, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15, color="black"), 
        axis.title.x = element_text(size=15, color="black"),
        axis.text.y = element_text(size=15, color="black"),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+ 
  #geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 6)+
  geom_text(data = aicw_emp, aes(x=model, y=y, label = paste("w =", mean), vjust=0.3)) +  
  ylim(0, 100)+
  theme_classic()+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of selected models",
       title ="a) Observed fruit sizes ")
p2_1

p2_2<-dd_Afr %>%
  filter(delta==0 & data == "simulated") %>%
  ggplot(aes( x=factor(model)))+
  geom_bar(width=.7, position="dodge2")+
  scale_x_discrete(limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  #geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
  geom_text(data = aicw_sim, aes(x=model, y=y, label = paste("w =", mean), vjust=0.3)) +  
  
  ylim(0, 100)+
  theme_classic()+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of selected models",
       title ="b) Simulated fruit sizes")
p2_2

## Theta dot-line plot:  ==============================


group.colors <- c(theta0 = "#636363", theta_1 = "blue")
afr_multi_theta$g_logTheta <- factor(afr_multi_theta$g_logTheta, levels = c("theta_1", "theta0"))
afr_multi_theta <- afr_multi_theta %>% filter(v_logTheta > 0.3)

p3 <- ggplot(data=afr_multi_theta, aes(y = v_logTheta, x = g_logTheta), show.legend=F)+
  geom_point(aes(color=g_logTheta), cex = 2)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   labels=c("theta0" = "Elsewhere", "theta_1" = "Africa"))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(theta~"log(Optimum fruit size) [cm]"))+
  scale_color_manual("Parameters", labels = c("Theta0", "Theta1"), values = group.colors)+
  geom_point(data = afr_multi_theta_MCC, aes(y = v_logTheta, x = g_logTheta), cex =2)+
  geom_line(data= afr_multi_theta_MCC, aes(group = id), col="black")+
  labs(subtitle = expression(bold(Observed)~n[OUM]~'= 80,'~n[OUMV]~'= 16,'~n[OUMVA]~'= 5'~~'\r\n'))+
  theme(legend.position="none")+
  theme( plot.subtitle= element_text(size=9),
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

p3
### =======================================
group.colors <- c(X1_sigma.sq = "blue",sigma0 = "#636363")
afr_multi_theta$g_sigma <- factor(afr_multi_theta$g_sigma, levels = c("X1_sigma.sq", "sigma0"))
table(afr_multi_theta$model)


p5 <- ggplot(data=afr_multi_sigma, aes(y = v_sigma, x = g_sigma), show.legend=F)+
  geom_point(aes(color=g_sigma), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_sigma.sq", "sigma0"),
                   labels=c("X1_sigma.sq" = "Africa", "sigma0" = "Elsewhere"))+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(sigma^2~"Evolutionary rate"))+
  scale_color_manual("Parameters", labels = c("Africa", "Elsewhere"), values = group.colors)+
  geom_point(data = afr_multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  geom_line(data= afr_multi_sigma_MCC, aes(group = id), col="black")+
  labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 5'~~'\r\n'))+
  theme(legend.position="none")+
  theme( plot.subtitle= element_text(size=9),
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


p5


########

## Models with convergence problems ##

setwd(dir = "~/GitHub/FruitsAfrica/Sensitivity_Extinctions/")

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
       title ="b) Simulated fruit sizes ")+
  theme_classic()

