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
setwd(dir = "Data/OUwie/out/")
# =======================
rm(list = ls())

# =======================

# 100 trees results:
Afr_emp <- read.table("merged_data_AVGL.txt", header=T)
Afr_BBM <- read.table("merged_data_BM.txt", header=T)

# MCC results:
MCC <- read.table("merged_MCC_all.txt", header=T)

# make data id column (empirical, simulated)
Afr_emp$data <- "empirical"
Afr_BBM$data <- "simulated"
MCC <- MCC %>% mutate(data = ifelse(id %in% c("ouwie_Afr_BM", "ouwie_md_BM"), "simulated", "empirical"))

unique(MCC$id)
# make regime id column (Africa/Elsewhere, Moist/Dry)
Afr_emp$regime <- "Megafauna"
Afr_BBM$regime <- "Megafauna"
MCC <- MCC %>% mutate(regime = ifelse(id %in% c("ouwie_Afr_BM", "ouwie_Afr_avgl"), "Megafauna", "Vegetation"))

# make MCC/set id column (MCC, set_trees)
Afr_emp$tree <- "set_trees"
Afr_BBM$tree <- "set_trees"
MCC$tree <- "MCC"

# merge different regimes seperately
MCC_afr<-MCC %>% filter(regime == "Megafauna")
Afr_merged <- merge(Afr_emp, Afr_BBM, all=T)
Afr_merged <- merge(Afr_merged, MCC_afr, all=T)

all_merged_wide <- Afr_merged %>% filter(regime == "Megafauna")

#
#write.table(all_merged_wide, file = "all_merged_wide.txt", sep = "\t", row.names = F, col.names = T)
#


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
MCC_best <- MCC %>% filter(delta == 0) 
all_best <- all_merged_wide %>% filter(delta==0)

multi_theta <- all_best %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
multi_sigma <- all_best %>% filter(model %in% c("OUMV", "OUMVA"))


## Fig 1: proportion of selected models ===============
dd_Afr <- Afr_merged  
group.colors <- c( OUM = "#cdcdcd", OUMV = "#cdcdcd", OUMVA = "#cdcdcd")

empty_models <-c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")

aicw <- dd_Afr %>%
  filter(delta==0 & data == "empirical" & tree == "set_trees") %>% group_by(model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw
aicw <- dd_Afr %>%
  filter(delta==0 & data == "simulated" & tree == "set_trees") %>% group_by(model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw

aicw <- all_merged_wide %>% filter(delta==0 & tree == "MCC") %>% group_by(data, regime, model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw

# Make summary table for ouwie: =====================

sum_tab <- all_merged_wide %>% filter(delta==0) %>% group_by(tree, data, model) %>% dplyr::summarize(
  n = n(),
  alpha0 = mean(X0_alpha) , 
  alpha1 = mean(X1_alpha),
  sigma0 = mean(sigma0),
  sigma1 = mean(X1_sigma.sq),
  theta0 = mean(theta0),
  theta1 = mean(theta_1),
  aicw = mean(w)
)
write.csv(sum_tab, "../../../Results/OUwie_SumTab_AfrOnly.csv")

# =============================================

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
  geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
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
  geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
  ylim(0, 100)+
  theme_classic()+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of selected models",
       title ="b) Simulated fruit sizes")
p2_2



pdf("Fig1_Count_selectedModels_grid.pdf",width=5, height=6)
gridExtra::grid.arrange(p2_1,p2_2)
dev.off()

## Theta dot-line plot:  ==============================


group.colors <- c(theta0 = "#636363", theta_1 = "blue")
afr_multi_theta$g_logTheta <- factor(afr_multi_theta$g_logTheta, levels = c("theta_1", "theta0"))

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
  labs(subtitle = expression(bold(Observed)~n[OUM]~'= 82,'~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
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

### =======================================
group.colors <- c(X1_sigma.sq = "blue",sigma0 = "#636363")
afr_multi_theta$g_sigma <- factor(afr_multi_theta$g_sigma, levels = c("X1_sigma.sq", "sigma0"))


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
  labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
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



gridExtra::grid.arrange(p3,p5)
