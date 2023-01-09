## New script: plots from OUwie models from Sensitivity Analysis with Extinctions (25, 50, 75)  =================================

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
setwd(dir = "~/GitHub/FruitsAfrica/Sensitivity_Extinctions2/Data")
# =======================
rm(list = ls())

multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}  

# =======================

# 100 trees results:
Afr_25 <- read.table("merged_data_raw25.txt", header=T)
Afr_50 <- read.table("merged_data_raw50.txt", header=T)
Afr_75 <- read.table("merged_data_raw75.txt", header=T)


# negative eigenvalues, odd estimates, negative AICc (convergence problems)
mods3 <- Afr_25[!is.na(Afr_25$eigval1) & Afr_25$eigval1 >= 0 & !is.na(Afr_25$eigval2) & Afr_25$eigval2 >= 0 & !is.na(Afr_25$eigval3) & Afr_25$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights

processed_list_25 <- list()
id_list <- unique(Afr_25$id)
for (i in(1:length(id_list))){
  print(i)
  mods5 <-mods4 %>% filter(id == i)
  fit <- geiger::aicw(mods5$AICc)
  mods6<-cbind(mods5, fit)
  processed_list_25[[i]] <- mods6
}

processed_25 <- multi_join(processed_list_25, full_join)
head(processed_25)
####

# negative eigenvalues, odd estimates, negative AICc (convergence problems)
mods3 <- Afr_50[!is.na(Afr_50$eigval1) & Afr_50$eigval1 >= 0 & !is.na(Afr_50$eigval2) & Afr_50$eigval2 >= 0 & !is.na(Afr_50$eigval3) & Afr_50$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights

processed_list_50 <- list()
id_list <- unique(Afr_50$id)
for (i in(1:length(id_list))){
  print(i)
  mods5 <-mods4 %>% filter(id == i)
  fit <- geiger::aicw(mods5$AICc)
  mods6<-cbind(mods5, fit)
  processed_list_50[[i]] <- mods6
}

processed_50 <- multi_join(processed_list_50, full_join)
head(processed_50)
####

# negative eigenvalues, odd estimates, negative AICc (convergence problems)
mods3 <- Afr_75[!is.na(Afr_75$eigval1) & Afr_75$eigval1 >= 0 & !is.na(Afr_75$eigval2) & Afr_75$eigval2 >= 0 & !is.na(Afr_75$eigval3) & Afr_75$eigval3 >= 0,]
mods4 <- mods3[mods3$theta0 >= 0 & mods3$theta_1 >=0,] 
mods4 <- mods4[mods4$AICc > 0,] 

# Model comparison: calculate AICc weights

processed_list_75 <- list()
id_list <- unique(Afr_75$id)
for (i in(1:length(id_list))){
  print(i)
  mods5 <-mods4 %>% filter(id == i)
  fit <- geiger::aicw(mods5$AICc)
  mods6<-cbind(mods5, fit)
  processed_list_75[[i]] <- mods6
}

processed_75 <- multi_join(processed_list_75, full_join)
head(processed_75)


# ============================================================================================================= #


# Convergence Problems:
unique(Afr_25$model)
unique(processed_25$model)

issues<-setdiff(Afr_25,processed_25[,1:18])
issues

write.table(issues, file = "issues_25_new.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)

# OUM  OUMA  OUMV OUMVA 
# 4   101     2    83 

sum(table(issues$model)) # 190 in total

# write to .txt
write.table(processed_25, file = "processed_25.txt", sep = "\t", row.names = F, col.names = T)

####

# Convergence Problems:
unique(Afr_50$model)
unique(processed_50$model)

issues<-setdiff(Afr_50,processed_50[,1:18])
issues

write.table(issues, file = "issues_50_new.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)

# OU1   OUM  OUMA  OUMV OUMVA 
# 4     3   101     5    86 

sum(table(issues$model)) # 199 in total

# write to .txt
write.table(processed_50, file = "processed_50.txt", sep = "\t", row.names = F, col.names = T)

# Convergence Problems:
unique(Afr_75$model)
unique(processed_75$model)

issues<-setdiff(Afr_75,processed_75[,1:18])
issues

write.table(issues, file = "issues_75_new.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)

# OU1   OUM  OUMA  OUMV OUMVA 
# 5     2   101     7    82 

sum(table(issues$model)) # 197 in total

# write to .txt
write.table(processed_75, file = "processed_75.txt", sep = "\t", row.names = F, col.names = T)

# ============================================================================================================= #

# Merge together:

processed_25$threshold <- "25"
processed_50$threshold <- "50"
processed_75$threshold <- "75"

Afr_merged <- merge(processed_25, processed_50, all=T)
Afr_merged <- merge(Afr_merged, processed_75, all=T)
str(Afr_merged)
write.table(Afr_merged, file = "Sensitivity_new2023_modelfit_merged_all.txt", sep = "\t", row.names = F, col.names = T)


# Best fit models:
BestFit <- Afr_merged %>% filter(delta == 0)
write.table(BestFit, file = "BestFit_models_all_Sensitivity2023.txt", sep = "\t", row.names = F, col.names = T)



# transform parameters to long-format
data_long <- gather(Afr_merged, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long)
head(data_long)

# round values for better visualization
data_long$value <- round(as.numeric(data_long$value),3)
data_long$w <- round(data_long$w, 2)
data_long$delta <- round(data_long$delta, 2)

# ================================================================================

# dataframes grouped by variables:
data_long_sep <- gather(Afr_merged, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_logThetaSE, v_logThetaSE, c(theta_se0, theta_se_1), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_sigma, v_sigma, c(sigma0, X1_sigma.sq), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_alpha, v_alpha, c(X0_alpha, X1_alpha), factor_key=TRUE)
head(data_long_sep)

dd_master_sep <- data_long_sep

# ================================================================================

dd_master_bestfit_sep <- unique(subset(dd_master_sep, delta==0))
head(dd_master_bestfit_sep)
multi_theta <- dd_master_bestfit_sep %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
multi_sigma <- dd_master_bestfit_sep %>% filter(model %in% c("OUMV", "OUMVA", "BMS"))

multi_theta_MCC <- multi_theta %>% filter(id == 101)
multi_sigma_MCC <- multi_sigma %>% filter(id == 101) # no models with multiple sigma had the best fit for the MCC trees

# =================================


# ================================================================================
MCC <- dd_master_bestfit_sep %>% filter(id == 101)
MCC_25 <- MCC %>% filter(threshold == "25")
MCC_50 <- MCC %>% filter(threshold == "50")
MCC_75 <- MCC %>% filter(threshold == "75")

MCC <- rbind(MCC_25, MCC_50)
MCC <- rbind(MCC, MCC_75)
BestFit2 <- BestFit %>% filter(id != 101 | delta==0)

all_best <- rbind(dd_master_bestfit_sep, MCC, all=T)

str(MCC)
str(dd_master_bestfit_sep)
multi_theta <- all_best %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
multi_sigma <- all_best %>% filter(model %in% c("OUMV", "OUMVA"))


## Fig 1: proportion of selected models ===============
dd_Afr <- BestFit  
group.colors <- c( OUM = "#cdcdcd", OUMV = "#cdcdcd", OUMVA = "#cdcdcd")

empty_models <-c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")

aicw_emp <- dd_Afr %>%
  filter(delta==0) %>% group_by(threshold, model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw_emp$mean <- round(aicw_emp$mean,2)

dd_Afr %>% filter(threshold == "25") %>% count(model)
dd_Afr %>% filter(threshold == "50") %>% count(model)
dd_Afr %>% filter(threshold == "75") %>% count(model)

aicw_emp$y <- c(3, 20, 81, 5, 35, 70, 3, 60, 46)
# =============================================
library("ggrepel")       
p2_1<-dd_Afr %>%
  filter(delta==0) %>%
  ggplot(aes( x=factor(model), fill=factor(threshold)))+
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
  scale_color_grey()+
  scale_fill_grey()+
  #geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 6)+
  geom_text(data = aicw_emp, aes(x=model, y=y, label = paste("w =", mean), vjust=0.3, color=threshold)) +  
  ylim(0, 100)+
  theme_classic()+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of selected models")
p2_1


## Theta dot-line plot:  ==============================


group.colors <- c(theta0 = "#636363", theta_1 = "blue")
multi_theta$g_logTheta <- factor(multi_theta$g_logTheta, levels = c("theta_1", "theta0"))
multi_theta <- multi_theta %>% filter(v_logTheta > 0.3)

p3 <- ggplot(data=multi_theta, aes(y = v_logTheta, x = g_logTheta), show.legend=F)+
  geom_point(aes(color=g_logTheta), cex = 2)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   labels=c("theta0" = "Elsewhere", "theta_1" = "Africa"))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(theta~"log(Optimum fruit size) [cm]"))+
  scale_color_manual("Parameters", labels = c("Theta0", "Theta1"), values = group.colors)+
  geom_point(data = multi_theta_MCC, aes(y = v_logTheta, x = g_logTheta), cex =2)+
  geom_line(data= multi_theta_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(25% largest:)~n[OUM]~'= 80,'~n[OUMV]~'= 16,'~n[OUMVA]~'= 5'~~'\r\n'))+
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
                              l = 12))+ # Left margin) #change font size of legend title 
  facet_wrap(threshold~.)
p3
### =======================================
group.colors <- c(X1_sigma.sq = "blue",sigma0 = "#636363")
multi_theta$g_sigma <- factor(multi_theta$g_sigma, levels = c("X1_sigma.sq", "sigma0"))
table(multi_theta$model)


p5 <- ggplot(data=multi_sigma, aes(y = v_sigma, x = g_sigma), show.legend=F)+
  geom_point(aes(color=g_sigma), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_sigma.sq", "sigma0"),
                   labels=c("X1_sigma.sq" = "Africa", "sigma0" = "Elsewhere"))+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(sigma^2~"Evolutionary rate"))+
  scale_color_manual("Parameters", labels = c("Africa", "Elsewhere"), values = group.colors)+
  geom_point(data = multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  geom_line(data= multi_sigma_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 5'~~'\r\n'))+
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
                              l = 12))+ # Left margin) #change font size of legend title 
  facet_wrap(threshold~.)

p5


########

## Models with convergence problems ##


issues_25 <- read.table("issues_25.txt", header=T)
issues_50 <- read.table("issues_50.txt", header=T)
issues_75 <- read.table("issues_75.txt", header=T)

issues_25$threshold <- factor("25")
issues_50$threshold <- factor("50")
issues_75$threshold <- factor("75")


issues <- merge(issues_25, issues_50, all=T)
issues <- merge(issues, issues_75, all=T)
str(issues)

###############################################################################
library(dplyr)
library(rstatix)
library(ggplot2)

issues_dat <- issues %>%
  group_by(threshold, model) %>%  summarise(count = n())
issues_dat$y <- c(6, 104, 5, 87, 6, 7, 104, 9, 90, 9, 5, 104, 11, 85)

issues %>%
  ggplot(aes(x = factor(model)))+ 
  geom_bar(width=.7, position="dodge2", aes(fill=threshold))+
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
  geom_text(data = issues_dat, aes(x=model, y=y, label = paste("n =", count), vjust=0.3, hjust = 0.5, nudge_y = 0.3, color=threshold)) + 
  ylim(0, 105)+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of convergence failures")+
  scale_color_grey()+
  scale_fill_grey()+
  theme_classic()

# ================================================================================================== #


dd_master_sep <- dd_master_bestfit_sep
str(dd_master_sep)
dd_master_sep$df

library(rstatix); library(dplyr); library(Rmisc); library(tidyr); library(qwraps2)

dd_master_sep %>% group_by(threshold) %>% get_summary_stats(type="full")

afr_25 <-dd_master_sep %>% filter(threshold == "25")
afr_50 <-dd_master_sep %>% filter(threshold == "50")
afr_75 <-dd_master_sep %>% filter(threshold == "75")


CI(afr_25$v_logTheta, ci=0.95)
CI(afr_50$v_logTheta, ci=0.95)
CI(afr_75$v_logTheta, ci=0.95)

library(tidyverse)

CI_tab <- dd_master_sep %>%
  group_by(threshold, g_logTheta) %>%
  dplyr::summarise(N=n(),
                   mean.ci = list(mean_ci(v_logTheta)),
                   "Percent"= n_perc(v_logTheta > 0)) %>% 
  unnest_wider(mean.ci)

sum_tab <- Afr_merged %>% filter(delta==0 & id != 101) %>% group_by(threshold, model) %>% dplyr::summarize(
  n = n(),
  alpha0 = mean(X0_alpha) , 
  alpha1 = mean(X1_alpha),
  sigma0 = mean(sigma0),
  sigma1 = mean(X1_sigma.sq),
  theta0 = mean(theta0),
  theta1 = mean(theta_1),
  aicw = mean(w)
)
sum_tab$trees <- "set_of_trees"


sum_tab_MCC <- Afr_merged %>% filter(delta==0 & id == 101) %>% group_by(threshold, model) %>% dplyr::summarize(
  n = n(),
  alpha0 = mean(X0_alpha) , 
  alpha1 = mean(X1_alpha),
  sigma0 = mean(sigma0),
  sigma1 = mean(X1_sigma.sq),
  theta0 = mean(theta0),
  theta1 = mean(theta_1),
  aicw = mean(w)
)
sum_tab_MCC$trees <- "MCC"

sum_tab2 <- merge(sum_tab, sum_tab_MCC, all=T)

write.csv(sum_tab2, "../Figures/OUwie_SumTab_25_50_75.csv")



