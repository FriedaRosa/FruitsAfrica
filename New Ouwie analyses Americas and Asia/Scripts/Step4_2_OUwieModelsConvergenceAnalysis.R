# First: Set working directory to the folder where the (101x) .rds files from the Cluster are located with the raw model output tables 
# For empirical and simulated data, the upper part script has to be run twice 
# where each time the right line of code is used to read in and save data

setwd("~/Andressa/output_AM_avgl/output_AM_avgl/")
setwd("../../output_AS_avgl/output_AS_avgl/")


setwd("~/Andressa/output_AM_BM/output_AM_BM/")
setwd("output_AS_BM")

# Function to merge data frames in list to single dataset which includes all model results and parameter estimates for all trees

multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}  


# Evaluate convergence problems ===================
out <- list()
out2 <- list()
input_files <-gsub("\\.rds$","", list.files(pattern=".*Raw.*\\.rds$", include.dirs = T))  
n_input<-(1:length(input_files))

for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  
  # Check convergence issues:
  # negative eigenvalues, odd estimates, negative AICc (convergence problems)
  mods4 <- output[!is.na(output$eigval1) & output$eigval1 >= 0 & !is.na(output$eigval2) & output$eigval2 >= 0 & !is.na(output$eigval3) & output$eigval3 >= 0,]
  mods4 <- mods4[mods4$theta0 >= 0.001 & mods4$theta_1 >=0,] 
  mods4 <- mods4[mods4$`0_alpha` >= 0 & mods4$`1_alpha` >=0,] 
  mods4 <- mods4[mods4$sigma0 >= 0.001 & mods4$`1_sigma.sq` >=0,]
  mods4 <- mods4[mods4$theta_se0 >= 0.001 & mods4$theta_se_1 >=0 & mods4$theta_se0 <= 10 & mods4$theta_se_1 <= 10,]
  #mods4 <- mods4[mods4$AICc >=  & mods4$theta_se_1 >=0,]
  
  # Model comparison: calculate AICc weights
  library(geiger)
  fit<- geiger::aicw(mods4$AICc)
  mods5<-cbind(mods4, fit)
  
  # select best model based on models without convergence problems
  
  out2[[i]] <- mods5
  out[[i]] <- output
}

merged_data_raw <- multi_join(out, full_join) # raw data
merged_data_processed <- multi_join(out2, full_join)


## only for Asia avgl: ======

# id82 <- merged_data_raw[merged_data_raw$id==82 & !merged_data_raw$model == "OUMVA"  & !merged_data_raw$model == "OUMA",]
# fit2<-geiger::aicw(merged_data_raw[merged_data_raw$id==82 & !merged_data_raw$model == "OUMVA"& !merged_data_raw$model == "OUMA",]$AICc)
# id82<- cbind(id82, fit2)
# 
# merged_data_processed2 <- merged_data_processed[!merged_data_processed$id == 82,]
# merged_data_processed3 <- rbind(merged_data_processed2, id82 )
# merged_data_processed <- merged_data_processed3
###

summary(merged_data_processed)
bestfit <- merged_data_processed[merged_data_processed$delta == 0,]


head(merged_data_processed)
str(merged_data_processed)
str(merged_data_raw)
summary(merged_data_processed)
n_input

## Convergence Problems:
unique(merged_data_raw$model)
unique(merged_data_processed$model) #

issues<-setdiff(merged_data_raw,merged_data_processed[,1:18])

## only for Asia avgl: =====
# x <- merged_data_raw[merged_data_raw$id==82 & merged_data_raw$model == "OUMVA",]
# issues2<- rbind(issues, x )
# issues <- issues2
####

table(issues$model)
table(bestfit$model)
sum(table(issues$model))


## write to .txt
# Americas
# write.table(merged_data_processed, file = "../merged_data_BM_AM.txt", sep = "\t", row.names = F, col.names = T)
# write.table(merged_data_processed, file = "../merged_data_AVGL_AM.txt", sep = "\t", row.names = F, col.names = T)

 # Asia
# write.table(merged_data_processed, file = "../merged_data_BM_AS.txt", sep = "\t", row.names = F, col.names = T)
# write.table(merged_data_processed, file = "../merged_data_AVGL_AS.txt", sep = "\t", row.names = F, col.names = T)
 
 
 # Americas
# write.table(issues, file = "../issues_AVGL_AM.txt", sep = "\t", row.names = F, col.names = T)
# write.table(issues, file = "../issues_BM_AM.txt", sep = "\t", row.names = F, col.names = T)

# Asia
# write.table(issues, file = "../issues_AVGL_AS.txt", sep = "\t", row.names = F, col.names = T)
# write.table(issues, file = "../issues_BM_AS.txt", sep = "\t", row.names = F, col.names = T)

## Models with convergence problems ##

setwd("~/ArbeitImUrlaub/GitHub/FruitsAfrica/Data/OUwie/out/")

# Africa
Afr_avgl <- read.table("issues_AVGL.txt", header=T)
Afr_BM <- read.table("issues_BM.txt", header=T)
Afr_avgl$regime <- "Africa"
Afr_avgl$dat <- "empirical"
Afr_BM$regime <- "Africa"
Afr_BM$dat <- "simulated"
Afr_avgl$group <- "Afr_emp"
Afr_BM$group <- "Afr_BM"

Africa <- merge(Afr_avgl, Afr_BM, all=T)

# Americas
setwd("~/ArbeitImUrlaub/Andressa/output_AM_avgl/")
AM_avgl <- read.table("issues_AVGL_AM.txt", header=T)
AM_BM <- read.table("../output_AM_BM/issues_BM_AM.txt", header=T)

AM_avgl$regime <- "Americas"
AM_avgl$dat <- "empirical"
AM_avgl$group <- "AM_avgl"

AM_BM$regime <- "Americas"
AM_BM$dat <- "simulated"
AM_BM$group <- "AM_BM"
Americas <- merge(AM_avgl, AM_BM, all=T)



# Asia
setwd("~/ArbeitImUrlaub/Andressa/output_AS_avgl/")
AS_avgl <- issues
AS_avgl <- read.table("issues_AVGL_AS.txt", header=T)
AS_avgl$regime <- "Asia"
AS_avgl$dat <- "empirical"
AS_avgl$group <- "AS_avgl"

AS_BM <- read.table("../output_AS_BM/issues_BM_AS.txt", header=T)
AS_BM$regime <- "Asia"
AS_BM$dat <- "simulated"
AS_BM$group <- "AS_BM"
Asia <- merge(AS_avgl, AS_BM, all=T)


all_avgl0 <- merge(Africa, Americas, all=T)
all_avgl<- merge(all_avgl0, Asia, all=T)


table(all_avgl$group)

str(all_avgl)
all_avgl$model <- as.factor(all_avgl$model)
all_avgl$regime <- as.factor(all_avgl$regime)
all_avgl$dat <- as.factor(all_avgl$dat)
all_avgl$group <- as.factor(all_avgl$group)

###############################################################################
library(dplyr)
library(rstatix)
library(ggplot2)


all_avgl %>% filter(dat == "empirical") %>%
  ggplot(aes(x = factor(model)))+ 
  facet_wrap(regime~.)+
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
       title ="a) Proportion of convergence failures (empirical) ")+
  theme_classic()

all_avgl %>% filter(dat == "simulated") %>%
  ggplot(aes(x = factor(model)))+ 
  facet_wrap(regime~.)+
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
       title ="b) Proportion of convergence failures (simulated) ")+
  theme_classic()



all_avgl %>% filter(dat == "empirical") %>%
  ggplot(aes(x = factor(model)))+ 
  facet_wrap(regime~.)+
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
       title ="a) Proportion of convergence failures (empirical) ")+
  theme_classic()


all_avgl %>% filter(dat == "empirical" &regime %in% c("Asia", "Americas")) %>%
  ggplot(aes(x = factor(model)))+
  facet_wrap(regime~.)+
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
       y = "Proportion of convergence failures ",
       title ="a) Proportion of convergence failures (empirical)")+
  theme_classic()

all_avgl %>% filter(dat == "simulated" &regime %in% c("Asia", "Americas")) %>%
  ggplot(aes(x = factor(model)))+
  facet_wrap(regime~.)+
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
       y = "Proportion of convergence failures ",
       title ="b) Proportion of convergence failures (simulated)")+
  theme_classic()






## Parameter estimate plots ==============


# Libraries  =======================
library(tidyr); library(dplyr); library(scales); library(viridis); library(viridisLite); library(ggplot2)
# =======================
rm(list = ls())
setwd("~/ArbeitImUrlaub/Andressa/output_AS_avgl/")
# =======================

# 100 trees results:
Asia_emp <- read.table("merged_data_AVGL_AS.txt", header=T)
Amer_emp <- read.table("../output_AM_avgl/merged_data_AVGL_AM.txt", header=T)
Afr_emp <- read.table("../../../GitHub/FruitsAfrica/Data/OUwie/out/merged_data_AVGL.txt", header=T)

Asia_sim <- read.table("../output_AS_BM/merged_data_BM_AS.txt", header=T)
Amer_sim <- read.table("../output_AM_BM/merged_data_BM_AM.txt", header=T)
Afr_sim <- read.table("../../../GitHub/FruitsAfrica/Data/OUwie/out/merged_data_BM.txt", header=T)



Asia_sim$regime <- "Asia"
Amer_sim$regime <- "Americas"
Afr_sim$regime <- "Africa"

Asia_emp$regime <- "Asia"
Amer_emp$regime <- "Americas"
Afr_emp$regime <- "Africa"

emp0 <- merge(Afr_emp, Asia_emp, all=T)
emp <- merge(emp0, Amer_emp, all=T)
emp$dat <- factor("empirical")

sim0 <- merge(Afr_sim, Asia_sim, all=T)
sim <- merge(sim0, Amer_sim, all=T)
sim$dat <- factor("simulated")

all_merged <- merge(emp, sim, all = T)
summary(all_merged)

# make data id column (empirical, simulated)
all_merged$data <- "empirical"

# make MCC/set id column (MCC, set_trees)
all_merged$tree <- "set_trees"

all_merged_wide <- all_merged

#
write.table(all_merged_wide, file = "../../all_merged_wide_Afr_Asi_Amer.txt", sep = "\t", row.names = F, col.names = T)

# transform parameters to long-format
## megafauna:
data_long <- gather(all_merged_wide, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long)
head(data_long)

# round values
data_long$value <- round(as.numeric(data_long$value),3)
data_long$w <- round(data_long$w, 2)
data_long$delta <- round(data_long$delta, 2)


# merge to masterfile ============================================================
dd_master <- data_long
head(dd_master)
summary(dd_master)

str(dd_master)
dd_master$data <- as.factor(dd_master$data)
dd_master$regime <- as.factor(dd_master$regime)
dd_master$tree <- as.factor(dd_master$tree)

write.table(dd_master, file = "../../all_merged_long_Afr_Asi_Amer.txt", sep = "\t", row.names = F, col.names = T)
# ================================================================================

# dataframes grouped by variables:
## megafauna
data_long_sep <- gather(all_merged_wide, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_logThetaSE, v_logThetaSE, c(theta_se0, theta_se_1), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_sigma, v_sigma, c(sigma0, X1_sigma.sq), factor_key=TRUE)
data_long_sep <- gather(data_long_sep, g_alpha, v_alpha, c(X0_alpha, X1_alpha), factor_key=TRUE)
head(data_long_sep)


# merge to parameter seperated masterfile ============================================================
dd_master_sep <- data_long_sep
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
multi_theta <- dd_master_bestfit_sep %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
multi_sigma <- dd_master_bestfit_sep %>% filter(model %in% c("OUMV", "OUMVA", "BMS"))
multi_alpha <- dd_master_bestfit_sep %>% filter(model %in% c("OUMA", "OUMVA"))



# =================================


# ================================================================================
all_best <- all_merged_wide %>% filter(delta==0)
#multi_theta <- all_best %>% filter(model %in% c("OUM", "OUMV", "OUMVA", "OUMA"))
#multi_sigma <- all_best %>% filter(model %in% c("OUMV", "OUMVA"))


## Fig 1: proportion of selected models ===============
dd <- all_best  
group.colors <- c( OUM = "#cdcdcd", OUMV = "#cdcdcd", OUMVA = "#cdcdcd")

empty_models <-c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")

aicw_emp <- dd %>%
  filter(delta==0 & data == "empirical" & tree == "set_trees") %>% group_by(regime,model) %>%  summarise_at(vars(w), list(mean = mean, sd = sd))
aicw_emp$mean <- round(aicw_emp$mean,2)
aicw_emp$y <- c(83, 17, 3, 
                17, 2, 46, 39, 
                7, 6, 39, 52)+0.5

aicw_emp$y <- c(83, 17, 7, 
                17, 46, 17, 
                6, 2,
                3, 50, 39)+0.5

write.csv(aicw_emp, "~/GitHub/FruitsAfrica/New Ouwie analyses Americas and Asia/Results/AIC_weights_Afr_Amer_Asi.csv")



dd2 <- dd %>% filter(regime =="Asia")
table(dd2$model)
dd1 <- dd %>% filter(regime == "Americas")
table(dd1$model)

dd0 <- dd %>% filter(regime =="Africa")
table(dd0$model)
# Make summary table for ouwie: =====================

sum_tab <- all_merged_wide %>% filter(delta==0) %>% group_by(regime, tree, data, model) %>% dplyr::summarize(
  n = n(),
  alpha0 = mean(X0_alpha) , 
  alpha1 = mean(X1_alpha),
  sigma0 = mean(sigma0),
  sigma1 = mean(X1_sigma.sq),
  theta0 = mean(theta0),
  theta1 = mean(theta_1),
  aicw = mean(w)
)
write.csv(sum_tab, "~/GitHub/FruitsAfrica/New Ouwie analyses Americas and Asia/Results/OUwie_SumTab_Americas_Asia_Africa.csv")

# =============================================
library("ggrepel")      

group.colors <- c(Africa = "#ecd70d",Americas = "#21609a", Asia = "#843c54")


p2_1<-dd %>%
  filter(delta==0 & data == "empirical") %>%
  ggplot(aes( x=factor(model), fill=regime))+
  geom_bar(width=.7, position="dodge2")+
  scale_x_discrete(limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
  scale_fill_manual("Regimes", labels = c("Africa", "Americas", "Asia"), values = group.colors)+
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
  #geom_text(data = aicw_emp, aes(x=model, y=y, label = paste("w =", mean), vjust=0.3, angle=45)) +  
  theme(text = element_text(angle = 180, hjust = 1))+
  ylim(0, 100)+
  theme_classic()+
  labs(x = "Evolutionary Trait Models", 
       y = "Proportion of selected models",
       title ="Proportion of selected models")
p2_1


## Theta dot-line plot:  ==============================


group.colors <- c(theta_1 = "blue", theta0 = "#636363",)
multi_theta$g_logTheta <- factor(multi_theta$g_logTheta, levels = c("theta_1", "theta0"))

p3 <- 
  multi_theta %>% filter(dat=="empirical" & 
                           delta == "0") %>%
  ggplot(aes(y = v_logTheta, x = g_logTheta), show.legend=T)+
  facet_wrap(~regime, scales="fixed")+
  geom_point(aes(color=g_logTheta), cex = 2)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, labels=NULL, limits=c("theta_1", "theta0"))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(title="a) Model estimated trait optima across realms (empirical)")+
  scale_y_continuous(name = expression(theta~"log(Optimum fruit size) [cm]"))+
  scale_color_manual("Parameters", labels = c("Theta0", "Theta1"), values = group.colors)+
  #geom_point(data = multi_theta_MCC, aes(y = v_logTheta, x = g_logTheta), cex =2)+
  #geom_line(data= multi_theta_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUM]~'= 82,'~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  #theme(legend.position="none")+
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


p4 <- 
  multi_theta %>% filter(dat=="simulated" & 
                           delta == "0") %>%
  ggplot(aes(y = v_logTheta, x = g_logTheta), show.legend=T)+
  facet_wrap(~regime, scales="fixed")+
  geom_point(aes(color=g_logTheta), cex = 2)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, labels = NULL, limits=c("theta_1", "theta0"))+
  theme_classic()+
  labs(title="b) Model estimated trait optima across realms (simulated)")+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(theta~"log(Optimum fruit size) [cm]"))+
  scale_color_manual("Parameters", labels = c("Theta0", "Theta1"), values = group.colors)+
  #geom_point(data = multi_theta_MCC, aes(y = v_logTheta, x = g_logTheta), cex =2)+
  #geom_line(data= multi_theta_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUM]~'= 82,'~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  #theme(legend.position="none")+
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

p4

### =======================================
group.colors <- c(X1_sigma.sq = "blue",sigma0 = "#636363")
multi_theta$g_sigma <- factor(multi_theta$g_sigma, levels = c("X1_sigma.sq", "sigma0"))


p5 <- multi_sigma %>% filter(dat=="empirical") %>%
  ggplot(aes(y = v_sigma, x = g_sigma), show.legend=F)+
  geom_point(aes(color=g_sigma), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_sigma.sq", "sigma0"),
                   labels=NULL)+
  facet_wrap(~regime, scales="fixed")+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(sigma^2~"Evolutionary rate"))+
  scale_color_manual("Parameters", labels = NULL, values = group.colors)+
  #geom_point(data = afr_multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  #geom_line(data= afr_multi_sigma_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  theme(legend.position="none")+
  labs(title="a) Model estimated trait evolutionary rate across realms (empirical)")+
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

## not possible for simulated
p6 <- multi_sigma %>% filter(dat=="simulated") %>%
  ggplot(aes(y = v_sigma, x = g_sigma), show.legend=F)+
  geom_point(aes(color=g_sigma), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_sigma.sq", "sigma0"),
                   labels=NULL)+
  facet_wrap(~regime, scales="fixed")+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(sigma^2~"Evolutionary rate"))+
  scale_color_manual("Parameters", labels = NULL, values = group.colors)+
  #geom_point(data = afr_multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  #geom_line(data= afr_multi_sigma_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  theme(legend.position="none")+
  labs(title="b) Model estimated trait evolutionary rate across realms (simulated)")+
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


p6

## ===============================================

group.colors <- c(X1_alpha = "blue", X0_alpha = "#636363")
multi_alpha$g_alpha <- factor(multi_alpha$g_alpha, levels = c("X1_alpha", "X0_alpha"))


p7 <- multi_alpha %>% filter(dat=="empirical") %>%
  ggplot(aes(y = v_alpha, x = g_alpha), show.legend=F)+
  geom_point(aes(color=g_alpha), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_alpha", "X0_alpha"),
                   labels=NULL)+
  facet_wrap(~regime, scales="fixed")+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(alpha~"Selective force"))+
  scale_color_manual("Parameters", labels = NULL, values = group.colors)+
  #geom_point(data = afr_multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  #geom_line(data= afr_multi_sigma_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  theme(legend.position="none")+
  labs(title="a) Model estimated selective force across realms (empirical)")+
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


p7


p8 <- multi_alpha %>% filter(dat=="simulated") %>%
  ggplot(aes(y = v_alpha, x = g_alpha), show.legend=F)+
  geom_point(aes(color=g_alpha), cex = 3)+
  geom_line(aes(group=id), col="light grey")+
  scale_x_discrete(name =NULL, 
                   limits = c("X1_alpha", "X0_alpha"),
                   labels=NULL)+
  facet_wrap(~regime, scales="fixed")+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_y_continuous(name =expression(alpha~"Selective force"))+
  scale_color_manual("Parameters", labels = NULL, values = group.colors)+
  #geom_point(data = afr_multi_sigma_MCC, aes(y = v_sigma, x = g_sigma), cex =3)+
  #geom_line(data= afr_multi_sigma_MCC, aes(group = id), col="black")+
  #labs(subtitle = expression(bold(Observed)~n[OUMV]~'= 16,'~n[OUMVA]~'= 2'~~'\r\n'))+
  theme(legend.position="none")+
  labs(title="a) Model estimated selective force across realms (simulated)")+
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


p8


