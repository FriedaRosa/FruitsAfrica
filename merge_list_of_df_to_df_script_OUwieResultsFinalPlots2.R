## Script to merge list of data frames to big data frame ##

library(dplyr)
library(scales)
library(tidyr)  
library(viridis)
library(viridisLite)
library(ggplot2)
library(rstatix)


# files to be merged to list
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_x_avgl/raw/")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_x_BM/raw/")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_x_avgl/raw/")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_x_BM/raw/")

setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_001/raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_1//raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_10//raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/Afr_8//raw")

setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_001/raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_1//raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_10//raw")
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/md_8//raw")



setwd("~/LinuxBackup/FruitSizeEvo_Project/ClusterOutput/AfricaWorld_empirical/raw")
setwd("~/LinuxBackup/FruitSizeEvo_Project/ClusterOutput/MoistDry_empirical/raw")

rm(list = ls())

# Function to merge data frames in list
multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}  

out <- list()


input_files <-gsub("\\.rds$","", list.files(pattern="\\.rds$"))  

n_input<-(1:length(input_files))

for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  output2 <- output[output$AICc > 0,]
  df2 <- output2
  
  dd_master_sep <- df2[!is.na(df2$eigval1) & df2$eigval1 >= 0 & !is.na(df2$eigval2) & df2$eigval2 >= 0 & !is.na(df2$eigval3) & df2$eigval3 >= 0,]
  df4 <- dd_master_sep %>% filter(theta0 > -7 & theta0 < 7 & theta_1 > -7 & theta_1 < 7) #exclude weird estimates
  #df4 <- df4[df4$expTheta0 >= 0 & df4$expTheta1 >=0 & df4$expTheta0 <= 60 & df4$expTheta1 <=60,] 
  
  fit<-geiger::aicw(df4$AICc)
  output2<-cbind(df4, fit)
  out[[i]] <- output2
}

out2 <- list()
for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  out2[[i]] <- output
}


    
# Merge data frames 
merged_data <- multi_join(out, full_join)
merged_data_raw <- multi_join(out2, full_join)

head(merged_data)
str(merged_data)
str(merged_data_raw)
summary(merged_data)
n_input


# Convergence Problems:
unique(merged_data_raw$model)
unique(merged_data$model)


issues<-setdiff(merged_data_raw,merged_data[,1:18])
issues

write.table(issues, file = "issues_md.txt", sep = "\t", row.names = F, col.names = T)
table(issues$model)
sum(table(issues$model))
n_input


# write to .txt
write.table(merged_data, file = "merged_data.txt", sep = "\t", row.names = F, col.names = T)


summary(merged_data)
ggplot()+
  geom_bar(data=issues, aes(x=model, fill=model))
ggplot()+
  geom_bar(data=merged_data_raw, aes(x=model))




################################################################################

#                 Merge root simulations to Masterfile                         #

################################################################################
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/")

md_001 <- read.table("md_001/raw/merged_data.txt", header=T)
md_1 <- read.table("md_1/raw/merged_data.txt", header=T)
md_10 <- read.table("md_10/raw/merged_data.txt", header=T)
md_8 <- read.table("md_8/raw/merged_data.txt", header=T)
md_4 <- read.table("md_x_BM/raw/merged_data.txt", header=T)
md_avgl <- read.table("md_x_avgl/raw/merged_data.txt", header=T)

Afr_001 <- read.table("Afr_001/raw/merged_data.txt", header=T)
Afr_1 <- read.table("Afr_1/raw/merged_data.txt", header=T)
Afr_10 <- read.table("Afr_10/raw/merged_data.txt", header=T)
Afr_8 <- read.table("Afr_8/raw/merged_data.txt", header=T)
Afr_4 <- read.table("Afr_x_BM/raw/merged_data.txt", header=T)
Afr_avgl <- read.table("Afr_x_avgl/raw/merged_data.txt", header=T)

# make root identifier column
md_001$root <- "0.01"
md_1$root <- "1"
md_10$root <- "10"
md_8$root <- "8"
md_4$root <- "4"
md_avgl$root <- "observed"

Afr_001$root <- "0.01"
Afr_1$root <- "1"
Afr_10$root <- "10"
Afr_8$root <- "8"
Afr_4$root <- "4"
Afr_avgl$root <- "observed"

# make sim/observed identifier column
md_avgl$df <- "observed"
md_001$df <- "simulated"
md_1$df <- "simulated"
md_10$df <- "simulated"
md_8$df <- "simulated"
md_4$df <- "simulated"

Afr_avgl$df <- "observed"
Afr_001$df <- "simulated"
Afr_1$df <- "simulated"
Afr_10$df <- "simulated"
Afr_8$df <- "simulated"
Afr_4$df <- "simulated"

# merge Africa and merge moist/dry together seperately
dd_md <- merge(md_001, md_1, all=T)
dd_md <- merge(dd_md, md_4, all=T)
dd_md <- merge(dd_md, md_8, all=T)
dd_md <- merge(dd_md, md_10, all=T)
dd_md <- merge(dd_md, md_avgl, all=T)

dd_md$expTheta0 <- exp(dd_md$theta0)
dd_md$expTheta1 <- exp(dd_md$theta_1)

head(dd_md)

dd_Afr <- merge(Afr_001, Afr_1, all=T)
dd_Afr <- merge(dd_Afr, Afr_4, all=T)
dd_Afr <- merge(dd_Afr, Afr_8, all=T)
dd_Afr <- merge(dd_Afr, Afr_10, all=T)
dd_Afr <- merge(dd_Afr, Afr_avgl, all=T)

dd_Afr$expTheta0 <- exp(dd_Afr$theta0)
dd_Afr$expTheta1 <- exp(dd_Afr$theta_1)

head(dd_Afr)

# make regime identifier column
dd_md$regime <- "Open/Closed"
dd_Afr$regime <- "Africa"

# transform parameters to long format
dd_Afr <- merged_data
dd_Afr$regime <- "Africa"

dd_md <- merged_data
dd_md$regime <- "Open/Closed"
library(tidyr)

#data_long_md <- gather(dd_md, parameter, value, c(expTheta0, expTheta1, theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
data_long_md <- gather(dd_md, parameter, value, c( theta0, theta_1, theta_se0, theta_se_1, sigma0, "1_sigma.sq", "0_alpha", "1_alpha"), factor_key=TRUE)
str(data_long_md)
head(data_long_md)

data_long_Afr <- gather(dd_Afr, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, "1_sigma.sq", "0_alpha", "1_alpha"), factor_key=TRUE)
#data_long_Afr <- gather(dd_Afr, parameter, value, c(expTheta0, expTheta1,theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long_Afr)
head(data_long_Afr)

# grouped by variable:
# data_long_Afr_sep <- gather(dd_Afr, g_expTheta, v_expTheta, c(expTheta0, expTheta1), factor_key=TRUE)
#data_long_Afr_sep <- gather(data_long_Afr_sep, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_Afr_sep <- gather(dd_Afr, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_logThetaSE, v_logThetaSE, c(theta_se0, theta_se_1), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_sigma, v_sigma, c(sigma0, "1_sigma.sq"), factor_key=TRUE)
data_long_Afr_sep <- gather(data_long_Afr_sep, g_alpha, v_alpha, c("0_alpha", "1_alpha"), factor_key=TRUE)
head(data_long_Afr_sep)

# data_long_md_sep <- gather(dd_md, g_expTheta, v_expTheta, c(expTheta0, expTheta1), factor_key=TRUE)
# data_long_md_sep <- gather(data_long_md_sep, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_md_sep <- gather(dd_md, g_logTheta, v_logTheta, c(theta0, theta_1), factor_key=TRUE)
data_long_md_sep <- gather(data_long_md_sep, g_logThetaSE, v_logThetaSE, c(theta_se0, theta_se_1), factor_key=TRUE)
data_long_md_sep <- gather(data_long_md_sep, g_sigma, v_sigma, c(sigma0, "1_sigma.sq"), factor_key=TRUE)
data_long_md_sep <- gather(data_long_md_sep, g_alpha, v_alpha, c("0_alpha", "1_alpha"), factor_key=TRUE)
head(data_long_md_sep)

# merge to masterfile
dd_master <- merge(data_long_md, data_long_Afr, all=T)
head(dd_master)
summary(dd_master)

str(dd_master)
dd_master$root <- as.factor(c("empirical"))
dd_master$regime <- as.factor(dd_master$regime)

#write.table(dd_master, file = "master_file_BM_rootstates.txt", sep = "\t", row.names = F, col.names = T)



# merge to parameter seperated masterfile
dd_master_sep <- merge(data_long_md_sep, data_long_Afr_sep, all=T)
head(dd_master_sep)
summary(dd_master_sep)

str(dd_master_sep)
dd_master_sep$root <- as.factor(c("empirical"))
dd_master_sep$regime <- as.factor(dd_master_sep$regime)

#write.table(dd_master_sep, file = "sep_master_file_BM_rootstates.txt", sep = "\t", row.names = F, col.names = T)



################################################################################

#                       Make Bestfit dataset                                   #

################################################################################
### read back in ==================
dd_master_sep <- read.table("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/sep_master_file_BM_rootstates.txt", header=T)
dd_master <- read.table("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/master_file_BM_rootstates.txt", header=T)

head(dd_master_sep)
head(dd_master)
library(dplyr)


dd_master$root <- factor(dd_master$root, levels=c('observed', '0.01', '1', '4', '8', '10'))
dd_master$df <- factor(dd_master$df, levels=c('simulated', 'observed'))
dd_master$parameter <- factor(dd_master$parameter, levels=c("expTheta0", "expTheta1", "theta0", "theta_1" , "theta_se0", "theta_se_1", "sigma0" , "X1_sigma.sq", "X0_alpha", "X1_alpha" ))
dd_master$parameter <- factor(dd_master$parameter, levels=c("theta0", "theta_1" , "theta_se0", "theta_se_1", "sigma0" , "1_sigma.sq", "0_alpha", "1_alpha" ))

dd_master$regime <- factor(dd_master$regime, levels=c("Open/Closed", "Africa"))


dd_master_sep$root <- factor(dd_master_sep$root, levels=c("observed", '0.01', '1', '4', '8', '10'))
dd_master_sep$df <- factor(dd_master_sep$df, levels=c('simulated', 'observed'))
dd_master_sep$regime <- factor(dd_master_sep$regime, levels=c("Open/Closed", "Africa"))
# parameters:
dd_master_sep$g_expTheta <- factor(dd_master_sep$g_expTheta, levels=c("expTheta0", "expTheta1"))
dd_master_sep$g_logTheta <- factor(dd_master_sep$g_logTheta, levels=c("theta0", "theta_1"))

dd_master_sep$g_logThetaSE <- factor(dd_master_sep$g_logThetaSE, levels=c("theta_se0", "theta_se_1"))
dd_master_sep$g_sigma <- factor(dd_master_sep$g_sigma, levels=c("sigma0", "1_sigma.sq"))
dd_master_sep$g_alpha <- factor(dd_master_sep$g_alpha, levels=c("0_alpha", "1_alpha"))





dd_master_bestfit <- subset(dd_master, delta==0)
head(dd_master_bestfit)
str(dd_master_bestfit)

levels(dd_master_bestfit$root)
unique(dd_master_bestfit$parameter)

dd_master_bestfit_sep <- subset(dd_master_sep, delta==0)
head(dd_master_bestfit_sep)
str(dd_master_bestfit_sep)
levels(dd_master_bestfit_sep$root)
unique(dd_master_bestfit_sep$root)

dd_master_bestfit_sep[is.na(dd_master_bestfit_sep$root),]

## AICcWeight ~ Models plot: ==========================
dd_Afr%>%
  ggplot(aes(y=w, x=factor(model)))+
  geom_boxplot(aes(fill=model), outlier.color = "white", width=1)+ 
  stat_summary(fun=mean, geom="point", size=2)+   #dot for the mean
  stat_boxplot(geom = "errorbar", width = 0.2) +  
  theme_light()+
  scale_fill_viridis(discrete=T)+
  scale_x_discrete(name ="OUwie Models", breaks=c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA"))+
  scale_y_continuous(name = "Akaike weight [w]")+
  labs(fill = "OUwie Models")



dd_Afr$root <- factor(c("empirical"))

dd_md$root <- factor(c("empirical"))

plot_data_Afr <- dd_Afr %>% 
  group_by(model) %>% 
  count(model) %>%
  mutate(percent = n/100)

### Best fit plots: (sorted by root)

p3<-dd_Afr %>%
  filter(delta==0) %>%
  ggplot(aes( x=factor(root), fill=factor(model)))+
  ggtitle("b) Regime: Megafaunal Community Stability")+
  geom_bar( width=1, position="dodge")+
  scale_fill_viridis(discrete=T)+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  labs(fill = "OUwie Models")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+ 
  ylim(0, 100)+
  theme_light()


p4<-dd_md %>%
  filter(delta==0) %>%
  ggplot(aes( x=factor(root), fill=factor(model)))+
  ggtitle("a) Regime: Open vs. Closed Vegetation")+
  geom_bar( width=1, position="dodge")+
  scale_fill_viridis(discrete=T)+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  labs(fill = "OUwie Models")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  ylim(0, 100)+
  theme_light()


pdf("Bestfit_grid.pdf",width=12, height=8)
gridExtra::grid.arrange(p4,p3)
dev.off()



## parameter estimates plots

head(dd_master_sep)


br <- c("theta0", "theta_1", "theta_se0", "theta_se_1", "sigma0", "1_sigma.sq", "0_alpha", "1_alpha")
lbs <- c("log(Theta World)", "log(Theta Africa)", "SE log(Theta World)", "SE log(Theta Africa)", "Sigma World", "Sigma Africa", "Alpha World", "Alpha Africa")

#br <- c("expTheta0", "expTheta1", "theta0", "theta_1", "theta_se0", "theta_se_1", "sigma0", "X1_sigma.sq", "X0_alpha", "X1_alpha")
lbs <- c("Theta World", "Theta Africa", "log(Theta World)", "log(Theta Africa)", "SE log(Theta World)", "SE log(Theta Africa)", "Sigma World", "Sigma Africa", "Alpha World", "Alpha Africa")


#lab1 <- expression(theta["World"])
#  lab2 <- expression(theta["Africa"]) 
  lab3 <- expression(paste("log(", theta["World"], ")"))
  lab4 <- expression(paste("log(", theta["Africa"], ")"))
  lab5 <- expression(paste("SE log(", theta["World"], ")"))
  lab6 <- expression(paste("SE log(", theta["Africa"], ")"))
  lab7 <- expression(sigma["World"])
  lab8 <- expression(sigma["Africa"]) 
  lab9 <- expression(alpha["World"])
  lab10 <-expression(alpha["Africa"]) 
  
#lbs <- c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8, lab9, lab10)
  lbs <- c(lab3, lab4, lab5, lab6, lab7, lab8, lab9, lab10)

png("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/Afr_obsFL_violin_all.tiff", height = 20, width=20, res=200, units="cm")
dd_master_bestfit %>%
  filter(root=="empirical" & value != "1") %>%
  ggplot(aes(parameter, value))+
  ggtitle("Observed Fruit length")+
  geom_violin(scale = "width", aes(fill=parameter), alpha=0.9) + 
  geom_jitter(height = 0, width = 0.2, alpha=0.1, col="light grey")+
  scale_fill_viridis(discrete=T, labels = lbs)+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  scale_x_discrete(name ="OUwie Parameters", breaks=br,
                   labels = lbs)+
  scale_y_continuous(name = "Estimated values")+
  labs(fill = "OUwie Parameter")+
  facet_grid(model~.)
dev.off()

###
png("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/Afr_obsFL_violin_Thetas.tiff", height = 20, width=20, res=200, units="cm")
dd_master_bestfit %>%
  filter(root=="empirical" & value != "1" & parameter %in% c("theta0", "theta_1")) %>%
  ggplot(aes(parameter, value))+
  ggtitle("Observed Fruit length")+
  geom_violin(scale = "width", aes(fill=parameter), alpha=0.9) + 
  geom_jitter(height = 0, width = 0.2, alpha=0.1, col="light grey")+
  scale_fill_viridis(discrete=T, labels = lbs)+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  scale_x_discrete(name ="OUwie Parameters", breaks=br,
                   labels = lbs)+
  scale_y_continuous(name = "Estimated values")+
  labs(fill = "OUwie Parameter")+
  facet_grid(model~.)
dev.off()

## Simulated


yname <- expression(paste("log(estimated optimum trait)"))
png("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/Afr_SimFL_violin_Thetas.tiff", height = 20, width=20, res=200, units="cm")
dd_master_bestfit %>%
  filter(value != "1" & parameter %in% c("theta0", "theta_1")) %>%
  ggplot(aes(parameter, value, fill=root))+
  ggtitle("Observed and Simulated Fruit v_logThetagth")+
  geom_violin(scale = "width", aes(fill=root), alpha=0.9, draw.quantiles=T) + 
  geom_jitter(height = 0, width = 0.2, alpha=0.1, col="light grey")+
  scale_fill_viridis(discrete=T)+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  scale_x_discrete(name ="Regime", breaks=c("theta0", "theta_1"), labels=c("World", "Africa"))+
  scale_y_continuous(name = yname, limits=c(-1, 3.5))+
  labs(fill = "Root State [cm]")+
  facet_grid(model~.)
dev.off()

## Sigma:

yname <- expression(paste("log(estimated evolutionary rate)"))
png("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/Afr_SimFL_violin_sigmas.tiff", height = 20, width=20, res=200, units="cm")
dd_master_bestfit %>%
  filter(value != "1" & parameter %in% c("sigma0", "sigma_1") & model %in% c("OUMV", "OUMVA")) %>%
  ggplot(aes(parameter, value, fill=root))+
  ggtitle("Observed and Simulated Fruit v_logThetagth")+
  geom_violin(scale = "width", aes(fill=root), alpha=0.9, draw.quantiles=T) + 
  geom_jitter(height = 0, width = 0.2, alpha=0.1, col="light grey", aes(fill=root))+
  scale_fill_viridis(discrete=T)+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  scale_x_discrete(name ="Regime", breaks=c("sigma0", "sigma_1"), labels=c("World", "Africa"))+
  scale_y_continuous(name = yname)+
  labs(fill = "Root State [cm]")+
  facet_grid(model~.)
dev.off()

#### WIP NEU #####
par(mar=c(1.5,1.5,1.5,1.5))
#1. open/Closed
p5<-dd_master_bestfit_sep %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(y=v_logTheta, x=root, fill=g_logTheta))+
  ggtitle("a) Regime: Open versus Closed Vegetation (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed ", "Open "))+
  scale_y_continuous(name="log(Estimated Optimum)", limits=c(0.5, 1.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(theta))



# Africa
p6 <-dd_master_bestfit_sep %>%
  filter(regime=="Africa") %>%
  ggplot(aes(y=v_logTheta, x=root, fill=g_logTheta))+
  ggtitle("b) Regime: Megafaunal Community Stability (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Africa", "World"))+
  scale_y_continuous(name="log(Estimated Optimum)", limits=c(0.5, 1.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(theta))

pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/bestfit_thetas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p5,p6)
dev.off()



### All models (not only best fit)
#1. open/Closed
p7<-dd_master_sep %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(y=v_logTheta, x=root, fill=g_logTheta))+
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed Vegetation", "Open Vegetation"))+
  scale_y_continuous(name="log(Estimated Optimum)", limits=c(0.5, 3.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(theta))


# Africa
p8 <-dd_master_sep %>%
  filter(regime=="Africa") %>%
  ggplot(aes(y=v_logTheta, x=root, fill=g_logTheta))+
  ggtitle("b) Regime: Megafaunal Community Stability (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Africa", "World"))+
  scale_y_continuous(name="log(Estimated Optimum)", limits=c(0.5, 3.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(theta))


pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/allmodels_thetas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p7,p8)
dev.off()


####### Sigmas
#1. open/Closed
p9<-dd_master_bestfit_sep %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(y=v_sigma, x=root, fill=g_sigma))+
  ggtitle("a) Regime: Open versus Closed Vegetation (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed ", "Open "))+
  scale_y_continuous(name="Evolutionary Rate", limits=c(-0.8, 4))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(sigma))



# Africa
p10 <-dd_master_bestfit_sep %>%
  filter(regime=="Africa") %>%
  ggplot(aes(y=v_sigma, x=root, fill=g_sigma))+
  ggtitle("b) Regime: Megafaunal Community Stability (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Africa", "World"))+
  scale_y_continuous(name="Evolutionary Rate", limits=c(-0.8, 4))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(sigma))

pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/bestfit_sigmas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p9,p10)
dev.off()



### All models (not only best fit)
#1. open/Closed
p11<-dd_master_sep %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(y=v_sigma, x=root, fill=g_sigma))+
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed Vegetation", "Open Vegetation"))+
  scale_y_continuous(name="Evolutionary Rate", limits=c(0,0.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(sigma))


# Africa
p12 <-dd_master_sep %>%
  filter(regime=="Africa") %>%
  ggplot(aes(y=v_sigma, x=root, fill=g_sigma))+
  ggtitle("b) Regime: Megafaunal Community Stability (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("World", "Africa"))+
  scale_y_continuous(name="Evolutionary Rate", limits=c(0,0.5))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(sigma))


pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/allmodels_sigmas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p11,p12)
dev.off()

#### Alphas:
####### Sigmas
#1. open/Closed



p13<-dd_master_bestfit_sep %>%  
  filter(regime=="Open/Closed", model %in% c("OUMA", "OUMVA")) %>%
  ggplot(aes(y=v_alpha, x=root, fill=g_alpha))+
  ggtitle("a) Regime: Open versus Closed Vegetation (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed ", "Open "))+
  scale_y_continuous(name="Selection Pressure", limits=c(0,0.3))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(alpha))



# Africa
p14 <-dd_master_bestfit_sep %>%
  filter(regime=="Africa", model %in% c("OUMA", "OUMVA")) %>%
  ggplot(aes(y=v_alpha, x=root, fill=g_alpha))+
  ggtitle("b) Regime: Megafaunal Community Stability (best fit models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Africa", "World"))+
  scale_y_continuous(name="Selection Pressure", limits=c(0,0.3))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(alpha))

pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/bestfit_alphas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p13,p14)
dev.off()



### All models (not only best fit)
#1. open/Closed
p15<-dd_master_sep %>%
  filter(regime=="Open/Closed", model %in% c("OUMA", "OUMVA")) %>%
  ggplot(aes(y=v_alpha, x=root, fill=g_alpha))+
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Closed Vegetation", "Open Vegetation"))+
  scale_y_continuous(name="Selective Force", limits=c(0,0.3))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(alpha))


# Africa
p16 <-dd_master_sep %>%
  filter(regime=="Africa", model %in% c("OUMA", "OUMVA"))  %>%
  ggplot(aes(y=v_alpha, x=root, fill=g_alpha))+
  ggtitle("b) Regime: Megafaunal Community Stability (all models)")+
  geom_boxplot(outlier.colour ="white")+
  scale_fill_viridis(discrete=T, labels=c("Africa", "World"))+
  scale_y_continuous(name="Selective Force", limits=c(0,0.3))+
  scale_x_discrete(name ="Data", breaks=c("0.01", "1", "4", "8", "10", "observed"), 
                   labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(fill = expression(alpha))


pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/allmodels_alphas_grid.pdf",width=10, height=8)
gridExtra::grid.arrange(p15,p16)
dev.off()

### Root size -> thetas relationship

dd_master_sep %>%
  filter(regime=="Africa") %>%
  ggplot(aes(group=g_logTheta, y=v_logTheta, x=root))+
  geom_jitter(aes(color=g_logTheta))+
  geom_line(aes(color=g_logTheta))


##########
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# best fit
dd_master_bestfit_sep2<-dd_master_bestfit_sep %>%
  filter(v_logTheta != "0")

dd_master_bestfit_sep3 <- data_summary(dd_master_bestfit_sep2, varname="v_logTheta", 
                    groupnames=c("regime","root","g_logTheta"))

dd_master_bestfit_sep3_sig <- data_summary(dd_master_bestfit_sep2, varname="v_sigma", 
                                       groupnames=c("regime","root","g_sigma", "model"))

dd_master_bestfit_sep3_alph <- data_summary(dd_master_bestfit_sep2, varname="v_alpha", 
                                       groupnames=c("regime","root","g_alpha", "model"))


head(dd_master_bestfit_sep3)

# all models
dd_master_sep2<-dd_master_sep %>%
  filter(v_logTheta != "0")

dd_master_sep3 <- data_summary(dd_master_sep2, varname="v_logTheta", 
                               groupnames=c("regime","root","g_logTheta"))

dd_master_sep3_sig <- data_summary(dd_master_sep2, varname="v_sigma", 
                               groupnames=c("regime","root","g_sigma", "model"))
dd_master_sep3_alph <- data_summary(dd_master_sep2, varname="v_alpha", 
                                   groupnames=c("regime","root","g_alpha", "model"))

head(dd_master_sep3)

dd_master_sep3$root <- factor(dd_master_sep3$root, levels=c("0.01", "1", "4", "8", "10", "observed"))
dd_master_bestfit_sep3$root <- factor(dd_master_bestfit_sep3$root, levels=c("0.01", "1", "4", "8", "10", "observed"))

# Use position_dodge to move overlapped errorbars horizontally
p18<-dd_master_sep3 %>%
  filter(regime=="Africa") %>%
  ggplot(aes(x=g_logTheta, y=v_logTheta, group=regime, cex=regime, color=regime)) + #
  #ggplot(aes(x=g_logTheta, y=v_logTheta)) + 
  ggtitle("b) Regime: Megafaunal Community Stability (all models)")+
  geom_errorbar(aes(ymin=v_logTheta-sd, ymax=v_logTheta+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  #scale_color_viridis(discrete=T, labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Megafauna Regime",breaks=c("theta0", "theta_1"), labels=c("World", "Africa"))+
  scale_y_continuous(name="log(Estimated Optimum)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")



p17<-dd_master_sep3 %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(x=g_logTheta, y=v_logTheta, group=regime, cex=regime, color=regime)) + 
  #ggplot(aes(x=g_logTheta, y=v_logTheta)) + 
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_errorbar(aes(ymin=v_logTheta-sd, ymax=v_logTheta+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  #scale_color_viridis(discrete=T, labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("theta0", "theta_1"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="log(Estimated Optimum)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")


pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/relationship_root_theta.pdf",width=10, height=8)
gridExtra::grid.arrange(p17,p18)
dev.off()


# sigma + alpha:

p17_2<-dd_master_sep3_sig %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(x=g_sigma, y=v_sigma, group=root, cex=root, color=root)) + 
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_errorbar(aes(ymin=v_sigma-sd, ymax=v_sigma+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("theta0", "theta_1"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="log(Estimated Sigma)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")

p17_3<-dd_master_sep3_alph %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(x=g_alpha, y=v_alpha, group=root, cex=root, color=root)) + 
  ggtitle("a) Regime: Open versus Closed Vegetation (all models)")+
  geom_errorbar(aes(ymin=v_alpha-sd, ymax=v_alpha+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("theta0", "theta_1"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="log(Estimated alpha)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")
gridExtra::grid.arrange(p17_2,p17_3)

pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/relationship_root_theta.pdf",width=10, height=8)
gridExtra::grid.arrange(p17,p18)
dev.off()




### Best fit models:

# Use position_dodge to move overlapped errorbars horizontally
p19<-dd_master_bestfit_sep3 %>%
  filter(regime=="Africa") %>%
  ggplot(aes(x=g_logTheta, y=v_logTheta, group=root, cex=root, color=root)) + 
  ggtitle("b) Regime: Megafaunal Community Stability (best fit models)")+
  geom_errorbar(aes(ymin=v_logTheta-sd, ymax=v_logTheta+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Megafauna Regime",breaks=c("theta0", "theta_1"), labels=c("World", "Africa"))+
  scale_y_continuous(name="log(Estimated Optimum)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")



p20<-dd_master_bestfit_sep3 %>%
  filter(regime=="Open/Closed") %>%
  ggplot(aes(x=g_logTheta, y=v_logTheta, group=root, cex=root, color=root)) + 
  ggtitle("a) Regime: Open versus Closed Vegetation (best fit models)")+
  geom_errorbar(aes(ymin=v_logTheta-sd, ymax=v_logTheta+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("theta0", "theta_1"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="log(Estimated Optimum)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")
####
# sigmas and alphas

p19_alph<-dd_master_sep3_alph %>%
  filter(regime=="Africa", model %in% c("OUMA", "OUMVA")) %>%
  ggplot(aes(x=g_alpha, y=v_alpha, group=root, cex=root, color=root)) + 
  ggtitle("d) Regime: Megafaunal Community Stability")+
  geom_errorbar(aes(ymin=v_alpha-sd, ymax=v_alpha+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Megafauna Regime",breaks=c("X0_alpha", "X1_alpha"), labels=c("World", "Africa"))+
  scale_y_continuous(name="Selective force")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")



p20_alph<-dd_master_sep3_alph %>%
  filter(regime=="Open/Closed", model %in% c("OUMA", "OUMVA")) %>%
  ggplot(aes(x=g_alpha, y=v_alpha, group=root, cex=root, color=root)) + 
  ggtitle("c) Regime: Open versus Closed Vegetation")+
  geom_errorbar(aes(ymin=v_alpha-sd, ymax=v_alpha+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("X0_alpha", "X1_alpha"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="Selective force")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")

gridExtra::grid.arrange(p20_alph,p19_alph)



p19_sig<-dd_master_bestfit_sep3_sig %>%
  filter(regime=="Africa", model %in% c("OUMV", "OUMVA")) %>%
  ggplot(aes(x=g_sigma, y=v_sigma, group=root, cex=root, color=root)) + 
  ggtitle("f) Regime: Megafaunal Community Stability")+
  geom_errorbar(aes(ymin=v_sigma-sd, ymax=v_sigma+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Megafauna Regime",breaks=c("sigma0", "X1_sigma.sq"), labels=c("World", "Africa"))+
  scale_y_continuous(name="Evolutionary Rates")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")



p20_sig<-dd_master_bestfit_sep3_sig %>%
  filter(regime=="Open/Closed", model %in% c("OUMV", "OUMVA")) %>%
  ggplot(aes(x=g_sigma, y=v_sigma, group=root, cex=root, color=root)) + 
  ggtitle("e) Regime: Open versus Closed Vegetation")+
  geom_errorbar(aes(ymin=v_sigma-sd, ymax=v_sigma+sd), width=.2, 
                position=position_dodge(0.40), lwd=1) +
  geom_line(lwd=1, position=position_dodge(0.40)) + 
  geom_point(position=position_dodge(0.4), cex=3.5)+
  scale_color_viridis(discrete=T, 
                      labels=c("Simulated 0.01cm", "Simulated 1cm", "Simulated 4cm", "Simulated 8cm", "Simulated 10cm", "Observed"))+
  theme_light()+
  scale_x_discrete(name="Vegetation Regime",breaks=c("sigma0", "X1_sigma.sq"), labels=c("Closed", "Open"))+
  scale_y_continuous(name="Evolutionary Rates)")+
  theme(title=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1, size=15), 
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))+
  labs(color = "Data")

gridExtra::grid.arrange(p20, p19, p20_alph,p19_alph, p20_sig, p19_sig)



pdf("~/Uni/Master/Masterthesis/MA_final_folder/MA_final_folder/Figures/relationship_root_theta_bestfit.pdf",width=10, height=8)
gridExtra::grid.arrange(p20,p19)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Quick summary plots                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggpubr); library(dplyr); library(scales)
ggboxplot(dd_master_bestfit, x = "parameter", y = "value", ylim=c(0,5))

ggplot(dd_master_bestfit, aes(parameter, value)) +
geom_boxplot(aes(fill=parameter)) +
  scale_y_continuous(limits=c(0,5))+
  facet_wrap(regime~root)+
  theme_bw()



dd_master_bestfit_sep %>% filter(regime=="Africa") %>%
  ggplot(aes(x=g_logTheta, y=v_logTheta, fill=root))+
  geom_boxplot()+
  facet_wrap(model~.)+
  scale_y_continuous(labels = comma)
