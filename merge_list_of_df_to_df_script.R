## Script to merge list of data frames to big data frame ##


rm(list = ls())


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
setwd("~/LinuxBackup/FruitSizeEvo_Project/ClusterOutput/AfricaWorld_empirical/bestfit")


out <- list()

input_files <- gsub("\\.rds$", "", list.files(pattern = "\\.rds$"))

for (i in input_files) {
  filepath <- file.path(paste(i, ".rds", sep = ""))
  
  output <- readRDS(filepath)
  str(output)
  output2 <- output[output$AICc > 0,]
  df2 <- output2
  
  df3 <- df2[!is.na(df2$eigval1) & df2$eigval1 >= 0 & !is.na(df2$eigval2) & df2$eigval2 >= 0 & !is.na(df2$eigval3) & df2$eigval3 >= 0,]
  df4 <- df3[df3$theta0 >= 0 & df3$theta_1 >=0,] #exclude weird estimates
  
  fit<-geiger::aicw(df4$AICc)
  output2<-cbind(df4, fit)
  out[[i]] <- output2
}

# Function to merge data frames in list
multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}    
    
# Merge data frames 
merged_data <- multi_join(out, full_join)

# write to .txt
write.table(merged_data, file = "merged_data.txt", sep = "\t", row.names = F, col.names = T)



################################################################################

#                 Merge root simulations to Masterfile                         #

################################################################################
setwd("~/Uni/Master/Masterthesis/EvolutionaryTraitModelsPalmFruits_git/Cluster_Download/woelke/")

md_001 <- read.table("md_001/raw/merged_data.txt", header=T)
md_1 <- read.table("md_1/raw/merged_data.txt", header=T)
md_10 <- read.table("md_10/raw/merged_data.txt", header=T)
md_8 <- read.table("md_8/raw/merged_data.txt", header=T)
md_4 <- read.table("md_x_BM/raw/merged_data.txt", header=T)

Afr_001 <- read.table("Afr_001/raw/merged_data.txt", header=T)
Afr_1 <- read.table("Afr_1/raw/merged_data.txt", header=T)
Afr_10 <- read.table("Afr_10/raw/merged_data.txt", header=T)
Afr_8 <- read.table("Afr_8/raw/merged_data.txt", header=T)
Afr_4 <- read.table("Afr_x_BM/raw/merged_data.txt", header=T)

# make root identifier column
md_001$root <- "0.01"
md_1$root <- "1"
md_10$root <- "10"
md_8$root <- "8"
md_4$root <- "4"

Afr_001$root <- "0.01"
Afr_1$root <- "1"
Afr_10$root <- "10"
Afr_8$root <- "8"
Afr_4$root <- "4"

# merge Africa and merge moist/dry together seperately
dd_md <- merge(md_001, md_1, all=T)
dd_md <- merge(dd_md, md_4, all=T)
dd_md <- merge(dd_md, md_8, all=T)
dd_md <- merge(dd_md, md_10, all=T)
head(dd_md)

dd_Afr <- merge(Afr_001, Afr_1, all=T)
dd_Afr <- merge(dd_Afr, Afr_4, all=T)
dd_Afr <- merge(dd_Afr, Afr_8, all=T)
dd_Afr <- merge(dd_Afr, Afr_10, all=T)
head(dd_Afr)

# make regime identifier column
dd_md$regime <- "moist_dry"
dd_Afr$regime <- "Africa"

# transform parameters to long format
library(tidyr)

data_long_md <- gather(dd_md, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long_md)
head(data_long_md)

data_long_Afr <- gather(dd_Afr, parameter, value, c(theta0, theta_1, theta_se0, theta_se_1, sigma0, X1_sigma.sq, X0_alpha, X1_alpha), factor_key=TRUE)
str(data_long_Afr)
head(data_long_Afr)

# merge to masterfile
dd_master <- merge(data_long_md, data_long_Afr, all=T)
head(dd_master)
summary(dd_master)

str(dd_master)
dd_master$root <- as.factor(dd_master$root)
dd_master$regime <- as.factor(dd_master$regime)

write.table(dd_master, file = "master_file_BM_rootstates.txt", sep = "\t", row.names = F, col.names = T)

################################################################################

#                       Make Bestfit dataset                                   #

################################################################################

dd_master_bestfit <- subset(dd_master, delta==0)
head(dd_master_bestfit)
str(dd_master_bestfit)
levels(dd_master_bestfit$root)
unique(dd_master_bestfit$root)

library(dplyr)
library(rstatix)

dd_master_bestfit %>%
  group_by(parameter) %>%
  get_summary_stats(value, type="mean_sd")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Quick summary plots                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggpubr)
ggboxplot(dd_master_bestfit, x = "parameter", y = "value", ylim=c(0,5))

ggplot(dd_master_bestfit, aes(parameter, value)) +
geom_boxplot(aes(fill=parameter)) +
  scale_y_continuous(limits=c(0,5))+
  facet_wrap(regime~root)+
  theme_bw()

