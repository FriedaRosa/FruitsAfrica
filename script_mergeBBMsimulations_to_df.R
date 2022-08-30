## 100 trees:
setwd("~/LinuxBackup/FruitSizeEvo_Project/ClusterOutput/BBM_Sim_Afr/")


# ======================
rm(list = ls())
BBM_list <- list()

# for all trees ====
input_files <-gsub("\\.rds$","", list.files(pattern=".*BBM.*\\.rds$"))  # only use raw data here

## ======

## extract data loop:

library(ape); library(geiger); library(phytools); library(tidyr)

n_input<-length(input_files)


for(i in seq_along(input_files)){
  
  ## for all trees =====
  filepath <- file.path(paste0(input_files[i], ".rds", sep=""))
  print(filepath)
  BBM_sim <- readRDS(filepath)
  
  BBM_df <- data.frame(BBM_sim = BBM_sim, tree=input_files[[i]])
  rownames(BBM_df) <- names(BBM_sim)
  # =========
  
  BBM_list[[i]] <- BBM_df 
  
  }

BBM_list

data <- data.frame(rownames = names(BBM_sim), BBM_simTEST = BBM_sim)

for (i in seq_along(input_files)){
  
  new <- as.data.frame(BBM_list[[i]])
  names(data) <- input_files[[i]]
  data[ , ncol(data) + 1] <- new
  
  }


              
BBM_df <- data.frame(BBM_sim)
