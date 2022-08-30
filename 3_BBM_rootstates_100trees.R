## 100 trees:
setwd("~/LinuxBackup/FruitSizeEvo_Project/ClusterOutput/BBM_Sim_Afr/")

# ======================
rm(list = ls())
BBM_list <- list()

# for all trees ====
input_files <-gsub("\\.RData$","", list.files(pattern=".*BM.*\\.RData$"))  # only use raw data here

## ======

## extract data loop:

library(ape); library(geiger); library(phytools); library(tidyr)

n_input<-length(input_files)


for(i in seq_along(input_files)){
  
  ## for all trees =====
  filepath <- file.path(paste0(input_files[i], ".RData", sep=""))
  print(filepath)
  load(filepath)

  # =========
  
  BBM_list[[i]] <- root_state 
  
}

BBM_list


# Function to merge data frames in list
multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}    

# Merge data frames 
merged_data <- multi_join(BBM_list, full_join)

merged_data$root_cm <- exp(merged_data$root)
merged_data$tree <- 1:100

range(merged_data$root_cm)

write.csv(merged_data, "FS_at_root_100sims.csv")





