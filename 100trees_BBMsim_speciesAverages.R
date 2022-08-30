setwd("ClusterOutput/BBM_Sim_Afr/")

BBM_list <- list()

input_files <-gsub("\\.rds$","", list.files(pattern="BBM_sim.*\\.rds$"))  

n_input<-(1:length(input_files))

for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  BBM_sim <- readRDS(filepath)
  BBM_sim_df <- as.data.frame(BBM_sim)
  BBM_list[[i]] <- BBM_sim_df}

library(data.table)
rbindlist(BBM_list)

BBM_df_100 <- do.call(cbind, BBM_list[1:100])

Sp_avg <- data.frame(rowMeans(BBM_df_100, na.rm=TRUE))
Sp_avg$SpecName <- rownames(Sp_avg)

write.csv2(Sp_avg, "../../mydata/BBMsim_Sp_mean_100trees.csv")

range(Sp_avg$rowMeans.BBM_df_100..na.rm...TRUE.)
range(BBM_sim)


