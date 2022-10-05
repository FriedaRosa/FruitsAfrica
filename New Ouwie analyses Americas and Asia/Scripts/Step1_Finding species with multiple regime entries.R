df <- read.csv("/home/acabral/Desktop/PhD iDiv/Co-supervision of master's thesis/Frieda/Manuscript New Phytologist/Version 3 - new analyses comparing Americas-elsewhere and Asia-elsewhere/Data/palms_final.csv")
sum(is.na(df$accRealm))

# Assigning species to unique regimes

# Finding multiple regimes per species ####
us <- unique(df$SpecName)
subset_sp<-list()
result<-list()
for (i in 1:length(us)){
  subset_sp[[i]]<-subset(df, SpecName==us[i])
  result[[i]] <- which(length(na.omit(unique(subset_sp[[i]]$accRealm)))>1)
}
position <- which(result == 1)
list_to_check <- us[position]
#List with names to check fruit type
list_to_check

#Including this data in the original matrix
index <- list()
for (i in 1:length(list_to_check)){
  index[[i]] <- which(df$SpecName == list_to_check[i]) 
}
new_index <- unlist(index)

df$regimes_to_check <- NA
df$regimes_to_check[new_index] <- "to_check"

write.csv(df,"Step1_duplicated_regimes_to_check")
