# First: Set working directory to the folder where the (101x) .rds files from the Cluster are located with the raw model output tables 
# For empirical and simulated data, the upper part script has to be run twice 
# where each time the right line of code is used to read in and save data

setwd("~/Afr_x_avgl/") 
setwd("~/Afr_x_BM/")

# Function to merge data frames in list to single dataset which includes all model results and parameter estimates for all trees

multi_join <- function(list_of_loaded_data, join_func, ...)
{
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}  


# Evaluate convergence problems ===================
out <- list()
input_files <-gsub("\\.rds$","", list.files(pattern=".*Raw.*\\.rds$"))  
n_input<-(1:length(input_files))

for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  out[[i]] <- output
}

out2 <- list()
input_files <-gsub("\\.rds$","", list.files(pattern=".*Processed.*\\.rds$"))  
for(i in seq_along(n_input)){
  filepath <- file.path(paste(input_files[i], ".rds", sep=""))
  print(filepath)
  output <- readRDS(filepath)
  out2[[i]] <- output
}

# Merge data frames 
merged_data <- multi_join(out, full_join) # raw
merged_data_raw <- multi_join(out2, full_join) # without convergence problems

head(merged_data)
str(merged_data)
str(merged_data_raw)
summary(merged_data)
n_input

## Convergence Problems:
unique(merged_data_raw$model)
unique(merged_data$model)

issues<-setdiff(merged_data_raw,merged_data[,1:18])
issues


table(issues$model)
sum(table(issues$model))


## write to .txt

# write.table(merged_data, file = "merged_data_BM.txt", sep = "\t", row.names = F, col.names = T)
# write.table(merged_data, file = "merged_data_AVGL.txt", sep = "\t", row.names = F, col.names = T)

# write.table(issues, file = "issues_AVGL.txt", sep = "\t", row.names = F, col.names = T)
# write.table(issues, file = "issues_BM.txt", sep = "\t", row.names = F, col.names = T)



## Models with convergence problems ##

setwd("~/GitHub/FruitsAfrica/Data/OUwie/out")

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
       title ="a) Simulated fruit sizes ")+
  theme_classic()




