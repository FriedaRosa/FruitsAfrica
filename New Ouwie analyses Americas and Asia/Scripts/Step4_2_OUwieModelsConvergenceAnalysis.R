# First: Set working directory to the folder where the (101x) .rds files from the Cluster are located with the raw model output tables 
# For empirical and simulated data, the upper part script has to be run twice 
# where each time the right line of code is used to read in and save data

setwd("~/Andressa/output_AM_avgl/output_AM_avgl/")
setwd("~/Andressa/output_AS_avgl/output_AS_avgl/")

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
  mods4 <- mods4[mods4$theta_se0 >= 0.001 & mods4$theta_se_1 >=0,]
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
issues


table(issues$model)
table(bestfit$model)
sum(table(issues$model))


## write to .txt
# Americas
# write.table(merged_data_processed, file = "../merged_data_BM_AM.txt", sep = "\t", row.names = F, col.names = T)
 write.table(merged_data_processed, file = "../merged_data_AVGL_AM.txt", sep = "\t", row.names = F, col.names = T)

 # Asia
# write.table(merged_data_processed, file = "../merged_data_BM_AS.txt", sep = "\t", row.names = F, col.names = T)
 write.table(merged_data_processed, file = "../merged_data_AVGL_AS.txt", sep = "\t", row.names = F, col.names = T)
 
 
 # Americas
write.table(issues, file = "../issues_AVGL_AM.txt", sep = "\t", row.names = F, col.names = T)
#write.table(issues, file = "../issues_BM_AM.txt", sep = "\t", row.names = F, col.names = T)

# Asia
write.table(issues, file = "../issues_AVGL_AS.txt", sep = "\t", row.names = F, col.names = T)
#write.table(issues, file = "../issues_BM_AS.txt", sep = "\t", row.names = F, col.names = T)



## Models with convergence problems ##

 setwd("~/Andressa/output_AM_avgl/")

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
setwd("~/Andressa/output_AM_avgl/")
AM_avgl <- read.table("issues_AVGL_AM.txt", header=T)
AM_BM <- read.table("../issues_BM_AM.txt", header=T)


AM_avgl <- issues
AM_avgl$regime <- "Americas"
AM_avgl$dat <- "empirical"
AM_avgl$group <- "AM_avgl"

AM_BM <- issues
AM_BM$regime <- "Americas"
AM_BM$dat <- "simulated"
AM_BM$group <- "AM_BM"

Americas <- merge(AM_avgl, AM_BM, all=T)



# Asia
AS_avgl <- issues
AS_avgl$regime <- "Asia"
AS_avgl$dat <- "empirical"
AS_avgl$group <- "AS_avgl"

AS_BM <- issues
AS_BM$regime <- "Asia"
AS_BM$dat <- "simulated"
AS_BM$group <- "AS_BM"

Asia <- merge(AS_avgl, AS_BM, all=T)

###############################################################################
library(dplyr)
library(rstatix)
library(ggplot2)

AM_avgl %>% filter(dat == "empirical") %>%
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
       title ="a) Americas (observed fruits) ")+
  theme_classic()


AS_avgl %>% filter(dat == "empirical") %>%
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
       title ="b) Asia (observed fruits) ")+
  theme_classic()


# Americas %>% filter(dat == "simulated") %>%
#   ggplot(aes(x = factor(model)))+ 
#   geom_bar(width=.7, position="dodge2")+
#   scale_x_discrete(NULL, limits = c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA"), drop = FALSE)+
#   theme(text= element_text(color="black"),
#         title=element_text(size=15, color = "black"),
#         axis.text=element_text(size=11.5, color = "black"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size=15, color="black"), 
#         axis.title.x = element_text(size=15, color="black"),
#         axis.text.y = element_text(size=15, color="black"),
#         axis.title.y = element_text(size=15),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=10))+ 
#   geom_text(aes(label = ifelse(..count.. > 0, paste0("n = ", ..count..), "")),  stat = "count", nudge_y = 3)+
#   ylim(0, 105)+
#   labs(x = "Evolutionary Trait Models", 
#        y = "Proportion of convergence failures",
#        title ="b) Americas (Simulated fruit sizes)")+
#   theme_classic()
# 



