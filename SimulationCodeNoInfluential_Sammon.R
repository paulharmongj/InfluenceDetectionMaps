## Simulation Script - No Outliers
## October 5, 2021

## Paul Harmon

#Description: This script simulates data and runs the different mapping techniques over a span of several different runs. 
# We will compare several scenarios. 



#### Preliminary code -----------------------------------------------------------#

# Library and source relevant functions
source("OutlierCompNew.R")
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(janitor)
library(magrittr)
library(knitr)
library(Rtsne)
library(robustHD)
library(GGally)
library(stringr)
library(plotly)
library(R.utils)
library(vegan)
library(mnormt)



#### Initialize lists and variables ---------------------------------------------#

n_sims = 100

#create a bunch of lists to store stuff later
plot_data_list <- list()
plot_tsne_results_list <- list()
plot_mds_results_list <- list()
plot_tsne_fstat_results_list <- list()
plot_mds_fstat_results_list <- list()
results_list <- list()
results_list_cmd <- list()
df3list = df2list = list()
perf_time <- list()
perf_time_cmd <- list()


#### Simulate the Data and Loop 

for(j in 1:n_sims){
  
  #### Simulate the Data with outliers and a 3 groups (simulates all from one group with mean 40)
  
  g1 <- rmnorm(n = 10, mean = rep(20, 20), diag(35, 20)) %>%
    t() %>%
    data.frame()
  g2 <- rmnorm(n = 10, mean = rep(40, 20), diag(35, 20)) %>%
    t() %>%
    data.frame()
  g3 <- rmnorm(n = 10, mean = rep(60, 20), diag(35, 20)) %>%
    t() %>%
    data.frame()
  
  
  yes_structure <- rbind(g1, g2, g3) %>%
    data.frame() %>%
    mutate(Group = c(rep("1", 20), rep("2", 20), rep("3", 20)))
  
  yes_structure_scale <- lapply(yes_structure[,1:10], scale) %>% as_tibble()
  
  ### No outliers - only 60 observations in this dataset 
  
  yes_structure_final <- yes_structure_scale
  yes_structure_final$Group <- c(yes_structure$Group)
  
  
  ##### Visualize the data and store in a list
  datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with No Outliers Added In") + guides(alpha = FALSE)
  plot_data_list[[j]] <- datplot
  
  
  
  #### Run the Permanova mapping (for t-SNE)
  tic = Sys.time()
  xTSNE <- MultiPermanova(yes_structure_final[,1:10], maptype = "Sammon")
  toc = Sys.time() 
  perf_time[[j]] <- toc - tic
  
  
  results_list[[j]] <- xTSNE
  
  df2 = tibble(Holdout = 1:length(xTSNE$modellist), PValues = sapply(xTSNE$modellist, pullPval), Fstat = sapply(xTSNE$modellist, pullFstat))
  df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
  df2list[[j]] = df2
  
  #### Sammon Plots########################################################
  #Visualize the results and store in a list
  sigplot <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: P-Values From Adonis")
  plot_tsne_results_list[[j]] <- sigplot
  
  #Visualize the f stats and store in a list
  sigplotf <- ggplot(df2, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: F Stats From Adonis")
  plot_tsne_fstat_results_list[[j]] <- sigplotf
  
  
  
  #### Run the Permanova mapping (for MDS)
  tic = Sys.time()
  xCMD <- MultiPermanova(yes_structure_final[,1:10], maptype = "mds")
  toc = Sys.time()
  toc - tic 
  results_list_cmd[[j]] <- xCMD
  perf_time_cmd[[j]] <- toc - tic
  
  ### CMD Plots #########################################################
  
  df3 = tibble(Holdout = 1:length(xCMD$modellist), PValues = sapply(xCMD$modellist, pullPval), Fstat = sapply(xCMD$modellist, pullFstat))
  df3$Colors <- ifelse(df3$PValues < 0.05, TRUE, FALSE)
  df3list[[j]] = df3
  
  
  #Visualize the results and store in a list
  sigplot2 <- ggplot(df3, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("CMD: P-Values From Adonis")
  plot_mds_results_list[[j]] <- sigplot2
  
  #Visualize the f stats and store in a list
  sigplot2f <- ggplot(df3, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("CMD: F Stats From Adonis")
  plot_mds_fstat_results_list[[j]] <- sigplot2f
  
  
}


#### Some methods to get to table of values to assess efficacy of detections -----#

# Correctly Identified Influential
#tmp = tail(df3list[[1]], 3)




#### Generate a table of values to assess efficacy of detections -----------------#
spuriousInfluence = function(result_df){
  det_rate <- length(which(result_df$PValues <= 0.05))/nrow(result_df)
  return(det_rate)}

#then the non-influential 
# meantest = function(df){
#   mean(df$PValues)
# }
#note - these will give values per run - so we need to take a look across the simulations

#t-SNE
ci1 <- sapply(df2list, spuriousInfluence)
cni1 <- 1-ci1


#MDS
ci2 <- sapply(df3list, spuriousInfluence)
cni2 <- 1 - ci2


#mean detection rates
mean(ci1)
mean(cni1)

mean(ci2)
mean(cni2)

## Code runs in about 2.5 hours with 100 simulations

## Which components do we save? 
##

#saves the output tables with p-values and F-stats (we can make plots from this)
saveRDS(df2list, paste0(Sys.Date(), "_simulationNoInfluence_Sammon_list.RDS"))
saveRDS(df3list, paste0(Sys.Date(), "_simulationNoInfluence_mds_list.RDS"))


### Some post-hoc analyisis 
df1 <- tibble(CMD = sapply(xCMD$modellist, pullPval), TSNE = sapply(xTSNE$modellist, pullPval)) %>% mutate(Index = 1:nrow(df1))
df1 %>% pivot_longer(1:2, "group") %>% ggplot(aes(x = value, y = group, group = group)) + geom_boxplot(aes(value, fill = group), alpha = 0.2) + geom_label(aes(label = Index, group = group, fill = group), position = "jitter", alpha = 0.9) + xlab("P-Value") + ylab("Group") + ggtitle("P-Values by Index") + scale_fill_viridis_d("Map Type", labels = c("t-SNE", "MDS"), begin = 0.5, end = 0.9) + theme_bw()


### Stack points on top of each other
par(mfrow = c(1,2))

xCMD$plot_map[[1]] %>% plot(pch = 20, col = 1, main = "MDS Holdout Plot")
for(i in 2:length(xCMD$plot_map)){
  points(xCMD$plot_map[[i]], pch = 20, col = i)
}


xTSNE$plot_map[[1]] %>% plot(pch = 20, col = 1, main = "TSNE Holdout Plot")
for(i in 2:length(xTSNE$plot_map)){
  points(xTSNE$plot_map[[i]], pch = 20, col = i)
}


