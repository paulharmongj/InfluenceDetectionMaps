## Simulation Script
## January 14, 2022

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
library(MASS)
start = Sys.time()


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
results_list_sammon <- list()
df4list = df3list = df2list = list()
perf_time <- list()
perf_time_cmd <- list()


#### Simulate the Data and Loop 

#outer loop - sets the number of elements that should differ from middle group

for(k in 1:5){
  
  #initializes the list of lists where I'm storing output for now
  plot_data_list[[k]] <- list()
  plot_tsne_results_list[[k]] <- list()
  plot_mds_results_list[[k]] <- list()
  plot_tsne_fstat_results_list[[k]] <- list()
  plot_mds_fstat_results_list[[k]] <- list()
  results_list[[k]] <- list()
  results_list_cmd[[k]] <- list()
  results_list_sammon[[k]] <- list()
  df4list[[k]] <-df3list[[k]] <- df2list[[k]] <- list()
  perf_time[[k]] <- list()
  perf_time_cmd[[k]] <- list()
  
  
  
  #Inner loop - simulates the data and runs the simulations
  
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
    
    ### adds outliers (forcing outlier point to be in upper group)
    add_outliers <- rbind(c(rnorm(k, 1.5, .3), rnorm(10-k, 0, .3))) %>% as_tibble()
    
    names(add_outliers) <- names(yes_structure_scale)
    
    yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
    yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))
    
    
    ##### Visualize the data and store in a list
    datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In") + guides(alpha = FALSE)
    plot_data_list[[k]][[j]] <- datplot
    
    
    
    #### Run the Permanova mapping (for t-SNE)
    
     xTSNE <- MultiPermanova(yes_structure_final[,1:10], maptype = "tSNE", perp_val = 20)
     results_list[[k]][[j]] <- xTSNE
    
     
    
     df2 = tibble(Holdout = 1:length(xTSNE$modellist), PValues = sapply(xTSNE$modellist, pullPval), Fstat = sapply(xTSNE$modellist, pullFstat))
     df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
     df2list[[k]][[j]] = df2
    
    
    
    #### Run the Permanova mapping (for MDS)
    
    xCMD <- MultiPermanova(yes_structure_final[,1:10], maptype = "mds")
    results_list_cmd[[k]][[j]] <- xCMD
    
    
    ### CMD Plots #########################################################
    
    df3 = tibble(Holdout = 1:length(xCMD$modellist), PValues = sapply(xCMD$modellist, pullPval), Fstat = sapply(xCMD$modellist, pullFstat))
    df3$Colors <- ifelse(df3$PValues < 0.05, TRUE, FALSE)
    df3list[[k]][[j]] = df3
    
    
    #### Run the Permanova mapping (for MDS)
    
    xSammon <- MultiPermanova(yes_structure_final[,1:10], maptype = "Sammon")
    results_list_sammon[[k]][[j]] <- xSammon
    
    ### Sammon ############################################################
    
    df4 = tibble(Holdout = 1:length(xSammon$modellist), PValues = sapply(xSammon$modellist, pullPval), Fstat = sapply(xSammon$modellist, pullFstat))
    df4$Colors <- ifelse(df4$PValues < 0.05, TRUE, FALSE)
    df4list[[k]][[j]] = df4
    
    
    
    
  }
}

#### Some methods to get to table of values to assess efficacy of detections -----#

# Correctly Identified Influential
#tmp = tail(df3list[[1]], 3)

correctInfluence = function(result_df){
  tmp <- tail(result_df,1)
  det_rate <- length(which(tmp$PValues <= 0.05))/nrow(tmp)
  return(det_rate)}

# Correctly Identified Non-Influential
#head(df3list[[1]], 60)

correctNonInfluence = function(result_df){
  tmp <- head(result_df,60)
  det_rate <- length(which(tmp$PValues > 0.05))/nrow(tmp)
  return(det_rate)}


#### Generate a table of values to assess efficacy of detections -----------------#

#note - these will give values per run - so we need to take a look across the simulations

# #t-SNE
# ci1 <- sapply(df2list, correctInfluence)
# cni1 <- sapply(df2list, correctNonInfluence)
# 
# 
# #MDS
# ci2 <- sapply(df3list, correctInfluence)
# cni2 <- sapply(df3list, correctNonInfluence)
# 
# 
# #mean detection rates
# mean(ci1)
# mean(cni1)
# 
# mean(ci2)
# mean(cni2)



## Which components do we save? 
##
end = Sys.time()

#saves the output tables with p-values and F-stats (we can make plots from this)
saveRDS(df2list, paste0(Sys.Date(), "_SparsitySimulation_tsne_list.RDS"))
saveRDS(df3list, paste0(Sys.Date(), "_SparsitySimulation_mds_list.RDS"))
saveRDS(df4list, paste0(Sys.Date(), "_SparsitySimulation_sammon_list.RDS"))



