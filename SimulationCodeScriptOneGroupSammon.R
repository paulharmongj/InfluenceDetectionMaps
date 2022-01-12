## Simulation Script
## October 5, 2021

## Paul Harmon

#Description: This script simulates data and runs the different mapping techniques over a span of several different runs. 
# We will compare several scenarios. 

#This code looks at a SINGLE GROUP



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
  
  g1 <- rmnorm(n = 10, mean = rep(40, 60), diag(35, 60)) %>%
    t() %>%
    data.frame()

  yes_structure <- g1 %>%
    data.frame() %>%
    mutate(Group = c(rep("1", 60)))
  
  yes_structure_scale <- lapply(yes_structure[,1:10], scale) %>% as_tibble()
  
  ### adds outliers - but not as a single group. Here we have an outlier that includes a really big one, a really small observation, and a very abnormal one
  add_outliers <- rbind(rnorm(10, 3, .5)) %>% as_tibble()
  names(add_outliers) <- names(yes_structure_scale)
  
  yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
  yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))
  
  
  ##### Visualize the data and store in a list
  datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In") + guides(alpha = FALSE)
  plot_data_list[[j]] <- datplot
  
  
  
  #### Run the Permanova mapping (for t-SNE)
  tic = Sys.time()
  xTSNE <- MultiPermanova(yes_structure_final[,1:10], maptype = "Sammon", perp_val = 20)
  toc = Sys.time() 
  perf_time[[j]] <- toc - tic
  
  
  results_list[[j]] <- xTSNE
  
  df2 = tibble(Holdout = 1:length(xTSNE$modellist), PValues = sapply(xTSNE$modellist, pullPval), Fstat = sapply(xTSNE$modellist, pullFstat))
  df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
  df2list[[j]] = df2
  
  #### TSNE Plots########################################################
  #Visualize the results and store in a list
  sigplot <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: P-Values From Adonis")
  plot_tsne_results_list[[j]] <- sigplot
  
  #Visualize the f stats and store in a list
  sigplotf <- ggplot(df2, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: F Stats From Adonis")
  plot_tsne_fstat_results_list[[j]] <- sigplotf
  
  
  
  #### Run the Permanova mapping (for MDS)
  # tic = Sys.time()
  # xCMD <- MultiPermanova(yes_structure_final[,1:10], maptype = "mds")
  # toc = Sys.time()
  # toc - tic 
  # results_list_cmd[[j]] <- xCMD
  # perf_time_cmd[[j]] <- toc - tic
  
  ### CMD Plots #########################################################
  
  # df3 = tibble(Holdout = 1:length(xCMD$modellist), PValues = sapply(xCMD$modellist, pullPval), Fstat = sapply(xCMD$modellist, pullFstat))
  # df3$Colors <- ifelse(df3$PValues < 0.05, TRUE, FALSE)
  # df3list[[j]] = df3
  # 
  # 
  # #Visualize the results and store in a list
  # sigplot2 <- ggplot(df3, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("CMD: P-Values From Adonis")
  # plot_mds_results_list[[j]] <- sigplot2
  # 
  # #Visualize the f stats and store in a list
  # sigplot2f <- ggplot(df3, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("CMD: F Stats From Adonis")
  # plot_mds_fstat_results_list[[j]] <- sigplot2f
  # 
  
}


#### Some methods to get to table of values to assess efficacy of detections -----#

# Here, we just assess which ones are significant because any detections are spurious

spuriousInfluence = function(result_df){
  tmp <- result_df
  det_rate <- length(which(result_df$PValues <= 0.05))/nrow(tmp)
  return(det_rate)}

#then the non-influential 
meantest = function(df){
  mean(df$PValues)
}

#### Generate a table of values to assess efficacy of detections -----------------#

#note - these will give values per run - so we need to take a look across the simulations

#sammon
ci1 <- sapply(df2list, spuriousInfluence)



#MDS
# ci2 <- sapply(df3list, spuriousInfluence)
# 1-ci2


#mean detection rates
mean(ci1)
mean(cni1)

# mean(ci2)
# mean(cni2)



## Code runs in about 2.5 hours with 100 simulations

## Which components do we save? 
##

#saves the output tables with p-values and F-stats (we can make plots from this)
saveRDS(df2list, paste0(Sys.Date(), "_simulation_Onegroup_sammon_list.RDS"))
#saveRDS(df3list, paste0(Sys.Date(), "_simulationNone_mds_list.RDS"))




#readRDS("2021-10-08_simulation1_mds_list.RDS")






