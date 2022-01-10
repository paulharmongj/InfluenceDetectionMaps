## Simulation Script
## October 28, 2021

## Paul Harmon

#Description: This script simulates data and runs the different mapping techniques over a span of several different runs. 
# We will compare several scenarios. 

#This code looks at a block-like correlation/covariance structure where the outlier breaks that structure



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
library(lavaan)

## function from stackexhange: https://stackoverflow.com/questions/33026183/r-make-symmetric-matrix-from-lower-diagonal


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
  
  ## Starts with 2-factor latent factor model
  myModel <- ' 
             # latent variable definitions 
               f1 =~ x1 + x2 + x3 + x4 + x5
               f2 =~ x6 + x7 + x8 + x9 + x10 
               
           '
  
  #Simulates data from the latent factor model
  nobs= 100
  temp1 = lavaan::simulateData(myModel, sample.nobs = nobs)
  
  
  #Show correlations
  corrplot(cor(temp1))
  #datplot <- temp1 %>% mutate(Id = 1:nrow(temp1)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, group = Id, color= factor(Id))) + geom_line(alpha = 0.5) + geom_point() + guides(color = FALSE) + ggtitle("Simulated Data From Factor Analysis") + guides(alpha = FALSE)
  #datplot
  
  
  add_outliers <- rnorm(10,0,1)
  
  names(add_outliers) <- names(temp1)
  
  yes_structure_final <- bind_rows(temp1, add_outliers) 
  yes_structure_final$Group <- c(rep("NonOutlier",nobs), rep("Outlier",1))
  
  
  ##### Visualize the data and store in a list
  datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In") + guides(alpha = FALSE)
  #datplot
  
  #cor(yes_structure_final[,1:10]) %>% corrplot("number")
  
  plot_data_list[[j]] <- datplot
 
  
  
  #### Run the Permanova mapping (for t-SNE)
  tic = Sys.time()
  xTSNE <- MultiPermanova(yes_structure_final[,1:10], maptype = "tSNE", perp_val = 20)
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

# Here, we just assess which ones are significant because any detections are spurious

correctInfluence = function(result_df){
  tmp <- tail(result_df,1)
  det_rate <- length(which(tmp$PValues <= 0.05))/nrow(tmp)
  return(det_rate)}

# Correctly Identified Non-Influential
#head(df3list[[1]], 60)

correctNonInfluence = function(result_df){
  tmp <- head(result_df,100)
  det_rate <- length(which(tmp$PValues > 0.05))/nrow(tmp)
  return(det_rate)}
#### Generate a table of values to assess efficacy of detections -----------------#

#note - these will give values per run - so we need to take a look across the simulations

#t-SNE
ci1 <- sapply(df2list, correctInfluence)
cni1 <- sapply(df2list, correctNonInfluence)



#MDS
ci2 <- sapply(df3list, correctInfluence)
cni2 <- sapply(df3list, correctNonInfluence)

#mean detection rates
mean(ci1)
mean(cni1)

mean(ci2)
mean(cni2)


#### SAVE OUTPUT ####

#saves the output tables with p-values and F-stats (we can make plots from this)
saveRDS(df2list, paste0(Sys.Date(), "_simulation2FM_tsne_list.RDS"))
saveRDS(df3list, paste0(Sys.Date(), "_simulation2FM_mds_list.RDS"))





