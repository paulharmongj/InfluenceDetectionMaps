## Simulations

#### Preliminary stuff
source("OutlierCompFunctions.R")
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

plot_data_list <- list()
plot_tsne_results_list <- list()
plot_mds_results_list <- list()
results_list <- list()

for(j in 1:3){
  
  
  #### Simulate the Data
  set.seed(34536)
  g1 <- rmnorm(n = 10, mean = rep(20, 10), diag(25, 10)) %>%
    t() %>%
    data.frame()
  g2 <- rmnorm(n = 10, mean = rep(50, 20), diag(25, 20)) %>%
    t() %>%
    data.frame()
  g3 <- rmnorm(n = 10, mean = rep(80, 30), diag(25, 30)) %>%
    t() %>%
    data.frame()
  
  
  yes_structure <- rbind(g1, g2, g3) %>%
    data.frame() %>%
    mutate(Group = c(rep("1", 10), rep("2", 20), rep("3", 30)))
  
  yes_structure_scale <- lapply(yes_structure[,1:10], scale) %>% as_tibble()
  
  ## adds outliers in the last 3 indices (simulating from normal based on 5 sd's from mean)
  add_outliers <- rbind(rnorm(10, 5, .5), rnorm(10, 5, .5), rnorm(10, 5, .5)) %>% as_tibble()
  names(add_outliers) <- names(yes_structure_scale)
  
  
  yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
  yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",3))
  
  
  #Visualize the data and store in a list
  datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In")
  plot_data_list[[j]] <- datplot
  
  
  
  #### Run the Permanova mapping
  
  
  x3 <- MultiPermanova(yes_structure_final[,1:10])
  
  results_list[[j]] <- x3
  
  df2 = tibble(Holdout = 1:length(x3$tsnelist), PValues = sapply(x3$tsnelist, pullPval))
  df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
  
  #Visualize the results and store in a list
  sigplot <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: P-Values From Adonis")
  plot_tsne_results_list[[j]] <- sigplot
  
  
  
  df3 = tibble(Holdout = 1:length(x3$cmdlist), PValues = sapply(x3$cmdlist, pullPval))
  df3$Colors <- ifelse(df3$PValues < 0.05, TRUE, FALSE)
  
  #Visualize the results and store in a list
  sigplot2 <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: P-Values From Adonis")
  plot_mds_results_list[[j]] <- sigplot2
  
  
}

#### Post-Simulations

## How many influential points were correctly classified? 
## How many spurious points were calculated?
## Comparison between t-SNE and MDS


### Identify the number of CORRECTLY Classified influential points
#extract the p-values from the last 3 from MDS
length(which(tail(df3, 3)$PValues < 0.05))/3
# extract the pvalues from the last 3 of the t-SNE
length(which(tail(df2, 3)$Pvalues < 0.05))/3


### Identify rate of CORRECTLY Classified non-influential points
#extract the p-values from the last 3 from MDS
length(which(head(df3, 60)$PValues >= 0.05))/60
# extract the pvalues from the last 3 of the t-SNE
length(which(tail(df2, 60)$PValues >= 0.05))/60





