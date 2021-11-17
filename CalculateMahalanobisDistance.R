## Power Curve Based on Mahalanobis Distance: 
# Goal: Calculate the Mahalanobis Distance 


mdistlist <- list()
datplotlist <- list()

results_list <- list()
results_list_cmd <- list()

df2list <- list()
df3list <- list()


sdvec <- c(0.5,1,1.5,3)
n_sims = 5


for(j in 1:length(sdvec)){
  
  
  #instantiate lists to store results of each of the models
  results_list[[j]] <- list()
  results_list_cmd[[j]] <- list()
  
  df2list[[j]][[i]] <- list()
  df3list[[j]][[i]] <- list()
  
  for(i in 1:n_sims){
  
  
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
  
  ## Adds the outliers
  
  add_outliers <- rbind(rnorm(10, 0, j)) %>% as_tibble()
  
  names(add_outliers) <- names(yes_structure_scale)
  
  yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
  yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))
  
  
  ##### Visualize the data and store in a list
  datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In") + guides(alpha = FALSE)
  datplotlist[[j]] <- datplot
  
  
  
  ### Calculate the Mahalanobis Distance
  mdist <- mahalanobis(add_outliers, rep(0,10), diag(10))
  mdistlist[[j]] <- mdist
  
  
  ### Now run the simulations and get the detection rate
  
  #### Run the Permanova mapping (for t-SNE)
  
  xTSNE <- MultiPermanova(yes_structure_final[,1:10], maptype = "tSNE", perp_val = 20)
  results_list[[j]][[i]] <- xTSNE
  
  df2 = tibble(Holdout = 1:length(xTSNE$modellist), PValues = sapply(xTSNE$modellist, pullPval), Fstat = sapply(xTSNE$modellist, pullFstat))
  df2list[[j]][[i]] = df2
 
  
  #### Run the Permanova mapping (for MDS)
 
  xCMD <- MultiPermanova(yes_structure_final[,1:10], maptype = "mds")
  results_list_cmd[[j]][[i]] <- xCMD
  
  df3 = tibble(Holdout = 1:length(xCMD$modellist), PValues = sapply(xCMD$modellist, pullPval), Fstat = sapply(xCMD$modellist, pullFstat))
  df3list[[j]][[i]] = df3
  
  }
}


# 
# sapply(results_list_cmd[[1]][[1]]$modellist, pullPval)
# sapply(results_list_cmd[[1]][[2]]$modellist, pullPval)
# 
# 
# 
# 
# 
# results_mat <- matrix(0, nrow = 61, ncol = length(sdvec)*n_sims)
# k_index <- rep(sdvec, each = n_sims)
# for(k in 1:length(sdvec)*n_sims){
# 
#   for(j in 1:length(results_list_cmd[[k]])){
#     k_index[k] <- q
#     results_mat[j,k] <- sapply(results_list_cmd[[q]][[j]]$modellist, pullPval)
#     
#   }
# }


