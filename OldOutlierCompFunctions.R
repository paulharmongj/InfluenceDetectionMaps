# preliminary functions: 

## Vectorize Dist: Function that takes a dataframe and returns a vectorized version of a distance matrix applied to that data matrix. 
## This function is used for generating vectorized data when we step through the Mantel test. 

vectorize_dist <- function(df) {
  tmp <- dist(df)
  vec <- c(tmp)
  return(vec)}

## Pull P-Val and Pull F stat return the p-values and f-statstics from the list of output MultiPermanova puts out. 

pullPval <- function(x){temp <- x$aov.tab$`Pr(>F`[1]; return(temp)}

pullFstat <- function(x){temp <- x$aov.tab$`F.Model`[1]; return(temp)}


### Old multipermanova function

Generate_Holdout_Maps <- function(data, j, ...){  
  
  #hold out a single row of the data (in the jth row)
  #note - have to remove, as t-SNE doesn't like NAs
  data_holdout <- data[-j,]
  
  #calculate distances FIRST
  data_dist <- dist(data_holdout)
  
  #### Creates Maps (based on several methods)
  
  #generate the classical MDS
  mds1 <- cmdscale(data_dist, k = 2)
  
  #generate the t-SNE map and store x,y in map list
  x <- Rtsne(data_dist, perplexity = perplexity, pca = pca_options, is_distance = TRUE)
  
  return(list(mds = mds1, tsne = x$Y))
}



MultiPermanova <- function(data, nperms = 5000, perplexity = 10, pca_options = FALSE){
  
  # initialize the size of the data and output lists
  n <- nrow(data)
  maplist <- list()
  cmdlist <- list()
  
  #generate maps on holdout data
  for(j in 1:n){
    
    temp <- Generate_Holdout_Maps(data, j)
    cmdlist[[j]] <- temp$mds
    maplist[[j]] <- temp$tsne
    
  }
  
  
  # Insert NAs into the holdout indices 
  maplist_new <- list()
  cmdlist_new <- list()
  
  #
  for(j in 1:length(maplist)){
    #create temporary vectors with NAs for each of the XY coordingates
    tempX <- R.utils::insert(maplist[[j]][,1], j, NA)
    tempY <- R.utils::insert(maplist[[j]][,2], j, NA)
    
    tempcmdX <- R.utils::insert(cmdlist[[j]][,1], j, NA)
    tempcmdY <- R.utils::insert(cmdlist[[j]][,2], j, NA)
    
    maplist_new[[j]] <- tibble(X = tempX, Y = tempY)
    cmdlist_new[[j]] <- tibble(X = tempcmdX, Y = tempcmdY)
  }
  
  #vectorizes the distances and outputs a list of vectorized distance matrices 
  dist_list <- lapply(maplist_new, vectorize_dist)
  dist_list_cmd <- lapply(cmdlist_new, vectorize_dist)
  
  
  
  ## takes the vectorized distance matrices and smashes them into a single data frame - to build correlation matrix
  
  dist_mat <- dist_list %>% as_tibble(.name_repair = "minimal")
  dist_mat_cmd <- dist_list_cmd %>% as_tibble(.name_repair = "minimal")
  
  #had to remove the NAs here to get this to run (as it removes when I try to use complete.obs later on)
  #dist_df_na <- sapply(dist_mat, na.omit) %>% as_tibble(.name_repair = "minimal")
  #since I'm using pairwise complete obs, I shouldn't need to do this
  
  cormat <- cor(dist_mat, use = "pairwise.complete.obs", method = "pearson")
  cormat %>% dim()
  
  
  cormat_cmd <- cor(dist_mat_cmd, use = "pairwise.complete.obs", method = "pearson")
  cormat_cmd %>% dim()
  
  
  ##generate distance on correlation matrix
  cordist <- sqrt(2*(1-cormat)) %>% as.dist()
  cordist_cmd <- sqrt(2*(1-cormat_cmd)) %>% as.dist()
  
  model_list <- list()
  model_list_cmd <- list()
  
  
  
  #Step through the permanova method: 
  
  for(j in 1: nrow(as.matrix(cordist))){
    #generate the indicator variable for each feature
    indicator_df <- rep(0, nrow(maplist_new[[1]])) %>% as.data.frame()
    indicator_df_cmd <- rep(0, nrow(cmdlist_new[[1]])) %>% as.data.frame()
    
    names(indicator_df) <- names(indicator_df_cmd) <- "Indicator"
    indicator_df$Indicator[j] = 1
    indicator_df_cmd$Indicator[j] = 1
    
    #each item in the list is an Adonis test using the indicator for the holdout variable
    model_list[[j]] <- adonis(cordist ~ Indicator, data = indicator_df, permutations = nperms, parallel = 2)
    model_list_cmd[[j]] <- adonis(cordist_cmd ~ Indicator, data = indicator_df_cmd, permutations = nperms)
    
    #returns the list of ADONIS outputs
    
    
  }
  
  return(list(tsnelist = model_list, cmdlist = model_list_cmd, plot_tsne = maplist, plot_mds = cmdlist))
}
