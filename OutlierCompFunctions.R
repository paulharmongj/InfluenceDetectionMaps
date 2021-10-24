## Script to functionalize this process


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



### Dependent Functions

#generate maps on holdout data (jth iteration is held out)
# this function should be able to operate on multiple mapping techniques or a single one 

## MultiPermanova: Function that takes an ordination method: Either t-SNE or MDS and produces a 

Generate_Holdout_Maps <- function(data, j, maptype = "MDS",...){  
  
  #hold out a single row of the data (in the jth row)
  #note - have to remove, as t-SNE doesn't like NAs
  data_holdout <- data[-j,]
  
  #calculate distances FIRST
  data_dist <- dist(data_holdout)
  
  #### Creates Maps (based on several methods: tSNE and MDS currently)
  
  if(maptype %in% c("MDS", "mds")){
    #generate the classical MDS
    map <- cmdscale(data_dist, k = 2)
  }
  else if(maptype == "tSNE"){
    #generate the t-SNE map and store x,y in map list
    map <- Rtsne(data_dist, perplexity = perplexity, pca = pca_options, is_distance = TRUE)$Y
  }
  else{print("No maptype selected.")}
  
  
  return(list(map = map))
}



## Generates the Permanova test 
# 1: listobject: list of maps (maplist, cmdlist) to evaluate
# 2: nperms - should inherit from above
# 3: j: index to be looped through 
# HDpermanova <- function(listobject, nperms, j){
#   #generate the indicator variable for each feature
#   indicator_df <- rep(0, nrow(listobject[[1]])) %>% as.data.frame()
#   
#   names(indicator_df) <- "Indicator"
#   indicator_df$Indicator[j] = 1
#   
#   #each item in the list is an Adonis test using the indicator for the holdout variable
#   #note that this is set to do parallel processing with 2 cores
#   model_output <- adonis(cordist ~ Indicator, data = indicator_df, permutations = nperms, parallel = 2)
#   
#   return(list(out = model_output)) 
# }


MultiPermanova <- function(data, nperms = 5000, perplexity = 10, pca_options = FALSE, maptype = "MDS"){
  
  # initialize the size of the data and output lists
  n <- nrow(data)
  maplist <- list()
  
  
  #generate maps on holdout data
  for(j in 1:n){
    
    temp <- Generate_Holdout_Maps(data, j, maptype = maptype)
    maplist[[j]] <- temp$map
    
  }
  
  
  # Insert NAs into the holdout indices 
  maplist_new <- list()
  
  #
  for(j in 1:length(maplist)){
    #create temporary vectors with NAs for each of the XY coordingates
    tempX <- R.utils::insert(maplist[[j]][,1], j, NA)
    tempY <- R.utils::insert(maplist[[j]][,2], j, NA)
    
    
    maplist_new[[j]] <- tibble(X = tempX, Y = tempY)
    
  }
  
  #vectorizes the distances and outputs a list of vectorized distance matrices 
  dist_list <- lapply(maplist_new, vectorize_dist)
  
  
  ## takes the vectorized distance matrices and smashes them into a single data frame - to build correlation matrix
  
  dist_mat <- dist_list %>% as_tibble(.name_repair = "minimal")
  
  #had to remove the NAs here to get this to run (as it removes when I try to use complete.obs later on)
  #dist_df_na <- sapply(dist_mat, na.omit) %>% as_tibble(.name_repair = "minimal")
  #since I'm using pairwise complete obs, I shouldn't need to do this
  
  cormat <- cor(dist_mat, use = "pairwise.complete.obs", method = "pearson")
  cormat %>% dim()
  
  
  ##generate distance on correlation matrix
  cordist <- sqrt(2*(1-cormat)) %>% as.dist()
  
  model_list <- list()
  
  
  #Step through the permanova method: 
  
  for(j in 1: nrow(as.matrix(cordist))){
    #generate the indicator variable for each feature
    indicator_df <- rep(0, nrow(maplist_new[[1]])) %>% as.data.frame()
    
    
    names(indicator_df) <- "Indicator"
    indicator_df$Indicator[j] = 1
    
    ### create the permutation matrix
    #initialize
    permat1 <- matrix(0, nrow = length(indicator_df$Indicator), ncol = length(indicator_df$Indicator))
    permat1[1,] <- 1:length(indicator_df$Indicator)
    
    for(i in 2:nrow(permat1)){
      permat1[i,] <- c(i:length(indicator_df$Indicator), 1:(i-1))
    }
    
    
    #each item in the list is an Adonis test using the indicator for the holdout variable
    model_list[[j]] <- adonis(cordist ~ Indicator, data = indicator_df, permutations = permat1, parallel = 2)
    
    #returns the list of ADONIS outputs
    
    
  }
  
  return(list(modellist = model_list, plot_map = maplist))
}





