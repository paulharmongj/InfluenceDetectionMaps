

## Slightly different idea: write function to replicate Carnegie Classifications 
## and then do the holdout step 

cc15 <- read.csv("https://raw.githubusercontent.com/paulharmongj/Carnegie_SEM/master/data/CC2015data.csv")
carndat = cc15


AGcc <- function(x){
  #rank the data
  ranked <- data.frame(x[,1:3],sapply(x[,-c(1:3)],minrank)) 
  #get pc's
  pca.ranked <- prcomp(ranked[,-c(1:4)], scale = TRUE)
  summary <- summary(pca.ranked)
  standard.score <- scale(pca.ranked$x[,1], scale = TRUE, center = TRUE)
  #needs to return the standardized scores
  return(list(scorez = standard.score, sum =summary))
}
#function for percap
PCcc <- function(x){
  #rank the data
  ranked.dat <- data.frame(sapply(x,minrank)) 
  #get pc's
  pc.ranked <- prcomp(ranked.dat, scale = TRUE)
  summary <- summary(pc.ranked)
  standard.score <- scale(pc.ranked$x[,1], scale = TRUE, center = TRUE)
  return(list(scorez = standard.score, sum = summary))
}


# make a slight addition to generate_holdout_maps

Generate_Holdout_Maps <- function(data, j, maptype = "MDS", perp_val = 30){  
  
  #hold out a single row of the data (in the jth row)
  #note - have to remove, as t-SNE doesn't like NAs
  data_holdout <- data[-j,]
  
  #calculate distances FIRST
  data_dist <- data_holdout
  
  #### Creates Maps (based on several methods: tSNE and MDS currently)
  
  if(maptype %in% c("MDS", "mds")){
    #generate the classical MDS in 2 dimensions
    map <- cmdscale(data_dist, k = 2)
  }
  else if (maptype %in% c("Carnegie","carnegie")){
    map <- ccmap(data_dist) #treats carnegie methodology as mapping technique - not converted to dist to start
  }
  else{print("No maptype selected.")}
  
  
  return(list(map = map))
}



ccmap <- function(carndat){
  
  #function for ranking the data
  minrank <- function(x){rank(x, ties.method = "min")}
  
  
  #dataset that we want to use
  cc2015Ps<-
    na.omit(carndat[,c("NAME","BASIC2010","BASIC2015","FACNUM","HUM_RSD","OTHER_RSD","SOCSC_RSD","STEM_RSD","PDNFRSTAFF","S.ER.D","NONS.ER.D")])
  
  #calculate the ranked data
  cc2015.r <- data.frame(cc2015Ps[,1:3],sapply(cc2015Ps[,-c(1:3)],minrank)) 
  
  cc2015percap <- cc2015Ps[,c("PDNFRSTAFF","S.ER.D","NONS.ER.D")]/cc2015Ps$FACNUM
  colnames(cc2015percap) <- c("PDNRSTAFF_PC", "S.ER.D_PC", "NONS.ER.D_PC")
  cc2015percap.r<-data.frame(sapply(cc2015percap,minrank))
  
  #percap
  pc1 <- -1*(PCcc(cc2015percap.r)$scorez)
  
  #agg
  ag1 <- -1*(AGcc(cc2015.r[,-c(1:4)])$scorez)
  
  
  #dfcc = tibble(ag = -ag1$scorez, pc = -pc1$scorez, col = cc2015.r$BASIC2015)
  #ggplot(dfcc, aes(ag, pc, color = factor(col))) + geom_point()
  ccmap = tibble(x = ag1, y = pc1)
  return(as.matrix(ccmap))
}





##### Try some testing

Generate_Holdout_Maps(carndat, maptype = "Carnegie", j = 2)


outcc = MultiPermanova(carndat, maptype = "Carnegie")










