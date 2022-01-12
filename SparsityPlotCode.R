df3list <- sparsesim

TF <- matrix(TRUE, nrow = length(df3list[[1]]), ncol = length(df3list))
templist <- list()

for (i in 1:(length(df3list)-1)){
  
  temp <- rep(TRUE, 100)
  ## Identification Rate at given K
  for(j in 1:length(df3list[[1]])){
    temp[j] <- tail(df3list[[i]][[j]]$PValues,1)<0.05
    
  }
  TF[,i] <- temp
  templist[[i]] <- sum(temp)/100
  
}


plot(x = c(1,2,3,4), y = unlist(templist), type = "l", xlab = "Number of Points In Different Mean Group", ylab = "Detection Rate")
points(x = c(1,2,3,4), y = unlist(templist), pch = 20)
title("Detection Rates By Sparsity of Outlier - MDS")
