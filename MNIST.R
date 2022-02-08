# library(readr)
# library(dplyr)
# mnist_raw <- read_csv("https://pjreddie.com/media/files/mnist_test.csv", col_names = FALSE)
# dim(mnist_raw)
# 
# head(mnist_raw)
# names(mnist_raw)

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


### INTRO: 

## This file loads information from MNIST testing dataset (wanted a smaller random sample)
## Higher-dimensional data. 
## We take a carefully chosen sample with three numbers (0, 3, 7) with only a single 7 ideally to identify as influential


#stratified sample with relatively equal representation of digits
#run with ~1000 observations


### From Kaggle - Custom dataset

mnist_test = read.csv("mnist_test.csv")
dim(mnist_test)
table(mnist_test$label)

#pick a group of observations based on 0, 3 and 7
set.seed(1234)
index <- sample(1:892, 100)
mnist_select <- mnist_test %>% dplyr::filter(label ==0) %>% slice(index)
mnist_select_3 <- mnist_test %>% dplyr::filter(label == 3)  %>% slice(index)
mnist_influence <- mnist_test %>%  dplyr::filter(label ==7)

mnist_final <- bind_rows(mnist_select, mnist_select_3, mnist_influence[1,])
#write.csv(mnist_final, "mnist_final2.csv")


### 
dim(mnist_final) #gives us a 1991X785 sample


### MDS

out1 = MultiPermanova(mnist_final[,-1], maptype = "MDS")
saveRDS(out1, "131MNIST_results_MDS.rds")



### t-SNE

#pick a perplexity
# par(mfrow = c(3,2))
# perplist = c(5,10,20,30,50,100)
# plotlist = list()
# for(i in 1:6){
#   df1 = Rtsne(na.omit(mnist_final[,-1]), perplexity = perplist[i])$Y %>% as_tibble()
#   names(df1) = c("V1","V2")
#   df1 <- df1 %>% add_column(label =na.omit(mnist_final[])$label)
#   plotlist[[i]] <- ggplot(df1, aes(V1, V2, color = factor(label))) + geom_point() + ggtitle(paste0("Perp: ", perplist[i])) + scale_color_viridis_d() + theme_bw()
#   
# }

#ggarrange(plotlist = plotlist, nrow = 2, ncol = 3)

out2 = MultiPermanova(mnist_final, maptype = "tSNE", perp_val = 20)
saveRDS(out2, "131MNIST_results_tSNE.rds")



### Quick Functions to analyze 

createMap <- function(inputmap, label){
  dat <- inputmap %>% as_tibble()
  names(dat) <- c("V1","V2")
  #dat$Label = factor(label)
  plot <- ggplot(dat, aes(V1,V2)) + geom_point() + ggtitle(paste0("Obs Held Out Map"))
  return(plot)
}


### Some analysis of results

dfnew = tibble(Holdout = 1:length(out1$modellist), PValues = sapply(out1$modellist, pullPval), Fstat = sapply(out1$modellist, pullFstat))
dfnew$Colors <- ifelse(dfnew$PValues < 0.05, TRUE, FALSE)
dfnew$Label <- factor(mnist_final$label)


mdsplot <- ggplot(dfnew, aes(Holdout, Fstat, color = Label, shape = Label)) + geom_point() + ggtitle("MDS Pseudo F-Stats: MNIST Sample") + theme_bw() + scale_color_viridis_d(begin = .2, end = 0.9, direction = -1) + ylim(-1,6)

which(dfnew$PValues < 0.05)

which.max(dfnew$Fstat)



plotlist <- lapply(out1$plot_map, createMap)
p1 <- plotlist[[1]] + ggtitle("MDS: 1st Obs (0) Held Out") + theme_bw()
p2 <- plotlist[[102]] + ggtitle("MDS: 102nd Obs (3) Held Out") + theme_bw()
p3 <- plotlist[[which.max(dfnew$Fstat)]] + ggtitle("MDS:Obs with Largest Pseudo F (0) Held Out") + theme_bw()

p4 <- plotlist[[201]] + ggtitle("MDS: 7 Held Out") + theme_bw()

library(ggpubr)
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


### Some analysis of results for tSNE

dftsne= tibble(Holdout = 1:length(out2$modellist), PValues = sapply(out2$modellist, pullPval), Fstat = sapply(out2$modellist, pullFstat))
dftsne$Colors <- ifelse(dftsne$PValues < 0.05, TRUE, FALSE)
dftsne$Label <- factor(mnist_final$label)

tsneplot <- ggplot(dftsne, aes(Holdout, Fstat, color = Label, shape = Label)) + geom_point() + ggtitle("tSNE Pseudo F-Stats: MNIST Sample") + theme_bw() + scale_color_viridis_d(begin = .2, end = 0.9, direction = -1) + ylim(-1,6) 

which.max(dftsne$Fstat)

#new tSNE plots
plotlist <- lapply(out2$plot_map, createMap)
p1 <- plotlist[[1]] + ggtitle("tSNE: 1st Obs (0) Held Out") + theme_bw()
p2 <- plotlist[[102]] + ggtitle("tSNE: 102nd Obs (3) Held Out") + theme_bw()
p3 <- plotlist[[which.max(dftsne$Fstat)]] + ggtitle("tSNE:Largest F Obs (0) Held Out") + theme_bw()

p4 <- plotlist[[201]] + ggtitle("tSNE: 7 Obs Held Out") + theme_bw()


ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


## plot of f stats
library(gridExtra)
ggarrange(mdsplot, tsneplot, nrow = 1, common.legend = TRUE, legend = "right" )






#### New dataset containing more observations of each value ####

mnist_sample <- mnist_test %>% slice_sample(n = 300)
table(mnist_sample$label)

















