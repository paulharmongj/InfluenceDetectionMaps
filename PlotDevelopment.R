## some live code development


#x3 is a list of outputs from multipermanova
out <- x3
out$tsnelist


library(tidyverse)
library(ggplot2)
library(ggridges)



histlist1 <- list()
histlist2 <- list()

for(j in 1:length(out$tsnelist)){
  histlist1[[j]] <- out$tsnelist[[j]]$f.perms %>% as_tibble() #%>% ggplot(aes(V1)) + geom_histogram(bins = 10)
  histlist2[[j]] <-  out$cmdlist[[j]]$f.perms %>% as_tibble() #%>% ggplot(aes(V1)) + geom_density() + scale_x_log10()
}


## Left Panel (stacks the distributions on top of each other)

df = histlist1 %>% bind_rows() %>% as_tibble()
names(df) = "Density"
df$Holdout = rep(1:length(histlist1),each = nrow(histlist1[[1]])) %>% as.factor()

permplot = ggplot(df, aes(Density, color = Holdout)) + geom_density(alpha = 0.7) + 
  scale_color_viridis_d(option = "B", end = 0.9) + guides(color = FALSE) + 
   coord_flip() + theme_bw() + xlab("") # + ggtitle("Permutation Distribution") 

permplot

## Right Panel - Original Plot (shows the f-statistics)
df2 = tibble(Holdout = 1:length(out$tsnelist), PValues = sapply(out$tsnelist, pullPval), Fstat = sapply(out$tsnelist, pullFstat))
df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
df2list[[j]] = df2

#### TSNE Plots########################################################
#Visualize the results and store in a list
sigplot <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() #+ ggtitle("t-SNE: P-Values From Adonis")
sigplot

#Visualize the f stats and store in a list
sigplotf <- ggplot(df2, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ylab("") + scale_color_viridis_d("Significant",option = "B", end = 0.8)#+ ggtitle("t-SNE: F Stats From Adonis")
sigplotf

library(ggpubr)

p1 <- ggarrange(permplot, sigplotf, label.y = "F-Statistic")
method = "T-SNE"
annotate_figure(p1, top = text_grob(paste0("F Statistics with Permutation Distribution: ", method), color = "black", face = "bold", size = 14),
                    left = text_grob("F-Statistic", size = 12, face = "bold", rot = 90))





#### Try for CMD ####################################################################################
df = histlist2 %>% bind_rows() %>% as_tibble()
names(df) = "Density"
df$Holdout = rep(1:length(histlist1),each = nrow(histlist1[[1]])) %>% as.factor()

permplot = ggplot(df, aes(Density+1, color = Holdout)) + geom_density(alpha = 0.7) + 
  scale_color_viridis_d(option = "A", end = 0.9) + guides(color = FALSE) + 
  coord_flip() + theme_bw() + xlab("") + scale_x_log10() 

permplot

## Right Panel - Original Plot (shows the f-statistics)
df2 = tibble(Holdout = 1:length(out$cmdlist), PValues = sapply(out$cmdlist, pullPval), Fstat = sapply(out$cmdlist, pullFstat))
df2$Colors <- ifelse(df2$PValues < 0.05, TRUE, FALSE)
df2list[[j]] = df2

#### CMD Plots########################################################
#Visualize the results and store in a list
sigplot <- ggplot(df2, aes(Holdout, PValues)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() #+ ggtitle("t-SNE: P-Values From Adonis")
sigplot

#Visualize the f stats and store in a list
sigplotf <- ggplot(df2, aes(Holdout, 1+Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ylab("") + scale_color_viridis_d("Significant",option = "A", end = 0.8) + scale_y_log10()
#+ ggtitle("t-SNE: F Stats From Adonis")
sigplotf

library(ggpubr)

p1 <- ggarrange(permplot, sigplotf, label.y = "F-Statistic")
method = "MDS"
annotate_figure(p1, top = text_grob(paste0("F Statistics with Permutation Distribution: ", method), color = "black", face = "bold", size = 14),
                left = text_grob("F-Statistic", size = 12, face = "bold", rot = 90))













#### Mark's code (for reference) ####

f1Fs_long %>% ggplot(mapping=aes(y=Number, x=value, fill = factor(stat(quantile)))) +
  stat_density_ridges(geom = "density_ridges_gradient",
                      calc_ecdf = TRUE,
                      quantiles = 0.95, 
                      alpha=0.3,
                      scale=2
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#A0A0A0A0", "#FF0000A0"),
    labels = c("(0, 0.95]", "(0.95, 1]")
  )  +
  scale_x_continuous(limits = c(-1, 3.5)) +
  geom_point(data=observed, mapping=aes(y=Number, x=FA), inherit.aes = F) +
  geom_line(data=observed, mapping=aes(y=Number, x=FA), inherit.aes=F) +
  coord_flip() + 
  theme_ridges() + 
  geom_vline(xintercept=0, col="grey")

