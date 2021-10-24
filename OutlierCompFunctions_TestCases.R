## Test Cases Four Outlier Comp Functions

x1 <- MultiPermanova(nfl_fin2)
x2 <- MultiPermanova(penguinsX2[,1:4])




## Simulated Data (with clear structure and 3 major ourliers)



# "Clear Structure": with 10, 20, and 30 observations per group
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



x3 <- MultiPermanova(yes_structure_final[,1:10])



## test: 

df2 = tibble(Holdout = 1:length(x3$tsnelist), PValues = sapply(x3$tsnelist, pullPval), Fstat = sapply(x3$tsnelist, pullFstat))
df2$Colors <- ifelse(df2$PValues < 0.05, "tomato2", "blue4")
ggplot(df2, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("t-SNE: F Stats From Adonis")


which(df2$PValues < 0.05)

df2 = tibble(Holdout = 1:length(x3$cmdlist), PValues = sapply(x3$cmdlist, pullPval), Fstat = sapply(x3$cmdlist, pullFstat))
df2$Colors <- ifelse(df2$PValues < 0.05, "tomato2", "blue4")
ggplot(df2, aes(Holdout, Fstat)) + geom_point(aes(color = Colors), size = 2) + geom_line(alpha = 0.5, color = "grey") + theme_bw() + ggtitle("MDS: F Stats From Adonis")


which(df2$PValues < 0.05)



## play with permutation matrix

#initialize
permat1 <- matrix(0, nrow = length(indicator_df$Indicator), ncol = length(indicator_df$Indicator))
permat1[1,] <- 1:length(indicator_df$Indicator)

for(j in 2:nrow(permat1)){
  permat1[j,] <- c(j:length(indicator_df$Indicator), 1:(j-1))
}



#rep(1:length(indicator_df$Indicator), length(indicator_df$Indicator))



