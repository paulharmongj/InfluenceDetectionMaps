#test isomap and sammon mapping
library(MASS)

#####  Test Data  ################################################################################
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

### adds outliers - but not as a single group. Here we have an outlier that includes a really big one, a really small observation, and a very abnormal one
#add_outliers <- rbind(rnorm(10, 3, .5), rnorm(10, -3, .5), rnorm(10, c(rep(-2,3),rep(0,3),rep(2,3),0), .5)) %>% as_tibble()

add_outliers <- rbind(rnorm(10, -3, .5)) %>% as_tibble()

names(add_outliers) <- names(yes_structure_scale)

yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))
#######################################################################################################

d = dist(yes_structure_final[,1:10])
sammon(d, y = cmdscale(d,2), k =2)


#isomap

isomap(dist = d, ndim = 2, k = )

#umap
custom.config = umap.defaults
custom.config$random_state = 123
d.mat = as.matrix(d)
umap.dist = umap(d.mat, config=custom.config, input="dist")
umap.dist






