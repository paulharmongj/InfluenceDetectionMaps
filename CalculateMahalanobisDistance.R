## Power Curve: 
# Goal: Calculate the Mahalanobis Distance 

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

add_outliers <- rbind(rnorm(10, 0, 3)) %>% as_tibble()

names(add_outliers) <- names(yes_structure_scale)

yes_structure_final <- bind_rows(yes_structure_scale, add_outliers) 
yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))


##### Visualize the data and store in a list
datplot <- yes_structure_final %>% mutate(Id = 1:nrow(yes_structure_final), Alpha = ifelse(Group %in% 'Outlier', 1,0.8)) %>%  pivot_longer(1:10) %>% ggplot(aes(name, value, color = Group, group = Id, alpha = Alpha)) + geom_line() + geom_point() + ggtitle("Simulated Data with Outliers Added In") + guides(alpha = FALSE)
datplot


### Calculate the Mahalanobis Distance
mdist <- mahalanobis(add_outliers, rep(0,10), diag(10))



