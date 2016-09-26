library(tidyverse)
theme_set(theme_bw())
library(stringr)
library(statnet)
el = read_csv("data/derived/edgelist.csv")
dyadAttr = read_csv("data/derived/dyadInfo.csv")
nodeAttr = read_csv("data/derived/nodeInfo.csv")

ggplot(dyadAttr, aes(x = relatedness, y = propOverlap, color = sexPair)) + 
  geom_point()

ggplot(nodeAttr, aes(x = age_months, y = weight_kg, color = sex)) + 
  geom_point()

n = 
  igraph::graph_from_data_frame(el, directed = TRUE, vertices = nodeAttr) %>%
  intergraph::asNetwork()
saveRDS(n, "data/derived/network.RDS")

vSize = 3 * n %v% "weight_kg" / max(n %v% "weight_kg", na.rm = TRUE)
vSize[is.na(vSize)] = .1

count(el, cat1, cat2) %>%
  ggplot(aes(n)) + 
  geom_density(adjust = .2, fill = "lightgray") +
  xlab("Number of kills shared (per directed dyad)")

plot(n
     , vertex.col = "sex"
     , vertex.cex = vSize
     , displaylabels = TRUE
     , edge.lwd = 2)