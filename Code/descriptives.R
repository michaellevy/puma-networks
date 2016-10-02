# From organizeData.R and makeNetworks.R come before this
library(tidyverse)
theme_set(theme_bw())
library(stringr)
library(statnet)
library(netUtils)
library(stargazer)
# el = read_csv("data/derived/edgelist.csv")
# dyadAttr = read_csv("data/derived/dyadInfo.csv")
# nodeAttr = read_csv("data/derived/nodeInfo.csv")
pumanet = readRDS("data/derived/network.RDS")
collapse_pumanet = readRDS("data/derived/simpleNet.RDS")

# ???
nets = list(multi = pumanet, simple = collapse_pumanet)

#basic network stats: density for valued, reciprocity, transitivity, clustering, modulatrity for community detection, sex homophily, triad census
sapply(nets, network.density)
sapply(nets, function(x) 
  summary(x ~ mutual + edges)[1] / summary(x ~ mutual + edges)[2])
sapply(nets, clusteringCoef)
lapply(nets, function(x) summary(x ~ nodematch("sex", diff = TRUE)))
lapply(nets, triad.census) %>% do.call(rbind, .)
# http://learning-raph.com/wp-content/uploads/2016/03/davis-leinhart-triad.png
