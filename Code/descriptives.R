library(statnet)
library(tidyverse)
theme_set(theme_bw())
# Note that this uses mean-imputed values in calculations

# Dyadic
dyadAttr = read_csv("data/derived/dyadInfo.csv")
el = read_csv("data/derived/edgelist.csv")
dyadAttr = count(el, cat1, cat2) %>%
  full_join(dyadAttr, .) 
dyadAttr$n[is.na(dyadAttr$n)] = 0
# Correlation of number of kills shared with relatedness and territory overlap:
round(cor(dyadAttr[, c(3, 4, 8)])[1:2, 3], 3)

# Nodal
net = readRDS("data/derived/network.RDS")
library(netUtils)
ndf = 
  sapply(c("indegree", "outdegree"), function(io) degree(net, gmode = "digraph", cmode = io)) %>%
  data.frame(vAttrDF(net), .) %>%
  mutate(male = as.integer(sex == "M"), female = as.integer(sex == "F")) %>%
  rename(shares = outdegree, receives = indegree, weight = weight_kg, age = age_months) %>%
  .[, sapply(., is.numeric)] %>%
  .[, c(3, 4, 1, 2, 5)]

cors = round(cor(ndf), 2)
prunedcors = cors[1:2, 2:5]
prunedcors[2, 1] = NA
prunedcors