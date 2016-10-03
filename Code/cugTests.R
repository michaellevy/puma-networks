library(statnet)
library(tidyverse)
source("Code/functions.R")
empNet = readRDS("data/derived/simpleNet.RDS")
sexes = empNet %v% "sex"

st = table(empNet %v% "sex")  # 9 Fs, 5 Ms. So, possible ties:
ft = c(st[1] * (st[1] - 1), rep(st[1] * st[2], 2), st[2] * (st[2] - 1))

calcStats(empNet)


# CUG Tests

nsim = 1e4
nodes = network.size(empNet)
edges = network.edgecount(empNet)
simNets = 
  replicate(nsim, {
    n = network(rgnm(n = 1, nv = nodes, m = edges))
    n %v% "sex" = sample(sexes)
    n
    }, simplify = FALSE)

system.time({
  simStats = lapply(simNets, calcStats)  
})

ss = lapply(simStats, unlist, recursive = FALSE) %>% do.call(rbind, .)
es = unlist(calcStats(empNet))

cugs = data.frame(
  statistic = lapply(calcStats(empNet), names) %>% do.call(c, .),
  empirical_value = es,
  p_sim_greater = sapply(seq_along(es), function(i) sum(es[i] < ss[, i])) / nrow(ss),
  row.names = NULL)
  
cugs = cugs[order(cugs$p_sim_greater), ]
saveRDS(cugs, "models/CUGtests.RDS")