library(statnet)
library(tidyverse)
library(stringr)
library(netUtils)
source("Code/functions.R")
net = readRDS("data/derived/network.RDS")

# For both tests it's about 10s per 100 nsim.
nsim = 1e4
sexes = net %v% "sex"
simNets = 
  replicate(nsim, {
    n = simMultiNet(network.size(net), network.edgecount(net))
    n %v% "sex" = sample(sexes)
    n
    }, simplify = FALSE)

system.time({
  vCUG = runGUGs(net, simNets)
})
saveRDS(vCUG, "models/valuedCUGs.RDS")
# The edges are much more distributed in these, of course. Maybe more useful to simulate binary networks.
par(mfrow = c(1, 2))
plot(collapseMultiNet(net)
     , arrowhead.cex = 2
     , vertex.col = "sex"
     , vertex.cex = 3
     , edge.lwd = "NumberTies"
     , usecurve = TRUE
     , edge.curve = .001
     , main = "Empirical"
)
plot(collapseMultiNet(simNets[[1]])
     , arrowhead.cex = 2
     , vertex.col = "sex"
     , vertex.cex = 3
     , edge.lwd = "NumberTies"
     , usecurve = TRUE
     , edge.curve = .001
     , main = "Simulated"
)

# Simulate binary networks
cNet = collapseMultiNet(net)
simNets = 
  replicate(nsim, {
    n = network(rgnm(n = 1, nv = network.size(cNet), m = network.edgecount(cNet)))
    n %v% "sex" = sample(cNet %v% "sex")
    n
  }, simplify = FALSE)
system.time({
  CUG = runGUGs(cNet, simNets)
})
saveRDS(CUG, "models/binaryCUGs.RDS")
