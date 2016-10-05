library(statnet)
library(tidyverse)
theme_set(theme_bw())
library(stringr)
library(netUtils)
library(stargazer)
source("Code/functions.R")
net = readRDS("data/derived/network.RDS")

if(FALSE) {  # Running CUG tests
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
    vCUG = runGUGs(simNets, net)
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
    CUG = runGUGs(simNets, cNet)
  })
  saveRDS(CUG, "models/binaryCUGs.RDS")
}


vstats = readRDS("models/valuedCUGs.RDS")
bstats = readRDS("models/binaryCUGs.RDS")

# ps = full_join(CUG, vCUG, by = "statistic")[, c(1, 3, 5)] %>%
#   rename(binary_p = p_sim_greater.x, valued_p = p_sim_greater.y) 
# 
# ggplot(ps, aes(binary_p, valued_p)) + 
#   geom_point() + 
#   geom_text(aes(label = statistic), hjust = "left", nudge_x = .01, size = 3)

# filter(ps, !str_detect(statistic, "(^tc)|(sex)"))
# The binary stats make sense for the binary-simulated networks, likewise for valued stats.
# Present them separately:

filter(vstats, str_detect(statistic, "^V|(Var)")) %>% 
  arrange(p_sim_greater) %>%
  mutate(statistic = str_replace(statistic, "^V", "")) %>%
  filter(!statistic %in% c("reciprocity", "transitivity")) %>%
  stargazer(type = "text", out = "results/valuedCUG.txt", summary = FALSE, rownames = FALSE)

filter(bstats, !str_detect(statistic, "(^tc)|(sex)")) %>% 
  arrange(p_sim_greater) %>%
  stargazer(type = "text", out = "results/binaryCUG.txt", summary = FALSE, rownames = FALSE)


#### Or, maybe we want the whole simulated distribution.
# Triad census versus valued is silly because it's effectively less dense. Running on binary
nsim = 1e3  #  ~ 1 min / 1k sims now
net = collapseMultiNet(net)
sexes = net %v% "sex"
simNets = 
  replicate(nsim, {
    n = simMultiNet(network.size(net), network.edgecount(net))
    n %v% "sex" = sample(sexes)
    n
  }, simplify = FALSE)

system.time({
  vCUG = runGUGs(simNets, returnAllSims = TRUE)
})
vCUG = as.data.frame(vCUG) 

emp = calcStats(net)[[3]]
emp = data.frame(triad = str_replace(names(emp), "tc-", ""), count = emp)

tc = vCUG[, str_detect(names(vCUG), "triad")]
names(tc) = str_replace(names(tc), "triadCensus\\.tc-", "")

gather(tc, triad, count) %>%
  ggplot(aes(triad, count)) + 
  geom_boxplot() +
  geom_point(data = emp, shape = 4, color = "red", size = 3) +
  scale_y_log10() + annotation_logticks(sides = "l") 
ggsave("results/binaryCUGtriads.png", width = 8, height = 4)