library(statnet)
library(tidyverse)
library(stringr)
library(netUtils)
library(stargazer)
net = readRDS("data/derived/network.RDS")
binaryNet = collapseMultiNet(net)
dyadAttr = read_csv("data/derived/dyadInfo.csv")
set.seed(5723)

related = spread(dyadAttr[, 1:3], cat1, relatedness)[, -1] %>% as.matrix()
ord = match(net %v% "vertex.names", colnames(related))
related = related[ord, ord]
overlap = spread(dyadAttr[, c(1, 2, 4)], cat1, propOverlap)[, -1] %>% as.matrix()
overlap = overlap[ord, ord]
diag(related) = diag(overlap) = 1


mBinary = ergm(binaryNet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + 
             gwidegree(.5, fixed = TRUE) + gwodegree(.5, fixed = TRUE) + 
             dgwesp(.25, TRUE, type = "OTP") + dgwesp(.25, TRUE, type = "ITP"))
bTerms = c("density", "reciprocity", "relatedness", "territory overlap",
           "male receiving", "male sharing", "even-distribution receiving", 
           "even-distribution sharing", "transitivity", "cyclicality")
stargazer(mBinary, type = "text", single.row = TRUE, covariate.labels = bTerms)

mValued = ergm(binaryNet ~ sum + nonzero + mutual(form = "geometric") +
                 nodeifactor("sex") + nodeofactor("sex") + 
                 edgecov(overlap) + edgecov(related) +
                 transitiveweights("geomean","max","geomean") +
                 cyclicalweights("geomean","max","geomean") ,
               reference = ~ Poisson,
               response = "NumberTies",
               control = control.ergm(MCMC.prop.weights="0inflated",
                                      MCMC.prop.args=list(p0=0.5)))
vTerms = c("density", "any shares", "reciprocity", "male receiving", "male sharing",
           "territory overlap", "relatedness", "transitivity", "cyclicality")
stargazer(mValued, type = "text", single.row = TRUE)  

stargazer(mBinary, mValued, type = "text", single.row = TRUE)

vSummary = 
  broom::tidy(mValued) %>%
  mutate(term = vTerms,
         sig = ifelse(p.value < .001, "***",
                      ifelse(p.value < .01, "**", 
                             ifelse(p.value < .05, "*",
                                    ifelse(p.value <= .1, ".", ""))))) %>%
  mutate(valued_edges = paste0(sprintf("%.2f", estimate), sig, " (", sprintf("%.2f", std.error), ")")) %>%
  select(term, valued_edges)

bSummary = 
  broom::tidy(mBinary) %>%
  mutate(term = bTerms,
         sig = ifelse(p.value < .001, "***",
                      ifelse(p.value < .01, "**", 
                             ifelse(p.value < .05, "*",
                                    ifelse(p.value <= .1, ".", ""))))) %>%
  mutate(binary_edges = paste0(sprintf("%.2f", estimate), sig, " (", sprintf("%.2f", std.error), ")")) %>%
  select(term, binary_edges)

full_join(bSummary, vSummary, by = "term")[c(1:6, 9:10, 7:8, 11), ] %>%
  write_csv("results/ergmTable.csv", na = ".")
  # stargazer(type = "text", summary = FALSE, rownames = FALSE,
  #           notes  = "*** p < .001, ** p < .01, * p < .05",
  #           out = "results/modelsTable.txt") 
