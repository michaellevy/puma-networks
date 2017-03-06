library(statnet)
library(tidyverse)
library(stringr)
library(netUtils)
library(stargazer)
net = readRDS("data/derived/network.RDS")
binaryNet = collapseMultiNet(net)
dyadAttr = read_csv("data/derived/dyadInfo.csv")

related = spread(dyadAttr[, 1:3], cat1, relatedness)[, -1] %>% as.matrix()
ord = match(net %v% "vertex.names", colnames(related))
related = related[ord, ord]
overlap = spread(dyadAttr[, c(1, 2, 4)], cat1, propOverlap)[, -1] %>% as.matrix()
overlap = overlap[ord, ord]
diag(related) = diag(overlap) = 1

mBinary = ergm(binaryNet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + 
             gwidegree(.5, fixed = TRUE) + gwodegree(.5, fixed = TRUE) + 
             dgwesp(.25, TRUE, type = "OTP") + dgwesp(.25, TRUE, type = "ITP"),
             control = control.ergm(seed = 39478))
bTerms = c("density", "reciprocity", "relatedness", "territory overlap",
           "male receiving", "male sharing", "even-distribution receiving", 
           "even-distribution sharing", "transitivity", "cyclicality")
stargazer(mBinary, type = "text", single.row = TRUE, covariate.labels = bTerms)

# In-text odds:
bests = broom::tidy(mBinary)
# Change in odds if there's a reciprocal tie being formed:
exp(bests$estimate[bests$term == "mutual"])
# Change in odds if there's a first transitive triangle formed:
exp(bests$estimate[bests$term == "gwesp.OTP.fixed.0.25"])


mValued = ergm(binaryNet ~ sum + nonzero + mutual(form = "geometric") +
                 nodeifactor("sex") + nodeofactor("sex") + 
                 edgecov(overlap) + edgecov(related) +
                 transitiveweights("geomean","max","geomean") +
                 cyclicalweights("geomean","max","geomean") ,
               reference = ~ Poisson,
               response = "NumberTies",
               control = control.ergm(MCMC.prop.weights="0inflated",
                                      MCMC.prop.args=list(p0=0.5),
                                      seed = 5026))
vTerms = c("density", "any shares", "reciprocity", "male receiving", "male sharing",
           "territory overlap", "relatedness", "transitivity", "cyclicality")
stargazer(mValued, type = "text", single.row = TRUE, covariate.labels = vTerms)  

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

et = full_join(bSummary, vSummary, by = "term")
et = rbind(
  data.frame(term = "Behavioral Ecology Hypotheses", binary_edges = "", valued_edges = ""),
  et[c(3, 4, 5, 6), ],
  data.frame(term = "Social Behavior Hypotheses", binary_edges = "", valued_edges = ""),
  et[c(2, 9, 10), ],
  data.frame(term = "Basic Network Attributes", binary_edges = "", valued_edges = ""),
  et[c(1, 11, 7, 8), ]
)
et$term = stringi::stri_trans_totitle(et$term)
write_csv(et, "results/ergmTable.csv", na = " ")
# stargazer(type = "text", summary = FALSE, rownames = FALSE,
#           notes  = "*** p < .001, ** p < .01, * p < .05",
#           out = "results/modelsTable.txt") 

# Goodness of fit
bgof = gof(mBinary)
pdf("results/binaryGOF.pdf")
par(mfrow = c(2, 2))
plot(bgof)
dev.off()
