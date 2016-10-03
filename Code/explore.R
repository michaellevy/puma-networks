library(tidyverse)
theme_set(theme_bw())
library(stringr)
library(statnet)
library(netUtils)
library(stargazer)
el = read_csv("data/derived/edgelist.csv")
dyadAttr = read_csv("data/derived/dyadInfo.csv")
nodeAttr = read_csv("data/derived/nodeInfo.csv")
pumanet = readRDS("data/derived/network.RDS")
collapse_pumanet = readRDS("data/derived/simpleNet.RDS")

ggplot(dyadAttr, aes(x = relatedness, y = propOverlap, color = sexPair)) + 
  geom_point()

ggplot(nodeAttr, aes(x = age_months, y = weight_kg, color = sex)) + 
  geom_point(size = 3)


vSize = 3 * pumanet %v% "weight_kg" / max(pumanet %v% "weight_kg", na.rm = TRUE)
vSize[is.na(vSize)] = .1

plot(pumanet
     , vertex.col = "sex"
     , vertex.cex = vSize
     , displaylabels = TRUE
     , edge.lwd = 2)

count(el, cat1, cat2) %>%
  ggplot(aes(x=n)) + 
  geom_histogram(bins = 11) +
  xlab("Number of kills shared (per directed dyad)")


num_shares=count(el, cat1, cat2)
dyadAttr=left_join(dyadAttr,num_shares)
dyadAttr$n[is.na(dyadAttr$n)]=0
names(dyadAttr)[8]="NumShares"

#correlation matrix to see if genetic relatedness and range overlap covary
round(cor(dyadAttr[,c("NumShares", "relatedness", "propOverlap")]), 2)

#regression number of shares in dyadic unit analysis to see basic relationships
share_regression<-lm(NumShares~propOverlap+relatedness+sexPair, dyadAttr)
summary(share_regression)

### Work with the individual edges with number of kills-shared as edge attr
wt = collapse_pumanet %v% "weight_kg"
wt[is.na(wt)] = mean(wt, na.rm = TRUE)
age = collapse_pumanet %v% "age_months"
age[is.na(age)] = mean(age, na.rm = TRUE)

plot(collapse_pumanet
     , vertex.col = "sex"
     , arrowhead.cex = 2
     , vertex.cex = age / 40
     , displaylabels = TRUE
     , label.cex = .85
     , edge.lwd = "n"
     , usecurve = TRUE
     , edge.curve = .02
     )

# Copy network and delete uncollared F's to see how the descriptives change
trim_pumanet = collapse_pumanet
missingF = c("UncF", "UncF2")
delete.vertices(trim_pumanet, which(trim_pumanet %v% "vertex.names" %in% missingF))
nets = list(full = collapse_pumanet, trimmed = trim_pumanet)

par(mfrow = c(1, 2))
lapply(nets, function(x)
  plot(x
       , vertex.col = "sex"
       , vertex.cex = 3
       , displaylabels = TRUE
       , edge.lwd = "n")
)

#basic network stats: density for valued, reciprocity, transitivity, clustering, modulatrity for community detection, sex homophily, triad census
sapply(nets, network.density)
sapply(nets, function(x) 
  summary(x ~ mutual + edges)[1] / summary(x ~ mutual + edges)[2])
sapply(nets, clusteringCoef)
lapply(nets, function(x) summary(x ~ nodematch("sex", diff = TRUE)))
lapply(nets, triad.census) %>% do.call(rbind, .)
# http://learning-raph.com/wp-content/uploads/2016/03/davis-leinhart-triad.png

#ERGM--range overlap, relateness, sex pairs, node age, node weight, reciprocity, transitivity, GWDEGREE(hierarchy?)
m1 = ergm(collapse_pumanet ~ edges + mutual)
m1t = ergm(trim_pumanet ~ edges + mutual)
stargazer(m1, m1t, type = "text")  # Good, little difference there at least.

related = spread(dyadAttr[, 1:3], cat1, relatedness)[, -1] %>% as.matrix()
overlap = spread(dyadAttr[, c(1, 2, 4)], cat1, propOverlap)[, -1] %>% as.matrix()
diag(related) = diag(overlap) = 1
# Need to order those by vertex.names
ord = match(trim_pumanet %v% "vertex.names", colnames(related))
related = related[ord, ord]
overlap = overlap[ord, ord]

m2 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap))
summary(m2)

m3 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
          nodeocov('age_months') + nodeicov('age_months'))
stargazer(m1t, m2, m3, type = "text")

m4 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeocov('weight_kg') + nodeicov('weight_kg'))
stargazer(m1t, m2, m4, type = "text")
# Interesting: Bigger cats give less and get more. How does sex play with that:

m5 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeifactor("sex") + nodeofactor("sex"))

m6 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeocov('weight_kg') + nodeicov('weight_kg') + 
            nodeifactor("sex") + nodeofactor("sex"))

stargazer(m1t, m2, m4, m5, m6, type = "text")
# Alternatively, males give less and take a little more. Highly collinear so cancel.

m7 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE))

# m8 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
#             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = TRUE))
# computationally singular -- insufficient data for this.
stargazer(m1t, m2, m5, m7, type = "text")

m8 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeocov('weight_kg') + nodeicov('weight_kg') + nodematch("sex", diff = FALSE))
stargazer(m7, m8, type = "text")  # Makes no difference

m9 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
            gwidegree(.25, fixed = TRUE))

m10 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
            nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
            gwidegree(.25, fixed = TRUE) + gwodegree(.25, fixed = TRUE))
stargazer(m10, type = "text", single.row = TRUE)  # Makes no difference

# m11 = ergm(trim_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
#              nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
#              ctriple)
# Both the triples don't mix well, even with nothing else in there:
# m11 = ergm(trim_pumanet ~ edges + mutual + ctriple)


####### Imputation was here some data are not going to be loaded here And I 
# did a comparison of trimmed vs imputed data. Now the mean-imputed values are there
# from the beginning.


related = spread(dyadAttr[, 1:3], cat1, relatedness)[, -1] %>% as.matrix()
ord = match(collapse_pumanet %v% "vertex.names", colnames(related))
related = related[ord, ord]
overlap = spread(dyadAttr[, c(1, 2, 4)], cat1, propOverlap)[, -1] %>% as.matrix()
overlap = overlap[ord, ord]
diag(related) = diag(overlap) = 1

m10Imp = ergm(collapse_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
             gwidegree(.25, fixed = TRUE) + gwodegree(.25, fixed = TRUE))
stargazer(m10, m10Imp, type = "text", single.row = TRUE)

m11 = ergm(collapse_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
             gwidegree(.25, fixed = TRUE) + gwodegree(.25, fixed = TRUE) + gwesp(.25, TRUE),
           control = control.ergm(init = c(coef(m10Imp), 0)))

m12 = ergm(collapse_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = FALSE) +
             gwidegree(.25, fixed = TRUE) + gwodegree(.25, fixed = TRUE) + 
             dgwesp(.25, TRUE, type = "OTP") + dgwesp(.25, TRUE, type = "ITP"),
           control = control.ergm(init = c(coef(m10Imp), 0, 0)))

# Try paring that down (and change bump decay)
m13 = ergm(collapse_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + 
             gwidegree(.5, fixed = TRUE) + gwodegree(.5, fixed = TRUE) + 
             dgwesp(.25, TRUE, type = "OTP") + dgwesp(.25, TRUE, type = "ITP"),
           control = control.ergm(init = coef(m12)[c(1:6, 8:11)]))

m14 = ergm(collapse_pumanet ~ edges + mutual + edgecov(related) + edgecov(overlap) +
             nodeifactor("sex") + nodeofactor("sex") + 
             gwidegree(.5, fixed = TRUE) + gwodegree(.5, fixed = TRUE) + 
             dgwesp(.25, TRUE, type = "OTP"),
           control = control.ergm(init = coef(m13)[-length(coef(m13))]))

stargazer(m10Imp, m11, m12, m13, m14, type = "text", single.row = TRUE)

bterms = c("density", "reciprocity", "relatedness", "territory overlap",
           "male receiving", "male sharing", "even-distribution receiving", 
           "even-distribution sharing", "transitivity", "cyclicality")
stargazer(m13, type = "text", single.row = TRUE, covariate.labels = bterms)

# Valued edges
as.matrix(collapse_pumanet, attrname = "n")
summary(collapse_pumanet ~ sum,
        reference = ~ Geometric,
        response = "n")
m = sum(collapse_pumanet %e% "n") / network.dyadcount(collapse_pumanet)
initSum = log(1 - 1 / (m + 1))
v1 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "geometric"),
     reference = ~ Geometric,
     response = "n",
     control = control.ergm(init = c(initSum, 0, 0)))

v2 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min"),
          reference = ~ Geometric,
          response = "n",
          control = control.ergm(init = c(initSum, 0, 0)))
stargazer(v1, v2, type = "text")

v3 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") + 
            nodeisqrtcovar,
         reference = ~ Geometric,
         response = "n",
         control = control.ergm(init = c(coef(v2), 0)))
stargazer(v1, v2, v3, type = "text")  # No in-popularity effect

# v4 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") + 
#             transitiveweights("min","max","min") +
#             cyclicalweights("min","max","min"),
#           reference = ~ Geometric,
#           response = "n",
#           control = control.ergm(init = c(coef(v2), 0, 0)))
# stargazer(v1, v2, v3, v4, type = "text")
# Computationally singular

# v4 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
#             transitiveweights("min","max","min"),
#           reference = ~ Geometric,
#           response = "n",
#           control = control.ergm(init = c(coef(v2), 0)))
# Error in eigen(crossprod(x1c), symmetric = TRUE) : 
#   infinite or missing values in 'x' 

# Got ahead of myself. Add in non-social effects
v5 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
            nodeifactor("sex") + nodeofactor("sex") + nodematch("sex"),
          reference = ~ Geometric,
          response = "n",
          control = control.ergm(init = c(coef(v2), 0, 0, 0)))
stargazer(v1, v2, v5, type = "text")  # Males give less, not clear they get more

# v6 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
#             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex", diff = TRUE),
#           reference = ~ Geometric,
#           response = "n",
#           control = control.ergm(init = c(coef(v2), 0, 0, 0, 0)))
# Singular

v6 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
            nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
            edgecov(overlap) + edgecov(related),
          reference = ~ Geometric,
          response = "n",
          control = control.ergm(init = c(coef(v5), 0, 0)))
stargazer(v1, v2, v5, v6, type = "text")
# Getting somewhere: Zero-inflation, mutuality, sex-bias, range-overlap

# Previous estimation of v7 gave trans term of .47 but hessian approx singular. 
# Starting with those coefs and bumping mcmc failed
# Try reestimating with controls jacked up.
system.time({
  v7 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
              nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
              edgecov(overlap) + edgecov(related) +
              transitiveweights("min","max","min") ,
            reference = ~ Geometric,
            response = "n",
            control = control.ergm(init = c(coef(v6), 0),
                                   MCMLE.trustregion = 100,  # Not sure if this is appropriate for geometric dist
                                   MCMC.samplesize = 5e3,
                                   MCMLE.maxit = 50
            ))
})
stargazer(v1, v2, v5, v6, v7, type = "text")
# Didn't converge after 50 iterations. Lots of transitivity but huge negative nonzero wo SE
# Try w/o nonzero term and with trustregion at default

system.time({  # 9 min
  v8 = ergm(collapse_pumanet ~ sum + mutual(form = "min") +
              nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
              edgecov(overlap) + edgecov(related) +
              transitiveweights("min","max","min") ,
            reference = ~ Geometric,
            response = "n",
            control = control.ergm(init = c(coef(v6)[-2], 0),
                                   MCMLE.trustregion = 100,  # Not sure if this is appropriate for geometric dist
                                   MCMC.samplesize = 5e3,
                                   MCMLE.maxit = 50
            ))
})
saveRDS(v8, "models/count-ergm-trans.RDS")
# Negative transativity! But it falls out when cyclic triads are included

stargazer(v1, v2, v5, v6, v8, type = "text")

terms = c("density", "reciprocity", "male-receiving", "male-sharing", "same-sex",
          "range-overlap", "relatedness", "transitivity")
stargazer(v8, type = "text", single.row = TRUE, covariate.labels = terms)


system.time({  # 90 sec
  v9 = ergm(collapse_pumanet ~ sum + mutual(form = "min") +
              nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
              edgecov(overlap) + edgecov(related) +
              transitiveweights("min","max","min") ,
            reference = ~ Geometric,
            response = "n",
            control = control.ergm(init = c(coef(v6)[-2], 0),
                                   # MCMLE.trustregion = 100,  # Not sure if this is appropriate for geometric dist
                                   # MCMC.samplesize = 5e3,
                                   # MCMLE.maxit = 50
            ))
})
stargazer(v8, v9, type= "text")  # Controls were unnecessary. The nonzero terms was just getting in the way

v10 = ergm(collapse_pumanet ~ sum + mutual(form = "min") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("min","max","min") +
             cyclicalweights("min","max","min") ,
           reference = ~ Geometric,
           response = "n",
           control = control.ergm(init = c(coef(v8), 0)))

v11 = ergm(collapse_pumanet ~ sum + mutual(form = "min") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Geometric,
           response = "n",
           control = control.ergm(init = c(coef(v8), 0)))

v12 = ergm(collapse_pumanet ~ sum + mutual(form = "geometric") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Geometric,
           response = "n",
           control = control.ergm(init = coef(v11)))
stargazer(v8, v10, v11, v12, type = "text")  

v13 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "geometric") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Geometric,
           response = "n",
           control = control.ergm(init = c(coef(v12)[1], -2, coef(v12)[2:length(coef(v12))])))
# Error in eigen(crossprod(x1c), symmetric = TRUE) : 
#   infinite or missing values in 'x' 
### Maybe should be using 0-infl pois
### Try first for model 12 which was converging

v12p = ergm(collapse_pumanet ~ sum + mutual(form = "geometric") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Poisson,
           response = "n",
           control = control.ergm(init = coef(v12)))
stargazer(v12, v12p, type = "text")  # Pois looks better resolved. Try nonzero in there

# For comparison
v6p = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related),
           reference = ~ Poisson,
           response = "n",
           control = control.ergm(init = coef(v6),
                                  MCMLE.trustregion=100,
                                  MCMC.prop.weights="0inflated",
                                  MCMC.prop.args=list(p0=0.5)))
# Super efficient sampling!
stargazer(v6, v6p, type = "text")  # Very similar story

v13 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "geometric") +
             nodeifactor("sex") + nodeofactor("sex") + nodematch("sex") +
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Poisson,
           response = "n",
           control = control.ergm(init = c(coef(v12)[1], -2, coef(v12)[2:length(coef(v12))])),
           MCMLE.trustregion=100,
           MCMC.prop.weights="0inflated",
           MCMC.prop.args=list(p0=0.5))
summary(v13)
# non-convergence. Tried longer chains it stalled out (likelihood didn't improve x 40)

# Try trimming it down a little.  Want transitivity and nonzero
v13 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "min") +
             nodeifactor("sex") + nodeofactor("sex") + 
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Poisson,
           response = "n",
           control = control.ergm(init = c(coef(v6p)[-6], 0, 0),
                                  # MCMLE.trustregion=100,
                                  MCMC.prop.weights="0inflated",
                                  MCMC.prop.args=list(p0=0.5)))
terms = c("density", "any shares", "reciprocity", "male receiving", "male sharing",
          "territory overlap", "relatedness", "transitivity", "cyclicality")
stargazer(v13, type = "text", single.row = TRUE, covariate.labels = terms)

v14 = ergm(collapse_pumanet ~ sum + nonzero + mutual(form = "geometric") +
             nodeifactor("sex") + nodeofactor("sex") + 
             edgecov(overlap) + edgecov(related) +
             transitiveweights("geomean","max","geomean") +
             cyclicalweights("geomean","max","geomean") ,
           reference = ~ Poisson,
           response = "n",
           control = control.ergm(init = coef(v13),
                                  MCMC.prop.weights="0inflated",
                                  MCMC.prop.args=list(p0=0.5)))
stargazer(v14, type = "text", single.row = TRUE, covariate.labels = terms)
# Winner! Right now, the two leading models are v14 and m13. Go ahead and save all.
ergms = Filter(function(x) inherits(get(x), "ergm"), ls())
save(list = ergms, file = "models/exploratoryERGMs.RDA")