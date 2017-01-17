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

# Add a column to dyadAttr with shares in the opposite direction for descriptive stats
back = data.frame(cat1 = dyadAttr$cat2, cat2 = dyadAttr$cat1, invShares = dyadAttr$n)
dyadAttr = left_join(dyadAttr, back) 

da = select(dyadAttr, relatedness, propOverlap, n, invShares)
dyadCor = Hmisc::rcorr(as.matrix(da), type = "pearson")
dcr = dyadCor$r[2:1, 4:3]
dcp = dyadCor$P[2:1, 4:3]

# columns 4:3 because overlap is fraction of cat1's territory also occupied by cat2
# Below it's receiving:sharing
# And first entry here (with the 4:3 switch) is how much cat1 gets from cat2
# plot to make sure that's right:
select(dyadAttr, propOverlap, n, invShares) %>%
  gather(type, num, -propOverlap) %>%
  ggplot(aes(x = propOverlap, y = num, color = type)) +
  geom_point() + 
  geom_smooth(method = "lm")
# As more of cat1's territory is also occupied by cat2, the more cat1 shares with cat2, 
# but it's not the case that cat1 also receives more from cat2.

# Correlation of number of kills shared with relatedness and territory overlap:
round(cor(dyadAttr[, c(3, 4, 8)])[1:2, 3], 3)

# Nodal
net = readRDS("data/derived/network.RDS")
library(netUtils)
ndf = 
  sapply(c("indegree", "outdegree"), function(io) degree(net, gmode = "digraph", cmode = io)) %>%
  data.frame(vAttrDF(net), .) %>%
  mutate(male = as.integer(sex == "M"), female = as.integer(sex == "F")) %>%
  rename(sharing = outdegree, receiving = indegree, weight = weight_kg, age = age_months) %>%
  .[, sapply(., is.numeric)] %>%
  .[, c(3, 4, 1, 2, 5)]

nodeCor = Hmisc::rcorr(as.matrix(ndf), type = "pearson")
ncr = nodeCor$r[1:2, 3:5]
ncp = nodeCor$P[1:2, 3:5]

makeTab = function(r, p, rdig = 3, pdig = 2) {
  rows = nrow(r)
  cols = ncol(r)
  m = matrix(nrow = rows, ncol = cols)
  for(i in 1:rows)
    for(j in 1:cols)
      m[i, j] = paste0(round(r[i, j], rdig), " (p = ", round(p[i, j], pdig), ")")
  attr(m, "dimnames") = attr(r, "dimnames") 
  return(m)
}
cortab = rbind(t(makeTab(ncr, ncp)), makeTab(dcr, dcp))
ctdf = broom::tidy(cortab)
ctdf$.rownames = c("Age (months)", "Weight (kg)", "Male (vs. female)",
                   "Spatial overlap", "Relatedness")
colnames(ctdf) = c(" ", "Receiving (In-Degree)", "Sharing (Out-Degree)")
write_csv(ctdf, "results/corr_table.csv")


# In-text statistics
## Fraction of dyads with any tolerance
### In either direction
is.directed(net)
m = net[,]
dyads = m[upper.tri(m)] + t(m)[upper.tri(m)]
sum(dyads > 0) / length(dyads)
### Fraction of directed-dyads with tolerance
sum(m) / (length(m) - network.size(net))  # minus because no loops possible

# Max tolerance for a dyad
x = intergraph::asIgraph(net) %>%
  igraph::as_adj(sparse = FALSE)
max(x + t(x))
