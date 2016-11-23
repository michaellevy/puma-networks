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

dyadCor = cor(select(dyadAttr, relatedness, propOverlap, n, invShares))[2:1, 4:3]
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

cors = cor(ndf)
prunedcors = cors[1:2, 2:5]
prunedcors[2, 1] = NA

(ct = round(rbind(t(prunedcors[, 2:4]), dyadCor), 3))
htmlTable::htmlTable(ct)

corrplot::corrplot(cors, method = "ellipse", type = "lower", 
                   diag = F, tl.col = "black", addCoef.col = "black")
