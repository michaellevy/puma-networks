library(statnet)
library(tidyverse)
library(netUtils)
net = readRDS("data/derived/network.RDS")

df = netUtils::vAttrDF(net)
df[grep("Unc", df$vertex.names), c(2, 4)] <- NA
df$outgoing_shares = degree(net, cmode = "outdegree")
df$incoming_shares = degree(net, cmode = "indegree")
df = 
  sapply(seq_len(network.size(net)), function(v) {
    es = get.edges(net, v, neighborhood = "combined")
    if(!length(es)) 
      return(c(0, 0))
    partners = 
      lapply(es, function(x) unlist(x[1:2])) %>%
      do.call(rbind, .)
    apply(partners, 2, function(x) length(unique(setdiff(x, v))))
  }) %>% 
  t() %>%
  cbind(df, .) %>%
  rename(cats_tolerated = inl, cats_tolerated_by = outl) %>%
  mutate(betweenness_centrality = betweenness(net),
         eigenvector_centrality = evcent(net))
write_csv(df, "results/cat_table.csv")
