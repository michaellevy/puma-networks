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

# Find communities and merge with above data frame:
g = intergraph::asIgraph(net)
# comm = cluster_walktrap(g)  # or, yielding identical communities:
comm = igraph::cluster_edge_betweenness(g)
comms = 
  purrr::map_df(names(igraph::communities(comm)), function(x) 
    data.frame(vertex.names = igraph::V(g)$vertex.names[igraph::communities(comm)[[x]]],
               community = x)) 

df = left_join(df, comms)

# Sort by community, sex, and weight, in that order, 
# Round the decimal values,
# Sort and rename columns, 
# and write to file:
df %>%
  arrange(community, 
          desc(stringr::str_extract(vertex.names, "M|F")),
          weight_kg) %>% 
  mutate(
    weight_kg = round(weight_kg, 0),
    betweenness_centrality = round(betweenness_centrality, 0),
    eigenvector_centrality = round(eigenvector_centrality, 2)
  ) %>%
  select(
    `Puma ID` = vertex.names,
    `Age (months)` = age_months,
    `Weight (kg)` = weight_kg,
    Community = community,
    `Outgoing tolerance` = outgoing_shares,
    `Pumas tolerated` = cats_tolerated,
    `Receiving tolerance` = incoming_shares,
    `Pumas that exhibited tolerance to` = cats_tolerated_by,
    `Betweenness centrality` = betweenness_centrality,
    `Eigenvector centrality` = eigenvector_centrality
  ) %>%
  write_csv("results/cat_table.csv")
