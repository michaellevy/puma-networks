runGUGs = function(empNet, simNets) {
  
  if(is.multiplex(empNet) != is.multiplex(simNets[[1]]))
    stop("Multiplexity of empirical and simulated networks don't match.")
  
  simStats = lapply(simNets, calcStats)  
  
  ss = lapply(simStats, unlist, recursive = FALSE) %>% do.call(rbind, .)
  es = unlist(calcStats(empNet))
  
  cugs = data.frame(
    statistic = lapply(calcStats(empNet), names) %>% do.call(c, .),
    empirical_value = es,
    p_sim_greater = sapply(seq_along(es), function(i) sum(es[i] < ss[, i])) / nrow(ss),
    row.names = NULL)
  
  return(cugs[order(cugs$p_sim_greater), ])
  
}

# Function to calculate descriptive statistics for their own sake and CUG tests
calcStats = function(net) {
  
  # A few stats that don't make sense for valued: tc, recip, trans, clust
  mp = is.multiplex(net)
  if(mp) {
    mNet = net
    net = collapseMultiNet(net) 
  }
  
  # Fraction of dyads that share at least once by combo of sexes
  st = table(net %v% "sex")  
  ft = c(st[1] * (st[1] - 1), rep(st[1] * st[2], 2), st[2] * (st[2] - 1))
  sexMix = summary(net ~ nodemix("sex"))
  sexFraction = structure(sexMix / ft, names = paste0("sex-", str_sub(names(sexMix), 9, 11)))
  
  tc = triad.census(net)
  tc = structure(as.vector(tc), names = paste0("tc-", colnames(tc)))
  # http://learning-raph.com/wp-content/uploads/2016/03/davis-leinhart-triad.png

  mod = igraph::cluster_walktrap(intergraph::asIgraph(net))
  
  desc = c(inCentralization = centralization(net, "degree", cmode = "indegree"),
           outCentralization = centralization(net, "degree", cmode = "outdegree"),
           reciprocity = unname(summary(net ~ mutual)) / network.edgecount(net),
           transitivity = gtrans(net),  # Weak form, unvalued. Other options for this
           modularity = igraph::modularity(intergraph::asIgraph(net), igraph::membership(mod)),
           clustering = clusteringCoef(net))  # Undirected, unvalued
  
  if(mp) {
    mod = igraph::cluster_walktrap(intergraph::asIgraph(mNet))
    
    vDesc = c(VinCentralization = centralization(mNet, "degree", cmode = "indegree"),
              VoutCentralization = centralization(mNet, "degree", cmode = "outdegree"),
              Vreciprocity = unname(summary(net ~ mutual(form = "geometric"), response = "NumberTies")),
              Vtransitivity = unname(summary(net ~ transitiveweights("geomean", "max", "geomean"), response = "NumberTies")),
              Vmodularity = igraph::modularity(intergraph::asIgraph(mNet), igraph::membership(mod)),
              edgeVar = var(c(net %e% "NumberTies", rep(0, network.size(net) * (network.size(net) - 1) - network.edgecount(net)))))
    desc = c(desc, vDesc)
  }
  
  return(
    list(descriptives = desc,
         sexMix = sexFraction,
         triadCensus = tc)
  )
}
