# Function to calculate descriptive statistics for their own sake and CUG tests
calcStats = function(net, possTies = ft) {
  
  mod = igraph::cluster_walktrap(intergraph::asIgraph(net))
  
  sexMix = summary(net ~ nodemix("sex"))
  sexFraction = structure(sexMix / possTies, names = paste0("sex-", str_sub(names(sexMix), 9, 11)))
  
  tc = triad.census(net)
  tc = structure(as.vector(tc), names = paste0("tc-", colnames(tc)))
  # http://learning-raph.com/wp-content/uploads/2016/03/davis-leinhart-triad.png
  
  desc = c(density = network.density(net),
           reciprocity = summary(net ~ mutual + edges)[[1]] / summary(net ~ mutual + edges)[[2]],
           transitivity = gtrans(net),  # Weak form, unvalued. Other options for this
           clustering = clusteringCoef(net),  # Undirected
           modularity = igraph::modularity(intergraph::asIgraph(net), igraph::membership(mod)))
  
  return(
    list(descriptives = desc,
         sexMix = sexFraction,
         triadCensus = tc)
  )
}
