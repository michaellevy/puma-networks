library(statnet)
library(readxl)
# install.packages(integraph)
d = read_excel("data/raw/Data_MLubell1_Clean2.xlsx", na = "NODATA")
d = d[, c(5:6, 1:4, 7:ncol(d))]  # First two columns must be the edgelist
n = intergraph::asNetwork(igraph::graph.data.frame(d)) #double colons call function from a package without loading it; intergraph 
# Make sure it looks right:
degree=degree(n, gmode="digraph")
communities<-igraph::cluster_edge_betweenness(intergraph::asIgraph(n)) #call igraph functions to get cluster by edge betweenness e.g. Girvan Newman
get.vertex.attribute(n,"vertex.names")[communities[[1]]] #allowsyou to examine a subset of vertices as identified by communities
colors<-rep("red", network.size(n)) #create vector of colors, 13 entries corresponding to number nodes in network
colors[communities[[2]]]<-"blue"
sex = substr(n %v% "vertex.names", 1, 1)  # Get the sex from the first letter of the vertex names
sex[sex == "U"] = "F"  # Make the UncF female
plot.network(n, displaylabels = TRUE, edge.label = "Genetic.relatedness", 
             vertex.cex=log(degree), vertex.col=colors, 
             vertex.sides = c(4, 50)[as.factor(sex)])



# Note M62 seems to be an isolate and needs to be added manually. This code is overkill, but leaving it in case it's useful
# Hmm, four cats appear in the second sheet but not the network, add them as isolates:
# Note also that statnet modifies objects "in place". We don't assign this back to n, in contrast to how R works. 
# add.vertices(n, nv = length(isos), vattr = lapply(isos, function(x) list("vertex.names" = x)))
# structure(degree(n), names = n %v% "vertex.names")  # Check that worked right
