library(statnet)
library(readxl)
library(dplyr)  # For convenience functions, esp the pipe "%>%" which takes the output of what's on the left and makes it the input to what's on the right/on the next line
# install.packages(integraph)
d = read_excel("Data_MLubell1_Clean2.xlsx", na = "NODATA")
d = d[, c(5:6, 1:4, 7:ncol(d))]  # First two columns must be the edgelist
n = intergraph::asNetwork(igraph::graph.data.frame(d)) #double colons call function from a package without loading it; intergraph 
# Make sure it looks right:
degree=degree(n, gmode="digraph")
communities<-igraph::cluster_edge_betweenness(intergraph::asIgraph(n)) #call igraph functions to get cluster by edge betweenness e.g. Girvan Newman
get.vertex.attribute(n,"vertex.names")[communities[[1]]] #allowsyou to examine a subset of vertices as identified by communities
colors<-rep("red", network.size(n)) #create vector of colors, 13 entries corresponding to number nodes in network
colors[communities[[2]]]<-"blue"
plot.network(n, displaylabels = TRUE, edge.label = "Genetic.relatedness", vertex.cex=log(degree), vertex.col=colors)

dyads = read_excel("Data_MLubell1_Clean2.xlsx", na = "x", sheet = "All possible combos")
dyads = dyads[-1, c(1, 2, 5, 6)]  # Remove the second row and keep only meaningful columns
dyads = dyads[!apply(dyads, 1, function(x) sum(is.na(x)) == length(x)), ]  # Remove rows at the bottom that are totally empty
nrow(dyads)
13 * 12 / 2  # Possible dyads. Hmm, why the extra row?
# Create columns for each cat:
dyads = 
  strsplit(dyads$`Cat Pair`, "_") %>%  # Split the pairs into component cats
  do.call(rbind, .) %>%  # Each row is returned as an item in a list. Bind them together as rows.
  data.frame(stringsAsFactors = FALSE) %>%  # Make it a data.frame, keeping from being factors
  structure(names = c("cat1", "cat2")) %>%   # Give it column names
  cbind(dyads, .)  # Column bind it to the dyads data.frame
# Which cats appear in second sheet but not the first?
isos = unique(c(dyads$cat1, dyads$cat2))[!unique(c(dyads$cat1, dyads$cat2)) %in% get.vertex.attribute(n, "vertex.names")]
# Hmm, four cats appear in the second sheet but not the network, add them as isolates:
# Note also that statnet modifies objects "in place". We don't assign this back to n, in contrast to how R works. 
add.vertices(n, nv = length(isos), vattr = lapply(isos, function(x) list("vertex.names" = x)))
structure(degree(n), names = n %v% "vertex.names")  # Check that worked right

# Okay, let's get all the dyads in an edgelist form, whether they're in the network or not:
e1 = d[, c(1, 2, 7, 10)]  # Just the four columns of interest from the first sheet
e1 = e1[!duplicated(e1), ]   # Just the not-duplicated rows
names(e1) = c("cat1", "cat2", "relatedness", "propOverlap")
# From second sheet want in the order given and in the other order:
e2a = structure(dyads[, c(5, 6, 2, 3)],  # Cat1 first and proportion of its overlap
                names = c("cat1", "cat2", "relatedness", "propOverlap"))
e2b = structure(dyads[, c(6, 5, 2, 4)],  # Cat2 first and proportion of its overlap
                names = c("cat1", "cat2", "relatedness", "propOverlap"))
allDyads = rbind(e1, e2a, e2b)  # Row-bind them together
sum(duplicated(allDyads[, 1:2]))  # Any duplicates?
sum(duplicated(allDyads[, ]))  # Do they all have the same info in columns 3 or 4? No, shit.
# Find duplicate dyads with different relatedness or overlap:
problemDyads = 
  group_by(allDyads, cat1, cat2) %>%
  summarise(numRelatednesses = n_distinct(relatedness),
            numOverlaps = n_distinct(propOverlap)) %>%
  filter(numRelatednesses > 1 | numOverlaps > 1) %>%
  semi_join(allDyads, .) %>%
  as.data.frame()

network.size(n) * (network.size(n) - 1)  # Number of possible dyads
nrow(allDyads[!duplicated(allDyads[, 1:2]), ])  # Number of dyads we have relatedness and range info for
# Accounting for the two unknown cats, there are still 210 possible dyads
# Find the dyads we don't have:
possibleDyads = 
  expand.grid(cat1 = n %v% "vertex.names", cat2 = n %v% "vertex.names", stringsAsFactors = FALSE) %>%
  filter(cat1 != cat2)

anti_join(possibleDyads, allDyads) %>%  # Filter possibleDyads to those not in allDyads
  filter(!cat1 %in% c("UncF", "UncF2") & !cat2 %in% c("UncF", "UncF2")) %>%  # Get rid of the unknowns
  rbind(problemDyads[!duplicated(problemDyads[, 1:2]), 1:2]) %>%  # Row bind these with the ones with conflicting info
  arrange(cat1, cat2)  %>%
  write.csv("unknownPairs.csv", row.names = FALSE)