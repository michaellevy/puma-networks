library(tidyverse)
library(stringr)
library(network)
# source(Code/organizeData.R) if not done already

nodeAttr = read_csv("data/raw/attributes.csv")
dyadAttr = read_csv("data/derived/dyadInfo.csv")
el = read_csv("data/derived/edgelist.csv")

# make a sex variable, and add unknown females
nodeAttr$sex = str_sub(nodeAttr$cat, 1, 1)
nodeAttr = 
  rbind(nodeAttr, data.frame(cat = c("UncF", "UncF2"),
                             weight_kg = NA,
                             age_months = NA,
                             sex = "F"))

########## Mean impute missing values
# Age use overall mean; weight use mean among females only.
nodeAttr$weight_kg[is.na(nodeAttr$weight_kg)] = mean(nodeAttr$weight_kg[nodeAttr$sex == "F"], na.rm = TRUE)
nodeAttr$age_months[is.na(nodeAttr$age_months)] = mean(nodeAttr$age_months, na.rm = TRUE)
write_csv(nodeAttr, "data/derived/nodeInfo.csv")


# Create all the possile pairings for the UncF's, in both directions
uncNames = c("UncF", "UncF2")
pairs = expand.grid(cat1 = uncNames,
                    cat2 = unique(c(dyadAttr$cat1, dyadAttr$cat2)),
                    stringsAsFactors = FALSE) %>%
  rbind(uncNames) %>%
  mutate(relatedness = mean(dyadAttr$relatedness),
         propOverlap = mean(dyadAttr$propOverlap),
         sex1 = "F",
         sex2 = str_sub(cat2, 1, 1),
         sex2 = ifelse(sex2 == "U", "F", sex2),
         sexPair = paste0(sex2, sex1))
dyadAttr = 
  rbind(pairs,
        select(pairs, cat1 = cat2, cat2 = cat1, relatedness, propOverlap, 
               sex1 = sex2, sex2 = sex1, sexPair),
        dyadAttr) 

# Save network objects
pumanet = 
  igraph::graph_from_data_frame(el, directed = TRUE, vertices = nodeAttr) %>%
  intergraph::asNetwork()
# Remove F59 since the data on her was spotty
delete.vertices(pumanet, which(pumanet %v% "vertex.names" == "F59"))
saveRDS(pumanet, "data/derived/network.RDS")
