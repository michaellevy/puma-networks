# puma-networks
Prey sharing networks of mountain lions

![map with network](https://github.com/michaellevy/puma-networks/blob/master/results/mapWithNetwork.png?raw=true)

Pipeline, as I use it as I remove M62:

1. organizeData.R:  
- Read edgelist in first few rows -> write `edgelist.csv`
- Read other dyads and Elbroch's response to email -> write `dyadInfo.csv`

1. makeNetworks.R: 
- Read dyadInfo.csv and edgelist.csv from organizeData; attributes.csv from raw
- Clean by removing F59, M62
- Impute, create all possible pairings
- Write nodeInfo.csv
- Write network.RDS

1. catTable.R, write cat_table.csv which goes into paper

1. descriptives.R, manually export correlation with in- and out-degree table to paper. Calculate statistics for use in-text.

1. catTable.R: write cat_table.csv for paper

...

