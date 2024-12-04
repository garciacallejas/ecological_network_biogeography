
# z-scores from observed and null networks

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")

observed.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

null.metrics.bipartite <- read.csv2("results/bipartite_network_metrics_null.csv")
null.metrics.unipartite <- read.csv2("results/unipartite_network_metrics_null.csv")

# dropping na's because null metrics are obtained for some networks
# that were discarded
# this is because null metrics are taken from the adjacency matrices
# in the folders "data/uni-bipartite_adjacency_matrices", and these folders
# may contain matrices that are actually discarded. 
# TODO I should clean that up
null.metrics <- bind_rows(null.metrics.bipartite,null.metrics.unipartite) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
  drop_na()

# -------------------------------------------------------------------------

null.means <- null.metrics %>%
  group_by(network_id,network_interaction,network_topology_type,metric) %>%
  summarise(null.mean = mean(value, na.rm = T),
            null.sd = sd(value,na.rm = T))
# 
z.scores.long <- left_join(observed.metrics,null.means) %>%
  # na.omit() %>%
  mutate(z_score = (value-null.mean)/null.sd) %>%
  dplyr::select(network_id,network_interaction,network_topology_type,
         metric,value,z_score) %>%
  drop_na()

# -------------------------------------------------------------------------

write.csv2(z.scores.long,"results/metrics_z_scores.csv",row.names = F)
