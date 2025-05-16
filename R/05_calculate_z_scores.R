
# z-scores from observed and null networks
# note that metrics for the three null models
# were stored differently and therefore, are read differently here
# -------------------------------------------------------------------------

library(tidyverse)
# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")

# metrics from observed networks

observed.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

# -------------------------------------------------------------------------
# metrics from connectance null model networks
null.connectance.metrics.bipartite <- list.files("results/connectance_null_bipartite_metrics/",full.names = T) %>%
  map_dfr(read.csv2) %>%
  bind_rows() %>%
  # when calculating the null connectance replicates I stored them with "_connectance" as a suffix to network_id. Remove it
  mutate(network_id = substr(network_id,1,nchar(network_id)-12))

null.connectance.metrics.unipartite <- list.files("results/connectance_null_unipartite_metrics/",full.names = T) %>%
  map_dfr(read.csv2) %>%
  bind_rows() %>%
  # when calculating the null connectance replicates I stored them with "_connectance" as a suffix to network_id. Remove it
  mutate(network_id = substr(network_id,1,nchar(network_id)-12))

# -------------------------------------------------------------------------
# metrics from degree distribution null model networks
# 1 - soft-constrained
null.degdist.soft.metrics.bipartite <- read.csv2("results/soft_degdist_null_bipartite_metrics/bipartite_network_metrics_null.csv") %>%
  mutate(null_model = "degdist_soft")
null.degdist.soft.metrics.unipartite <- read.csv2("results/soft_degdist_null_unipartite_metrics/unipartite_network_metrics_null.csv") %>%
  mutate(null_model = "degdist_soft")

# 2 - hard-constrained
null.degdist.hard.metrics.bipartite <- list.files("results/hard_degdist_null_bipartite_metrics/",full.names = T) %>%
  map_dfr(read.csv2) %>%
  bind_rows() %>%
  # when calculating the null degdist replicates I stored them with "_degdist" as a suffix to network_id. Remove it
  mutate(network_id = substr(network_id,1,nchar(network_id)-8),
         null_model = "degdist_hard")

null.degdist.hard.metrics.unipartite <- list.files("results/hard_degdist_null_unipartite_metrics/",full.names = T) %>%
  map_dfr(read.csv2) %>%
  bind_rows() %>%
  # when calculating the null degdist replicates I stored them with "_degdist" as a suffix to network_id. Remove it
  mutate(network_id = substr(network_id,1,nchar(network_id)-8),
         null_model = "degdist_hard")
# -------------------------------------------------------------------------

# dropping na's because null metrics are obtained for some networks
# that were discarded
# this is because null metrics are taken from the adjacency matrices
# in the folders "data/uni-bipartite_adjacency_matrices", and these folders
# may contain matrices that are actually discarded. 
# it works, but I should clean that up
null.degdist.soft.metrics <- bind_rows(null.degdist.soft.metrics.bipartite,
                                       null.degdist.soft.metrics.unipartite) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
  drop_na()

null.degdist.hard.metrics <- bind_rows(null.degdist.hard.metrics.bipartite,
                                  null.degdist.hard.metrics.unipartite) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
  drop_na()

null.connectance.metrics <- bind_rows(null.connectance.metrics.bipartite,
                                      null.connectance.metrics.unipartite) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
  drop_na()

null.metrics <- bind_rows(null.connectance.metrics,null.degdist.soft.metrics,null.degdist.hard.metrics)

# -------------------------------------------------------------------------

null.means <- null.metrics %>%
  group_by(null_model,network_id,network_interaction,network_topology_type,metric) %>%
  summarise(null.mean = mean(value, na.rm = T),
            null.sd = sd(value,na.rm = T))
# 
z.scores.long <- left_join(observed.metrics,null.means) %>%
  # na.omit() %>%
  mutate(z_score = (value-null.mean)/null.sd) 
z.scores.long$z_score_full <- ifelse(z.scores.long$null.sd == 0,
                                     (z.scores.long$value-z.scores.long$null.mean),
                                     z.scores.long$z_score)

z.scores.long <- z.scores.long %>%
  dplyr::select(null_model,network_id,network_interaction,network_topology_type,
         metric,value,z_score_full) %>%
  drop_na()

# -------------------------------------------------------------------------

write.csv2(z.scores.long,"results/metrics_z_scores.csv",row.names = F)
