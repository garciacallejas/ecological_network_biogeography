
library(tidyverse)
library(colorblindr)
library(patchwork)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
# -------------------------------------------------------------------------
# comment dataframes that I do not currently use

net.collection <- read_csv2("results/network_collection.csv")

# metrics from observed networks

network.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

# -------------------------------------------------------------------------
# metrics from degree distribution null model networks
null.metrics.bipartite <- read.csv2("results/bipartite_network_metrics_null.csv")
# null.metrics.unipartite <- read.csv2("results/unipartite_network_metrics_null.csv")

# -------------------------------------------------------------------------
# metrics from connectance null model networks
# null.connectance.metrics.bipartite <- list.files("results/metrics_bipartite_connectance_null/",full.names = T) %>%
#   map_dfr(read.csv2) %>%
#   bind_rows() %>%
#   # when calculating the null connectance replicates I stored them with "_connectance" as a suffix to network_id. Remove it
#   mutate(network_id = substr(network_id,1,nchar(network_id)-12))

# null.connectance.metrics.unipartite <- list.files("results/metrics_unipartite_connectance_null/",full.names = T) %>%
#   map_dfr(read.csv2) %>%
#   bind_rows() %>%
#   # when calculating the null connectance replicates I stored them with "_connectance" as a suffix to network_id. Remove it
#   mutate(network_id = substr(network_id,1,nchar(network_id)-12))

# -------------------------------------------------------------------------

# null.degdist.metrics <- bind_rows(null.metrics.bipartite,null.metrics.unipartite) %>%
#   left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
#   drop_na()
# 
# null.connectance.metrics <- bind_rows(null.connectance.metrics.bipartite,
#                                       null.connectance.metrics.unipartite) %>%
#   left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
#   drop_na()
# 
# null.metrics <- bind_rows(null.degdist.metrics,null.connectance.metrics)

# network.metrics.long <- read.csv2("results/metrics_observed_null_long.csv")

# TEST------------------------
# metrics.zscores <- read.csv2("results/metrics_z_scores.csv")
metrics.zscores <- read.csv2("results/metrics_z_scores_TEST.csv")
# TEST------------------------

metrics.zscores$study_common <- net.collection$study_common[match(metrics.zscores$network_id,net.collection$network_id)]

# metrics.zscores.clean <- subset(metrics.zscores, !is.infinite(z_score) & abs(z_score) < 10)
# metrics.zscores.clean <- subset(metrics.zscores, !is.infinite(z_score))

# -------------------------------------------------------------------------
# I could explore drawing observed values and distribution of nulls, but
# there are 80 combinations of 10 metrics * 2 topologies * 4 interactions,
# and each of those is not even a single network, but potentially many
# so it starts to be very complex.

# For now, plot the z-scores for these combinations

# -------------------------------------------------------------------------
# zscores.plot.full.data <- subset(metrics.zscores.clean, !(metric %in% c("richness","degree.shannon","degree.mean",
#                                                                          "degree.sd","degree.skewness",
#                                                                   "degree.kurtosis","connectance","link.density")))
zscores.plot.full.data <- subset(metrics.zscores, !(metric %in% c("richness","link.density","centrality.degree")))
# zscores.plot.full.data <- metrics.zscores

zscores.plot.full.data$metric[zscores.plot.full.data$metric == "modularity.infomap"] <- "infomap \nmodularity"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
# zscores.plot.full.data$metric[zscores.plot.full.data$metric == "link.density"] <- "link \ndensity"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "interaction.overlap"] <- "interaction \noverlap"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.eigen"] <- "eigenvector \ncentralization"
# zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.degree"] <- "degree \ncentralization"
# zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.betweenness"] <- "betweenness \ncentrality"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "nestedness"] <- "NODF nestedness"

# zscores.plot.full.data$metric <- factor(zscores.plot.full.data$metric, levels = c(#"connectance","link \ndensity",
#                                                                                   "interaction \noverlap",
#                                                                                   "NODF nestedness",
#                                                                                   "betweenness \nmodularity",
#                                                                                   "infomap \nmodularity",
#                                                                                   # "betweenness \ncentrality",
#                                                                                   "degree \ncentralization",
#                                                                                   "eigenvector \ncentralization"
#                                                                                   ))

zscores.plot.full.data$metric <- factor(zscores.plot.full.data$metric, levels = c(#"connectance","link \ndensity",
  #"richness",
  "connectance","degree.mean","degree.sd","degree.skewness","degree.kurtosis",
  "degree.shannon",#"link.density",
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  #"degree \ncentralization",
  "eigenvector \ncentralization"
))
zscores.plot.full.data$z_score <- zscores.plot.full.data$z_score_full
degdist.metrics <- c("richness","connectance","degree.mean","degree.sd","degree.skewness","degree.kurtosis",
                     "degree.shannon","link.density")

median.values <- zscores.plot.full.data %>% group_by(null_model,metric,network_interaction) %>%
  summarise(median.value = median(z_score_full,na.rm = T))

# -------------------------------------------------------------------------

zscores.connectance.plot <- ggplot(subset(zscores.plot.full.data, null_model == "connectance"),
                                   # zscores.connectance.plot <- ggplot(subset(zscores.plot.full.data, null_model == "connectance" & metric %in% degdist.metrics), 
                                   aes(x = network_interaction, y = z_score)) + 
  geom_boxplot(aes(fill = network_interaction),
               outlier.colour = "grey30",outlier.size = .7) + 
  # geom_text(data = pair.comp.labels, aes(label = group, y = ypos, x = network_interaction),
  #           position = position_dodge(2),
  #           color = "grey",
  #           show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black") +
  facet_grid(cols = vars(metric)) +
  # facet_grid(cols = vars(metric), rows = vars(null_model)) +
  theme_bw() +
  # coord_cartesian(ylim = c(-6.5,10)) +
  coord_cartesian(ylim = c(-15,15)) +
  # ylim(-6.5,6) +
  # coord_flip() +
  # scale_x_discrete(limits = rev(levels(zscores.plot.data$metric))) +
  labs(x = "", y = "z-score") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),legend.position = "bottom")+
  # theme(legend.justification=c(1,-0.01), legend.position=c(0.99,0)) +
  ggtitle("Hard-constrained connectance") +
  NULL
# zscores.connectance.plot

zscores.degdist.soft.plot <- ggplot(subset(zscores.plot.full.data, null_model == "degdist_soft"),
                               aes(x = network_interaction, y = z_score)) +
  # zscores.degdist.soft.plot <- ggplot(subset(zscores.plot.full.data, null_model == "degdist_soft" & metric %in% degdist.metrics), 
  #                                     aes(x = network_interaction, y = z_score)) + 
  geom_boxplot(aes(fill = network_interaction),
               outlier.colour = "grey30",outlier.size = .7) + 
  # geom_text(data = pair.comp.labels, aes(label = group, y = ypos, x = network_interaction),
  #           position = position_dodge(2),
  #           color = "grey",
  #           show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black") +
  facet_grid(cols = vars(metric)) +
  # facet_grid(cols = vars(metric), rows = vars(null_model)) +
  theme_bw() +
  # coord_cartesian(ylim = c(-6.5,10)) +
  coord_cartesian(ylim = c(-5,5)) +
  # coord_flip() +
  # scale_x_discrete(limits = rev(levels(zscores.plot.data$metric))) +
  labs(x = "", y = "z-score") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),legend.position = "bottom")+
  # theme(legend.justification=c(1,-0.01), legend.position=c(0.99,0)) +
  ggtitle("Soft-constrained degree distribution") +
  NULL
# zscores.degdist.soft.plot

zscores.degdist.hard.plot <- ggplot(subset(zscores.plot.full.data, null_model == "degdist_hard"),
# zscores.degdist.hard.plot <- ggplot(subset(zscores.plot.full.data, null_model == "degdist_hard" & metric %in% degdist.metrics), 
                                    aes(x = network_interaction, y = z_score)) + 
  geom_boxplot(aes(fill = network_interaction),
               outlier.colour = "grey30",outlier.size = .7) + 
  # geom_text(data = pair.comp.labels, aes(label = group, y = ypos, x = network_interaction),
  #           position = position_dodge(2),
  #           color = "grey",
  #           show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black") +
  facet_grid(cols = vars(metric)) +
  # facet_grid(cols = vars(metric), rows = vars(null_model)) +
  theme_bw() +
  # coord_cartesian(ylim = c(-6.5,10)) +
  coord_cartesian(ylim = c(-5,5)) +
  # coord_flip() +
  # scale_x_discrete(limits = rev(levels(zscores.plot.data$metric))) +
  labs(x = "", y = "z-score") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),legend.position = "bottom")+
  # theme(legend.justification=c(1,-0.01), legend.position=c(0.99,0)) +
  ggtitle("Hard-constrained degree distribution") +
  NULL
# zscores.degdist.hard.plot

# -------------------------------------------------------------------------
full_zscores_plot <- zscores.connectance.plot/zscores.degdist.soft.plot/zscores.degdist.hard.plot + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
# -------------------------------------------------------------------------

# ggsave("results/images/aux_zscores_full_plot.pdf",full_zscores_plot,
#        device = cairo_pdf,
#        width = 12,height = 10,dpi = 300)
