
library(tidyverse)
library(colorblindr)
library(patchwork)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")
network.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

null.metrics.bipartite <- read.csv2("results/bipartite_network_metrics_null.csv")
null.metrics.unipartite <- read.csv2("results/unipartite_network_metrics_null.csv")

null.metrics <- bind_rows(null.metrics.bipartite,null.metrics.unipartite) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")]) %>%
  drop_na()

# network.metrics.long <- read.csv2("results/metrics_observed_null_long.csv")
metrics.zscores <- read.csv2("results/metrics_z_scores.csv")
metrics.zscores$study_common <- net.collection$study_common[match(metrics.zscores$network_id,net.collection$network_id)]

metrics.zscores.clean <- subset(metrics.zscores, !is.infinite(z_score) & abs(z_score) < 10)

# -------------------------------------------------------------------------
# I could explore drawing observed values and distribution of nulls, but
# there are 80 combinations of 10 metrics * 2 topologies * 4 interactions,
# and each of those is not even a single network, but potentially many
# so it starts to be very complex.

# For now, plot the z-scores for these combinations

# -------------------------------------------------------------------------
zscores.plot.full.data <- subset(metrics.zscores, !(metric %in% c("richness","degree.shannon","degree.mean",
                                                                         "degree.sd","degree.skewness",
                                                                  "degree.kurtosis","connectance","link.density")))

zscores.plot.full.data$metric[zscores.plot.full.data$metric == "modularity.infomap"] <- "infomap \nmodularity"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
# zscores.plot.full.data$metric[zscores.plot.full.data$metric == "link.density"] <- "link \ndensity"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "interaction.overlap"] <- "interaction \noverlap"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.eigen"] <- "eigenvector \ncentralization"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.degree"] <- "degree \ncentralization"
# zscores.plot.full.data$metric[zscores.plot.full.data$metric == "centrality.betweenness"] <- "betweenness \ncentrality"
zscores.plot.full.data$metric[zscores.plot.full.data$metric == "nestedness"] <- "NODF nestedness"

zscores.plot.full.data$metric <- factor(zscores.plot.full.data$metric, levels = c(#"connectance","link \ndensity",
                                                                                  "interaction \noverlap",
                                                                                  "NODF nestedness",
                                                                                  "betweenness \nmodularity",
                                                                                  "infomap \nmodularity",
                                                                                  # "betweenness \ncentrality",
                                                                                  "degree \ncentralization",
                                                                                  "eigenvector \ncentralization"
                                                                                  ))

zscores.plot.data <- subset(zscores.plot.full.data, 
                              !(network_interaction %in% c("unspecified","competition","other mutualism")) &
                               network_topology_type == "bipartite")

# sample size?
zscores.num <- zscores.plot.data %>% 
  dplyr::select(network_id, network_interaction) %>%
  unique() %>%
  group_by(network_interaction) %>%
  summarise(num = n())
metrics.num <- network.metrics %>%
  dplyr::select(network_id, network_interaction, network_topology_type) %>%
  unique() %>%
  group_by(network_interaction,network_topology_type) %>%
  summarise(num = n())

# -------------------------------------------------------------------------
# plot a single example of a network, for a given metric.

example_network_id <- "AI22_1128"
example_metric <- "nestedness"

my.zscore <- metrics.zscores$z_score[metrics.zscores$network_id == example_network_id &
                                       metrics.zscores$metric == example_metric]

my.null.data <- subset(null.metrics.bipartite, 
                       network_id == example_network_id &
                         metric == example_metric)
my.obs.data <- network.metrics %>%
  filter(metric == example_metric & network_id == example_network_id)

example.metric.plot <- ggplot(my.null.data, aes(x = value)) + 
  geom_density(alpha = .4) + 
  geom_point(data = my.obs.data, aes(x = value, 
                                    y = 0),
             size = 4) +
  theme_bw() +
  # xlim(6,22) +
  theme(legend.position="none")+
  labs(x = example_metric,y = "density estimate") +
  NULL
# example.metric.plot

# -------------------------------------------------------------------------
# statistical analyses
# are these patterns (zscores.full.plot) significant?
# one model for each metric

my.metrics <- c("centrality.degree","centrality.eigen","interaction.overlap",
                "modularity.betweenness","modularity.infomap","nestedness")
pair.comparisons <- list()
for(i.metric in 1:length(my.metrics)){
  my.zscore.data <- subset(metrics.zscores.clean, metric == my.metrics[i.metric])
  my.model <- glmmTMB(z_score ~ network_interaction + (1|study_common), data = my.zscore.data)
emmeans.res <- emmeans(my.model, ~network_interaction)
my.pairs <- as.data.frame(pairs(emmeans.res, adjust = "tukey"))[c("contrast","p.value")]

pair.comparisons[[length(pair.comparisons)+1]] <- tibble(metric = my.metrics[i.metric],my.pairs)
}

pair.comp.df <- bind_rows(pair.comparisons)

# look at these pairwise comparisons and manually add labels
pair.comp.labels <- expand_grid(metric = c(
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  "degree \ncentralization",
  "eigenvector \ncentralization"
),network_interaction = unique(metrics.zscores$network_interaction))

pair.comp.labels$ypos <- -6.5
pair.comp.labels$group <- "a"
pair.comp.labels$group[pair.comp.labels$metric == "interaction \noverlap" & pair.comp.labels$network_interaction == "frugivory"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "interaction \noverlap" & pair.comp.labels$network_interaction == "herbivory"] <- "c"
pair.comp.labels$group[pair.comp.labels$metric == "interaction \noverlap" & pair.comp.labels$network_interaction == "parasitism"] <- "bc"
pair.comp.labels$group[pair.comp.labels$metric == "interaction \noverlap" & pair.comp.labels$network_interaction == "pollination"] <- "b"

pair.comp.labels$group[pair.comp.labels$metric == "NODF nestedness" & pair.comp.labels$network_interaction == "frugivory"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "NODF nestedness" & pair.comp.labels$network_interaction == "herbivory"] <- "c"
pair.comp.labels$group[pair.comp.labels$metric == "NODF nestedness" & pair.comp.labels$network_interaction == "parasitism"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "NODF nestedness" & pair.comp.labels$network_interaction == "pollination"] <- "b"

pair.comp.labels$group[pair.comp.labels$metric == "betweenness \nmodularity" & pair.comp.labels$network_interaction == "herbivory"] <- "c"
pair.comp.labels$group[pair.comp.labels$metric == "betweenness \nmodularity" & pair.comp.labels$network_interaction == "frugivory"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "betweenness \nmodularity" & pair.comp.labels$network_interaction == "parasitism"] <- "abc"
pair.comp.labels$group[pair.comp.labels$metric == "betweenness \nmodularity" & pair.comp.labels$network_interaction == "pollination"] <- "ac"

pair.comp.labels$group[pair.comp.labels$metric == "infomap \nmodularity" & pair.comp.labels$network_interaction == "herbivory"] <- "c"
pair.comp.labels$group[pair.comp.labels$metric == "infomap \nmodularity" & pair.comp.labels$network_interaction == "frugivory"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "infomap \nmodularity" & pair.comp.labels$network_interaction == "parasitism"] <- "cd"
pair.comp.labels$group[pair.comp.labels$metric == "infomap \nmodularity" & pair.comp.labels$network_interaction == "pollination"] <- "d"

pair.comp.labels$group[pair.comp.labels$metric == "degree \ncentralization" & pair.comp.labels$network_interaction == "herbivory"] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "degree \ncentralization" & pair.comp.labels$network_interaction == "frugivory"] <- "ab"
pair.comp.labels$group[pair.comp.labels$metric == "degree \ncentralization" & pair.comp.labels$network_interaction == "parasitism"] <- "ab"

pair.comp.labels$group[pair.comp.labels$metric == "eigenvector \ncentralization" & pair.comp.labels$network_interaction == "herbivory"] <- "b"

# -------------------------------------------------------------------------
# after doing the analyses, plot the boxplots with the groups that are 
# significantly different

zscores.plot <- ggplot(zscores.plot.full.data, aes(x = network_interaction, y = z_score)) + 
  geom_boxplot(aes(fill = network_interaction),
               outlier.colour = "grey30",outlier.size = .7) + 
  geom_text(data = pair.comp.labels, aes(label = group, y = ypos, x = network_interaction), 
            position = position_dodge(2),
            color = "grey",
            show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black") +
  facet_grid(cols = vars(metric)) +
  theme_bw() +
  ylim(-6.5,6) +
  # coord_flip() +
  # scale_x_discrete(limits = rev(levels(zscores.plot.data$metric))) +
  labs(x = "", y = "z-score") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),legend.position = "bottom")+
  # theme(legend.justification=c(1,-0.01), legend.position=c(0.99,0)) +
  NULL
# zscores.plot

# zscores.unip.data <- subset(zscores.plot.full.data,
#                             !(network_interaction %in% c("unspecified","competition",
#                                                          "unspecified mutualism")) &
#                               network_topology_type == "unipartite")
# 
# zscores.uniplot <- ggplot(zscores.unip.data, aes(x = metric, y = z_score)) + 
#   geom_boxplot(aes(fill = network_interaction),
#                outlier.colour = "grey30",outlier.size = .7) + 
#   scale_fill_OkabeIto() +
#   geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
#   geom_hline(yintercept = -2, linetype = "dashed", color = "black") +
#   # facet_grid(.~network_topology_type) +
#   theme_bw() +
#   ylim(-6,6) +
#   scale_x_discrete(limits = rev(levels(zscores.unip.data$metric))) +
#   coord_flip() +
#   labs(x = "", y = "") +
#   theme(legend.title=element_blank())+
#   # theme(legend.justification=c(1,-0.05), legend.position=c(0.99,0)) +
#   # theme(#axis.ticks.y = element_blank(), 
#   #       axis.text.y = element_blank()) +
#   NULL
# # zscores.uniplot
# 
# zscores.full.plot <- zscores.plot + zscores.uniplot

# -------------------------------------------------------------------------

# ggsave("results/images/z_score_example_plot.pdf",example.metric.plot,
#        device = cairo_pdf,
#        width = 6,height = 4,dpi = 300)
# 
ggsave("results/images/z_scores_plot_v2.pdf",zscores.plot,
       device = cairo_pdf,
       width = 12,height = 6,dpi = 300)

