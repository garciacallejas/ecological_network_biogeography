
library(tidyverse)
library(colorblindr)
library(patchwork)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
library(gghalves)
library(ggbeeswarm) 
library(corrplot)
# -------------------------------------------------------------------------
# comment dataframes that I do not currently use

net.collection <- read_csv2("results/network_collection.csv")

# metrics from observed networks

network.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

metrics.zscores <- read.csv2("results/metrics_z_scores.csv")
metrics.zscores$study_common <- net.collection$study_common[match(metrics.zscores$network_id,net.collection$network_id)]

# pair comparisons among groups
observed_comparisons_labels <- read.csv2("results/observed_metrics_pairwise_comparisons.csv")
z_comparisons_labels <- read.csv2("results/z_scores_pairwise_comparisons_hard_degdist.csv")

# -------------------------------------------------------------------------
# prepare data for plots

observed_metrics_clean <- network.metrics
# in the main text I will plot the hard-constrained degdist null model
zscores_clean <- metrics.zscores#subset(metrics.zscores, null_model == "degdist_hard")

my.metrics <- c("interaction \noverlap",
                "NODF nestedness",
                "betweenness \nmodularity",
                "infomap \nmodularity",
                "eigenvector \ncentralisation")

my.metrics.conn <- c("connectance","interaction \noverlap",
                     "NODF nestedness",
                     "betweenness \nmodularity",
                     "infomap \nmodularity",
                     "eigenvector \ncentralisation")
my.metrics.conn.degmean <- c("connectance","degree.mean","interaction \noverlap",
                             "NODF nestedness",
                             "betweenness \nmodularity",
                             "infomap \nmodularity",
                             "eigenvector \ncentralisation")
# -------------------------------------------------------------------------
observed_metrics_clean <- subset(observed_metrics_clean, !(metric %in% c("richness","degree.shannon","degree.mean",
                                                                         "degree.sd","degree.skewness",
                                                                  "degree.kurtosis","connectance","link.density","centrality.degree")))

observed_metrics_clean$metric[observed_metrics_clean$metric == "modularity.infomap"] <- "infomap \nmodularity"
observed_metrics_clean$metric[observed_metrics_clean$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
observed_metrics_clean$metric[observed_metrics_clean$metric == "interaction.overlap"] <- "interaction \noverlap"
observed_metrics_clean$metric[observed_metrics_clean$metric == "centrality.eigen"] <- "eigenvector \ncentralisation"
observed_metrics_clean$metric[observed_metrics_clean$metric == "nestedness"] <- "NODF nestedness"

observed_metrics_clean$metric <- factor(observed_metrics_clean$metric, levels = my.metrics)
observed_comparisons_labels$metric <- factor(observed_comparisons_labels$metric, levels = my.metrics)

# -------------------------------------------------------------------------
zscores_clean <- subset(zscores_clean, !(metric %in% c("richness","degree.shannon","degree.mean",
                                                                         "degree.sd","degree.skewness",
                                                                         "degree.kurtosis","connectance","link.density","centrality.degree")))

zscores_clean$metric[zscores_clean$metric == "modularity.infomap"] <- "infomap \nmodularity"
zscores_clean$metric[zscores_clean$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
# zscores_clean$metric[zscores_clean$metric == "link.density"] <- "link \ndensity"
zscores_clean$metric[zscores_clean$metric == "interaction.overlap"] <- "interaction \noverlap"
zscores_clean$metric[zscores_clean$metric == "centrality.eigen"] <- "eigenvector \ncentralisation"
# zscores_clean$metric[zscores_clean$metric == "centrality.degree"] <- "degree \ncentralisation"
# zscores_clean$metric[zscores_clean$metric == "centrality.betweenness"] <- "betweenness \ncentrality"
zscores_clean$metric[zscores_clean$metric == "nestedness"] <- "NODF nestedness"

zscores_clean$metric <- factor(zscores_clean$metric, levels = my.metrics)
z_comparisons_labels$metric <- factor(z_comparisons_labels$metric, levels = my.metrics)

# -------------------------------------------------------------------------
# observed metrics plot, Fig. 2 of main text

observed_metrics_conn_clean <- subset(network.metrics, !(metric %in% c("richness","degree.shannon",
                                                                         "degree.sd","degree.skewness",
                                                                         "degree.kurtosis","link.density","centrality.degree")))

observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "modularity.infomap"] <- "infomap \nmodularity"
observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "interaction.overlap"] <- "interaction \noverlap"
observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "centrality.eigen"] <- "eigenvector \ncentralisation"
observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "nestedness"] <- "NODF nestedness"
# observed_metrics_conn_clean$metric[observed_metrics_conn_clean$metric == "degree.mean"] <- "deg. dist. average"

observed_metrics_conn_clean$metric <- factor(observed_metrics_conn_clean$metric, levels = my.metrics.conn.degmean)

# 1. Connectance and network metrics
conn_metrics_data <- observed_metrics_conn_clean %>% pivot_wider(names_from = "metric", values_from = "value") %>%
  pivot_longer(cols = all_of(my.metrics),names_to = "metric") %>%
  select(network_id,network_topology_type,network_interaction,connectance,degree.mean,metric,value)
conn_metrics_data$metric <- factor(conn_metrics_data$metric, levels = my.metrics)

# -------------------------------------------------------------------------
# spearman correlations
my.int <- sort(unique(observed_metrics_conn_clean$network_interaction))
conn_metrics_wide <- observed_metrics_conn_clean %>% pivot_wider(names_from = "metric", values_from = "value")
conn.corr <- list()
for(i.metric in 1:length(my.metrics)){
  for(i.type in 1:length(my.int)){
    my.conn <- as.numeric(conn_metrics_wide$connectance[conn_metrics_wide$network_interaction == my.int[i.type]])
    my.metric.data <- pull(conn_metrics_wide[conn_metrics_wide$network_interaction == my.int[i.type],
                                                   which(names(conn_metrics_wide) == my.metrics[i.metric])])
    my.test <- cor.test(x = my.conn,y = my.metric.data,conf.level = 0.95,method = "spearman")
    my.cor <- my.test$estimate
    my.pvalue <- my.test$p.value
    conn.corr[[length(conn.corr)+1]] <- tibble(network_interaction = my.int[i.type],
                                               metric_1 = "connectance",
                                               metric_2 = my.metrics[i.metric],
                                               spearman_rho = my.cor,
                                               p_value = my.pvalue)
    
  }
}

conn_corr_df <- bind_rows(conn.corr)

deg.corr <- list()
for(i.metric in 1:length(my.metrics)){
  for(i.type in 1:length(my.int)){
    my.deg <- as.numeric(conn_metrics_wide$degree.mean[conn_metrics_wide$network_interaction == my.int[i.type]])
    my.metric.data <- pull(conn_metrics_wide[conn_metrics_wide$network_interaction == my.int[i.type],
                                             which(names(conn_metrics_wide) == my.metrics[i.metric])])
    my.test <- cor.test(x = my.deg,y = my.metric.data,conf.level = 0.95,method = "spearman",exact = F)
    my.cor <- my.test$estimate
    my.pvalue <- my.test$p.value
    deg.corr[[length(deg.corr)+1]] <- tibble(network_interaction = my.int[i.type],
                                               metric_1 = "degree.mean",
                                               metric_2 = my.metrics[i.metric],
                                               spearman_rho = my.cor,
                                             p_value = my.pvalue)
    
  }
}

deg_corr_df <- bind_rows(deg.corr)

full_corr_df <- bind_rows(conn_corr_df,deg_corr_df)

pd <- 0.2
corr_metrics_plot <- ggplot(full_corr_df, aes(x = spearman_rho, y = metric_2)) + 
  geom_point(aes(fill = network_interaction), 
             shape = 21,size = 1.5, 
             position = position_dodge(pd)) +
  facet_grid(cols = vars(metric_1)) + 
  scale_fill_OkabeIto(order = c(1:8)) +
  labs(x = expression(rho), y = "") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size = 2))) +
  NULL

# Supplementary figure
# ggsave("results/images/observed_metrics_correlations.pdf",corr_metrics_plot,
#        device = cairo_pdf,
#        width = 6,height = 3,dpi = 300)

# -------------------------------------------------------------------------
# plot
conn_metrics_plot <- ggplot(conn_metrics_data, aes(x = connectance, y = value)) + 
  geom_point(aes(fill = network_interaction),shape = 21, alpha = 0.7) + 
  facet_wrap(facets = vars(metric),nrow = 1,scales = "free_y") + 
  scale_fill_OkabeIto(order = c(1:8),guide = "none") +
  ylab("metric value") +
  theme_bw() +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(size=12),
        strip.text.x = element_blank(),
        legend.position = "none")+
  NULL
# conn_metrics_plot

deg_metrics_plot <- ggplot(conn_metrics_data, aes(x = degree.mean, y = value)) + 
  geom_point(aes(fill = network_interaction),shape = 21, alpha = 0.7) + 
  facet_wrap(facets = vars(metric),nrow = 1,scales = "free_y") + 
  scale_fill_OkabeIto(order = c(1:8)) +
  ylab("metric value") +
  xlab("degree distribution mean") +
  theme_bw() +
  theme(legend.title=element_blank())+
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(size=12),
        strip.text.x = element_blank(),
        legend.position = "bottom")+
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  NULL
# deg_metrics_plot

observed_comparisons_labels$ypos <- -0.05
observed_metrics_plot <- ggplot(observed_metrics_clean, 
                               aes(x = network_interaction, y = value)) + 
  
  # geom_half_point(aes(color = network_interaction),
  #                 transformation = position_quasirandom(width = 0.1),
  #                 side = "l", size = 0.5, alpha = 0.5) +
  # geom_half_boxplot(aes(fill = network_interaction), side = "r",outlier.size = 0.7,outlier.colour = "grey30") +
  
  geom_boxplot(aes(fill = network_interaction),
               outlier.colour = "grey30",outlier.size = .7) +
  geom_text(data = observed_comparisons_labels, aes(label = group, y = ypos, x = network_interaction),
            position = position_dodge(3),
            color = "grey",
            show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8), guide ="none") +
  scale_color_OkabeIto(order = c(1:8), guide = "none")  +
  
  facet_wrap(facets = vars(metric),nrow = 1,scales = "free_y",axis.labels = "margins") +
  # facet_grid(cols = vars(metric), rows = vars(null_model)) +
  theme_bw() +
  # coord_cartesian(ylim = c(-6.5,10)) +
  # coord_flip() +
  # scale_x_discrete(limits = rev(levels(zscores.plot.data$metric))) +
  labs(x = "", y = "metric value") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.title=element_blank())+
  # theme(strip.background = element_blank(),legend.position = "none")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")+
  # theme(legend.justification=c(1,-0.01), legend.position=c(0.99,0)) +
  NULL
# observed_metrics_plot

obs_metrics_combined_plot <- observed_metrics_plot/conn_metrics_plot/deg_metrics_plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# ggsave("results/images/observed_metrics_combined_plot.pdf",obs_metrics_combined_plot,
#        device = cairo_pdf,
#        width = 12,height = 8,dpi = 300)

# -------------------------------------------------------------------------
# zscore plot
zscores_clean$null_model[zscores_clean$null_model == "degdist_soft"] <- "soft\ndeg. dist."
zscores_clean$null_model[zscores_clean$null_model == "degdist_hard"] <- "hard\ndeg. dist."
zscores_clean$null_model <- factor(zscores_clean$null_model, levels = c("connectance",
                                                                        "soft\ndeg. dist.",
                                                                        "hard\ndeg. dist."))

zscores_aggr <- zscores_clean %>% 
  filter(z_score_full < 1e4) %>%
  group_by(null_model,network_interaction,metric) %>%
  summarise(avg_zscore = mean(z_score_full),
            sd_zscore = sd(z_score_full))
zscores_aggr$avg_zscore[zscores_aggr$avg_zscore>41] <- NA
zscores_aggr$sd_zscore[is.na(zscores_aggr$avg_zscore)] <- NA

pd <- 0.3

zscores_line_plot <- ggplot(zscores_aggr, aes(x = null_model, y = avg_zscore)) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "darkgrey") +
  geom_linerange(aes(ymin = avg_zscore - sd_zscore,
                      ymax = avg_zscore + sd_zscore,
                      color = network_interaction),position = position_dodge(pd)) +
  geom_line(aes(color = network_interaction, group = network_interaction),position = position_dodge(pd)) +
  geom_point(aes(fill = network_interaction), shape = 21, position = position_dodge(pd)) +
  scale_fill_OkabeIto(order = c(1:8),guide="none") +
  scale_color_OkabeIto(order = c(1:8),guide="none") +
  coord_cartesian(ylim = c(-5,41)) +
  facet_wrap(facets = vars(metric),nrow = 1) +
  theme_bw() +
  labs(y = "z-score",x="") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "none")+
  NULL
# zscores_line_plot

z_comparisons_labels$ypos <- -3.7
zscores_degdisthard_line_plot <- ggplot(subset(zscores_aggr, null_model == "hard\ndeg. dist."), 
                                        aes(x = network_interaction, y = avg_zscore)) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "darkgrey") +
  geom_linerange(aes(ymin = avg_zscore - sd_zscore,
                     ymax = avg_zscore + sd_zscore,
                     color = network_interaction),linewidth = 1.3) +
  # geom_line(aes(color = network_interaction, group = network_interaction),position = position_dodge(pd)) +
  geom_point(aes(fill = network_interaction), shape = 21, 
             size = 3) +
  geom_text(data = z_comparisons_labels, aes(label = group, y = ypos, x = network_interaction),
            position = position_dodge(2),
            color = "grey",
            show.legend = FALSE ) +
  scale_fill_OkabeIto(order = c(1:8)) +
  scale_color_OkabeIto(order = c(1:8)) +
  coord_cartesian(ylim = c(-4,4)) +
  facet_wrap(facets = vars(metric),nrow = 1) +
  theme_bw() +
  labs(y = "z-score",x="") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")+
  NULL
# zscores_degdisthard_line_plot

zscores_plot <- zscores_line_plot/zscores_degdisthard_line_plot + plot_layout(guides = "collect")  & theme(legend.position = "bottom")
# zscores_plot

# ggsave("results/images/z_scores_line_plot.pdf",zscores_plot,
#        device = cairo_pdf,
#        width = 13,height = 7,dpi = 300)

