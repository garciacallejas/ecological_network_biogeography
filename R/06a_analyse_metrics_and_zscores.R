
library(tidyverse)
library(colorblindr)
library(patchwork)
library(glmmTMB)
library(broom.mixed)
library(DHARMa)
library(emmeans)
# -------------------------------------------------------------------------
# comment dataframes that I do not currently use

net.collection <- read_csv2("results/network_collection.csv")

# metrics from observed networks

network.metrics <- list.files("results/metrics",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction")])

metrics.zscores <- read.csv2("results/metrics_z_scores.csv")

metrics.zscores$study_common <- net.collection$study_common[match(metrics.zscores$network_id,net.collection$network_id)]

# -------------------------------------------------------------------------
# prepare data for plot

observed_metrics_clean <- network.metrics
# in the main text I will plot the hard-constrained degdist null model
zscores_clean <- subset(metrics.zscores, null_model == "degdist_hard")

# -------------------------------------------------------------------------
observed_metrics_clean <- subset(observed_metrics_clean, !(metric %in% c("richness","degree.shannon","degree.mean",
                                                                         "degree.sd","degree.skewness",
                                                                         "degree.kurtosis","connectance","link.density","centrality.degree")))

observed_metrics_clean$metric[observed_metrics_clean$metric == "modularity.infomap"] <- "infomap \nmodularity"
observed_metrics_clean$metric[observed_metrics_clean$metric == "modularity.betweenness"] <- "betweenness \nmodularity"
# observed_metrics_clean$metric[observed_metrics_clean$metric == "link.density"] <- "link \ndensity"
observed_metrics_clean$metric[observed_metrics_clean$metric == "interaction.overlap"] <- "interaction \noverlap"
observed_metrics_clean$metric[observed_metrics_clean$metric == "centrality.eigen"] <- "eigenvector \ncentralisation"
# observed_metrics_clean$metric[observed_metrics_clean$metric == "centrality.degree"] <- "degree \ncentralisation"
# observed_metrics_clean$metric[observed_metrics_clean$metric == "centrality.betweenness"] <- "betweenness \ncentrality"
observed_metrics_clean$metric[observed_metrics_clean$metric == "nestedness"] <- "NODF nestedness"

observed_metrics_clean$metric <- factor(observed_metrics_clean$metric, levels = c(#"connectance","link \ndensity",
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  # "degree \ncentralisation",
  "eigenvector \ncentralisation"
))
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

zscores_clean$metric <- factor(zscores_clean$metric, levels = c(#"connectance","link \ndensity",
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  # "degree \ncentralisation",
  "eigenvector \ncentralisation"
))

# -------------------------------------------------------------------------
# statistical analyses for differences across groups in

# 1) observed metric values

my.metrics <- c("interaction \noverlap","NODF nestedness","betweenness \nmodularity",
                "infomap \nmodularity","eigenvector \ncentralisation")
observed.pair.comparisons <- list()

for(i.metric in 1:length(my.metrics)){
  my.metric.data <- subset(observed_metrics_clean, metric == my.metrics[i.metric])
  my.model <- glmmTMB(value ~ network_interaction + (1|network_id), data = my.metric.data)
emmeans.res <- emmeans(my.model, ~network_interaction)
my.pairs <- as.data.frame(pairs(emmeans.res, adjust = "tukey"))[c("contrast","p.value")]

observed.pair.comparisons[[length(observed.pair.comparisons)+1]] <- tibble(metric = my.metrics[i.metric],my.pairs)
}

observed_comparisons_df <- bind_rows(observed.pair.comparisons)

# look at these pairwise comparisons and manually add labels
observed_comparisons_labels <- expand_grid(metric = c(
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  # "degree \ncentralisation",
  "eigenvector \ncentralisation"
),network_interaction = unique(observed_metrics_clean$network_interaction))

observed_comparisons_labels$ypos <- -0.5
observed_comparisons_labels$group <- "a"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "interaction \noverlap" & observed_comparisons_labels$network_interaction == "frugivory"] <- "a"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "interaction \noverlap" & observed_comparisons_labels$network_interaction == "herbivory"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "interaction \noverlap" & observed_comparisons_labels$network_interaction == "parasitism"] <- "a"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "interaction \noverlap" & observed_comparisons_labels$network_interaction == "pollination"] <- "c"

observed_comparisons_labels$group[observed_comparisons_labels$metric == "NODF nestedness" & observed_comparisons_labels$network_interaction == "frugivory"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "NODF nestedness" & observed_comparisons_labels$network_interaction == "herbivory"] <- "c"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "NODF nestedness" & observed_comparisons_labels$network_interaction == "parasitism"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "NODF nestedness" & observed_comparisons_labels$network_interaction == "pollination"] <- "d"

observed_comparisons_labels$group[observed_comparisons_labels$metric == "betweenness \nmodularity" & observed_comparisons_labels$network_interaction == "frugivory"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "betweenness \nmodularity" & observed_comparisons_labels$network_interaction == "herbivory"] <- "c"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "betweenness \nmodularity" & observed_comparisons_labels$network_interaction == "parasitism"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "betweenness \nmodularity" & observed_comparisons_labels$network_interaction == "pollination"] <- "d"

observed_comparisons_labels$group[observed_comparisons_labels$metric == "infomap \nmodularity" & observed_comparisons_labels$network_interaction == "frugivory"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "infomap \nmodularity" & observed_comparisons_labels$network_interaction == "herbivory"] <- "c"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "infomap \nmodularity" & observed_comparisons_labels$network_interaction == "parasitism"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "infomap \nmodularity" & observed_comparisons_labels$network_interaction == "pollination"] <- "d"

observed_comparisons_labels$group[observed_comparisons_labels$metric == "eigenvector \ncentralisation" & observed_comparisons_labels$network_interaction == "frugivory"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "eigenvector \ncentralisation" & observed_comparisons_labels$network_interaction == "herbivory"] <- "c"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "eigenvector \ncentralisation" & observed_comparisons_labels$network_interaction == "parasitism"] <- "b"
observed_comparisons_labels$group[observed_comparisons_labels$metric == "eigenvector \ncentralisation" & observed_comparisons_labels$network_interaction == "pollination"] <- "d"

# -------------------------------------------------------------------------
# 2 - z-scores

z.pair.comparisons <- list()

for(i.metric in 1:length(my.metrics)){
  my.metric.data <- subset(zscores_clean, metric == my.metrics[i.metric])
  my.model <- glmmTMB(z_score_full ~ network_interaction + (1|study_common), data = my.metric.data)
  emmeans.res <- emmeans(my.model, ~network_interaction)
  my.pairs <- as.data.frame(pairs(emmeans.res, adjust = "tukey"))[c("contrast","p.value")]
  
  z.pair.comparisons[[length(z.pair.comparisons)+1]] <- tibble(metric = my.metrics[i.metric],my.pairs)
}

z_comparisons_df <- bind_rows(z.pair.comparisons)

# look at these pairwise comparisons and manually add labels
z_comparisons_labels <- expand_grid(metric = c(
  "interaction \noverlap",
  "NODF nestedness",
  "betweenness \nmodularity",
  "infomap \nmodularity",
  # "betweenness \ncentrality",
  # "degree \ncentralisation",
  "eigenvector \ncentralisation"
),network_interaction = unique(observed_metrics_clean$network_interaction))

z_comparisons_labels$ypos <- -7.5
z_comparisons_labels$group <- "a"
z_comparisons_labels$group[z_comparisons_labels$metric == "interaction \noverlap" & z_comparisons_labels$network_interaction == "frugivory"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "interaction \noverlap" & z_comparisons_labels$network_interaction == "herbivory"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "interaction \noverlap" & z_comparisons_labels$network_interaction == "parasitism"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "interaction \noverlap" & z_comparisons_labels$network_interaction == "pollination"] <- "b"

z_comparisons_labels$group[z_comparisons_labels$metric == "NODF nestedness" & z_comparisons_labels$network_interaction == "frugivory"] <- "ab"
z_comparisons_labels$group[z_comparisons_labels$metric == "NODF nestedness" & z_comparisons_labels$network_interaction == "herbivory"] <- "ab"
z_comparisons_labels$group[z_comparisons_labels$metric == "NODF nestedness" & z_comparisons_labels$network_interaction == "parasitism"] <- "ab"
z_comparisons_labels$group[z_comparisons_labels$metric == "NODF nestedness" & z_comparisons_labels$network_interaction == "pollination"] <- "b"

z_comparisons_labels$group[z_comparisons_labels$metric == "betweenness \nmodularity" & z_comparisons_labels$network_interaction == "herbivory"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "betweenness \nmodularity" & z_comparisons_labels$network_interaction == "frugivory"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "betweenness \nmodularity" & z_comparisons_labels$network_interaction == "parasitism"] <- "b"
z_comparisons_labels$group[z_comparisons_labels$metric == "betweenness \nmodularity" & z_comparisons_labels$network_interaction == "pollination"] <- "b"

z_comparisons_labels$group[z_comparisons_labels$metric == "infomap \nmodularity" & z_comparisons_labels$network_interaction == "herbivory"] <- "a"
z_comparisons_labels$group[z_comparisons_labels$metric == "infomap \nmodularity" & z_comparisons_labels$network_interaction == "frugivory"] <- "a"
z_comparisons_labels$group[z_comparisons_labels$metric == "infomap \nmodularity" & z_comparisons_labels$network_interaction == "parasitism"] <- "a"
z_comparisons_labels$group[z_comparisons_labels$metric == "infomap \nmodularity" & z_comparisons_labels$network_interaction == "pollination"] <- "a"

z_comparisons_labels$group[z_comparisons_labels$metric == "eigenvector \ncentralisation" & z_comparisons_labels$network_interaction == "herbivory"] <- "ab"
z_comparisons_labels$group[z_comparisons_labels$metric == "eigenvector \ncentralisation" & z_comparisons_labels$network_interaction == "frugivory"] <- "ab"
z_comparisons_labels$group[z_comparisons_labels$metric == "eigenvector \ncentralisation" & z_comparisons_labels$network_interaction == "parasitism"] <- "a"
z_comparisons_labels$group[z_comparisons_labels$metric == "eigenvector \ncentralisation" & z_comparisons_labels$network_interaction == "pollination"] <- "b"

# -------------------------------------------------------------------------

write.csv2(observed_comparisons_labels,"results/observed_metrics_pairwise_comparisons.csv",row.names = F)
write.csv2(z_comparisons_labels,"results/z_scores_pairwise_comparisons_hard_degdist.csv",row.names = F)
