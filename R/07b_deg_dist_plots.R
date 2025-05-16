# -----------------------------------------------------------------------

library(tidyverse)
# library(glmmTMB)
library(DHARMa)
library(performance)
library(effects)
library(corrplot)
library(GGally)
library(colorblindr)
library(patchwork)
library(gghalves)
library(ggbeeswarm) 
# library(jtools)
library(scales)
# library(huxtable)
library(broom)
library(broom.mixed)
library(xtable)
library(ggeffects)
library(glmmTMB)
library(emmeans)
# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")

observed.metrics <- list.files("results/metrics/",full.names = T) %>%
  map_dfr(read.csv2) %>%
  left_join(net.collection[,c("network_id","network_topology_type","network_interaction","study_common")])

degree.metrics <- observed.metrics %>% filter(metric %in% c("degree.mean",
                                                            "degree.sd",
                                                            "degree.skewness",
                                                            "degree.kurtosis") &
                                                network_interaction %in% c("pollination",
                                                                           "food web",
                                                                           "herbivory",
                                                                           "frugivory",
                                                                           "parasitism")) 
degree.metrics.wide <- degree.metrics %>%
  pivot_wider(names_from = metric)

degree.metrics.for.plots <- degree.metrics
degree.metrics.for.plots$metric[degree.metrics.for.plots$metric == "degree.mean"] <- "degree distribution\n average"
degree.metrics.for.plots$metric[degree.metrics.for.plots$metric == "degree.sd"] <- "degree distribution\n standard deviation"
degree.metrics.for.plots$metric[degree.metrics.for.plots$metric == "degree.skewness"] <- "degree distribution\n skewness"
degree.metrics.for.plots$metric[degree.metrics.for.plots$metric == "degree.kurtosis"] <- "degree distribution\n kurtosis"

degree.metrics.for.plots$metric <- factor(degree.metrics.for.plots$metric,levels = c("degree distribution\n average",
                                                                                     "degree distribution\n standard deviation",
                                                                                     "degree distribution\n skewness",
                                                                                     "degree distribution\n kurtosis"))

# -------------------------------------------------------------------------
# first plot: richness and degree distribution metrics
# this is supplementary figure S1
rich.deg.data <- subset(observed.metrics, metric %in% c("richness","degree.mean",
                                                        "degree.sd",
                                                        "degree.skewness",
                                                        "degree.kurtosis")) %>%
  select(-study_common) %>%
  pivot_wider(names_from = metric,values_from = value)
rich.deg.plot <- ggpairs(rich.deg.data,columns = 4:8,mapping = ggplot2::aes(colour=network_interaction))
for(i in 1:rich.deg.plot$nrow) {
  for(j in 1:rich.deg.plot$ncol){
    rich.deg.plot[i,j] <- rich.deg.plot[i,j] +
      scale_fill_OkabeIto(darken = .2) +
      scale_color_OkabeIto(darken = .2)
  }
}
# rich.deg.plot

# ggsave("results/images/degree_metrics_correlation.pdf",rich.deg.plot,
#        device = cairo_pdf,
#        width = 11,height = 9,dpi = 300)

# -------------------------------------------------------------------------



