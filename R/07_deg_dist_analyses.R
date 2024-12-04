
# statistical analyses and plots on network metrics and degree distributions

# INPUTS
# - network collection file
# - network metric data from folder results/metrics

# OUTPUTS
# - Figures 4,5 of manuscript
# - associated models and tables

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
# first plot: richness against degree distribution metrics
rich.deg.data <- subset(observed.metrics, metric %in% c("richness","degree.mean",
                                                        "degree.sd",
                                                        "degree.skewness",
                                                        "degree.kurtosis")) %>%
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
# analyse the differences between metrics across interaction types
my.metrics <- unique(degree.metrics.for.plots$metric)

pair.comparisons <- list()
for(i.metric in 1:length(my.metrics)){
  my.data <- subset(degree.metrics.for.plots, metric == my.metrics[i.metric])
  my.model <- glmmTMB(value ~ network_interaction + (1|study_common), data = my.data)
  emmeans.res <- emmeans(my.model, ~network_interaction)
  my.pairs <- as.data.frame(pairs(emmeans.res, adjust = "tukey"))[c("contrast","p.value")]
  
  pair.comparisons[[length(pair.comparisons)+1]] <- tibble(metric = my.metrics[i.metric],my.pairs)
}

pair.comp.df <- bind_rows(pair.comparisons)

# look at these pairwise comparisons and manually add labels
pair.comp.labels <- expand_grid(metric = my.metrics,
                                network_interaction = unique(degree.metrics.for.plots$network_interaction))

pair.comp.labels$ypos <- 0
pair.comp.labels$group <- "a"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n average" & 
                         pair.comp.labels$network_interaction %in% c("frugivory","herbivory","parasitism","pollination")] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n standard deviation" & 
                         pair.comp.labels$network_interaction %in% c("frugivory","herbivory","parasitism","pollination")] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n skewness" & 
                         pair.comp.labels$network_interaction %in% c("herbivory","pollination")] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n skewness" & 
                         pair.comp.labels$network_interaction %in% c("parasitism")] <- "ab"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n kurtosis" & 
                         pair.comp.labels$network_interaction %in% c("herbivory","pollination")] <- "b"
pair.comp.labels$group[pair.comp.labels$metric == "degree distribution\n kurtosis" & 
                         pair.comp.labels$network_interaction %in% c("frugivory","parasitism")] <- "ab"
# -------------------------------------------------------------------------

# plot metrics of the degree distribution
# I do not plot all four metrics together because I log-transform mean, sd, and kurtosis

metric.plot.list <- list()
for(i.metric in 1:length(my.metrics)){
  my.data <- subset(degree.metrics.for.plots, metric == my.metrics[i.metric])
  my.label.data <- subset(pair.comp.labels,metric == my.metrics[i.metric])
  
  if(i.metric != 3){
    
    my.label.data$ypos <- min(log(my.data$value)) - .5
    
    my.plot <- ggplot(my.data, aes(x = network_interaction, 
                                                    y = log(value))) +
      geom_half_point(aes(color = network_interaction),
                      transformation = position_quasirandom(width = 0.1),
                      side = "l", size = 0.5, alpha = 0.5) +
      geom_half_boxplot(aes(fill = network_interaction), side = "r",outlier.size = 0.8) +
      geom_text(data = my.label.data, aes(label = group, y = ypos, x = network_interaction), 
                position = position_dodge(2),
                color = "darkgrey",
                show.legend = FALSE ) +
      
      # geom_boxplot(aes(fill = network_topology_type)) +
      # facet_wrap(~metric, scales = "free") +
      scale_color_OkabeIto() +
      scale_fill_OkabeIto() +
      labs(x = "", y = paste("log(",my.metrics[i.metric],")",sep="")) +
      theme_bw() +
      theme(legend.position="none") +
      theme(strip.background = element_blank())+
      # scale_y_discrete(breaks=NULL, limits=rev) +
      # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
      # guides(color=guide_legend(override.aes = list(shape=21))) + 
      # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank()) +
      NULL
  }else{
    
    my.label.data$ypos <- min(my.data$value) - .5
    
    my.plot <- ggplot(my.data, aes(x = network_interaction, 
                                   y = value)) +
      geom_half_point(aes(color = network_interaction),
                      transformation = position_quasirandom(width = 0.1),
                      side = "l", size = 0.5, alpha = 0.5) +
      geom_half_boxplot(aes(fill = network_interaction), side = "r",outlier.size = 0.8) +
      geom_text(data = my.label.data, aes(label = group, y = ypos, x = network_interaction), 
                position = position_dodge(2),
                color = "darkgrey",
                show.legend = FALSE ) +
      # geom_boxplot(aes(fill = network_topology_type)) +
      # facet_wrap(~metric, scales = "free") +
      scale_color_OkabeIto() +
      scale_fill_OkabeIto() +
      labs(x = "", y = my.metrics[i.metric]) +
      theme_bw() +
      theme(legend.position="none") +
      theme(strip.background = element_blank())+
      # scale_y_discrete(breaks=NULL, limits=rev) +
      # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
      # guides(color=guide_legend(override.aes = list(shape=21))) + 
      # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank()) +
      NULL
  }
  metric.plot.list[[length(metric.plot.list)+1]] <- my.plot
}

degree.metrics.plot <- wrap_plots(metric.plot.list,ncol = 2)
# degree.metrics.plot

# ggsave("results/images/degree_dist_metrics.pdf",degree.metrics.plot,
#        device = cairo_pdf,
#        width = 9,height = 5,dpi = 300)

# -------------------------------------------------------------------------
# correlations for bipartite networks

# deg.met.bi <- subset(degree.metrics.wide, network_topology_type == "bipartite")
# 
# degree.metrics.pairs.plot <- ggpairs(deg.met.bi, columns = 4:7, 
#         ggplot2::aes(colour=network_interaction))
# for(i in 1:degree.metrics.pairs.plot$nrow) {
#   for(j in 1:degree.metrics.pairs.plot$ncol){
#     degree.metrics.pairs.plot[i,j] <- degree.metrics.pairs.plot[i,j] + 
#       scale_fill_OkabeIto(darken = .2) +
#       scale_color_OkabeIto(darken = .2)  
#   }
# }

# -------------------------------------------------------------------------
# other plots for bipartite networks

# mean and sd

# mean.sd.plot <- ggplot(deg.met.bi, aes(x = degree.mean, y = degree.sd)) + 
#   geom_point(aes(color = network_interaction)) +
#   theme_bw() +
#   NULL
# mean.sd.plot

# mean and skewness
# mean.sk.plot <- ggplot(deg.met.bi, aes(x = degree.mean, y = degree.skewness)) + 
#   geom_point(aes(color = network_interaction)) +
#   theme_bw() +
#   NULL
# mean.sk.plot

# mean and kurtosis
# mean.ku.plot <- ggplot(deg.met.bi, aes(x = degree.mean, y = log(degree.kurtosis))) + 
#   geom_point(aes(color = network_interaction)) +
#   theme_bw() +
#   NULL
# mean.ku.plot

# sd and skewness
# sd.sk.plot <- ggplot(deg.met.bi, aes(x = degree.sd, y = degree.skewness)) + 
#   geom_point(aes(color = network_interaction)) +
#   theme_bw() +
#   NULL
# sd.sk.plot

# sd and kurtosis
# sd.ku.plot <- ggplot(deg.met.bi, aes(x = degree.sd, y = log(degree.kurtosis))) + 
#   geom_point(aes(color = network_interaction)) +
#   theme_bw() +
#   NULL
# sd.ku.plot

# skewness and kurtosis
# sk.ku.plot <- ggplot(deg.met.bi, aes(x = degree.skewness, y = degree.kurtosis)) + 
#   geom_point(aes(color = network_interaction)) +
#   facet_wrap(~network_interaction, scales = "free")+
#   theme_bw() +
#   NULL
# sk.ku.plot

# -------------------------------------------------------------------------
# correlations
deg.matrix <- as.matrix(degree.metrics.wide[,c("degree.mean",
                                               "degree.sd","degree.skewness",
                                               "degree.kurtosis")])
deg.cor.test <- cor.mtest(degree.metrics.wide[,c("degree.mean",
                                                 "degree.sd","degree.skewness",
                                                 "degree.kurtosis")], conf.level = 0.95)
deg.cor <- cor(deg.matrix, method = "spearman")

# -------------------------------------------------------------------------
# try regressions
# for degree.mean, degree.sd, and degree.skewness
# -> kurtosis is very strongly correlated with skewness

# type of interaction will be a covariate, 
# to allow us to test if there are differences among them
deg.met.data <- degree.metrics.wide %>%
  left_join(net.collection) %>%
  dplyr::select(network_id, network_interaction, 
                study_common,
         network_topology_type, 
         # network_main_animal_taxa, 
         # habitat_type, 
         degree.mean,
         degree.sd,
         degree.skewness,
         human_footprint_2009, 
         network_sampling_intensity,
         bio_01, bio_04, bio_12, bio_15) %>%
  drop_na()

names(deg.met.data)[which(names(deg.met.data) == "bio_01")] <- "annual_mean_temp"
names(deg.met.data)[which(names(deg.met.data) == "bio_04")] <- "temp_seasonality"
names(deg.met.data)[which(names(deg.met.data) == "bio_12")] <- "annual_mean_prec"
names(deg.met.data)[which(names(deg.met.data) == "bio_15")] <- "prec_seasonality"

# table(deg.met.data$network_interaction)

# animal taxa is very much related to interaction type,so don't use it
# table(deg.met.data$network_main_animal_taxa,deg.met.data$network_interaction)

# unsure what to do with habitat type - some types are barely represented
# and there are many, surely wrong to model them as random factor
# table(deg.met.data$network_interaction, deg.met.data$habitat_type)

# -------------------------------------------------------------------------
# scale variables!
deg.met.scaled <- deg.met.data
deg.met.scaled$human_footprint_2009 <- scales::rescale(deg.met.scaled$human_footprint_2009)
deg.met.scaled$network_sampling_intensity <- scales::rescale(deg.met.scaled$network_sampling_intensity)
deg.met.scaled$annual_mean_temp <- scales::rescale(deg.met.scaled$annual_mean_temp)
deg.met.scaled$temp_seasonality <- scales::rescale(deg.met.scaled$temp_seasonality)
deg.met.scaled$annual_mean_prec <- scales::rescale(deg.met.scaled$annual_mean_prec)
deg.met.scaled$prec_seasonality <- scales::rescale(deg.met.scaled$prec_seasonality)

# just to be sure i don't use any other variable 
deg.met.scaled <- deg.met.scaled %>% dplyr::select(degree.mean,
                                                   degree.sd,
                                                   degree.skewness,
                                                   network_interaction,
                                                   study_common,
                                                   human_footprint_2009, 
                                                   network_sampling_intensity,
                                                   annual_mean_temp, 
                                                   temp_seasonality, 
                                                   annual_mean_prec, 
                                                   prec_seasonality)

predictors <- c("human_footprint_2009","network_sampling_intensity",
                  "annual_mean_temp", "temp_seasonality", "annual_mean_prec", "prec_seasonality")
predictors.names <- c("human footprint","sampling intensity","annual mean temperature",
                      "temperature seasonality","annual mean precipitation",
                      "precipitation seasonality")
interaction.names <- sort(unique(deg.met.scaled$network_interaction))
# -------------------------------------------------------------------------
# mean degree

# DGC - 15 may - deg.log.mean.model seems to work well, and gives qualitatively similar results to the specific analyses
# DGC - Nov 2024 - add study as random factor, it improves AIC in all cases, so keep it

deg.log.mean.model <- lm(log(degree.mean) ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                                                                                 annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality),
                      data = deg.met.scaled)

deg.log.mean.model.randstudy <- glmmTMB::glmmTMB(log(degree.mean) ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                                                                   annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality) + 
                                                   (1|study_common),
                         data = deg.met.scaled)

# AIC(deg.log.mean.model,deg.log.mean.model.randstudy)

tidy.mean.coefs <- broom.mixed::tidy(deg.log.mean.model.randstudy)
# print(xtable::xtable(as.data.frame(broom.mixed::tidy(deg.log.mean.model.randstudy,exponentiate = F)),floating=FALSE,
#       digits = 2,
#       latex.environments=NULL,
#       booktabs=FALSE))
mean.plots <- list()
avg.pred.data.list <- list()
for(i.pred in 1:length(predictors)){
  
  # my.pred <- ggeffects::predict_response(deg.log.mean.model, c(predictors[i.pred],"network_interaction"))
  # my.pred <- as.data.frame(ggeffects::predict_response(deg.log.mean.model, c(predictors[i.pred],"network_interaction")))
  my.pred <- as.data.frame(ggeffects::predict_response(deg.log.mean.model.randstudy, c(predictors[i.pred],"network_interaction")))
  
  my.pred$group <- factor(my.pred$group, levels = c("food web","frugivory","herbivory","parasitism","pollination"))
  my.pred$significant <- NA
  
  my.coefs <- tidy.mean.coefs[grepl(pattern = predictors[i.pred],x = tidy.mean.coefs$term),]
  
  for(i.int in 1:length(interaction.names)){
    my.signif <- ifelse(my.coefs$p.value[grepl(interaction.names[i.int],my.coefs$term)] < 0.05,T,F)
    my.pred$significant[my.pred$group == interaction.names[i.int]] <- my.signif
  }
  
  my.pred$predictor <- predictors.names[i.pred]
  my.pred$response <- "predicted deg. dist. average"
  avg.pred.data.list[[length(avg.pred.data.list)+1]] <- my.pred
  
  my.plot <- ggplot(my.pred,aes(y = predicted, x = x)) + 
    # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
    geom_line(aes(color = group, linetype = significant), linewidth = 1.3) +
    scale_color_OkabeIto(name = "interaction type",guide = "none")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
    scale_linetype_manual(values = c("dashed","solid"), guide = "none") +
    labs(x = paste0(predictors.names[i.pred]," (scaled)"), y = "predicted deg. dist. average") +
    theme_bw() +
    # theme(legend.justification=c(1,-0.01), legend.position=c(0.27,0.65)) +
    NULL
  
  mean.plots[[length(mean.plots)+1]] <- my.plot
}
avg.pred.data <- bind_rows(avg.pred.data.list)
avg.effects.plot <- wrap_plots(mean.plots,ncol = 6) 

# model checks
# remember that to get interpretable coefficients, I need to exponentiante
# and these are interpreted as multiplicative. See:
# https://stats.stackexchange.com/questions/161216/backtransform-coefficients-of-a-gamma-log-glmm

# summary(deg.log.mean.model)
# performance
# performance::check_model(deg.log.mean.model)
# dharma
# mean.simulationOutput <- simulateResiduals(fittedModel = deg.log.mean.model.randstudy)
# plotResiduals(mean.simulationOutput)
# testResiduals(mean.simulationOutput,plot = T)
# testOverdispersion(mean.simulationOutput)

# -------------------------------------------------------------------------
# standard deviation of degree

deg.log.sd.model <- lm(log(degree.sd) ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                                                              annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality),
                           data = deg.met.scaled)

deg.log.sd.model.randstudy <- glmmTMB::glmmTMB(log(degree.sd) ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                                                               annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality) + 
                         (1|study_common),
                       data = deg.met.scaled)

# AIC(deg.log.sd.model,deg.log.sd.model.randstudy)

tidy.sd.coefs <- broom.mixed::tidy(deg.log.sd.model.randstudy)
# print(xtable::xtable(as.data.frame(tidy(deg.log.sd.model.randstudy,exponentiate = F)),floating=FALSE,
#                      digits = 2,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

sd.plots <- list()
sd.pred.data.list <- list()
for(i.pred in 1:length(predictors)){
  
  # my.pred <- ggeffects::predict_response(deg.log.mean.model, c(predictors[i.pred],"network_interaction"))
  my.pred <- as.data.frame(ggeffects::predict_response(deg.log.sd.model.randstudy, c(predictors[i.pred],"network_interaction")))
  my.pred$group <- factor(my.pred$group, levels = c("food web","frugivory","herbivory","parasitism","pollination"))
  my.pred$significant <- NA
  
  my.coefs <- tidy.sd.coefs[grepl(pattern = predictors[i.pred],x = tidy.sd.coefs$term),]
  
  for(i.int in 1:length(interaction.names)){
    my.signif <- ifelse(my.coefs$p.value[grepl(interaction.names[i.int],my.coefs$term)] < 0.05,T,F)
    my.pred$significant[my.pred$group == interaction.names[i.int]] <- my.signif
  }
  
  my.pred$predictor <- predictors.names[i.pred]
  my.pred$response <- "predicted deg. dist. standard deviation"
  sd.pred.data.list[[length(sd.pred.data.list)+1]] <- my.pred
  
  my.plot <- ggplot(my.pred,aes(y = predicted, x = x)) + 
    # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
    geom_line(aes(color = group, linetype = significant), linewidth = 1.3) +
    scale_color_OkabeIto(name = "interaction type",guide = "none")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
    scale_linetype_manual(values = c("dashed","solid"), guide = "none") +
    labs(x = paste0(predictors.names[i.pred]," (scaled)"), y = "predicted deg. dist. standard deviation") +
    theme_bw() +
    # theme(legend.justification=c(1,-0.01), legend.position=c(0.27,0.65)) +
    NULL
  
  sd.plots[[length(sd.plots)+1]] <- my.plot
}

sd.pred.data <- bind_rows(sd.pred.data.list)
sd.effects.plot <- wrap_plots(sd.plots,ncol = 6) 

# model checks

# summary(deg.log.sd.model.randstudy)
# performance
# performance::check_model(deg.log.sd.model.randstudy)
# dharma
# sd.simulationOutput <- simulateResiduals(fittedModel = deg.log.sd.model.randstudy)
# plotResiduals(sd.simulationOutput)
# testResiduals(sd.simulationOutput,plot = T)
# testOverdispersion(sd.simulationOutput)

# plot effects
# plot(allEffects(deg.sd.model))

# note that I do not exponentiate the coefficients 
# print(xtable::xtable(as.data.frame(tidy(deg.sd.model,exponentiate = F)),floating=FALSE,
#                      digits = 2,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# skewness of degree

deg.sk.model <- lm(degree.skewness ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                     annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality),
                   data = deg.met.scaled)

deg.sk.model.randstudy <- glmmTMB::glmmTMB(degree.skewness ~ network_interaction:(human_footprint_2009 + network_sampling_intensity +
                                                            annual_mean_temp + temp_seasonality + annual_mean_prec + prec_seasonality) + 
                               (1|study_common),
                   data = deg.met.scaled)

# AIC(deg.sk.model,deg.sk.model.randstudy)

tidy.sk.coefs <- broom.mixed::tidy(deg.sk.model.randstudy)
# print(xtable::xtable(as.data.frame(tidy(deg.sk.model.randstudy)),floating=FALSE,
#                      digits = 2,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

sk.plots <- list()
sk.pred.data.list <- list()
for(i.pred in 1:length(predictors)){
  
  # my.pred <- ggeffects::predict_response(deg.log.mean.model, c(predictors[i.pred],"network_interaction"))
  my.pred <- as.data.frame(ggeffects::predict_response(deg.sk.model.randstudy, c(predictors[i.pred],"network_interaction")))
  my.pred$group <- factor(my.pred$group, levels = c("food web","frugivory","herbivory","parasitism","pollination"))
  my.pred$significant <- NA
  
  my.coefs <- tidy.sk.coefs[grepl(pattern = predictors[i.pred],x = tidy.sk.coefs$term),]
  
  for(i.int in 1:length(interaction.names)){
    my.signif <- ifelse(my.coefs$p.value[grepl(interaction.names[i.int],my.coefs$term)] < 0.05,T,F)
    my.pred$significant[my.pred$group == interaction.names[i.int]] <- my.signif
  }
  
  my.pred$predictor <- predictors.names[i.pred]
  my.pred$response <- "predicted deg. dist. skewness"
  sk.pred.data.list[[length(sk.pred.data.list)+1]] <- my.pred
  
  my.plot <- ggplot(my.pred,aes(y = predicted, x = x)) + 
    # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
    geom_line(aes(color = group, linetype = significant), linewidth = 1.3) +
    scale_color_OkabeIto(name = "interaction type",guide = "none")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
    scale_linetype_manual(values = c("dashed","solid"), guide = "none") +
    labs(x = paste0(predictors.names[i.pred]," (scaled)"), y = "predicted deg. dist. skewness") +
    theme_bw() +
    # theme(legend.justification=c(1,-0.01), legend.position=c(0.27,0.65)) +
    NULL
  
  sk.plots[[length(sk.plots)+1]] <- my.plot
}

sk.effects.plot <- wrap_plots(sk.plots,ncol = 6) 
sk.pred.data <- bind_rows(sk.pred.data.list)

# model checks

# summary(deg.sk.model.randstudy)
# performance
# performance::check_model(deg.sk.model)
# dharma
# sk.simulationOutput <- simulateResiduals(fittedModel = deg.sk.model.randstudy)
# plotResiduals(sk.simulationOutput)
# testResiduals(sk.simulationOutput,plot = T)
# testOverdispersion(sk.simulationOutput)

# plot effects
# plot(allEffects(deg.sk.model.randstudy))

# print(xtable::xtable(as.data.frame(tidy(deg.sk.model)),floating=FALSE,
#                      digits = 2,
#                      latex.environments=NULL,
#                      booktabs=FALSE))


# -------------------------------------------------------------------------
pred.data <- bind_rows(list(avg.pred.data,sd.pred.data,sk.pred.data))
pred.data$predictor <- factor(pred.data$predictor,levels = predictors.names)
pred.data$response <- factor(pred.data$response, levels = c("predicted deg. dist. average",
                                                            "predicted deg. dist. standard deviation",
                                                            "predicted deg. dist. skewness"))

all.effects.plot <- ggplot(pred.data,aes(y = predicted, x = x)) + 
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
  geom_line(aes(color = group, linetype = significant), linewidth = 1.1) +
  scale_color_OkabeIto()+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  scale_linetype_manual(values = c("dashed","solid"), guide = "none") +
  facet_grid(rows = vars(response),cols = vars(predictor), scales = "free_y") +
  labs(x = "predictor scaled value", y = "") +
  theme_bw() +
  theme(strip.background = element_blank())+
  guides(colour = guide_legend(position = "bottom", title = NULL)) + # place legend inside plot
  # theme(legend.text = element_text(size = 7),
  #       legend.key.spacing.y = unit(.2, "pt"),
  #       legend.box.margin = unit(0,"pt"),
  #       legend.justification.inside = c(.995, .995)) + # top right
  NULL
# all.effects.plot
# all.effects.plot <- avg.effects.plot/sd.effects.plot/sk.effects.plot

# -------------------------------------------------------------------------

# ggsave("results/images/model_effects.pdf",all.effects.plot,
#        device = cairo_pdf,
#        width = 12,height = 7.5,dpi = 300)



