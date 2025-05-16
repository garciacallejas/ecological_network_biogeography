
# network metrics for null connectance adjacency matrices

# -------------------------------------------------------------------------

library(igraph) # for centralization
library(maxnodf) # for nestedness
library(tidyverse)
library(purrr)
# library(MultitrophicFun)
source("R/aux_metric_functions.R")
source("R/aux_matrix_from_edge_list.R")
library(stringr)

set.seed(6373)

# -------------------------------------------------------------------------
library(foreach)
library(doParallel)
#
workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------
#uni- and/or bipartite networks?
calc.unipartite <- F
calc.bipartite <- T

# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")
link.collection <- read_csv2("results/network_links_collection.csv")
net.ids <- na.omit(unique(link.collection$network_id))

unipartite.conn.null.path <- "data/null_adjacency_matrices/connectance_null_unipartite_adjacency_matrices/"
bipartite.conn.null.path <- "data/null_adjacency_matrices/connectance_null_bipartite_adjacency_matrices/"

# unipartite.exceptions <- c("gateway_109","gateway_116","gateway_155","mangal_5013")
# bipartite.exceptions <- c("P22_132","WOL_227")

# -------------------------------------------------------------------------

if(calc.unipartite){
  
  uni.conn.null.matrices.list <- list.files(path = unipartite.conn.null.path, pattern = ".csv")
  
  my.conn.null.metrics <- foreach(i.null = 1:length(uni.conn.null.matrices.list), 
                                  # .combine=comb.fun, 
                                  .packages = c('tidyverse',
                                                "MultitrophicFun",
                                                "igraph",
                                                "maxnodf"
                                                # "infomapecology"
                                                # "intsegration",
                                                # "reshape2"
                                  )) %dopar% 
    {
      
      source("R/aux_metric_functions.R")
      
      # for(i.null in 1:length(uni.null.matrices.list)){
      my.null.net.name.full <- substr(uni.conn.null.matrices.list[i.null],start = 1,
                                      stop = (nchar(uni.conn.null.matrices.list[i.null])-4))
      my.name.parts <- str_split(my.null.net.name.full, "_", simplify=T)
      my.rep <- as.numeric(my.name.parts[length(my.name.parts)])
      my.null.net.name <- paste(my.name.parts[c(1:(length(my.name.parts)-2))],collapse="_")
      
      # NOTE: I saved connectance null replicates with write.csv2
      my.null.net <- read.csv2(paste(unipartite.conn.null.path,uni.conn.null.matrices.list[i.null],sep=""))
      my.null.net.data <- my.null.net[,2:ncol(my.null.net)]
      my.null.matrix <- as.matrix(my.null.net.data)
      rownames(my.null.matrix) <- my.null.net$X
      
      # metrics
      rich <- richness.matrix(my.null.matrix)
      con <- connectance.matrix(my.null.matrix,diago = F)
      ld <- link.density.matrix(my.null.matrix)
      
      deg.shannon <- NA
      try(deg.shannon <- degree.matrix(my.null.matrix,directed = F))
      
      deg.metrics <- NA
      try(deg.metrics <- degree.metrics.matrix(my.null.matrix,directed = F))
      if(inherits(deg.metrics, "list")){
        deg.mean <- deg.metrics[[1]]
        deg.sd <- deg.metrics[[2]]
        deg.skewness <- deg.metrics[[3]]
        deg.kurtosis <- deg.metrics[[4]]
      }else{
        deg.mean <- NA
        deg.sd <- NA
        deg.skewness <- NA
        deg.kurtosis <- NA
      }
      
      mod.betweenness <- NA
      # these take an unreasonable amount of time
      # if(!any(grepl(paste(unipartite.exceptions,collapse = "|"),my.null.net.name))){
      if(nrow(my.null.matrix)<300){
        try(mod.betweenness <- modularity.matrix(my.null.matrix,modularity.type = "edge_betweenness"))
      }
      
      mod.infomap <- NA
      try(mod.infomap <- modularity.matrix(my.null.matrix,modularity.type = "infomap"))
      
      nest <- 0
      if(con < 1){
        nest <- NA
        try(nest <- nestedness.matrix(my.null.matrix))
      }
      
      cent.degree <- NA
      try(cent.degree <- centralization.matrix(my.null.matrix,centrality.type = "degree"))
      
      # cent.betweenness <- NA
      # try(cent.betweenness <- centralization.matrix(my.null.matrix, centrality.type = "betweenness"))
      # 
      cent.eigen <- NA
      try(cent.eigen <- centralization.matrix(my.null.matrix, centrality.type = "eigen"))
      
      avg.overlap <- NA
      try(avg.overlap <- interaction.overlap.matrix(my.null.matrix,method = "horn"))
      if(inherits(avg.overlap,"data.frame")){
        avg.overlap <- mean(avg.overlap$overlap,na.rm = T)
      }
      
      # my.null.metrics[[length(my.null.metrics)+1]] <- data.frame(network_id = my.null.net.name,
      my.met <- data.frame(network_id = my.null.net.name,
                 null_model = "connectance",
                 replicate = my.rep,
                 richness = rich,
                 connectance = con,
                 degree.shannon = deg.shannon,
                 degree.mean = deg.mean,
                 degree.sd = deg.sd,
                 degree.skewness = deg.skewness,
                 degree.kurtosis = deg.kurtosis,
                 modularity.betweenness = mod.betweenness,
                 modularity.infomap = mod.infomap,
                 nestedness = nest,
                 centrality.degree = cent.degree,
                 # centrality.betweenness = cent.betweenness,
                 centrality.eigen = cent.eigen,
                 interaction.overlap = avg.overlap,
                 link.density = ld
      )
      
      my.met.long <- my.met %>% pivot_longer(richness:link.density,
                                             names_to = "metric",
                                             values_to = "value")
      write.csv2(my.met.long,paste("results/connectance_null_unipartite_metrics/network_metrics_",my.null.net.name.full,".csv",sep=""),
                 row.names = F)
      
      return(NULL)
      
    }
  # 
  # my.conn.null.long <- bind_rows(my.conn.null.metrics) %>%
  #   pivot_longer(richness:link.density,names_to = "metric",values_to = "value")
  # write.csv2(my.conn.null.long,paste("results/unipartite_network_metrics_connectance_null.csv",sep=""),row.names = F)
  
  
}

if(calc.bipartite){
  
  bi.conn.null.matrices.list <- list.files(path = bipartite.conn.null.path, pattern = ".csv")
  
  my.conn.null.metrics <- foreach(i.null = 1:length(bi.conn.null.matrices.list), 
                                  # .combine=comb.fun, 
                                  .packages = c('tidyverse',
                                                "MultitrophicFun",
                                                "igraph",
                                                "maxnodf"
                                                # "infomapecology"
                                                # "intsegration",
                                                # "reshape2"
                                  )) %dopar% 
    {
      
      source("R/aux_metric_functions.R")
      
      # for(i.null in 1:length(uni.null.matrices.list)){
      my.null.net.name.full <- substr(bi.conn.null.matrices.list[i.null],start = 1,stop = (nchar(bi.conn.null.matrices.list[i.null])-4))
      my.name.parts <- str_split(my.null.net.name.full, "_", simplify=T)
      my.rep <- as.numeric(my.name.parts[length(my.name.parts)])
      my.null.net.name <- paste(my.name.parts[c(1:(length(my.name.parts)-2))],collapse="_")
      
      # NOTE: I saved connectance null replicates with write.csv2
      my.null.net <- read.csv2(paste(bipartite.conn.null.path,bi.conn.null.matrices.list[i.null],sep=""))
      my.null.net.data <- my.null.net[,2:ncol(my.null.net)]
      my.null.matrix <- as.matrix(my.null.net.data)
      rownames(my.null.matrix) <- my.null.net$X
      
      # metrics
      rich <- richness.matrix(my.null.matrix)
      con <- connectance.matrix(my.null.matrix,diago = F)
      ld <- link.density.matrix(my.null.matrix)
      
      deg.shannon <- NA
      try(deg.shannon <- degree.matrix(my.null.matrix,directed = F))
      
      deg.metrics <- NA
      try(deg.metrics <- degree.metrics.matrix(my.null.matrix,directed = F))
      if(inherits(deg.metrics, "list")){
        deg.mean <- deg.metrics[[1]]
        deg.sd <- deg.metrics[[2]]
        deg.skewness <- deg.metrics[[3]]
        deg.kurtosis <- deg.metrics[[4]]
      }else{
        deg.mean <- NA
        deg.sd <- NA
        deg.skewness <- NA
        deg.kurtosis <- NA
      }
      
      mod.betweenness <- NA
      
      # some networks take exceptionally long for this metric
      # if(!any(grepl(paste(bipartite.exceptions,collapse = "|"),my.null.net.name))){
      if(nrow(my.null.matrix)<300){
        try(mod.betweenness <- modularity.matrix(my.null.matrix,modularity.type = "edge_betweenness"))
      }
      
      mod.infomap <- NA
      try(mod.infomap <- modularity.matrix(my.null.matrix,modularity.type = "infomap"))
      
      nest <- 0
      if(con < 1){
        nest <- NA
        try(nest <- nestedness.matrix(my.null.matrix))
      }
      
      cent.degree <- NA
      try(cent.degree <- centralization.matrix(my.null.matrix,centrality.type = "degree"))
      
      # cent.betweenness <- NA
      # try(cent.betweenness <- centralization.matrix(my.null.matrix, centrality.type = "betweenness"))
      # 
      cent.eigen <- NA
      try(cent.eigen <- centralization.matrix(my.null.matrix, centrality.type = "eigen"))
      
      avg.overlap <- NA
      try(avg.overlap <- interaction.overlap.matrix(my.null.matrix,method = "horn"))
      if(inherits(avg.overlap,"data.frame")){
        avg.overlap <- mean(avg.overlap$overlap,na.rm = T)
      }
      
      # my.null.metrics[[length(my.null.metrics)+1]] <- data.frame(network_id = my.null.net.name,
      my.met <- data.frame(network_id = my.null.net.name,
                 null_model = "connectance",
                 replicate = my.rep,
                 richness = rich,
                 connectance = con,
                 degree.shannon = deg.shannon,
                 degree.mean = deg.mean,
                 degree.sd = deg.sd,
                 degree.skewness = deg.skewness,
                 degree.kurtosis = deg.kurtosis,
                 modularity.betweenness = mod.betweenness,
                 modularity.infomap = mod.infomap,
                 nestedness = nest,
                 centrality.degree = cent.degree,
                 # centrality.betweenness = cent.betweenness,
                 centrality.eigen = cent.eigen,
                 interaction.overlap = avg.overlap,
                 link.density = ld
      )
      
      my.met.long <- my.met %>% pivot_longer(richness:link.density,
                                             names_to = "metric",
                                             values_to = "value")
      write.csv2(my.met.long,paste("results/connectance_null_bipartite_metrics/network_metrics_",my.null.net.name.full,".csv",sep=""),
                 row.names = F)
      
      return(NULL)
      
    }
  
  # my.conn.null.long <- bind_rows(my.conn.null.metrics) %>%
  #   pivot_longer(richness:link.density,names_to = "metric",values_to = "value")
  # write.csv2(my.null.long,paste("results/bipartite_network_metrics_connectance_null.csv",sep=""),row.names = F)
  
}# if bipartite

stopCluster(cl)
