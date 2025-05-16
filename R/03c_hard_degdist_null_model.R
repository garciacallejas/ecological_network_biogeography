
# create replicates from a null model fixing node degrees. In this case, 
# the curveball algorithm. This is a hard-constrained null model.

library(tidyverse)
library(vegan)
library(bipartite)

# -------------------------------------------------------------------------

unipartite.path <- "data/unipartite_adjacency_matrices/"
bipartite.path <- "data/bipartite_adjacency_matrices/"
uni.matrices.list <- list.files(path = unipartite.path, pattern = ".csv")
bi.matrices.list <- list.files(path = bipartite.path, pattern = ".csv")

num.replicates <- 100

# -------------------------------------------------------------------------

# unipartite
for(i.net in 1:length(uni.matrices.list)){

  my.net.name <- substr(uni.matrices.list[i.net],start = 1,stop = (nchar(uni.matrices.list[i.net])-4))

  my.net <- read.csv(paste(unipartite.path,uni.matrices.list[i.net],sep=""))
  my.net.data <- my.net[,2:ncol(my.net)]
  my.matrix <- as.matrix(my.net.data)
  rownames(my.matrix) <- my.net$X

  cat("unipartite network",i.net,"-",my.net.name,"...")
  
  uninm <- vegan::nullmodel(my.matrix, "curveball")
  null_array <- simulate(uninm, nsim = num.replicates)
  null_networks <- lapply(1:dim(null_array)[3], function(i) null_array[,,i])
  
  # null_networks <- permatfull(my.matrix,
  #                             times = num.replicates,
  #                             fixedmar = "both",
  #                             shuffle = "both",
  #                             mtype = "prab")$perm
  cat("done\n")
  
  # store each null separately
  for(i.null in 1:length(null_networks)){
    write.csv2(x = null_networks[[i.null]],file = paste("data/null_adjacency_matrices/hard_degdist_null_unipartite_adjacency_matrices/",
                                        my.net.name,"_degdist_sample_",
                                        sprintf("%03d", i.null),".csv",sep=""))
  }

}# for each unipartite


for(i.net in 1:length(bi.matrices.list)){
  
  my.net.name <- substr(bi.matrices.list[i.net],start = 1,stop = (nchar(bi.matrices.list[i.net])-4))
  
  my.net <- read.csv(paste(bipartite.path,bi.matrices.list[i.net],sep=""))
  my.net.data <- my.net[,2:ncol(my.net)]
  my.matrix <- as.matrix(my.net.data)
  rownames(my.matrix) <- my.net$X
  
  cat("bipartite network",i.net,"-",my.net.name,"...")
  
  # null_networks <- bipartite::vaznull(N = num.replicates,web = my.matrix) 
  # null_networks <- bipartite::swap.web(N = num.replicates,web = my.matrix) 
  btnm <- vegan::nullmodel(my.matrix, "curveball")
  null_array <- simulate(btnm, nsim = num.replicates)
  null_networks <- lapply(1:dim(null_array)[3], function(i) null_array[,,i])

  cat("done\n")
  
  # store each null separately
  for(i.null in 1:length(null_networks)){
    write.csv2(x = null_networks[[i.null]],file = paste("data/null_adjacency_matrices/hard_degdist_null_bipartite_adjacency_matrices/",
                                                        my.net.name,"_degdist_sample_",
                                                        sprintf("%03d", i.null),".csv",sep=""))
  }
  
}# for each bipartite


