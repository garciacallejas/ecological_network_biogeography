
# create replicates from a null model only fixing connectance

library(tidyverse)
library(vegan)
library(bipartite)

# -------------------------------------------------------------------------

unipartite.path <- "data/unipartite_adjacency_matrices/"
bipartite.path <- "data/bipartite_adjacency_matrices/"
uni.matrices.list <- list.files(path = unipartite.path, pattern = ".csv")
bi.matrices.list <- list.files(path = bipartite.path, pattern = ".csv")

# to match the deg.dist null model, that is more computationally demanding
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
  
  # this works equally well with unipartite matrices, so use it
  # for compatibility
  null_networks <- bipartite::shuffle.web(my.matrix,N = num.replicates,legacy = F)
  
  # null_networks <- permatfull(my.matrix,
  #                             times = num.replicates,
  #                             fixedmar = "none",
  #                             shuffle = "both",
  #                             mtype = "prab")$perm
  cat("done\n")
  
  # store each null separately
  for(i.null in 1:length(null_networks)){
    
    rownames(null_networks[[i.null]]) <- rownames(my.matrix)
    colnames(null_networks[[i.null]]) <- colnames(my.matrix)
    
    write.csv2(x = null_networks[[i.null]],file = paste("data/null_adjacency_matrices/connectance_null_unipartite_adjacency_matrices/",
                                        my.net.name,"_connectance_sample_",
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
  
  # null_networks <- bipartite::nullmodel(my.matrix, 
  #                                       method = "vaznull", 
  #                                       N = num.replicates) 
  null_networks <- bipartite::shuffle.web(my.matrix,N = num.replicates,legacy = F)
  
  cat("done\n")
  
  # store each null separately
  for(i.null in 1:length(null_networks)){
    
    rownames(null_networks[[i.null]]) <- rownames(my.matrix)
    colnames(null_networks[[i.null]]) <- colnames(my.matrix)
    
    write.csv2(x = null_networks[[i.null]],file = paste("data/null_adjacency_matrices/connectance_null_bipartite_adjacency_matrices/",
                                                        my.net.name,"_connectance_sample_",
                                                        sprintf("%03d", i.null),".csv",sep=""))
  }
  
}# for each bipartite


