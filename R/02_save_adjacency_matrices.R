
# store adjacency matrices of each network
# look at the options for unipartite matrices, 
# in line 22. They can be directed/undirected, with/out diagonal self-links.

library(tidyverse)
source("R/aux_matrix_from_edge_list.R")

# -------------------------------------------------------------------------

net.collection <- read_csv2("results/network_collection.csv")
link.collection <- read_csv2("results/network_links_collection.csv")

net.ids <- na.omit(unique(link.collection$network_id))

for(i.id in 1:length(net.ids)){
  partit <- net.collection$network_topology_type[which(net.collection$network_id == net.ids[i.id])]
  my.el <- subset(link.collection,network_id == net.ids[i.id])
  if(partit == "bipartite"){
    my.adjmat <- matrix_from_edge_list(my.el,unipartite = F)
    my.file.name <- paste("data/bipartite_adjacency_matrices/",net.ids[i.id],".csv",sep="")
  }else{
    my.adjmat <- matrix_from_edge_list(edge.list = my.el,unipartite = T,directed = F,diagonal = F,weighted = F)
    my.file.name <- paste("data/unipartite_adjacency_matrices/",net.ids[i.id],".csv",sep="")
  }
  
  write.csv(my.adjmat,file = my.file.name)
  
}

