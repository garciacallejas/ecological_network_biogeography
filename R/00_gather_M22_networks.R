
# curate datasets from Martins et al. 2022
# https://www.nature.com/articles/s41467-022-34355-w

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------
load("data/raw_interaction_datasets/Martins_2022/metadata_networks.Rdata")

# read all networks

net.path <- "data/raw_interaction_datasets/Martins_2022/networks/"
net.list <- list.files(path = net.path,pattern = "*.csv",full.names = T) %>%
  map(read.csv2)

# -------------------------------------------------------------------------

net.names <- sort(unique(metadata_networks$net_id))
net.ids <- paste("M22_",1:length(net.names),sep="")#na.omit(unique(link.collection$network_id))

# networks
m.nets <- data.frame(network_id = net.ids, 
                     original_id = net.names,
                     network_db = "M22",
                     network_interaction = "frugivory",
                     main_animal_taxa = "Birds",
                     network_year = metadata_networks$sampling_years,
                     network_lat = metadata_networks$lat,
                     network_lon = metadata_networks$lon,
                     habitat_type = NA,
                     network_topology_type = "bipartite",
                     network_spatial_type = "Point")
m.nets$network_year <- substr(m.nets$network_year,1,4)

link.list <- list()
for(i.net in 1:length(net.list)){
  tt <- net.list[[i.net]]
  tl <- pivot_longer(tt,cols = c(-X),names_to = "node_from",values_to = "freq") %>%
    filter(freq >0)
  names(tl)[1] <- "node_to" 
  tl$network_id <- m.nets$network_id[i.net]
  tl$link_type <- "frugivory"
  tl$original_id <- NA
  tl <- tl[,c("network_id","node_from","node_to","link_type","original_id")]
  link.list[[length(link.list)+1]] <- tl
}

m.links <- bind_rows(link.list)
m.links$node_from <- str_replace(m.links$node_from,"\\.","_")
m.links$node_to <- str_replace(m.links$node_to," ","_")

all.sp <- sort(unique(c(m.links$node_from,m.links$node_to)))
node.ids <- paste("M22_node_",1:length(all.sp),sep="")
m.nodes <- data.frame(node_id = node.ids,taxonomy_name = all.sp,taxonomy_rank = "species")

m.links$node_from.new <- m.nodes$node_id[match(m.links$node_from,m.nodes$taxonomy_name)]
m.links$node_to.new <- m.nodes$node_id[match(m.links$node_to,m.nodes$taxonomy_name)]

m.links$node_from <- m.links$node_from.new
m.links$node_to <- m.links$node_to.new

m.links$node_from.new <- m.links$node_to.new <- NULL

# -------------------------------------------------------------------------
write.csv2(m.nets,"data/M22_networks.csv",row.names = F)
write.csv2(m.links,"data/M22_links.csv",row.names = F)
write.csv2(m.nodes,"data/M22_nodes.csv",row.names = F)


