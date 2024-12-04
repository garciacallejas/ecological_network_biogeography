
# curate data from Parra et al. 2022
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13762

# -------------------------------------------------------------------------

library(tidyverse)
source("R/aux_metric_functions.R")
# -------------------------------------------------------------------------
p22.nets <- read.delim("data/raw_interaction_datasets/Parra_2022/networks_full_pollination_dataset.txt")
p22.net.metadata <- read.csv("data/raw_interaction_datasets/Parra_2022/metadata_full_pollination_dataset.csv")

# I checked manually the pdf "data-sources" and noted which ones come from web of life
nets.in.wol <- c(26:31,51,52,71,86,107,108,115:118)
p22.nets.nowol <- subset(p22.nets, !(id_network_aggreg %in% nets.in.wol))
p22.net.metanowol <- subset(p22.net.metadata, webcode %in% unique(p22.nets.nowol$webcode))
p22.net.metanowol$id_network <- p22.nets.nowol$id_network[match(p22.net.metanowol$webcode,p22.nets.nowol$webcode)]
net.ids <- unique(p22.net.metanowol$id_network)

p22.nets.clean <- data.frame(network_id = paste("P22_",formatC(1:length(net.ids), width = 3, format = "d", flag = "0"),sep=""),
                             original_id = net.ids,
                             network_db = "P22",
                             network_interaction = "pollination",
                             main_animal_taxa = NA,
                             network_year = NA,
                             network_lat = as.numeric(p22.net.metanowol$Latitude_dec),
                             network_lon = as.numeric(p22.net.metanowol$Longitude_dec),
                             habitat_type = NA,
                             network_topology_type = "bipartite",
                             network_spatial_type = "point")

for(i.net in 1:nrow(p22.nets.clean)){
  my.net <- subset(p22.nets.nowol,id_network == p22.nets.clean$original_id[i.net])
  my.orders <- unique(my.net$insectorder)
  if(length(my.orders) == 1){
    p22.nets.clean$main_animal_taxa[i.net] <- my.orders
  }else{
    p22.nets.clean$main_animal_taxa[i.net] <- "Insects"
  }
}

p22.pol <- unique(paste(p22.nets.nowol$insectgenus,p22.nets.nowol$insectspecies,sep="_"))
p22.plant <- unique(paste(p22.nets.nowol$plantgenus,p22.nets.nowol$plantspecies,sep="_"))

p22.nodes <- data.frame(node_id = NA,taxonomy_name = c(p22.plant,p22.pol),taxonomy_rank = "species")
p22.nodes$node_id <- paste("P22_node_",1:nrow(p22.nodes),sep="")
order.level <- grepl("Unidentified",x = p22.nodes$taxonomy_name)
order.level.2 <- grepl("n.i.",x = p22.nodes$taxonomy_name,fixed = T)
p22.nodes$taxonomy_rank[order.level | order.level.2] <- "order"

genus.level <- grepl("_sp",x = p22.nodes$taxonomy_name)
genus.level.2 <- grepl("sp.",x = p22.nodes$taxonomy_name, fixed = T)
genus.level.3 <- grepl("_n.i.",x = p22.nodes$taxonomy_name, fixed = T)
p22.nodes$taxonomy_rank[(genus.level | genus.level.2 | genus.level.3) & !order.level] <- "genus"

p22.nets.nowol$node_from <- paste(p22.nets.nowol$insectgenus,p22.nets.nowol$insectspecies,sep="_")
p22.nets.nowol$node_to <- paste(p22.nets.nowol$plantgenus,p22.nets.nowol$plantspecies,sep="_")

p22.links <- p22.nets.nowol[,c("id_network","node_from","node_to")]
p22.links$network_id <- p22.nets.clean$network_id[match(p22.links$id_network,p22.nets.clean$original_id)]

names(p22.links)[1] <- "original_id"
p22.links <- p22.links %>% dplyr::select(network_id,node_from,node_to) %>%
  mutate(link_type = "pollination",
         link_direction = NA,
         link_value = NA,
         link_units = NA)

p22.links$node_from.new <- p22.nodes$node_id[match(p22.links$node_from,p22.nodes$taxonomy_name)]
p22.links$node_to.new <- p22.nodes$node_id[match(p22.links$node_to,p22.nodes$taxonomy_name)]

p22.links$node_from <- p22.links$node_from.new
p22.links$node_to <- p22.links$node_to.new

p22.links$node_from.new <- p22.links$node_to.new <- NULL

# -------------------------------------------------------------------------
write.csv2(p22.nets.clean,"data/p22_networks.csv",row.names = F)
write.csv2(p22.links,"data/p22_links.csv",row.names = F)
write.csv2(p22.nodes,"data/p22_nodes.csv",row.names = F)


