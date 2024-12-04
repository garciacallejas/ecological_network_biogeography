
# curate data from Fricke and Svenning 2020
# https://www.nature.com/articles/s41586-020-2640-y

# -------------------------------------------------------------------------

library(tidyverse)

# gather networks from Fricke & Svenning 2020 
# 10.1038/s41586-020-2640-y

load(file = "data/raw_interaction_datasets/Fricke_2020/homogenization.RData")

# A 'long' version of the dataset (1 row per local network by species 
# combination) is called net.long0
# net.long is the subset of this with just interactions among species (not unidentified taxa)

# note: net.long keeps taxa identified at any level, not just species.
# net.long0 adds unidentified nodes.
# so, use net.long

fs <- net.long

# fs$network_id <- paste("FS20_",1:nrow(fs),sep="")
fs$network_db <- "FS20"

fs <- fs[,c("net.id","year","network_db","latitude","longitude","animal.accepted.species","plant.accepted.species",
            "target.bats","target.birds","target.primates","target.mamm.carns","target.mamm.herbs","target.mamm.others",
            "target.fish","target.herps")]

fs <- fs[complete.cases(fs),]

# before selecting columns, add the animal guild
fs$main_animal_taxa <- NA

for(i.net in 1:nrow(fs)){
  if(fs$target.bats[i.net] != "no"){
    fs$main_animal_taxa[i.net] <- "Bats"
  }else if(fs$target.birds[i.net] != "no"){
    fs$main_animal_taxa[i.net] <- "Birds"
  }else if(fs$target.fish[i.net] != "no"){
    fs$main_animal_taxa[i.net] <- "Fish"
  }else if(fs$target.herps[i.net] != "no"){
    fs$main_animal_taxa[i.net] <- "Reptiles"
  }else if(fs$target.mamm.carns[i.net] != "no" | fs$target.mamm.herbs[i.net] != "no" | 
           fs$target.mamm.others[i.net] != "no" | fs$target.primates[i.net] != "no"){
    fs$main_animal_taxa[i.net] <- "Mammals"
  }
}

# networks
fs.nets <- fs[,c("net.id","network_db","year","latitude","longitude","main_animal_taxa")]
names(fs.nets) <- c("original_id","network_db","network_year","network_lat","network_lon","main_animal_taxa")
fs.nets$network_topology_type <- "bipartite"
fs.nets$network_spatial_type <- "Point"
fs.nets$network_interaction <- "frugivory"
fs.nets$habitat_type <- NA

fs.nets <- arrange(unique(fs.nets),original_id)
fs.nets$network_id <- paste("FS20_",formatC(1:nrow(fs.nets), width = 3, format = "d", flag = "0"),sep="")
fs.nets <- fs.nets[,c("network_id","original_id","network_db","network_interaction",
                      "main_animal_taxa",
                      "network_year","network_lat","network_lon",
                      "habitat_type",
                      "network_topology_type","network_spatial_type")]

# nodes
all.nodes <- sort(unique(c(fs$animal.accepted.species,fs$plant.accepted.species)))
fs.nodes <- data.frame(node_id = paste("FS20_node_",1:length(all.nodes),sep=""),
                       taxonomy_name = all.nodes,taxonomy_rank = "species")

# links
fs.links <- fs
fs.links$network_id <- fs.nets$network_id[match(fs.links$net.id,fs.nets$original_id)]
fs.links <- fs.links[,c("network_id","animal.accepted.species","plant.accepted.species")]
names(fs.links)[2:3] <- c("node_from","node_to")
fs.links$link_type <- "frugivory"
fs.links$link_direction <- "directed"
fs.links$link_value <- NA_real_
fs.links$link_units <- NA_character_

# rename links
fs.links$node_from.new <- fs.nodes$node_id[match(fs.links$node_from,fs.nodes$taxonomy_name)]
fs.links$node_to.new <- fs.nodes$node_id[match(fs.links$node_to,fs.nodes$taxonomy_name)]

fs.links$node_from <- fs.links$node_from.new
fs.links$node_to <- fs.links$node_to.new

fs.links$node_from.new <- fs.links$node_to.new <- NULL

# save to disk
write.csv2(fs.nets,file = "data/fs_networks.csv",row.names = FALSE)
write.csv2(fs.nodes,file = "data/fs_nodes.csv",row.names = FALSE)
write.csv2(fs.links,file = "data/fs_links.csv",row.names = FALSE)


