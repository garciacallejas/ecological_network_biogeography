
# curate herbivory datasets from Elisa Thebault and Louise Laux

library(tidyverse)

# -------------------------------------------------------------------------
# net.coords <- read.csv2("data/thebault_datasets/thebault_herbivory_datasets/metadataherbi_22_02_2023_clean.csv")
net.coords <- read.delim("data/raw_interaction_datasets/Thebault_unpublished/metadata_herbivory.txt",
                         header = T,sep = "\t")
net.coords$latitude <- as.numeric(net.coords$latitude)
net.coords$longitude <- as.numeric(net.coords$longitude)

link.list <- read.csv2("data/raw_interaction_datasets/Thebault_unpublished/data_herbivory.csv")

# clean up
link.list$insectorder <- stringr::str_to_title(link.list$insectorder)
link.list$insectorder[link.list$insectorder == "Lepidpotera"] <- "Lepidoptera"
link.list$insectorder[link.list$insectorder == "Trysanoptera"] <- "Thysanoptera"
link.list$insectorder[link.list$insectorder %in% c("N.i","Ni")] <- "Insects"

# -------------------------------------------------------------------------
net.coords <- arrange(net.coords,web)
net.names <- net.coords$web

net.ids <- paste("TI",1:length(net.names),sep="")#na.omit(unique(link.collection$network_id))

# main_animal_taxa <- stringr::str_to_title(net.coords$feedingguild1)
# main_animal_taxa_2 <- stringr::str_to_title(net.coords$feedingguild2)

# networks
ti.nets <- data.frame(network_id = net.ids, 
                      original_id = net.names,
                      network_db = "TI22",
                      network_interaction = "herbivory",
                      main_animal_taxa = NA,
                      network_year = net.coords$samplingyear,
                      network_lat = net.coords$latitude,
                      network_lon = net.coords$longitude,
                      habitat_type = NA,
                      network_topology_type = "bipartite",
                      network_spatial_type = "Point")

for(i.net in 1:nrow(ti.nets)){
  my.net <- subset(link.list,web == ti.nets$original_id[i.net])
  orthoptera.only <- net.coords$feedingguild2[net.coords$web == ti.nets$original_id[i.net]]
  
  # there is a nested categorization here. First division is orthoptera/other,
  # and some of these may not have "insectorder" field recorded
  # if the network is not labelled as "Orthoptera",
  # take a look at that "insectorder" field
  
  if(orthoptera.only == "Orthoptera"){
    ti.nets$main_animal_taxa[i.net] <- "Orthoptera"
  }else{
    my.orders <- unique(my.net$insectorder)
    if(length(my.orders) == 1){
      ti.nets$main_animal_taxa[i.net] <- my.orders
    }else{
      ti.nets$main_animal_taxa[i.net] <- "Insects"
    }
  }# if not only orthoptera
}

# some are still blank - no orthoptera as per "feedingguild2", and no "insectorder"
# I assume these are "insects"
ti.nets$main_animal_taxa[ti.nets$main_animal_taxa == ""] <- "Insects"

# nodes
insects <- data.frame(taxonomy_name = paste(link.list$insectgenus,"_",
                                            link.list$insectspecies,sep=""),
                      genus = link.list$insectgenus,
                      species = link.list$insectspecies)

insects <- unique(insects)
toMatch <- c("sp."," sp","_sp")
insects$taxonomy_rank <- ifelse(grepl(paste(toMatch,collapse="|"),insects$taxonomy_name),"genus","species")
insects$node_id <- paste("TI_HERB_",1:nrow(insects),sep="")

plants <- data.frame(taxonomy_name = paste(link.list$plantgenus,"_",
                                           link.list$plantspecies,sep=""),
                     genus = link.list$plantgenus,
                     species = link.list$plantspecies)

plants <- unique(plants)
toMatch <- c("sp."," sp","_sp")
plants$taxonomy_rank <- ifelse(grepl(paste(toMatch,collapse="|"),plants$taxonomy_name),"genus","species")
plants$node_id <- paste("TI_PLANT_",1:nrow(plants),sep="")

ti.nodes <- rbind(plants,insects)
ti.nodes <- ti.nodes[,c("node_id","taxonomy_name","taxonomy_rank")]

# links
ti.links.list <- list()

for(i.net in 1:nrow(ti.nets)){
  my.list <- subset(link.list,web == ti.nets$original_id[i.net])
  
  my.list.2 <- left_join(my.list,insects,by = c("insectgenus" = "genus","insectspecies" = "species"))
  names(my.list.2)[which(names(my.list.2) == "node_id")] <- "node_from"
  my.list.2$taxonomy_name <- NULL
  my.list.2$taxonomy_rank <- NULL
  
  my.list.3 <- left_join(my.list.2,plants,by = c("plantgenus" = "genus","plantspecies" = "species"))
  names(my.list.3)[which(names(my.list.3) == "node_id")] <- "node_to"
  my.list.3$taxonomy_name <- NULL
  my.list.3$taxonomy_rank <- NULL
  
  my.links.clean <- data.frame(network_id = ti.nets$network_id[i.net],
                               node_from = my.list.3$node_from,
                               node_to = my.list.3$node_to,
                               link_type = "herbivory",
                               original_id = ti.nets$original_id[i.net])
  
  # tf.nodes.list[[i.net]] <- my.nodes
  ti.links.list[[i.net]] <- my.links.clean
}

ti.links <- bind_rows(ti.links.list)

# -------------------------------------------------------------------------
write.csv2(ti.nets,"data/ti22_networks.csv",row.names = F)
write.csv2(ti.links,"data/ti22_links.csv",row.names = F)
write.csv2(ti.nodes,"data/ti22_nodes.csv",row.names = F)


