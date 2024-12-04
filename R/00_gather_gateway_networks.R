
# curate data from gateway dataset
# https://www.nature.com/articles/s41559-019-0899-x

# -------------------------------------------------------------------------

library(tidyverse)
source("R/aux_metric_functions.R")

# gather networks from the GATEWAY database, downloaded in a local folder
rd <- read.csv("data/raw_interaction_datasets/gateway/283_3_283_2_FoodWebDataBase_2018_12_10.csv")

# nodes
gw.con.nodes <- rd[,c("con.taxonomy","con.taxonomy.level")]
gw.res.nodes <- rd[,c("res.taxonomy","res.taxonomy.level")]

names(gw.con.nodes) <- c("taxonomy_name","taxonomy_rank")
names(gw.res.nodes) <- names(gw.con.nodes)

gw.nodes <- bind_rows(gw.con.nodes,gw.res.nodes)
gw.nodes <-  arrange(unique(gw.nodes),taxonomy_name)
gw.nodes$node_id <- paste("gateway_node_",1:nrow(gw.nodes),sep="")
gw.nodes <- gw.nodes[,c("node_id","taxonomy_name","taxonomy_rank")]


# -------------------------------------------------------------------------
# these networks may have more than one type of link, between those listed in the data
# unique(rd$interaction.type)
# gw.link.types <-  rd %>% group_by(foodweb.name, interaction.type) %>%
#   summarise(links = n())
# for now, mark them as "food web" - eventually, I should refine these, since
# I divide mutualistic ones into plant-pollinator and plant-frugivore networks.

gw.nets <- unique(rd[,c("foodweb.name","sampling.start.year","sampling.end.year","longitude","latitude")])

# for now, I take year as the starting year of the sampling. Most networks have sampling periods <4 years:
# table(gw.nets$sampling.end.year - gw.nets$sampling.start.year)

names(gw.nets)[4:5] <- c("network_lon","network_lat")
gw.nets$network_id <- paste("gateway_",formatC(1:nrow(gw.nets), width = 3, format = "d", flag = "0"),sep="")
gw.nets$network_db <- "gateway"
gw.nets$network_spatial_type <- "Point"
gw.nets$network_interaction <- "food web"
gw.nets$main_animal_taxa <- NA
gw.nets$habitat_type <- NA
gw.nets$network_topology_type <- "unipartite"

gw.nets <- gw.nets[,c("network_id","foodweb.name","network_db","network_interaction",
                      "main_animal_taxa",
                      "sampling.start.year",
                      "network_lat","network_lon",
                      "habitat_type",
                      "network_topology_type",
                      "network_spatial_type")]
names(gw.nets)[2] <- "original_id"
names(gw.nets)[6] <- "network_year"

# main animal taxa
# gather the "metabolic type" of consumers, 
# and simply assign the "main animal taxa" column of the network
# to the one with most observations
# in parallel, assign habitat type - particularly for aquatic communities
for(i.net in 1:nrow(gw.nets)){
  
  my.eco.type <- table(rd$ecosystem.type[which(rd$foodweb.name == gw.nets$original_id[i.net])])
  if(length(my.eco.type) == 1){
    my.type <- names(my.eco.type)
    my.type <- stringr::str_to_title(my.type)
  }else{
    my.type <- names(my.eco.type[which(my.eco.type == max(my.eco.type))])
    my.type <- stringr::str_to_title(my.type)
  }# if-else one ecosystem type
  
  # I will add habitat information when combining datasets - so discard the 
  # "aboveground" type and override it later, 
  # but the other types are actually relevant
  if(my.type == "Terrestrial Aboveground"){
    my.type <- NA
  }
  gw.nets$habitat_type[i.net] <- my.type

# -------------------------------------------------------------------------
  
  my.animals <- table(rd$con.metabolic.type[which(rd$foodweb.name == gw.nets$original_id[i.net])])
  my.animal.type <- names(which(my.animals == max(my.animals)))
  my.animal.type <- stringr::str_to_title(my.animal.type)
  my.animal.type <- paste(my.animal.type,"s",sep="")
  
  gw.nets$main_animal_taxa[i.net] <- my.animal.type
}

# links

gw.links <- rd[,c("con.taxonomy","res.taxonomy","interaction.type","foodweb.name")]
gw.links$link_direction <- "directed"
gw.links$link_value <- NA_real_
gw.links$link_units <- NA_real_

gw.links$node_from <- gw.nodes$node_id[match(gw.links$con.taxonomy,gw.nodes$taxonomy_name)]
gw.links$node_to <- gw.nodes$node_id[match(gw.links$res.taxonomy,gw.nodes$taxonomy_name)]
gw.links$network_id <- gw.nets$network_id[match(gw.links$foodweb.name,gw.nets$original_id)]

gw.links <- gw.links[,c("network_id","node_from","node_to","interaction.type","link_direction","link_value","link_units")]
names(gw.links)[4] <- "link_type"

# subset replicates of gateway networks: 
# there are networks with experiments or treatments, these are to be discarded
# keep only one replicate when there are several in the same spatial coordinates
# that way, we also avoid temporal replicates

# this must be done semi-manually
gw.nets$spatial.coord <- paste(gw.nets$network_lat,gw.nets$network_lon,sep="_")
gw.nets$remove <- FALSE
# this selection removes experimental treatments labeled as such in the original id,
# and replicates with the same spatial coordinates, leaving only one web per coordinates
# (the first one in the dataset is the one kept)
to.remove <- c("001","002","004","005","006","008","020","038",
               "106","111","137","138","146","148","163")
# to search for all these patterns at once
to.remove.pattern <- paste(to.remove,collapse = "|")
gw.nets$remove[grep(to.remove.pattern,gw.nets$network_id)] <- TRUE
gw.nets.2 <- subset(gw.nets,remove == FALSE)
gw.nets.2$remove <- NULL
gw.nets.2$spatial.coord <- NULL

gw.links.2 <- subset(gw.links,network_id %in% gw.nets.2$network_id)

# -------------------------------------------------------------------------
# check topology of each network
for(i.net in 1:nrow(gw.nets.2)){
  my.el <- subset(gw.links, network_id == gw.nets.2$network_id[i.net])
  bip <- bipartiteness.edge.list(my.el)
  if(bip){
    gw.nets.2$network_topology_type[i.net] <- "bipartite"
  }
}

# save to disk
write.csv2(gw.nets.2,file = "data/gateway_networks.csv",row.names = FALSE)
write.csv2(gw.nodes,file = "data/gateway_nodes.csv",row.names = FALSE)
write.csv2(gw.links.2,file = "data/gateway_links.csv",row.names = FALSE)




