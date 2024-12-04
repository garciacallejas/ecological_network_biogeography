
# curate WEB OF LIFE datasets
# https://www.web-of-life.es/

# NOTE: this needs an active internet connection to retrieve data

library(tidyverse)
library(rjson)
library(jsonlite)
source("R/aux_metric_functions.R")

# -------------------------------------------------------------------------

base_url <- "https://www.web-of-life.es/"    

# -------------------------------------------------------------------------

# METADATA of all networks in the DB
wol.nets <- read.csv(paste0(base_url,"get_network_info.php"))

wol.nets$network_db <- "WOL"
wol.nets$main_animal_taxa <- NA
wol.nets$network_year <- as.numeric(substr(wol.nets$author,(nchar(wol.nets$author) - 3),nchar(wol.nets$author)))
wol.nets$network_id <- paste("WOL_",formatC(1:nrow(wol.nets), width = 3, format = "d", flag = "0"),sep="")
# is this true?
wol.nets$network_topology_type <- "bipartite"
wol.nets$network_spatial_type <- "point"
wol.nets$habitat_type <- NA


wol.nets <- wol.nets[,c("network_id",
                      "network_name",
                      "network_db",
                      "network_type",
                      "main_animal_taxa",
                      "network_year",
                      "latitude",
                      "longitude",
                      "habitat_type",
                      "network_topology_type",
                      "network_spatial_type")]

names(wol.nets) <- c("network_id","original_id",
                    "network_db","network_interaction",
                    "main_animal_taxa","network_year",       
                    "network_lat","network_lon",
                    "habitat_type",
                    "network_topology_type","network_spatial_type")
wol.nets$network_interaction <- recode(wol.nets$network_interaction, 
                                      Pollination = "pollination",
                                      `Anemone-Fish` = "other mutualism",
                                      `Food Webs` = "food web",
                                      `Host-Parasite` = "parasitism",
                                      `Plant-Ant` = "other mutualism",
                                      `Plant-Herbivore` = "herbivory",
                                      `Seed Dispersal` = "frugivory")

# -------------------------------------------------------------------------
# sets of links and nodes

json_url <- paste0(base_url,"get_networks.php") 
wol.links <- jsonlite::fromJSON(json_url)

wol.links$network_id <- wol.nets$network_id[match(wol.links$network_name,
                                                  wol.nets$original_id)]
wol.links$network_name <- NULL
wol.links <- wol.links[,c("network_id","species1","species2","connection_strength")]
names(wol.links) <- c("network_id","node_from","node_to","link_value")
wol.links$link_type <- NA
wol.links$link_direction <- NA
wol.links$link_units <- NA

for(i.id in 1:nrow(wol.nets)){
  
  my.id <- wol.nets$network_id[i.id]
  my.type <- wol.nets$network_interaction[i.id]
  wol.links$link_type[wol.links$network_id == my.id] <- my.type
  
  # check bipartiteness
  my.el <- subset(wol.links,network_id == my.id)
  bip <- bipartiteness.edge.list(my.el)
  if(!bip){
    wol.nets$network_topology_type[i.id] <- "unipartite"
  }
}

wol.links <- wol.links[,c("network_id","node_from","node_to",
                          "link_type","link_direction","link_value","link_units")]

all.nodes <- sort(unique(c(wol.links$node_from,wol.links$node_to)))
wol.nodes <- data.frame(node_id = paste("WOL_node_",1:length(all.nodes),sep=""),
                        taxonomy_name = all.nodes,
                        taxonomy_rank = NA)

# rename links
wol.links$node_from.new <- wol.nodes$node_id[match(wol.links$node_from,wol.nodes$taxonomy_name)]
wol.links$node_to.new <- wol.nodes$node_id[match(wol.links$node_to,wol.nodes$taxonomy_name)]

wol.links$node_from <- wol.links$node_from.new
wol.links$node_to <- wol.links$node_to.new

wol.links$node_from.new <- wol.links$node_to.new <- NULL

# -------------------------------------------------------------------------

write.csv2(wol.nets,file = "data/wol_networks.csv",row.names = FALSE)
write.csv2(wol.nodes,file = "data/wol_nodes.csv",row.names = FALSE)
write.csv2(wol.links,file = "data/wol_links.csv",row.names = FALSE)

