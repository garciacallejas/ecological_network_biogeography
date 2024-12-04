
# curate data from Boscolo et al. 2023
# https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3900

library(tidyverse)

# -------------------------------------------------------------------------

ai.full <- read.csv("data/raw_interaction_datasets/AI22/AtlanticForestInvertFloInteractionData_2022-07.csv")

ai.full.clean <- ai.full %>% dplyr::select(latitude_y,longitude_x,site_name_id,
                                           plant_species_complete_name,
                                           invertebrate_species_complete_name,
                                           campain_year_finish,
                                           invertebrate_behavior,invertebrate_order,
                                           invertebrate_family,
                                           invertebrate_species_complete_name) %>%
  filter(!is.na(plant_species_complete_name) & !is.na(invertebrate_species_complete_name))

# -------------------------------------------------------------------------
# observations are in many cases deaggregated for specific host plants
# I need a way to aggregate them into "networks", or local communities.
# A basic way is simply to aggregate them spatially, by lat-lon coincidences
# two decimals should be at most ~1km, around the equator
# three, at most ~100m

ai.full.clean$lat <- round(ai.full.clean$latitude_y,3)
ai.full.clean$lon <- round(ai.full.clean$longitude_x,3)

ai.full.clean$coords <- paste(ai.full.clean$lat,"_",ai.full.clean$lon,sep="")
ai.full.clean$network_year <- ai.full.clean$campain_year_finish

# furthermore, there are both mutualist and herbivore interactions
# these need to be differentiated
# a simple way is by looking at the order of the invertebrate partner:
# diptera/lepidoptera, are pollinators, as well as apoidea families: 
# Melittidae, Apidae, Megachillidae, Andrenidae, 
# Halictidae, Colletidae, Strenotritidae

ai.full.clean$network_interaction <- "herbivory"
families.pol <- c("Melittidae", "Apidae", "Megachillidae", "Andrenidae",
                  "Halictidae", "Colletidae", "Strenotritidae")
orders.pol <- c("Diptera","Lepidoptera")

pol.condition <- which(ai.full.clean$invertebrate_family %in% families.pol |
                         ai.full.clean$invertebrate_order %in% orders.pol)

ai.full.clean$network_interaction[pol.condition] <- "pollination"

ai.full.clean.2 <- ai.full.clean %>% 
  group_by(coords,network_interaction, network_year) %>%
  mutate(network_id = cur_group_id())
ai.full.clean.2$network_id <- paste("AI22_",ai.full.clean.2$network_id,sep="")

# -------------------------------------------------------------------------
# filter by year
# when more than one network in a given location, select the latest
ai.full.clean.3 <- unique(ai.full.clean.2[,c("network_id","lat","lon","network_year","network_interaction")])

ai.full.last.year <- ai.full.clean.3 %>% 
  group_by(lat,lon,network_interaction) %>%
  mutate(max.year = max(network_year))
ai.full.max.year <- subset(ai.full.last.year, network_year == max.year)

# -------------------------------------------------------------------------
ai.nets <- ai.full.max.year[,c("network_id","lat","lon","network_year","network_interaction")]
names(ai.nets)[2:4] <- c("network_lat","network_lon","network_year")

ai.nets$original_id <- ai.nets$network_id
ai.nets$network_db <- "AI22"
ai.nets$network_topology_type <- "bipartite"

# no need to fill this here, I will overlap the habitat type map in the
# "combine_data" script
ai.nets$habitat_type <- NA

ai.nets$network_spatial_type = "Polygon"
ai.nets$main_animal_taxa <- "Invertebrates"

# in the combine_data script I filter by number of interactions, etc
ai.nets <-  unique(ai.nets[,c("network_id","original_id","network_db","network_interaction",
                              "main_animal_taxa",
                       "network_year","network_lat","network_lon",
                       "habitat_type",
                       "network_topology_type","network_spatial_type")])

# -------------------------------------------------------------------------
# links

ai.links <- ai.full.clean.2[which(ai.full.clean.2$network_id %in% unique(ai.nets$network_id)),
                            c("network_id","invertebrate_species_complete_name",
                             "plant_species_complete_name","network_interaction")]
names(ai.links)[2:4] <- c("node_from","node_to","link_type")
ai.links$link_direction <- NA
ai.links$link_value <- NA
ai.links$link_units <- NA

# -------------------------------------------------------------------------
# nodes

all.nodes <- sort(unique(c(ai.links$node_from,ai.links$node_to)))

ai.nodes <- data.frame(node_id = paste("AI22_node_",1:length(all.nodes),sep=""),
                       taxonomy_name = all.nodes,taxonomy_rank = "species")

# -------------------------------------------------------------------------
# rename links
ai.links$node_from.new <- ai.nodes$node_id[match(ai.links$node_from,ai.nodes$taxonomy_name)]
ai.links$node_to.new <- ai.nodes$node_id[match(ai.links$node_to,ai.nodes$taxonomy_name)]

ai.links$node_from <- ai.links$node_from.new
ai.links$node_to <- ai.links$node_to.new

ai.links$node_from.new <- ai.links$node_to.new <- NULL

# -------------------------------------------------------------------------
write.csv2(ai.nets,"data/AI22_networks.csv",row.names = F)
write.csv2(ai.links,"data/AI22_links.csv",row.names = F)
write.csv2(ai.nodes,"data/AI22_nodes.csv",row.names = F)


