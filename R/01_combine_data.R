
# combine datasets from different sources
# and add environmental information

# IMPORTANT NOTE 
# this script needs additional data not included in the repo because of its size
# in particular, 

# 1 - the global habitat data types from Jung et al. 2020 (see around line 120)
# 2 - worldclim environmental data (see around line 545)
# 3 - human footprint global data from Venter et al. 2016 (see around line 590)

# download these data (see the relevant lines in the script) to folders of your choice and modify the following variables
# which are defined below:
# habitat.path
# env.data.path
# hf.data.path

# -------------------------------------------------------------------------

library(tidyverse)
library(raster)
library(taxadb)
library(taxize)
library(igraph)
library(geodata)
library(terra)
# library(rgdal)
library(sf)
library(stars)
library(stringi)
source("R/aux_matrix_from_edge_list.R")
source("R/aux_metric_functions.R")

# -------------------------------------------------------------------------
# filters

min.richness <- 20

# -------------------------------------------------------------------------

# link data
fs.links <- read.csv2("data/fs_links.csv")
gateway.links <- read.csv2("data/gateway_links.csv")
mangal.links <- read.csv2("data/mangal_links.csv")
wol.links <- read.csv2("data/wol_links.csv")
p22.links <- read.csv2("data/p22_links.csv")
ti22.links <- read.csv2("data/ti22_links.csv")
m22.links <- read.csv2("data/M22_links.csv")
ai22.links <- read.csv2("data/AI22_links.csv")
d21.links <- read.csv2("data/D21_links.csv")
b22.links <- read.csv2("data/B22_links.csv")

ai22.links$network_id <- as.character(ai22.links$network_id)

all.links <- bind_rows(fs.links[,1:4],
                       gateway.links[,1:4],
                       mangal.links[,1:4],
                       wol.links[,1:4],
                       p22.links[1:4],
                       ti22.links[1:4],
                       m22.links[1:4],
                       ai22.links[1:4],
                       d21.links[1:4],
                       b22.links[1:4])

# -------------------------------------------------------------------------
# node data
fs.nodes <- read.csv2("data/fs_nodes.csv")
gateway.nodes <- read.csv2("data/gateway_nodes.csv")
mangal.nodes <- read.csv2("data/mangal_nodes.csv")
wol.nodes <- read.csv2("data/wol_nodes.csv")
p22.nodes <- read.csv2("data/p22_nodes.csv")
ti22.nodes <- read.csv2("data/ti22_nodes.csv")
m22.nodes <- read.csv2("data/M22_nodes.csv")
ai22.nodes <- read.csv2("data/AI22_nodes.csv")
d21.nodes <- read.csv2("data/D21_nodes.csv")
b22.nodes <- read.csv2("data/B22_nodes.csv")

all.nodes <- bind_rows(fs.nodes,gateway.nodes,mangal.nodes,wol.nodes,p22.nodes,
                       ti22.nodes,m22.nodes,ai22.nodes,d21.nodes,b22.nodes)

# -------------------------------------------------------------------------
# network data

fs.nets <- read.csv2("data/fs_networks.csv")
gw.nets <- read.csv2("data/gateway_networks.csv")
mg.nets <- read.csv2("data/mangal_networks.csv")
wol.nets <- read.csv2("data/wol_networks.csv")
p22.nets <- read.csv2("data/p22_networks.csv")
ti22.nets <- read.csv2("data/ti22_networks.csv")
m22.nets <- read.csv2("data/M22_networks.csv")
ai22.nets <- read.csv2("data/AI22_networks.csv")
d21.nets <- read.csv2("data/D21_networks.csv")
b22.nets <- read.csv2("data/B22_networks.csv")

fs.nets$original_id <- as.character(fs.nets$original_id)
gw.nets$original_id <- as.character(gw.nets$original_id)
mg.nets$original_id <- as.character(mg.nets$original_id)
wol.nets$original_id <- as.character(wol.nets$original_id)
p22.nets$original_id <- as.character(p22.nets$original_id)
ti22.nets$original_id <- as.character(ti22.nets$original_id)
m22.nets$original_id <- as.character(m22.nets$original_id)
ai22.nets$original_id <- as.character(ai22.nets$original_id)
d21.nets$original_id <- as.character(d21.nets$original_id)
b22.nets$original_id <- as.character(b22.nets$original_id)

ai22.nets$network_id <- as.character(ai22.nets$network_id)

all.nets <- bind_rows(list(fs.nets,gw.nets,mg.nets,
                           wol.nets,p22.nets,ti22.nets,
                           m22.nets,ai22.nets,d21.nets,b22.nets))

# -------------------------------------------------------------------------
# include standardised citation
citation_list <- read.csv2("data/dataset_references/citation_list.csv")

all.nets.2 <- left_join(all.nets,citation_list[,c("network_db","original_id","study_common")])
# tt <- subset(all.nets.2, is.na(study_common))
# -------------------------------------------------------------------------
# include habitat type
# NOTE: this data is not included in the github repo
# it can be obtained from Jung et al. 2020
# https://zenodo.org/records/4058819

# DOWNLOAD the tif file "iucn_habitatclassification_composite_lvl1_ver004.tif"
# and the associated codes "habitat_codes.csv"

# to a folder of your choice, specified in "habitat.path"
habitat.path <- "/media/david/Elements SE/Work/datasets/global_habitat_types/"

r <- raster(paste0(habitat.path,"iucn_habitatclassification_composite_lvl1_ver004.tif"))
codes <- read.csv2(paste0(habitat.path,"habitat_codes.csv"))

coords <- all.nets[,c("network_lat","network_lon")]
coords <- coords[which(!is.na(coords$network_lat) & !is.na(coords$network_lon)),]
coords <- SpatialPoints(coords[, c('network_lon', 'network_lat')], 
                        proj4string=CRS('+proj=longlat +datum=WGS84'))

hab <- data.frame(coordinates(coords),extract(r, coords))
names(hab) <- c("network_lon","network_lat","habitat_code")
hab$habitat_type2 <- codes$IUCNLevel[match(hab$habitat_code,codes$NewCode)]
hab <- unique(hab[,c("network_lon","network_lat","habitat_type2")])
nets.hab <- left_join(all.nets.2,hab)
nets.hab$habitat_type2[which(!is.na(nets.hab$habitat_type))] <- nets.hab$habitat_type[which(!is.na(nets.hab$habitat_type))]
nets.hab$habitat_type <- nets.hab$habitat_type2
nets.hab$habitat_type2 <- NULL

# -------------------------------------------------------------------------
# reclassify pollinator and frugivory networks sampled in islands whose coordinates 
# show marine habitat type

pol.is <- which(nets.hab$network_interaction == "pollination" & 
                  (nets.hab$habitat_type == "11 Marine Deep Ocean Floor (Benthic and Demersal)"|
                  nets.hab$habitat_type == "9. Marine Neritic"))
nets.hab$habitat_type[pol.is] <- "Islands"

fru.is <- which(nets.hab$network_interaction == "frugivory" & 
                  (nets.hab$habitat_type == "11 Marine Deep Ocean Floor (Benthic and Demersal)"|
                     nets.hab$habitat_type == "9. Marine Neritic"))
nets.hab$habitat_type[fru.is] <- "Islands"

# -------------------------------------------------------------------------
# remove duplicated records, mostly between mangal and wol, and records without
# coordinates

wol.remove.codes <- as.character(c(110,115,131,126,122,099,133,039,029,135,100,109,121,137,132,106,119,101,102,009,141,
                      105,013,014,015,175:222,067,023,022,229,034,040,044,082,
                      049,060,033,012,011,005,021,006,024,061,074,045,083,080,027,084,144,
                      043,081,116,117,134,108,129,095,128,097,096,140,089,143,124,114,123,112,127,107,113,
                      139,125,025,136,142,118,111,044,082,104,103,138,130,098,052,094))
mangal.remove.codes <- as.character(c(13,946,12,931,944,88))
to.remove.wol <- paste("WOL_",wol.remove.codes,sep="")
to.remove.mangal <- paste("mangal_",mangal.remove.codes,sep="")

to.remove.coords <- nets.hab$network_id[which(is.na(nets.hab$network_lat) |
                                                is.na(nets.hab$network_lon))]

to.remove <- unique(c(to.remove.wol,to.remove.mangal,to.remove.coords))

# tt <- subset(nets.hab,network_id %in% to.remove)
nets.clean <- subset(nets.hab,!network_id %in% to.remove)
nets.clean <- unique(nets.clean)
nets.clean$habitat_type <- recode(nets.clean$habitat_type,
                                  "1. Forest" = "Forest",
                                  "11 Marine Deep Ocean Floor (Benthic and Demersal)" = "Marine Benthic",
                                  "12 Marine Intertidal" = "Marine Intertidal",
                                  "14 Artificial - Terrestrial" = "Anthropic habitats",
                                  "2. Savanna" = "Savanna",
                                  "3. Shrubland" = "Shrubland",
                                  "4. Grassland" = "Grassland",
                                  "5. Wetlands (inland)" = "Wetland",
                                  "6. Rocky Areas (e.g., inland cliffs, mountain peaks)" = "Rocky area",
                                  "8. Desert" = "Desert",
                                  "9. Marine Neritic" = "Marine Neritic",
                                  "Marine" = "Marine unspecified"
                                  )
# also, in general, remove records with duplicate coordinates
# these are mostly samplings from different years or treatments in given coordinates
nets.clean <- dplyr::arrange(nets.clean,network_db,desc(network_year))
nets.clean$dup.coords <- duplicated(nets.clean[,c("network_lat","network_lon")])
nets.clean.2 <- subset(nets.clean,!dup.coords)
nets.clean.2$dup.coords <- NULL

# further checks of matching habitat and animal guilds
# there are a couple of fish networks in savanna or shrubland
# somehow arbitrarily, but assign them to streams
nets.clean.2$habitat_type[which(nets.clean.2$main_animal_taxa == "Fish" & 
                                  !(nets.clean.2$habitat_type %in% c("Lakes","Marine Benthic",
                                                                     "Marine Intertidal","Marine Neritic",
                                                                     "Marine unspecified","Streams","Wetland")))] <- "Streams"

# dup.coords <- nets.clean[which(duplicated(nets.clean[,c("network_lat","network_lon")])),]
nets.clean.2$main_animal_taxa_2 <- nets.clean.2$main_animal_taxa

insects <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")
vertebrates <- c("Ectotherm Vertebrates","Endotherm Vertebrates")

nets.clean.2$main_animal_taxa_2[nets.clean.2$main_animal_taxa %in% insects] <- "Insects"
nets.clean.2$main_animal_taxa_2[nets.clean.2$main_animal_taxa %in% vertebrates] <- "Vertebrates"

# -------------------------------------------------------------------------
# a few (16) unipartite networks are listed as "herbivory" ones. Assign them to food webs.
nets.clean.2$network_interaction[nets.clean.2$network_interaction == "herbivory" & 
                                   nets.clean.2$network_topology_type == "unipartite"] <- "food web"

# -------------------------------------------------------------------------
# OTHER INFO

# richness and number of components

# here, specifically perform the following checks:
# 1 - remove fully connected networks
# 2 - in networks with more than one component, keep the largest one
# 3 - richness higher than minimum

# 1 - connectance

nets.clean.2$connectance <- NA

for(i.net in 1:nrow(nets.clean.2)){
  
  # select the network edge list
  my.net <- subset(all.links, network_id == nets.clean.2$network_id[i.net])
  
  # is it uni or bipartite?
  part <- nets.clean.2$network_topology_type[i.net]
  
  if(nrow(my.net)>0){
    
    if(part == "bipartite"){
      my.mat <- matrix_from_edge_list(my.net,weighted = F,unipartite = F)
    }else{
      my.mat <- matrix_from_edge_list(my.net,weighted = F,unipartite = T, 
                                      directed = F,diagonal = T)
    }
    
    nets.clean.2$connectance[i.net] <- connectance.matrix(my.mat)
  }
  
}

# fullcon <- subset(nets.clean.2, connectance == 1)

nets.clean.3 <- subset(nets.clean.2, !is.na(connectance) & connectance < 1)

# -------------------------------------------------------------------------
# 2 - richness

nets.clean.3$richness <- NA
nets.clean.3$bipartite.richness.smallest.layer <- NA

for(i.net in 1:nrow(nets.clean.3)){
  
  # select the network edge list
  my.net <- subset(all.links, network_id == nets.clean.3$network_id[i.net])
  
  # is it uni or bipartite?
  part <- nets.clean.3$network_topology_type[i.net]
  
  num.nodes.from <- length(unique(my.net$node_from))
  num.nodes.to <- length(unique(my.net$node_to))
  nets.clean.3$richness[i.net] <- length(unique(c(my.net$node_from,my.net$node_to)))
  
  if(part == "bipartite"){
    nets.clean.3$bipartite.richness.smallest.layer[i.net] <- min(num.nodes.from,num.nodes.to)
  }
}

# richness > min, and (either unipartite, or bipartite and richness of the smallest layer > min.richness/2)
nets.clean.4 <- subset(nets.clean.3, richness >= min.richness & (network_topology_type == "unipartite" | 
                                                                   (network_topology_type == "bipartite" & 
                                                                      bipartite.richness.smallest.layer >= (min.richness/2))
  ))

# -------------------------------------------------------------------------
# UPDATE: keep all components in all networks, as suggested by Elisa Thebault
# I will subsequently update the metrics
nets.clean.5 <- nets.clean.4
all.links.2 <- all.links
# 3 - for networks with more than one component, keep only the principal one
# 3.1 ensure that the principal component satisfies richness/connectance constraints

# first, get the number of components of all networks
# nets.clean.4$num_components <- NA
# 
# for(i.net in 1:nrow(nets.clean.4)){
#   my.net <- subset(all.links, network_id == nets.clean.4$network_id[i.net])
# 
#   my.graph <- igraph::graph_from_data_frame(my.net[,c("node_from","node_to")])
#   nets.clean.4$num_components[i.net] <- igraph::count_components(my.graph,mode = "weak")
# }
# 
# # now, only for those with >1 component, prune them
# multicomponent.nets <- subset(nets.clean.4, num_components > 1)
# 
# # flags for the tests of connectance and richness
# multicomponent.nets$principal.connectance <- NA
# multicomponent.nets$principal.richness <- NA
# multicomponent.nets$bipartite.principal.richness.smallest <- NA
# 
# # this is the edge list of the principal components
# multi.nets.to.keep <- list()
# 
# for(i.mult in 1:nrow(multicomponent.nets)){
#   
#   my.net <- subset(all.links, network_id == multicomponent.nets$network_id[i.mult])
#   my.graph <- igraph::graph_from_data_frame(my.net[,c("node_from","node_to")])
#   
#   # get components
#   my.components <- igraph::clusters(my.graph, mode="weak")
#   
#   # which is the biggest component?
#   biggest_component <- which.max(my.components$csize)
#   
#   # Extract the membership vector for nodes
#   membership <- my.components$membership
#   
#   # Create a data frame with node IDs and their corresponding component IDs
#   node_comp <- data.frame(
#     node = V(my.graph)$name,
#     component = membership
#   )
#   
#   # get the names of the nodes in the principal component
#   nodes.to.keep <- subset(node_comp, component == biggest_component)
#   node.names.to.keep <- sort(unique(nodes.to.keep$node))
#   
#   # this is the principal component in edge list format
#   my.principal.net <- my.net[which(my.net$node_from %in% node.names.to.keep),]
#   
#   multi.nets.to.keep[[length(multi.nets.to.keep)+1]] <- my.principal.net
#   
#   # now, check connectance and richness of it
#   # is it uni or bipartite?
#   part <- multicomponent.nets$network_topology_type[i.mult]
#   
#   if(nrow(my.principal.net)>0){
#     
#     if(part == "bipartite"){
#       my.mat <- matrix_from_edge_list(my.principal.net,weighted = F,unipartite = F)
#     }else{
#       my.mat <- matrix_from_edge_list(my.principal.net,weighted = F,unipartite = T, 
#                                       directed = F,diagonal = T)
#     }
#     
#     multicomponent.nets$principal.connectance[i.mult] <- connectance.matrix(my.mat)
#   }
#   
#   num.nodes.from <- length(unique(my.principal.net$node_from))
#   num.nodes.to <- length(unique(my.principal.net$node_to))
#   multicomponent.nets$principal.richness[i.mult] <- length(unique(c(my.principal.net$node_from,my.principal.net$node_to)))
#   
#   if(part == "bipartite"){
#     multicomponent.nets$bipartite.principal.richness.smallest[i.mult] <- min(num.nodes.from,num.nodes.to)
#   }
#   
# }
# 
# # subset - same conditions as above
# multicomp.nets.to.keep <- subset(multicomponent.nets, principal.connectance < 1 & 
#                                    principal.richness >= min.richness & 
#                                    (network_topology_type == "unipartite" | 
#                                       (network_topology_type == "bipartite" & 
#                                          bipartite.principal.richness.smallest >= (min.richness/2)))
#                                    ) %>%
#   dplyr::select(names(nets.clean.4))
# 
# # edge list of all the principal components
# multicomp.edgelist <- bind_rows(multi.nets.to.keep)
# multicomp.edgelist.to.keep <- subset(multicomp.edgelist, network_id %in% multicomp.nets.to.keep$network_id)
# 
# # now, substitute the original edge lists of the nets with multiple components
# # by the pruned versions
# 
# # first, remove all links from multicomponent networks
# pruned.links <- subset(all.links, !(network_id %in% multicomponent.nets$network_id))
# 
# # then, put back the links of the principal components kept
# all.links.2 <- bind_rows(pruned.links, multicomp.edgelist.to.keep)
# 
# # also, update the network dataframe
# nets.clean.4.2 <- subset(nets.clean.4, !(network_id %in% multicomponent.nets$network_id))
# nets.clean.5 <- bind_rows(nets.clean.4.2, multicomp.nets.to.keep)
# 
# # -------------------------------------------------------------------------
# # recheck everything to be sure
# 
# nets.clean.5$num_components <- NA
# 
# for(i.net in 1:nrow(nets.clean.5)){
#   my.net <- subset(all.links.2, network_id == nets.clean.5$network_id[i.net])
#   
#   my.graph <- igraph::graph_from_data_frame(my.net[,c("node_from","node_to")])
#   nets.clean.5$num_components[i.net] <- igraph::count_components(my.graph,mode = "weak")
# }

# summary(nets.clean.5$num_components)

# -------------------------------------------------------------------------

nets.clean.5$connectance <- NA

for(i.net in 1:nrow(nets.clean.5)){
  # select the network edge list
  my.net <- subset(all.links.2, network_id == nets.clean.5$network_id[i.net])
  # is it uni or bipartite?
  part <- nets.clean.5$network_topology_type[i.net]
  
  if(nrow(my.net)>0){
    
    if(part == "bipartite"){
      my.mat <- matrix_from_edge_list(my.net,weighted = F,unipartite = F)
    }else{
      my.mat <- matrix_from_edge_list(my.net,weighted = F,unipartite = T, 
                                      directed = F,diagonal = T)
    }
    nets.clean.5$connectance[i.net] <- connectance.matrix(my.mat)
  }
}

# summary(nets.clean.5$connectance)

# -------------------------------------------------------------------------
nets.clean.5$richness <- NA
nets.clean.5$bipartite.richness.smallest.layer <- NA

for(i.net in 1:nrow(nets.clean.5)){
  
  # select the network edge list
  my.net <- subset(all.links.2, network_id == nets.clean.5$network_id[i.net])
  
  # is it uni or bipartite?
  part <- nets.clean.5$network_topology_type[i.net]
  
  num.nodes.from <- length(unique(my.net$node_from))
  num.nodes.to <- length(unique(my.net$node_to))
  nets.clean.5$richness[i.net] <- length(unique(c(my.net$node_from,my.net$node_to)))
  
  if(part == "bipartite"){
    nets.clean.5$bipartite.richness.smallest.layer[i.net] <- min(num.nodes.from,num.nodes.to)
  }
}

# summary(nets.clean.5$richness)
# summary(nets.clean.5$bipartite.richness.smallest.layer)

# -------------------------------------------------------------------------
# other filters

# some mangal ids are repeated
# stick with a single one
# I don't have any criteria for choosing, so stick with the first one

rep.ids <- duplicated(nets.clean.5$network_id)
nets.clean.6 <- nets.clean.5[!rep.ids,]

# -------------------------------------------------------------------------
# links
links.clean <- unique(all.links.2)
links.clean.2 <- subset(links.clean, network_id %in% nets.clean.6$network_id)

# harmonise link types - but NOTE that interaction type is also recorded
# in the network dataset, and is taken from there in the subsequent analyses
# so, this is simply for completeness
links.clean.2$link_type[links.clean.2$link_type %in% c("bacterivorous","food web",
                                                       "predacious","predation","scavenger")] <- "food web"

links.clean.2$link_type[links.clean.2$link_type %in% c("detritivore","detritivorous")] <- "detritivory"
links.clean.2$link_type[links.clean.2$link_type %in% c("frugivory","frugivorous")] <- "frugivory"
links.clean.2$link_type[links.clean.2$link_type %in% c("herbivorous","herbivory","fungivorous")] <- "herbivory"
links.clean.2$link_type[links.clean.2$link_type %in% c("parasitic","parasitism","parasitoid")] <- "parasitism"

# nodes
nodes.clean <- unique(all.nodes)
nodes.clean.2 <- subset(nodes.clean, node_id %in% unique(c(links.clean.2$node_from,
                                                         links.clean.2$node_to)))
# -------------------------------------------------------------------------
# clean taxonomic names
# NOTE: this takes some time, a few hours at least

# test
# fs::dir_delete(taxadb:::taxadb_dir())

# nodes.clean.2$resolved.sp <- FALSE
# 
# for(i.node in 1:nrow(nodes.clean.2)){
#   my.clean.name <- str_to_sentence(taxadb::clean_names(nodes.clean.2$taxonomy_name[i.node],
#                                        lowercase = F))
#   my.resolved.name <- taxize::gnr_resolve(sci = my.clean.name,
#                                           best_match_only = T,
#                                           canonical = T,
#                                           with_canonical_ranks = F)
#   if(nrow(my.resolved.name) == 1){
#     epithet <- word(my.resolved.name$matched_name2,2)
#     if(!is.na(epithet)){
#       nodes.clean.2$resolved.sp[i.node] <- TRUE
#     }# if proper species name
#   }# if resolved name
# }

# -------------------------------------------------------------------------
# add environmental covariates
net.collection <- nets.clean.6

net.collection$lat.int <- round(net.collection$network_lat,2)
net.collection$lon.int <- round(net.collection$network_lon,2)
net.collection$int.coords <- paste(net.collection$lat.int,"_",net.collection$lon.int,sep="")

# -------------------------------------------------------------------------
# to extract point information from a raster, it can be done from the exact point
# given by the coordinates, or interpolating from the 4 closest values
# if done with the exact point, it returns 163 NAs from the worldclim vars, 
# and 172 from the human footprint index
# if done with the latter, 73 NAs from worldclim, 102 from HF

# method <- "simple"
method <- "bilinear" 

# -------------------------------------------------------------------------
# load and append environmental data

# NOTE: similar to habitat types, this data is not available in the repo
# gather it from worldclim 
# https://www.worldclim.org/data/worldclim21.html

# and store them in a folder of your choice
env.data.path <- "/media/david/Elements SE/Work/datasets/worldclim/"

vars <- c("tavg","tmax","tmin","prec")

paths <- paste(env.data.path,"wc2.1_30s_",vars,"/",sep="")
months <- stri_pad_left(1:12, 2, 0)

env.df <- net.collection[,c("network_id","lon.int","lat.int")]

for(i.path in 1:length(paths)){
  var.df <- env.df
  for(i.month in 1:length(months)){
    my.tif.name <- paste(paths[i.path],"wc2.1_30s_",vars[i.path],"_",months[i.month],".tif",sep="")
    my.tif <- terra::rast(my.tif.name)
    
    my.values <- terra::extract(my.tif,env.df[,c("lon.int","lat.int")],method = method)
    var.df <- cbind(var.df,data.frame(my.values[,2,drop = F]))
  }# i.month
  annual.var.df <- var.df %>% pivot_longer(-(network_id:lat.int),names_to = "month") %>%
    group_by(network_id,lon.int,lat.int) %>%
    summarise(var = mean(value))
  names(annual.var.df)[which(names(annual.var.df) == "var")] <- vars[i.path]
  env.df <- left_join(env.df,annual.var.df)
}# i.path

# bioclimatic variables are in another format
# bio.df <- env.df
bio.vars <- stri_pad_left(1:19, 2, 0)
for(i.var in 1:length(bio.vars)){
  my.tif.name <- paste(env.data.path,"/wc2.1_30s_bio/wc2.1_30s_bio_",bio.vars[i.var],".tif",sep="")
  my.tif <- terra::rast(my.tif.name)
  
  my.values <- terra::extract(my.tif,env.df[,c("lon.int","lat.int")],method = method)
  names(my.values)[2] <- paste("bio_",bio.vars[i.var],sep="")
  env.df <- cbind(env.df,data.frame(my.values[,2,drop = F]))
}

# -------------------------------------------------------------------------
# human footprint rasters

# NOTE: these datasets are not in the repo, can be obtained from 
# https://www.nature.com/articles/sdata201667

# download them to a folder of your choice
hf.data.path <- "/media/david/Elements SE/Work/datasets/HFP_maps/"

hpf93 <- terra::rast(paste0(hf.data.path,"HFP1993_WGS84.tif"))
hpf09 <- terra::rast(paste0(hf.data.path,"HFP2009_WGS84.tif"))

hf93.values <- terra::extract(hpf93,env.df[,c("lon.int","lat.int")],method = method)
names(hf93.values)[2] <- "human_footprint_1993"
hf09.values <- terra::extract(hpf09,env.df[,c("lon.int","lat.int")],method = method)
names(hf09.values)[2] <- "human_footprint_2009"

hf.df <- cbind(env.df,hf93.values[,2,drop = F])
hf.df <- cbind(hf.df,hf09.values[,2,drop = F])
hf.df$lon.int <- hf.df$lat.int <- NULL

# -------------------------------------------------------------------------
# final cleaning

net.collection$lat.int <- net.collection$lon.int <- net.collection$int.coords <- NULL
net.collection.all <- left_join(net.collection,hf.df) %>%
  dplyr::select(network_id,original_id,study_common,network_db,network_year,network_lat, #num_components,
         network_lon,network_interaction,main_animal_taxa,main_animal_taxa_2,richness,
         network_topology_type,habitat_type,tavg:human_footprint_2009)

names(net.collection.all)[which(names(net.collection.all) == "richness")] <- "network_richness"
names(net.collection.all)[which(names(net.collection.all) == "main_animal_taxa")] <- "network_main_animal_taxa"
names(net.collection.all)[which(names(net.collection.all) == "main_animal_taxa_2")] <- "network_main_animal_taxa_2"

# -------------------------------------------------------------------------
# network sampling intensity from Brimacombe et al 2022

sampling_intensity_net <- function(A){
  n.interactions <- sum(A)
  n.sp <- nrow(A)*ncol(A)
  
  si <- sqrt(n.interactions)/sqrt(n.sp)
  return(si)
}

net.collection.all$network_sampling_intensity <- NULL

for(i in 1:nrow(net.collection.all)){
  my.el <- subset(links.clean.2,network_id == net.collection.all$network_id[i])
  my.matrix <- matrix_from_edge_list(my.el)
  
  net.collection.all$network_sampling_intensity[i] <- sampling_intensity_net(my.matrix)
  
}# i.id

# -------------------------------------------------------------------------
# yet another final cleaning...
# there are only 8 bipartite food webs. Discard them
# likewise, there are only 3 unipartite parasitism networks. Discard them

bip.fw <- net.collection.all %>% filter(network_topology_type == "bipartite" &
                                            network_interaction == "food web")
unip.par <- net.collection.all %>% filter(network_topology_type == "unipartite" &
                                            network_interaction == "parasitism")

net.collection.clean <- net.collection.all %>% filter(!(network_topology_type == "bipartite" &
                                                          network_interaction == "food web") & 
                                                        !(network_topology_type == "unipartite" &
                                                            network_interaction == "parasitism"))

links.clean.3 <- subset(links.clean.2, network_id %in% net.collection.clean$network_id)
nodes.clean.3 <- subset(nodes.clean.2, node_id %in% unique(c(links.clean.3$node_from,
                                                             links.clean.3$node_to)))

# -------------------------------------------------------------------------
# I might as well clean up interaction types here - potentially 
# this workflow accepts any interaction specified in "network_interaction",
# but as of 12/2023 I "only" analyse food web, frug,herbi,parasitism,pollination

net.collection.clean.2 <- net.collection.clean %>% filter(network_interaction %in% 
                                                            c("food web","parasitism",
                                                              "pollination","herbivory","frugivory"))
links.clean.4 <- subset(links.clean.3, network_id %in% net.collection.clean.2$network_id)
nodes.clean.4 <- subset(nodes.clean.3, node_id %in% unique(c(links.clean.4$node_from,
                                                             links.clean.4$node_to)))
# -------------------------------------------------------------------------
# writing - network metadata, links, nodes, and adjacency matrices
# NOTE
# the null model (configuration model, see scripts 03_*) does not allow fully connected rows
# or columns
# there are some networks that have these fully connected nodes. 
# Here is where I decide what to do with them. options are:
# 1) discard these networks
# 2) remove those nodes (no!)
# 3) in these nodes, remove one single link at random

# option 1 is the most conservative, but depending on how many networks i discard,
# i might lose relevant data.
# option 3 is also not ideal, specially if there are many fully connected nodes in a net

# after checking, most nets discarded are food webs, so go with it, as I have enough food webs

# generate adjacency matrices
ids.to.discard <- NULL
for(i.id in 1:nrow(net.collection.clean.2)){
  partit <- net.collection.clean.2$network_topology_type[i.id]
  my.el <- subset(links.clean.4,network_id == net.collection.clean.2$network_id[i.id])
  if(partit == "bipartite"){
    my.adjmat <- matrix_from_edge_list(my.el,unipartite = F)
    my.file.name <- paste("data/bipartite_adjacency_matrices/",net.collection.clean.2$network_id[i.id],".csv",sep="")
  }else{
    my.adjmat <- matrix_from_edge_list(edge.list = my.el,unipartite = T,directed = F,diagonal = F,weighted = F)
    my.file.name <- paste("data/unipartite_adjacency_matrices/",net.collection.clean.2$network_id[i.id],".csv",sep="")
  }

  fully.connected.rows <- which(rowSums(my.adjmat) == ncol(my.adjmat))
  fully.connected.cols <- which(colSums(my.adjmat) == nrow(my.adjmat))
  
  if(length(fully.connected.rows) > 0){
    cat(net.collection.clean.2$network_id[i.id],"- fully connected rows:",fully.connected.rows,"\n")
  }

  if(length(fully.connected.cols) > 0){
    cat(net.collection.clean.2$network_id[i.id],"- fully connected columns:",fully.connected.cols,"\n")
  }
  
  if(length(fully.connected.cols) == 0 & length(fully.connected.rows) == 0){
    write.csv(my.adjmat,file = my.file.name)
  }else{
    ids.to.discard[length(ids.to.discard)+1] <- net.collection.clean.2$network_id[i.id]
  }
  
}

discarded.nets <- net.collection.clean.2[net.collection.clean.2$network_id %in% ids.to.discard,]

net.collection.clean.3 <- net.collection.clean.2[!(net.collection.clean.2$network_id %in% ids.to.discard),]
links.clean.5 <- subset(links.clean.4, network_id %in% net.collection.clean.3$network_id)
nodes.clean.5 <- subset(nodes.clean.4, node_id %in% unique(c(links.clean.5$node_from,
                                                             links.clean.5$node_to)))

# final_definitive_pleasenomoreofthis check: discard networks from same study, same type, same richness
# and coordinates similar up to one decimal point (e.g. 70.1X)
# as these are most likely duplicates as well.

net.collection.clean.3$round.lat <- round(net.collection.clean.3$network_lat,1)
net.collection.clean.3$round.lon <- round(net.collection.clean.3$network_lon,1)
net.collection.clean.3 <- net.collection.clean.3 %>%
  group_by(study_common,network_interaction,network_richness,round.lat,round.lon) %>%
  mutate(is_duplicate = row_number() > 1) %>%
  ungroup()

net.collection.clean.4 <- net.collection.clean.3 %>% 
  distinct(study_common,network_interaction,network_richness,round.lat,round.lon,.keep_all = TRUE)
net.collection.clean.4$round.lon <- NULL
net.collection.clean.4$round.lat <- NULL
net.collection.clean.4$is_duplicate <- NULL

links.clean.6 <- subset(links.clean.5, network_id %in% net.collection.clean.4$network_id)
nodes.clean.6 <- subset(nodes.clean.5, node_id %in% unique(c(links.clean.6$node_from,
                                                             links.clean.6$node_to)))

# -------------------------------------------------------------------------

write.csv2(links.clean.6,"results/network_links_collection.csv",row.names = F)
write.csv2(nodes.clean.6,"results/network_nodes_collection.csv",row.names = F)
write.csv2(net.collection.clean.4,"results/network_collection.csv",row.names = F)

# -------------------------------------------------------------------------
# discarded taxonomy cleaning

# 
# node.names <- sort(unique(nodes.clean$taxonomy_name))
# # test.names <- node.names[1:100]
# 
# names.clean <- bdc_clean_names(sci_names = node.names, save_outputs = FALSE)
# names.clean <- names.clean[,c("scientificName","names_clean")]
# names(names.clean)[2] <- "bdc"
# names.clean$taxadb <- taxadb::clean_names(names.clean$scientificName,lowercase = F)
# names.clean$bdc <- str_to_sentence(names.clean$bdc)
# names.clean$taxadb <- str_to_sentence(names.clean$taxadb)
# 
# tt <- unique(names.clean$taxadb)
# 
# # for some reason gnr_resolve crashes with the whole dataset
# resolved.list <- list()
# names.seq <- seq(1,length(tt),by = 100)
# for(i.seq in 1:(length(names.seq)-1)){
#   
#   if(i.seq != (length(names.seq)-1)){
#     my.names <- tt[names.seq[i.seq]:names.seq[i.seq+1]]
#   }else{
#     my.names <- tt[names.seq[i.seq]:length(tt)]
#   }
#   
#   resolved.list[[length(resolved.list)+1]] <- taxize::gnr_resolve(sci = my.names,
#                                                                   best_match_only = T,
#                                                                   canonical = T,
#                                                                   with_canonical_ranks = F)
#   
# }
# 
# resolved.names <- bind_rows(resolved.list)
# resolved.names$epithet <- word(resolved.names$matched_name,2)
# 
# 
# # resolved.names <- taxize::gnr_resolve(sci = tt[1:5000],#unique(names.clean$taxadb),
# #                                       best_match_only = T)
# # resolved.names$final.name <- str_to_sentence(taxadb::clean_names(resolved.names$matched_name))
# 
# # test.names.clean.df <- taxize::gnr_resolve(sci = test.names, best_match_only = TRUE)
# 
# library(bdc)
# 
# bdc_query_names_taxadb(
#   sci_name            = test.names,
#   replace_synonyms    = TRUE, # replace synonyms by accepted names?
#   suggest_names       = TRUE, # try to found a candidate name for misspelled names?
#   suggestion_distance = 0.9, # distance between the searched and suggested names
#   db                  = "all", # taxonomic database
#   # rank_name           = "Plantae", # a taxonomic rank
#   # rank                = "kingdom", # name of the taxonomic rank
#   parallel            = FALSE, # should parallel processing be used?
#   ncores              = 2, # number of cores to be used in the parallelization process
#   export_accepted     = FALSE # save names linked to multiple accepted names
# )
# 
# # -------------------------------------------------------------------------
# 
# 
# td_create(provider = "itis",overwrite = T)
# 
# node.names <- sort(unique(nodes.clean$taxonomy_name))
# node.names.clean <- sort(unique(clean_names(node.names,lowercase = F)))
# 
# # nodes.clean.2 <- nodes.clean
# # nodes.clean.2$taxonomy_name_clean <- clean_names(nodes.clean.2$taxonomy_name,lowercase = F)
# # 
# # taxadb::filter_name("lavandula",provider = "col")
# # 
# # taxadb::get_ids("lavandula stoechas",provider = "col")
# # 
# # taxize::get_ids("Carabus ovalis",db = "eol")
# # 
# # nodes.clean.2$db_id <- get_ids(nodes.clean.2$taxonomy_name_clean,db = "itis")
# 
# node.names.df <- data.frame(orig.name = node.names.clean) 
# node.names.df <- node.names.df %>%
#   mutate(id_col = get_ids(orig.name, db = "col"),
#          id_itis = get_ids(orig.name, db = "itis"))
# 
# sum(is.na(node.names.df$id_col))
# sum(is.na(node.names.df$id_itis))
# 
# # 
# # nodes.clean.2 <- nodes.clean.2 %>%
# #   mutate(accepted_name = get_names(id, "col"))
#   


