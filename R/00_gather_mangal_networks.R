
# curate datasets from mangal
# https://mangal.io/


# -------------------------------------------------------------------------

# mangal to link list
library(tidyverse)
library(tidygraph)

source("R/aux_metric_functions.R")
# library(rmangal)

# a few general notes:
# Divide antagonistic networks in micro (insect) vs macro (mammals, etc)
# I guess there will be differences between them

# download network objects from mangal and store them
# mgs <- search_datasets("", verbose = FALSE)
# mgn <- get_collection(mgs)
# save(mgn,file = "data/mangal_networks_08_21.Rdata")

# once they are downloaded, I can redo the cleaning with the local copy
load("data/raw_interaction_datasets/mangal/mangal_networks_02_22.Rdata")

num.net <- length(mgn)

# end up with the following dataframes:
# nodes: node_id - taxonomy.name - taxonomy.rank
# links: network_id - node_from - node_to - link_type - link-direction - link_value - link_units
# networks: network_id - database - lat - lon - spatial_type 

net.list <- list()
node.list <- list()
link.list <- list()

node.column.names <- c("node_id","taxonomy.name","taxonomy.rank")
link.column.names <- c("node_from","node_to","type","direction","value", "attribute.unit")

# check the interaction types that appear in each network
# or the network descriptions
# sink('data/mangal-descriptions.txt')
# for(i.net in 1:num.net){
# 
#   my.id <- paste("mangal_",mgn[[i.net]]$network$network_id,sep="")
# 
#   links.orig <- mgn[[i.net]]$interactions
# 
#   net.description <- mgn[[i.net]]$network$description
#   cat(my.id,"-",net.description,"\n")
#   
#   # net.type <-  names(table(links.orig$type))
#   # cat(my.id,"-",net.type,"\n")
# }
# sink()

# -------------------------------------------------------------------------
# I want to differentiate frugivory from pollination networks
# and mangal only has "mutualism" as link info
# so, when I find a "mutualism" network, I will check the network description
# and look for keywords - including some misspellings that I identified
poll.keywords <- c("insect","insects",
                   "Insect","Insects",
                   "pollinator","pollinators","pollination", "pollinating","polinator",
                   "Pollinator","Pollinators","Pollnator","Pollination", "Pollinating",
                   "flower","flowers",
                   "Flower","Flowers")

frug.keywords <- c("bird","birds",
                   "Bird","Birds",
                   "fruit","fruits",
                   "Fruit","Fruits",
                   "frugivory","frugivore",
                   "Frugivory","Frugivore",
                   "seed","seeds",
                   "Seed","Seeds")

# herbivory is also a category
herb.interactions <- c("herbivory")

# furthermore, all other consumptive interactions are part of the network category "food web"
food.web.interactions <- c("detritivore","scavenger","predation")

# I want to identify the main animal guild sampled, if any,
# at least for terrestrial systems. For aquatic ones, it's... too much
# so I have skimmed through the descriptions to select a few keywords to choose from

invertebrate.guild.keywords <- c("Insect","insect","Insects","insects","Pollnator",
                                 "sawflies","pollinator","bee","Ahpid","logs",
                                 "Pollination","plant-flower visitor","arthropod",
                                 "polinator","Flea","wildflower","gall","galls",
                                 "carrion","soybean","grassland")

stream.keywords <- c("Ria","fishery","Current","stream","River","spring",
                      "stream","springbrook","streams","rapids",
                      "freshwater")

lake.keywords <- c("Lake","lake","Loch","lakes","lagoon","Laguna")

marine.keywords <- c("anemons","bay","estuary","anemone","fishery","Bay",
                     "plankton","Sea","intertidal","shore",
                     "atoll","ice","Ice","pelagic","sand-bottom","Strait",
                     "Gulf","Benguela","Kelp","kelp","benthic","seagrass",
                     "reefs","antarctic","backwater","marine","seas","mangrove","Mussel",
                     "Mudflat","Aleutian","Shallow","marsh","mudflat","Shelf","shelf",
                     "oligohaline","Antarctic","reservoir","fish","fishes",
                     "Creek","offshore")

fish.keywords <- c("fish","Fish","fishes","Fishes","fishery")

hummingbird.keywords <- c("Hummingbirds","bromeliad")

bird.keywords <- c("bird","Bird","birds","Birds","owls","frugivore")

mammal.keywords <- c("Himalayas","primate")

# -------------------------------------------------------------------------

# test
# descriptions <- character()
# id <- character()
# for(i.net in 1:num.net){
#   id[length(id)+1] <- mgn[[i.net]]$network$network_id
#   descriptions[length(descriptions)+1] <- mgn[[i.net]]$network$description
# }
# dd <- data.frame(network_id = id, des = descriptions)

for(i.net in 1:num.net){
  
  my.id <- paste("mangal_",mgn[[i.net]]$network$network_id,sep="")
  
  nodes.orig <- mgn[[i.net]]$nodes
  links.orig <- mgn[[i.net]]$interactions
  
  orig.type <-  names(table(links.orig$type))[1]
  
  if(is.null(orig.type)){
    net.type <- "unspecified"
  }else if(orig.type == "mutualism"){
    net.description <- mgn[[i.net]]$network$description
    if(grepl(paste(poll.keywords,collapse="|"),net.description)){
      net.type <- "pollination"
    }else if(grepl(paste(frug.keywords,collapse="|"),net.description)){
      net.type <- "frugivory"
    }else{
      net.type <- "unspecified mutualism"
    }
  }else if(orig.type %in% herb.interactions){
    net.type <- "herbivory"
  }else if(orig.type %in% food.web.interactions){
    net.type <- "food web"
  }else{
    net.type <- orig.type
  }
  
  if(all(node.column.names %in% names(nodes.orig)) &
     all(link.column.names %in% names(links.orig))){
    
    nodes <- mgn[[i.net]]$nodes[,node.column.names]
    names(nodes) <- c("node_id","taxonomy_name","taxonomy_rank")
    nodes$node_id <- paste("mangal_node_",nodes$node_id,sep="")
    
    links <- mgn[[i.net]]$interactions[,link.column.names]
    links$network_id <- my.id
    names(links) <- c("node_from","node_to","link_type","link_direction","link_value","link_units","network_id")
    links <- links[,c("network_id","node_from","node_to","link_type","link_direction","link_value","link_units")]
    
    bip <- bipartiteness.edge.list(links)
    
    desc <- mgn[[i.net]]$network$description
    
    my.animal.guild <- NA
    if(grepl(paste0('\\b', invertebrate.guild.keywords, '\\b', collapse = '|'),desc)){
      my.animal.guild <- "Insects"
    }else if(grepl(paste0('\\b', hummingbird.keywords, '\\b', collapse = '|'),desc)){
      my.animal.guild <- "Hummingbirds"
    }else if(grepl(paste0('\\b', bird.keywords, '\\b', collapse = '|'),desc)){
      my.animal.guild <- "Birds"
    }else if(grepl(paste0('\\b', mammal.keywords, '\\b', collapse = '|'),desc)){
      my.animal.guild <- "Mammals"
    }else if(grepl(paste0('\\b', fish.keywords, '\\b', collapse = '|'),desc)){
      my.animal.guild <- "Fish"
    }else{
      my.animal.guild <- "Unspecified"
    }
    
    # aquatic habitats
    my.habitat <- NA
    if(grepl(paste0('\\b', stream.keywords, '\\b', collapse = '|'),desc)){
      my.habitat <- "Streams"
    }else if(grepl(paste0('\\b', marine.keywords, '\\b', collapse = '|'),desc)){
      my.habitat <- "Marine"
    }else if(grepl(paste0('\\b', lake.keywords, '\\b', collapse = '|'),desc)){
      my.habitat <- "Lakes"
    }
    
    network <- data.frame(network_id = my.id,
                          original_id = as.character(mgn[[i.net]]$network$network_id),
                          network_db = "mangal",
                          network_interaction = net.type,
                          main_animal_taxa = my.animal.guild,
                          network_year = as.numeric(substr(mgn[[i.net]]$network$date,1,4)),
                          network_lat = mgn[[i.net]]$network$geom_lat[[1]],
                          network_lon = mgn[[i.net]]$network$geom_lon[[1]],
                          habitat_type = my.habitat,
                          network_topology_type = ifelse(bip,"bipartite","unipartite"),
                          network_spatial_type = mgn[[i.net]]$network$geom_type)
    
    net.list[[length(net.list)+1]] <- network
    node.list[[length(node.list)+1]] <- nodes
    link.list[[length(link.list)+1]] <- links
    
  }# if valid nodes and links
}

mg.networks <- bind_rows(net.list)
mg.nodes <- bind_rows(node.list)
mg.links <- bind_rows(link.list)

mg.nodes.2 <- data.frame(node_id = NA,taxonomy_name = unique(mg.nodes$taxonomy_name))
mg.nodes.2$node_id <- paste("mangal_node_",formatC(1:nrow(mg.nodes.2), width = 4, format = "d", flag = "0"),sep="")
mg.nodes.2$taxonomy_rank <- mg.nodes$taxonomy_rank[match(mg.nodes.2$taxonomy_name,mg.nodes$taxonomy_name)]

# links
mg.nodes.temp <- mg.nodes.2
names(mg.nodes.temp)[1] <- "node_id_new"
mg.nodes.temp <- left_join(mg.nodes,mg.nodes.temp)

mg.links.2 <- mg.links
mg.links.2$node_from <- paste("mangal_node_",mg.links.2$node_from,sep="")
mg.links.2$node_to <- paste("mangal_node_",mg.links.2$node_to,sep="")

mg.links.2$node_from_new <- mg.nodes.temp$node_id_new[match(mg.links.2$node_from,mg.nodes.temp$node_id)]
mg.links.2$node_to_new <- mg.nodes.temp$node_id_new[match(mg.links.2$node_to,mg.nodes.temp$node_id)]

mg.links.3 <- mg.links.2[,c(1,8,9,4,5,6,7)]
names(mg.links.3)[c(2,3)] <- c("node_from","node_to")

# -------------------------------------------------------------------------

# save to disk
write.csv2(mg.networks,file = "data/mangal_networks.csv",row.names = FALSE)
write.csv2(mg.nodes.2,file = "data/mangal_nodes.csv",row.names = FALSE)
write.csv2(mg.links.3,file = "data/mangal_links.csv",row.names = FALSE)

