
# curate data from Dalsgaard et al. 2021
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2435.13784

library(tidyverse)

# -------------------------------------------------------------------------
# read metadata
net.metadata <- readxl::read_xlsx("data/raw_interaction_datasets/Dalsgaard_2021/000_Networks_metadata.xlsx")

# read all networks
net.path <- "data/raw_interaction_datasets/Dalsgaard_2021/"
net.list <- list.files(path = net.path,pattern = "*.csv",full.names = T) %>%
  map(read.csv2)
net.names <- list.files(path = net.path,pattern = "*.csv")
net.names <- substr(net.names,1,nchar(net.names)-4)

# -------------------------------------------------------------------------

net.ids <- paste("D21_",1:length(net.names),sep="")#na.omit(unique(link.collection$network_id))
net.metadata <- arrange(net.metadata,`Network ID`)

# networks
d21.nets <- data.frame(network_id = net.ids, 
                     original_id = net.names,
                     network_db = "D21",
                     network_interaction = "pollination",
                     main_animal_taxa = "Hummingbirds",
                     network_year = NA,
                     network_lat = net.metadata$Latitude,
                     network_lon = net.metadata$Longitude,
                     habitat_type = NA,
                     network_topology_type = "bipartite",
                     network_spatial_type = "Point")

link.list <- list()
for(i.net in 1:length(net.list)){
  tt <- net.list[[i.net]]
  tl <- pivot_longer(tt,cols = c(-X),names_to = "node_from",values_to = "freq") %>%
    filter(freq > 0)
  names(tl)[1] <- "node_to" 
  tl$network_id <- d21.nets$network_id[i.net]
  tl$link_type <- "pollination"
  tl$original_id <- d21.nets$original_id[i.net]
  tl <- tl[,c("network_id","node_from","node_to","link_type","original_id")]
  link.list[[length(link.list)+1]] <- tl
}

d21.links <- bind_rows(link.list)
d21.links$node_from <- str_replace(d21.links$node_from,"\\.","_")
d21.links$node_to <- str_replace(d21.links$node_to," ","_")

all.sp <- sort(unique(c(d21.links$node_from,d21.links$node_to)))
node.ids <- paste("D21_node_",1:length(all.sp),sep="")
d21.nodes <- data.frame(node_id = node.ids,taxonomy_name = all.sp,taxonomy_rank = "species")

# rename links
d21.links$node_from.new <- d21.nodes$node_id[match(d21.links$node_from,d21.nodes$taxonomy_name)]
d21.links$node_to.new <- d21.nodes$node_id[match(d21.links$node_to,d21.nodes$taxonomy_name)]

d21.links$node_from <- d21.links$node_from.new
d21.links$node_to <- d21.links$node_to.new

d21.links$node_from.new <- d21.links$node_to.new <- NULL

# -------------------------------------------------------------------------
write.csv2(d21.nets,"data/D21_networks.csv",row.names = F)
write.csv2(d21.links,"data/D21_links.csv",row.names = F)
write.csv2(d21.nodes,"data/D21_nodes.csv",row.names = F)


