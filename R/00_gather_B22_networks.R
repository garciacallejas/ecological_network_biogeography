
# curate datasets from Brimacombe et al. 2022
# https://onlinelibrary.wiley.com/doi/10.1111/geb.13593

# Here I follow a different approach than with other datasets
# I take all previous ones, and filter in advance for duplicated records
# since in a visual inspection a few networks seem to come from the same databases

# so, if running everything from scratch, THIS MUST BE THE LAST 00 script to run

library(tidyverse)
library(geosphere)

# -------------------------------------------------------------------------

# Brimacombe et al. 2022 metadata
df1 <- read.csv("data/raw_interaction_datasets/Brimacombe_2022/information_binary_networks.csv")

# collection of networks from the previous datasets
# to avoid processing everything and then filtering

fs.nets <- read.csv2("data/fs_networks.csv")
gw.nets <- read.csv2("data/gateway_networks.csv")
mg.nets <- read.csv2("data/mangal_networks.csv")
wol.nets <- read.csv2("data/wol_networks.csv")
p22.nets <- read.csv2("data/p22_networks.csv")
ti22.nets <- read.csv2("data/ti22_networks.csv")
m22.nets <- read.csv2("data/M22_networks.csv")
ai22.nets <- read.csv2("data/AI22_networks.csv")
d21.nets <- read.csv2("data/D21_networks.csv")

fs.nets$original_id <- as.character(fs.nets$original_id)
gw.nets$original_id <- as.character(gw.nets$original_id)
mg.nets$original_id <- as.character(mg.nets$original_id)
wol.nets$original_id <- as.character(wol.nets$original_id)
p22.nets$original_id <- as.character(p22.nets$original_id)
ti22.nets$original_id <- as.character(ti22.nets$original_id)
m22.nets$original_id <- as.character(m22.nets$original_id)
ai22.nets$original_id <- as.character(ai22.nets$original_id)
d21.nets$original_id <- as.character(d21.nets$original_id)

ai22.nets$network_id <- as.character(ai22.nets$network_id)

all.nets <- bind_rows(list(fs.nets,gw.nets,mg.nets,
                           wol.nets,p22.nets,ti22.nets,
                           m22.nets,ai22.nets,d21.nets))

# -------------------------------------------------------------------------
# Function to check if a row from df1 is included in df2
# is_included <- function(row, df) {
#   distances <- distVincentySphere(df[, c("Longitude", "Latitude")], row[, c("Longitude", "Latitude")])
#   min_distance <- min(distances)
#   return(min_distance < 100)  # Assuming 0 distance means exact match
# }
df1$min.dist <- NA

for(i1 in 1:nrow(df1)){
  my.distances <- distVincentyEllipsoid(c(df1$long[i1],df1$lat[i1]),as.matrix(all.nets[,c("network_lon","network_lat")]))
  df1$min.dist[i1] <- min(my.distances,na.rm = T)
}

# set an "overlap" distance wide enough - be conservative
df1.not.included <- subset(df1,min.dist > 10000)

# keep the first entries of similar coordinates
df1.unique <- df1.not.included[!duplicated(df1.not.included[,c("lat","long")]),]

# -------------------------------------------------------------------------
# now, go through each network

# networks
net.ids <- paste("B22_",1:nrow(df1.unique),sep="")
b22.nets <- data.frame(network_id = net.ids, 
                       original_id = df1.unique$Name,
                       network_db = "B22",
                       network_interaction = NA,
                       main_animal_taxa = NA,
                       network_year = NA,
                       network_lat = df1.unique$lat,
                       network_lon = df1.unique$long,
                       habitat_type = NA,
                       network_topology_type = "bipartite",
                       network_spatial_type = "Point")

link.list <- list()
for(i.net in 1:nrow(b22.nets)){
  
  # name is in lowercase...
  if(grepl("Belay",df1.unique$Name[i.net])){
    df1.unique$Name[i.net] <- paste("A_HP_belay",substr(df1.unique$Name[i.net],11,nchar(df1.unique$Name[i.net])),sep="")
  }else if(grepl("Lewis",df1.unique$Name[i.net])){
    df1.unique$Name[i.net] <- paste("A_PH_lewis",substr(df1.unique$Name[i.net],11,nchar(df1.unique$Name[i.net])),sep="")
  }else if(grepl("Novotny",df1.unique$Name[i.net])){
    df1.unique$Name[i.net] <- paste("A_PH_novotny",substr(df1.unique$Name[i.net],13,nchar(df1.unique$Name[i.net])),sep="")
  }else if(grepl("Tavakilian",df1.unique$Name[i.net])){
    df1.unique$Name[i.net] <- paste("A_PH_tavakilian",substr(df1.unique$Name[i.net],16,nchar(df1.unique$Name[i.net])),sep="")
  }else if(grepl("Ueckert",df1.unique$Name[i.net])){
    df1.unique$Name[i.net] <- paste("A_PH_ueckert",substr(df1.unique$Name[i.net],13,nchar(df1.unique$Name[i.net])),sep="")
  }
  
  my.net.name <- as.character(paste("data/Brimacombe_2022/binary_networks/",df1.unique$Name[i.net],".csv",sep=""))
  my.matrix <- read.csv(my.net.name)
  names(my.matrix) <- paste("B22net",i.net,"_col_",1:ncol(my.matrix),sep="")
  my.matrix$row.sp <- paste("B22net",i.net,"_row_",1:nrow(my.matrix),sep="")
  
  my.name.parts <- str_split(df1.unique$Name[i.net],pattern = "_")
  my.type.orig <- my.name.parts[[1]][2]
  
  if(my.type.orig == "HP"){
    my.type <- "parasitism"
  }else if(my.type.orig == "PH"){
    my.type <- "herbivory"
  }else if(my.type.orig == "PL"){
    my.type <- "pollination"
  }else if(my.type.orig == "SD"){
    my.type <- "frugivory"
  }
  
  b22.nets$network_interaction[i.net] <- my.type
  
  # links
  tl <- pivot_longer(my.matrix,cols = -row.sp,names_to = "node_from",values_to = "freq") %>%
    filter(freq > 0)
  names(tl)[1] <- "node_to" 
  tl$network_id <- b22.nets$network_id[i.net]
  tl$link_type <- "my.type"
  tl$original_id <- b22.nets$original_id[i.net]
  tl <- tl[,c("network_id","node_from","node_to","link_type","original_id")]
  link.list[[length(link.list)+1]] <- tl
  
}

b22.links <- bind_rows(link.list)
b22.links$node_from <- str_replace(b22.links$node_from,"\\.","_")
b22.links$node_to <- str_replace(b22.links$node_to," ","_")

all.sp <- sort(unique(c(b22.links$node_from,b22.links$node_to)))
node.ids <- paste("B22_node_",1:length(all.sp),sep="")
b22.nodes <- data.frame(node_id = node.ids,taxonomy_name = all.sp,taxonomy_rank = "species")

# -------------------------------------------------------------------------
write.csv2(b22.nets,"data/B22_networks.csv",row.names = F)
write.csv2(b22.links,"data/B22_links.csv",row.names = F)
write.csv2(b22.nodes,"data/B22_nodes.csv",row.names = F)




