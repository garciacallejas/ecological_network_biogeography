
# obtain basic descriptors of the network dataset

library(tidyverse)

# -------------------------------------------------------------------------
# raw datasets

fs.links <- read.csv2("data/fs_links.csv")
gateway.links <- read.csv2("data/gateway_links.csv")
mangal.links <- read.csv2("data/mangal_links.csv")
wol.links <- read.csv2("data/wol_links.csv")
p22.links <- read.csv2("data/p22_links.csv")
ti22.links <- read.csv2("data/ti22_links.csv")
m22.links <- read.csv2("data/M22_links.csv")
ai22.links <- read.csv2("data/AI22_links.csv")
d21.links <- read.csv2("data/D21_links.csv")

ai22.links$network_id <- as.character(ai22.links$network_id)

raw.links <- bind_rows(fs.links[,1:4],
                       gateway.links[,1:4],
                       mangal.links[,1:4],
                       wol.links[,1:4],
                       p22.links[1:4],
                       ti22.links[1:4],
                       m22.links[1:4],
                       ai22.links[1:4],
                       d21.links[1:4])

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

raw.nodes <- bind_rows(fs.nodes,gateway.nodes,mangal.nodes,wol.nodes,p22.nodes,
                       ti22.nodes,m22.nodes,ai22.nodes,d21.nodes)

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

raw.nets <- bind_rows(list(fs.nets,gw.nets,mg.nets,
                           wol.nets,p22.nets,ti22.nets,
                           m22.nets,ai22.nets,d21.nets))

# -------------------------------------------------------------------------
# filtered datasets
net.collection <- read_csv2("results/network_collection.csv")
all.links <- read_csv2("results/network_links_collection.csv")
all.nodes <- read_csv2("results/network_nodes_collection.csv")

# -------------------------------------------------------------------------
# some numbers

valid.network.interactions <- c("food web","frugivory","herbivory","parasitism","pollination")

raw.nets.2 <- subset(raw.nets,network_interaction %in% valid.network.interactions)
valid.nets.2 <- subset(net.collection, network_interaction %in% valid.network.interactions)

# how many networks?
length(unique(raw.nets.2$network_id))
length(unique(valid.nets.2$network_id))

# of each type?
table(raw.nets.2$network_interaction)
table(valid.nets.2$network_interaction)

# -------------------------------------------------------------------------

study.count <- valid.nets.2 %>% group_by(study_common) %>%
  summarise(count = n())

sum(study.count$count == 1)
sum(study.count$count > 1)

