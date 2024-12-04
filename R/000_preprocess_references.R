
# this script preprocesses and harmonises references across datasets
# so that I can identify which networks are potentially from the same study

# I go through each dataset, obtain the study associated to each network id
# and store them

# In a second step, manually, I check the syntax of every reference
# and combine those that are likely to be coming from the same study
# e.g. in some cases one can find arroyo82, Arroyo_1982, etc.

# lastly, from line 310 onwards, I read this curated list again and link it
# to the network IDs

library(tidyverse)

# -------------------------------------------------------------------------
db_refs_list <- list()
# -------------------------------------------------------------------------
# gateway

gateway_data <- read.csv("data_clean/raw_interaction_datasets/gateway/283_3_283_2_FoodWebDataBase_2018_12_10.csv")

# 
gateway_refs <- unique(data.frame(network_db = "gateway", 
                           network_id = gateway_data$foodweb.name,
                           reference = NA,
                           citation = gateway_data$link.citation))

db_refs_list[[length(db_refs_list)+1]] <- gateway_refs
write.csv2(gateway_refs,file = "data/dataset_references/gateway_references.csv")
# -------------------------------------------------------------------------
# mangal

load("data_clean/raw_interaction_datasets/mangal/mangal_networks_02_22.Rdata")
num.net <- length(mgn)

ref.list <- list()
for(i.net in 1:num.net){
  my.id <- mgn[[i.net]]$network$network_id
  my.ref <- paste(mgn[[i.net]]$reference$first_author,mgn[[i.net]]$reference$year,sep="")
  my.doi <- mgn[[i.net]]$reference$doi
  
  ref.list[[length(ref.list)+1]] <- data.frame(network_db = "mangal", 
                                               network_id = as.character(my.id),
                                               reference = NA,
                                               citation = my.ref,
                                               doi = my.doi)
  
}
mangal_refs <- unique(bind_rows(ref.list))

db_refs_list[[length(db_refs_list)+1]] <- mangal_refs
write.csv2(mangal_refs,file = "data/dataset_references/mangal_references.csv")

# -------------------------------------------------------------------------
# wol

base_url <- "https://www.web-of-life.es/"    
wol.nets <- read.csv(paste0(base_url,"get_network_info.php"))

wol_refs <- unique(data.frame(network_db = "WOL",
                       network_id = wol.nets$network_name,
                       reference = wol.nets$reference,
                       citation = wol.nets$author))

db_refs_list[[length(db_refs_list)+1]] <- wol_refs
write.csv2(wol_refs,file = "data/dataset_references/wol_references.csv")

# -------------------------------------------------------------------------
# p22

p22.net.metadata <- read.csv("data_clean/raw_interaction_datasets/Parra_2022/metadata_full_pollination_dataset.csv")
p22.nets <- read.delim("data_clean/raw_interaction_datasets/Parra_2022/networks_full_pollination_dataset.txt")

id_network <- p22.nets$id_network[match(p22.net.metadata$webcode,p22.nets$webcode)]

# net.ids <- unique(p22.net.metanowol$id_network)

p22_refs <- unique(data.frame(network_db = "P22",
                       network_id = as.character(id_network),#p22.net.metadata$webcode,
                       reference = NA,
                       citation = p22.net.metadata$Name_paper))

db_refs_list[[length(db_refs_list)+1]] <- p22_refs
write.csv2(p22_refs,file = "data/dataset_references/P22_references.csv")

# -------------------------------------------------------------------------
# TI22

net.coords <- read.delim("data_clean/raw_interaction_datasets/Thebault_unpublished/metadata_herbivory.txt",header = T,sep = "\t")

ti22_refs <- unique(data.frame(network_db = "TI22",
                        network_id = net.coords$web,
                        reference = NA,
                        citation = net.coords$web))

db_refs_list[[length(db_refs_list)+1]] <- ti22_refs
write.csv2(ti22_refs,file = "data/dataset_references/TI22_references.csv")

# -------------------------------------------------------------------------
# M22

load("data_clean/raw_interaction_datasets/Martins_2022/metadata_networks.Rdata")

# Extract author and year
m22.matches <- str_match(metadata_networks$reference, "^([A-Za-zÀ-ÖØ-öø-ÿ'’]+).*\\((\\d{4})\\)")

# Combine to form "Author_year"
m22.author_year <- paste0(m22.matches[,2], "_", m22.matches[,3])

M22_refs <- unique(data.frame(network_db = "M22",
                       network_id = metadata_networks$net_id,
                       reference = metadata_networks$reference,
                       citation = m22.author_year))

M22_refs$citation[which(M22_refs$citation == "NA_NA")] <- paste0("M22_",M22_refs$network_id[which(M22_refs$citation == "NA_NA")])

db_refs_list[[length(db_refs_list)+1]] <- M22_refs
write.csv2(M22_refs,file = "data/dataset_references/M22_references.csv")

# -------------------------------------------------------------------------
# AI22 - this is special in that I filter networks by location,
# so there is no network_id in the original data. I need to redo the
# whole process and keep the references for each id I generate

ai.full <- read.csv("data_clean/raw_interaction_datasets/AI22/AtlanticForestInvertFloInteractionData_2022-07.csv")

# code from gather_AI22_networks modified to keep references
ai.full.clean <- ai.full %>% dplyr::select(latitude_y,longitude_x,site_name_id,
                                           plant_species_complete_name,
                                           invertebrate_species_complete_name,
                                           campain_year_finish,
                                           invertebrate_behavior,invertebrate_order,
                                           invertebrate_family,
                                           invertebrate_species_complete_name,
                                           reference_citation) %>%
  filter(!is.na(plant_species_complete_name) & !is.na(invertebrate_species_complete_name))

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

# filter by year
# when more than one network in a given location, select the latest
ai.full.clean.3 <- unique(ai.full.clean.2[,c("network_id","lat","lon","network_year","network_interaction","reference_citation")])

ai.full.last.year <- ai.full.clean.3 %>% 
  group_by(lat,lon,network_interaction) %>%
  mutate(max.year = max(network_year))
ai.full.max.year <- subset(ai.full.last.year, network_year == max.year)

ai22.matches <- str_match(ai.full.max.year$reference_citation, "^([A-Za-zÀ-ÖØ-öø-ÿ'’]+).*?(\\d{4})")
ai22.author_year <- paste0(ai22.matches[,2], "_", ai22.matches[,3])

AI22_refs <- unique(data.frame(network_db = "AI22",
                      network_id = ai.full.max.year$network_id,
                      reference = ai.full.max.year$reference_citation,
                      citation = ai22.author_year))

# update citation when NA_NA
# assume that each of these is an independent study
AI22_refs$citation[which(AI22_refs$citation == "NA_NA")] <- AI22_refs$network_id[which(AI22_refs$citation == "NA_NA")]

db_refs_list[[length(db_refs_list)+1]] <- AI22_refs
write.csv2(AI22_refs,file = "data/dataset_references/AI22_references.csv")

# -------------------------------------------------------------------------
# d21

net.metadata <- readxl::read_xlsx("data_clean/raw_interaction_datasets/Dalsgaard_2021/000_Networks_metadata.xlsx")

net.ids <- paste("D21_",1:max(net.metadata$`Network ID`),sep="")#na.omit(unique(link.collection$network_id))
net.metadata$network_id <- paste(net.metadata$`Network ID`,"_",net.metadata$Country,sep="")
# it's ok, needs to be done twice because of "Trinidad & Tobago"
net.metadata$network_id <- str_replace(net.metadata$network_id," ","_")
net.metadata$network_id <- str_replace(net.metadata$network_id," ","_")

d21.matches <- str_match(net.metadata$Reference, "^([A-Za-zÀ-ÖØ-öø-ÿ'’]+).*?(\\(?\\d{4}\\)?)")
year <- gsub("[()]", "", d21.matches[,3])
d21.citations <- paste0(d21.matches[,2], "_", year)

D21_refs <- unique(data.frame(network_db = "D21",
                        network_id = net.metadata$network_id,
                        reference = net.metadata$Reference,
                        citation = d21.citations))

# update citation when NA_NA
# assume that all of these are from the same "study" as they are all
# their own unpublished data
D21_refs$citation[which(D21_refs$citation == "NA_NA")] <- "D21_unpublished_data"

db_refs_list[[length(db_refs_list)+1]] <- D21_refs
write.csv2(D21_refs,file = "data/dataset_references/D21_references.csv")

# -------------------------------------------------------------------------
# B22

df1 <- read.csv("data_clean/raw_interaction_datasets/Brimacombe_2022/information_binary_networks.csv")

orig.names <- df1$Name

clean.b22.names <- str_replace(orig.names, "^[^_]+_[^_]+_([A-Za-zÀ-ÖØ-öø-ÿ'’_]+)", "\\1")
# no real references, but unique names for each study
B22_refs <- unique(data.frame(network_db = "B22",
                       network_id = df1$Name,
                       reference = NA,
                       citation = clean.b22.names))

db_refs_list[[length(db_refs_list)+1]] <- B22_refs
write.csv2(B22_refs,file = "data/dataset_references/B22_references.csv")

# -------------------------------------------------------------------------
# FS20

load(file = "data_clean/raw_interaction_datasets/Fricke_2020/homogenization.RData")
fs <- net.long
fs20.orig.refs <- read.csv("data_clean/raw_interaction_datasets/Fricke_2020/network references.csv")

fs$network_db <- "FS20"

fs <- fs[,c("net.id","study.id","year","network_db","latitude","longitude","animal.accepted.species","plant.accepted.species",
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
fs.nets <- fs[,c("net.id","study.id","network_db","year","latitude","longitude","main_animal_taxa")]
names(fs.nets) <- c("original_id","study.id","network_db","network_year","network_lat","network_lon","main_animal_taxa")
fs.nets$network_topology_type <- "bipartite"
fs.nets$network_spatial_type <- "Point"
fs.nets$network_interaction <- "frugivory"
fs.nets$habitat_type <- NA

fs.nets <- arrange(unique(fs.nets),original_id)
fs.nets$network_id <- paste("FS20_",formatC(1:nrow(fs.nets), width = 3, format = "d", flag = "0"),sep="")
fs.nets <- fs.nets[,c("network_id","study.id","original_id","network_db","network_interaction",
                      "main_animal_taxa",
                      "network_year","network_lat","network_lon",
                      "habitat_type",
                      "network_topology_type","network_spatial_type")]

fs.nets$ref <- fs20.orig.refs$Reference[match(fs.nets$study.id,fs20.orig.refs$Study.ID)]

FS20_refs <- unique(data.frame(network_db = "FS20",
                        network_id = fs.nets$original_id,
                        reference = fs.nets$ref,
                        citation = fs.nets$study.id))

db_refs_list[[length(db_refs_list)+1]] <- FS20_refs
write.csv2(FS20_refs,file = "data/dataset_references/FS20_references.csv")

# -------------------------------------------------------------------------

# put everything together
db_refs_df <- bind_rows(db_refs_list)
# unique.refs <- sort(unique(db_refs_df$citation))

# From this, curate by hand the citations to mark those that come 
# from the same study, and harmonize them
# write.csv2(db_refs_df,file = "data/dataset_references/db_references.csv")
# the curated citations are in this file:

curated_refs <- read.csv2("data/dataset_references/curated_duplicated_references.csv")

# then, append the harmonised "citation" to the network_collection table
db_refs_full <- unique(left_join(db_refs_df,curated_refs, join_by(citation == study_orig)))

names(db_refs_full)[which(names(db_refs_full) == "network_id")] <- "original_id"
write.csv2(db_refs_full,file = "data/dataset_references/citation_list.csv",row.names = F)
