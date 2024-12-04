
# plot network coordinates, Fig. 1 of the manuscript

# -------------------------------------------------------------------------
library(tidyverse)
library(sf)
library("rnaturalearth") 
library("rnaturalearthdata")
library(colorblindr)

# -------------------------------------------------------------------------
all.nets <- read.csv2("results/network_collection.csv")

# -------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") 

sub.nets <- subset(all.nets,network_interaction %in% c("frugivory","parasitism","food web","pollination","herbivory"))
# tt <- subset(sub.nets, network_id == "TI1")
# herb.nets <- subset(sub.nets, network_interaction %in% c("herbivory"))

pw <- ggplot(data = world) + 
  geom_sf(color = "grey80", fill = "grey80") +
  # geom_point(data = herb.nets, aes(x = network_lon, y = network_lat)) +
  # geom_text(data = herb.nets, aes(x = network_lon, y = network_lat, label = network_id))
  geom_point(data = sub.nets, aes(x = network_lon, y = network_lat, color = network_interaction)) +
  labs(x="",y="") +
  scale_color_OkabeIto(name = "Network type") +
  theme_bw() +
  coord_sf(datum = NA) +
  NULL

# pw

# ggsave(filename = "results/images/network_locations.pdf",plot = pw,
#        device = cairo_pdf,
#        width = 12,height = 5,dpi = 300)


# -------------------------------------------------------------------------

# tt <- subset(sub.nets, network_id %in% c("TI1","TI129","mangal_88","TI7"))
