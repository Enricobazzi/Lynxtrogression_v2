# map plot
library(tidyverse)
library(raster)
library(rgdal)
library(sf)

# load data
coords <- read.table("data/samples_coordinates.csv", sep = ",", header = T)
wel_coords <- coords[coords$pop=="wel",]
sel_coords <- coords[coords$pop=="sel",]
lpa_coords <- coords[coords$pop=="lpa",]

world_map <- rgdal::readOGR("plots/ms_figures/ne_50m_land/ne_50m_land.shp")
extent <- c(-10, 64, 25, 75)
world_map.crop <- crop(world_map, extent)

ll_distr <- rgdal::readOGR("plots/ms_figures/redlist_species_data_ll/data_0.shp")
ll_distr.crop <- crop(ll_distr, extent)
cau <- c(25, 64, 25, 45)
ll_distr.cau <- crop(ll_distr, cau)
lp_distr <- rgdal::readOGR("plots/ms_figures/redlist_species_data_lp/data_0.shp")
lp_distr.crop <- crop(lp_distr, extent)

plot_map <- function(){
  plot(world_map.crop, col="lightgrey", axes = TRUE)
  plot(ll_distr.crop, col="#4876c1", border = NA, add=T)
  points(wel_coords$longitude, wel_coords$latitude,
         pch=21, lwd=1, cex=1.3, bg="#003c82")
  plot(ll_distr.cau, col="#7ac34b", border = NA, add=T)
  points(sel_coords$longitude, sel_coords$latitude,
         pch=21, lwd=1, cex=1.3, bg="#1b9773")
  plot(lp_distr.crop, col="#674697", border = NA, add=T)
  points(lpa_coords$longitude, lpa_coords$latitude,
         pch=21, lwd=1, cex=1.3, bg="#8d56a4")
}

pdf("plots/ms_figures/distribution_map_v2.pdf", width = 8, height = 8)
plot_map()
dev.off()
