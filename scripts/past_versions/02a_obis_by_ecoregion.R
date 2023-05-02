library(here)
library(sf)
library(tidyverse)
library(robis)

# World map
library(rnaturalearth)
world_map <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))

#meow <- sf::read_sf("rawdata/spatial_data/MEOW/meow_ecos.shp")
meow <- sf::read_sf("rawdata/spatial_data/MEOW/meow_ecos_simple.shp")
ppow <- sf::read_sf("rawdata/spatial_data/MEOW/ppow_simple_NEP.shp")

m1 <- meow %>%
  dplyr::select(c('ECOREGION', 'geometry')) %>%
  rename(region = ECOREGION) %>%
  mutate(type = 'meow')
m2 <- ppow %>%
  dplyr::select(c('PROVINC', 'geometry')) %>%
  rename(region = PROVINC) %>%
  mutate(type = 'ppow')

geo <- rbind(m1,m2)

#simplified in QGIS using Simplify function: method = Distance (Douglas-Peuker), Tolerance = 1
#st_simplify doesn't work great in R across dateline
#     meow1 <- st_simplify(meow)
#     plot(meow1)
#this isn't perfect. gaps and overlaps between ecoregions exist. The gaps are probably more
#concerning. e.g. Bamfield is in a gap. 

# Base map
kk <- ggplot() +
  geom_sf(data = world_map, size = .2, fill = "gray80", col = "gray90") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))
#Ecoregions
kk+  
  geom_sf(data = meow, aes(fill = ECOREGION), size = .2, col = 0, alpha=.3)+
  ggtitle("MEOW REALM")+
  ggtitle(paste("Marine Ecoregions of the World  (MEOW)(Spalding et al., 2007) - ", length(unique(meow$ECOREGION)),"Ecoregions"))+
  theme(legend.position = "none") +
  coord_sf(expand = FALSE)

#meow
meow$WKT <- st_as_text(meow$geometry) #creat well known text of ecoregions
st_is_valid(meow$geometry, reason = TRUE)
meow_NEP <- meow %>% #filter for regions of interest
  filter(ECOREGION %in% c("Aleutian Islands", "Gulf of Alaska", "North American Pacific Fijordland", 
                        "Puget Trough/Georgia Basin", "Oregon, Washington, Vancouver Coast and Shelf",
                        "Northern California", "Southern California Bight"))

meow_speclist <- data.frame(matrix(ncol = 11, nrow = 0))

for (i in 1:7) {
  region <- meow_NEP[i,]
  taxa <- checklist(geometry = region$WKT,
                    startdate = '1945-01-01') %>%
    select(c("scientificName", "taxonID", "taxonRank", "kingdom", "phylum", "class", "order", 
             "family", "genus", "species"))
  taxa$region <- paste(region$ECOREGION)
  meow_speclist <- rbind(meow_speclist, taxa) 
}

write_csv(meow_speclist, "./processeddata/species_lists/region/20230410_obis_by_ecoregion.csv")
#unique(meow_speclist$ecoregion) # got them all


#ppow - pelagic provinces
# selected 4 provinces covering Aleutian islands to California in a box
# I had to clip these provinces because of geometry errors - maybe there is a work around here
#     - provinces span entire oceans (east-west) 

ppow$WKT <- st_as_text(ppow$geometry) #creat well known text of ecoregions
st_is_valid(ppow$geometry, reason = TRUE)

ppow_NEP <- ppow %>% #filter for regions of interest
  filter(PROVINC %in% c("California Current", "Subarctic Pacific", "North Pacific Transitional", 
                          "North Central Pacific Gyre"))

ppow_speclist <- data.frame(matrix(ncol = 11, nrow = 0))

for (i in 1:4) {
  region <- ppow_NEP[i,]
  taxa <- checklist(geometry = region$WKT,
                    startdate = '1945-01-01') %>%
    select(c("scientificName", "taxonID", "taxonRank", "kingdom", "phylum", "class", "order", 
             "family", "genus", "species", ))
  taxa$region <- paste(region$PROVINC)
  ppow_speclist <- rbind(ppow_speclist, taxa) 
}
#unique(ppow_speclist$province) # got them all
write_csv(ppow_speclist, "./processeddata/species_lists/region/20230410_obis_by_pelagic_province.csv")

speclists <- rbind(ppow_speclist, meow_speclist)

write_csv(speclists, "./processeddata/species_lists/region/20230410_obis_by_meow_ppow_NEP.csv")





