library(here)


top500 <- read.delim("rawdata/top500_20230405/blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out",
                   h=TRUE,
                   fill = TRUE) %>%
  `colnames<-`(c("ASV", "subject", "accesion_num", "taxa_ID", "perc_ID", "coverage", "evalue", "bitscore", "source", "taxonomy")) %>%
  as.data.frame() %>%
  na.exclude() %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = " / ") 

top500_fish <- top500 %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes") 







require(sf)

# World map
library(rnaturalearth)
world_map <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))


# Base map
kk <- ggplot() +
  geom_sf(data = world_map, size = .2, fill = "gray80", col = "gray90") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))

meow <- sf::read_sf("rawdata/spatial_data/MEOW/meow_ecos.shp")
kk+  
  geom_sf(data = meow, aes(fill = ECOREGION), size = .2, col = 0, alpha=.3)+
  ggtitle("MEOW REALM")+
  ggtitle(paste("Marine Ecoregions of the World  (MEOW)(Spalding et al., 2007) - ", length(unique(meow$ECOREGION)),"realms"))+
  theme(legend.position = "none") +
  coord_sf(expand = FALSE)
