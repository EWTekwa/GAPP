library(here)
library(tidyverse)
library(rfishbase)


speclists <- read_csv("./processeddata/species_lists/region/20230410_obis_by_meow_ppow_NEP.csv")

spec_fish <- speclists %>%
  filter(!is.na(species)) %>%
  filter(class %in% c('Actinopteri', 'Elasmobranchii', 'Holocephali',
                      'Petromyzonti', 'Myxini') |
           order %in% c('Scorpaeniformes'))

#get distribution
dist_fish <- distribution(species_list = spec_fish$species, server = 'fishbase')
unique(dist_fish$FAO)

#filter for species in region
a1 <- dist_fish %>%
  filter(FAO %in% 'Pacific, Northwest' |
         FAO %in% 'Pacific, Eastern Central' |
         FAO %in% 'America, North - Inland waters') %>%
  select(c('Species')) %>%
  distinct()

#
a2 <- spec_fish$species[(!spec_fish$species %in% a1$Species)] %>%
  as.data.frame()








spec_inv <- speclists %>%
  filter(!is.na(species)) %>%
  filter(kingdom %in% 'Animalia' & 
           !class %in% c('Insecta', 'Aves', 'Chilopoda',) &
           !order %in% c('Araneae', 'Pseudoscorpionida'))

dist_inv <- distribution(species_list = spec_inv$species, server = 'sealifebase')

