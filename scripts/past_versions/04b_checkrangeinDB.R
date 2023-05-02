library(here)
library(tidyverse)
library(rfishbase)


speclists <- read_csv("./processeddata/species_lists/region/20230412_obis-fish_NEP.csv")

#fish by region, add fishbase names
fish <- speclists %>%
  filter(!is.na(species)) %>%
  filter(class %in% c('Actinopteri', 'Elasmobranchii', 'Holocephali',
                      'Petromyzonti', 'Myxini') |
           order %in% c('Scorpaeniformes')) %>%
  mutate(fishbase_name = validate_names(species_list = .$accepted_name_usage))

#fishlist, add fishbase names
spec_fish <- fish %>%
  dplyr::select(-c('region')) %>%
  distinct() %>%
  mutate(fishbase_name = validate_names(species_list = .$accepted_name_usage))

#get fishbase names
spec_fish$fishbase_name <- validate_names(species_list = spec_fish$accepted_name_usage)

#get distribution for fishbase names
dist_fish <- distribution(species_list = spec_fish$fishbase_name, server = 'fishbase')

#filter for species in region
fishbase_in_FAO <- dist_fish %>%
  filter(FAO %in% 'Pacific, Northeast' |
           FAO %in% 'Pacific, Eastern Central' |
           FAO %in% 'America, North - Inland waters') %>%
  select(c('Species')) %>%
  distinct()

#filter dataset (fish by region) for fish confirmed by fishbase
fish_by_region <- fish %>%
  filter(fishbase_name %in% fishbase_in_FAO$Species)
length(unique(fish_by_region$species)) #retained 2360 species

#what fish are not in the region (according to fishbase) and where are they from?
fishbase_not_in_FAO <- dist_fish %>%
  filter(!Species %in% fishbase_in_FAO$Species)
length(unique(fishbase_not_in_FAO$Species)) #removed 379 species

write.csv(fish_by_region, "./processeddata/species_lists/region/20230412_obis-fish_fishbase_NEP.csv")

inverts <- read.csv('./processeddata/species_lists/region/20230412_obis-invertebrates_nep.csv')

spec_inv <- speclists %>%
  filter(!is.na(species)) %>%
  filter(kingdom %in% 'Animalia' & 
           !class %in% c('Insecta', 'Aves', 'Chilopoda',) &
           !order %in% c('Araneae', 'Pseudoscorpionida'))

write.csv(spec_inv, "./processeddata/species_lists/region/20230412_obis-inverts_NEP.csv")

dist_inv <- distribution(species_list = spec_inv$species, server = 'sealifebase')

