#taxize
library(here)
library(tidyverse)
library(taxize)

#data ####
top500 <- read.delim("rawdata/top500_20230405/blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out",
                     h=TRUE,
                     fill = TRUE) %>%
  `colnames<-`(c("ASV", "subject", "accesion_num", "taxa_ID", "perc_ID", "coverage", "evalue", "bitscore", "source", "taxonomy")) %>%
  as.data.frame() %>%
  na.exclude() %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = " / ") 


#Obis fish - try might have to do something different for invertebrates because of multiple names
speclists <- read_csv("./processeddata/species_lists/region/20230410_obis_by_meow_ppow_NEP.csv")
sp_spec <- speclists %>% #just obs with species assignment
  filter(!is.na(species)) %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes")%>%
  select(species) %>%
  distinct() 
taxize_obis <- tol_resolve(sp_spec$species) #get names
gbifIDs_obis <- get_gbifid(sp_spec$species, rows = 1) #get gbif IDs
a1 <- taxize_obis %>% #get higher taxonomy
  mutate(gbifIDs = gbifIDs_obis)
classification <- classification(gbifIDs_obis)
 c1 <- cbind(classification)
d1 <- merge(a1, c1, by.x = "gbifIDs", by.y = "query") %>%
  select(c("gbifIDs", "search_string", "approximate_match", "number_matches", "kingdom", "phylum", "class", "order", 
           "family", "genus", "species")) %>%
  mutate(search_string = str_to_sentence(search_string))
e1 <- speclists %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes") %>%
  filter(taxonRank == "Species")  %>%
  select(c("species", "region")) %>%
  merge(., d1, by.x = "species", by.y = "search_string") %>%
  select(-c("species")) %>%
  rename(species = species.y)
write_csv(e1, "./processeddata/species_lists/region/20230410_fish_by_region.csv")


sp_spec_BL <- top500 %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes")  %>% 
  select(c("species")) %>%
  distinct() 
taxize_obis_BL <- tol_resolve(sp_spec_BL$species) #get names from Open tree of life
gbifIDs_obis_BL <- get_gbifid(sp_spec_BL$species, rows = 1) #get gbif IDs
#gbifIDs_obis_BL_df <- dplyr::bind_rows(gbifIDs_obis_BL, .id = "query")             ###use with get_gbifid_()
wormsIDs_obis_BL <- get_wormsid_(sp_spec_BL$species, marine_only = F, accepted = F) #get worms IDs

#filter for accepted
dd <- gbifIDs_obis_BL_df %>%
  filter(status == "ACCEPTED")

a2 <- taxize_obis_BL %>% #get higher taxonomy
  mutate(gbifIDs = gbifIDs_obis_BL)
classification_BL <- classification(gbifIDs_obis_BL)
c2 <- cbind(classification)                                     #THIS ISN'T WORKING##########################################
c2 <- dplyr::bind_rows(classification, .id = classification)
c2 <- unlist(classification)

q1 <- dplyr::bind_rows(classification[1])
q2 <- classification[1]
classification[1]





d2 <- merge(a2, c2, by.x = "gbifIDs", by.y = "query") %>%
  select(c("gbifIDs", "search_string", "approximate_match", "number_matches", "kingdom", "phylum", "class", "order", 
           "family", "genus", "species")) %>%
  mutate(search_string = str_to_sentence(search_string))
e2 <- speclists %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes") %>%
  filter(taxonRank == "Species")  %>%
  select(c("species", "region")) %>%
  merge(., d2, by.x = "species", by.y = "search_string") %>%
  select(-c("species")) %>%
  rename(species = species.y)
write_csv(e2, "./processeddata/species_lists/region/20230410_fish_top500.csv")


