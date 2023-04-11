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
gbifIDs_obis <- get_gbifid(sp_spec$species) #get gbif IDs
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
taxize_obis_BL <- tol_resolve(sp_spec_BL$species) #get names
gbifIDs_obis_BL <- get_gbifid_(sp_spec_BL$species) #get gbif IDs

gbifIDs_obis_BL_df <- data.frame(matrix(ncol = 22, nrow = 0))

for (i in 1:1000) {
  temp1 <- as.data.frame(gbifIDs_obis_BL[i])
  colnms <- colnames(gbifIDs_obis_BL[[i]])
  temp2 <- temp1 %>%
    `colnames<-`(colnms) %>%
    select(c("usagekey", "scientificname", "rank", "status", "matchtype", "canonicalname",   
             "confidence", "kingdom", "phylum", "order", "family", "genus",           
             "species", "kingdomkey", "phylumkey", "orderkey", "familykey", "genuskey",        
             "specieskey", "synonym")) %>%
    mutate(query = names(gbifIDs_obis_BL)[i])
  gbifIDs_obis_BL_df <- rbind(gbifIDs_obis_BL_df, temp2)
  }


length(gbifIDs_obis_BL)
str(gbifIDs_obis_BL)
aa1 <- cbind(gbifIDs_obis_BL) 
#try get_gbifid_() to get all names



a2 <- taxize_obis_BL %>% #get higher taxonomy
  mutate(gbifIDs = gbifIDs_obis)
classification <- classification(gbifIDs_obis_BL)
c2 <- cbind(classification)                                     #THIS ISN'T WORKING##########################################
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
write_csv(e2, "./processeddata/species_lists/region/20230410_fish_top500_taxized.csv")


