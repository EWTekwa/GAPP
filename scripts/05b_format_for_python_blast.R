#format for Python Blast search

#


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER", force = TRUE)
BiocManager::install("msa", force = TRUE)
library(here)
library(tidyverse)
library(taxonomizr)
library(Biostrings)
library(msa)
library(phytools)

#fish
fish <- read.csv("./processeddata/species_lists/region/20230412_obis-fish_fishbase_NEP.csv") %>%
#  select(c("ncbi_id")) %>%
#  distinct() %>% 
  filter(is_marine == "TRUE" | is_brackish == "TRUE" & is_freshwater == "FALSE",
         region == "North American Pacific Fijordland")

string <- unlist(as.vector(fish))

#merge top500 (after updating taxonomy) and fasta, group by family, make allignment, make consensus sequence - make consensus for all

top500_tax <- read_csv('rawdata/top500_20230405/20230413_top500-worrms.csv') 

ASVs <- read.delim("rawdata/top500_20230405/blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out", 
                          h = TRUE, fill = TRUE) %>% 
  # standardize names and remove overhanging x
  janitor::clean_names() %>%
  rename_with(., ~str_remove_all(string = .,
                                 pattern = "x_")) %>%
  separate_wider_delim(taxonomy, delim = ' / ', 
                       names = c('kingdom', 'phylum', 'class', 'order', 
                                 'family', 'genus', 'species')) %>%
  # column to match existing scripts
  dplyr::rename(asv = query_id) %>%
  filter(kingdom == "Eukaryota") %>%
  filter(identity_percentage > 98)

fastaFile <- readDNAStringSet("rawdata/top500_20230405/12S_ASV_sequences.length_var.fasta")
asv = names(fastaFile)
sequence = paste(fastaFile)
asv_sequences <- data.frame(asv, sequence)

#use NCBI_species
t1 <- merge(ASVs[c("asv", "species")],top500_tax, by.x = "species", by.y = "ncbi_species")

t2 <- merge(t1, asv_sequences, by = "asv", all.x = T)
unique(t2$order)

unique(fish$order)

#msaConvert()

Clupeiformes <- t2 %>%
  filter(order == "Clupeiformes") %>%
  dplyr::select(c("sequence")) %>%
  distinct()
Clupeiformes_set <- DNAStringSet(Clupeiformes$sequence)
Clupeiformes_aliC <- msa(Clupeiformes_set, method = "ClustalW") 
Clupeiformes_aliC
consensus <- msaConsensusSequence(Clupeiformes_aliM, type = "Biostrings")
consensus


mitofish <- readDNAStringSet("rawdata/MitoFish/mito-all.fasta")
mito_align <- msa(mitofish, method = "ClustalW") 




#all data aligned and a consensus provided
alignment <- msa(fastaFile, method = "ClustalW")
alignment_muscle <- msa(fastaFile, method = "Muscle")

consensus <- msaConsensusSequence(alignment, type = "Biostrings")
consensus

#writing a fasta file
Xfasta <- character(nrow(X) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", X$column1)
Xfasta[c(FALSE, TRUE)] <- X$column2

writeLines(Xfasta, "filename.fasta")

