#figuring out how to open blast files

library(XML)
library(tidyverse)

results <- xmlParse("scripts/python/blast.xml")
df <- xmlToDataFrame(results)
list<- xmlToList(results)
results

#columns are:
  # Entrez query: taxa ID, and other filters
  # Query ID (ASV###)
  # GI and accession numbers
  # Percent identity
  # alignment length (sequence overlap)
  # mismatch (number of mismatches)
  # gapopen (number of gap openings)
  # qstart (start of allignment in query)
  # qend (end of allignment in query)
  # sstart (start of alignment in subject)
  # send (end of alignment in subject)
  # e-value
  # Bit Score?
  
  
t1 <- read_delim("scripts/Python/blastn_all_sub.tab",
                 col_names = FALSE, delim = "\t")  %>%
  `colnames<-`(c("query", "consensus_group", "ID", "percent_identity", "align_length", "mismatch", 
                 "gapopen", "qstart", "qend", "sstart", "send", "e_value", "bit_score")) %>%
  separate(col = ID, into = c("del", "gi", "del2", "accession", "del3"), sep = "\\|") %>%
  dplyr::select(-c("del", "del2", "del3"))
