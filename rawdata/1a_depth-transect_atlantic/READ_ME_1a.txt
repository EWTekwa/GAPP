Part 1a - Depth transect - Atlantic


[A] LIBRARY PREP:
12S: 1x MiSeq V3 run using MiFish-U 12S primers
COI: 1x MiSeq V3 run using Leray COI primers


[B] DATA FILE OVERVIEW:
1. METADATA
Vector_metadata_MAL_Feb3,2023.csv -> Metadata for all samples and controls included in the sequencing run

2. 12S
DATA PENDING

3. COI
sequence_table.CO1.merged.w_ASV_names.txt -> 12S ASV table (raw, not filtered, not QCed)
CO1_ASV_sequences.no_blast_hit.fasta -> List of all COI ASV sequences
taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA+best_hit.txt -> Annotation of COI ASVs using hybrid approach of LCA + best hit
taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA_only.txt -> Annotation of COI ASVs using LCA only
CO1_ASV_sequences.blast_90.out -> Full blast output from all COI ASVs (minimum identity = 90%)
CO1_ASV_sequences.blast_96.out ->Full blast output from all COI ASVs (minimum identity = 96%)
eastcoast_Run20230127.CO1.tar.gz -> Full, raw, bioinformatics output from Evan
