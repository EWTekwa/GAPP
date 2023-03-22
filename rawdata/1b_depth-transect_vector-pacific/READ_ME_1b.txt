Part 1b - Depth transect - Vector Pacific


[A] LIBRARY PREP:
12S: 1x MiSeq V3 run using MiFish-U 12S primers
COI: 1x MiSeq V3 run using Leray COI primers


[B] DATA FILE OVERVIEW:
1. METADATA
Vector_metadata_MAL_Feb3,2023.csv -> Metadata for all samples and controls included in the sequencing run

2. 12S
sequence_table.12S.merged.w_ASV_names.length_var.txt -> 12S ASV table (not filtered, not QCed, with all ASVs including non-fish)
12S_ASV_sequences.length_var.fasta -> List of all 12S ASV sequences
taxonomy_table.12S.NCBI_NT.96sim.LCA_only.txt -> Annotation of 12S ASVs using LCA  (this is a good place to start
taxonomy_table.12S.NCBI_NT.96sim.LCA+besthit.txt -> Annotation of 12S ASVs using hybrid approach of LCA + best hit
12S_ASV_sequences.length_var.blast.out -> Full blast output from all 12S ASVs
vector2021_Run20230121.12S.tar.gz -> Full, raw, bioinformatics output from Evan

3. COI
sequence_table.CO1.merged.w_ASV_names.txt -> 12S ASV table (raw, not filtered, not QCed)
CO1_ASV_sequences.fasta -> List of all COI ASV sequences
taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.txt -> Annotation of COI ASVs using hybrid approach of LCA + best hit
CO1_ASV_sequences.blast_87.out -> Full blast output from all COI ASVs (minimum identity = 87%)
CO1_ASV_sequences.blast_96.out ->Full blast output from all COI ASVs (minimum identity = 96%)
vectorPacific_Re-Run20220916.CO1.tar.gz -> Full, raw, bioinformatics output from Evan
