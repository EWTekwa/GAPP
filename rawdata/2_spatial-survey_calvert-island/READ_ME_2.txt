Part 2 - Spatial Survey - Calvert Island

[A] LIBRARY PREP:
12S: 3x MiSeq V3 runs using MiFish-U 12S primers
COI: 6x MiSeq V3 runs using Leray COI primers


[B] DATA FILE OVERVIEW:
1. METADATA
Vector_metadata_MAL_Feb3,2023.csv -> Metadata for all samples and controls included in the sequencing run

2. 12S
DATA PENDING

3. COI
sequence_table.calvert2021_CO1.libraries1_6.combined.w_ASV_names.txt -> 12S ASV table (raw, not filtered, not QCed)
calvert2021_CO1.libraries1_6.combined.ASV_sequences.fasta -> List of all COI ASV sequences
taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.txt -> Annotation of COI ASVs using hybrid approach of LCA + best hit
CO1_ASV_sequences.blast_90.out -> Full blast output from all COI ASVs (minimum identity = 90%)
CO1_ASV_sequences.blast_96.out ->Full blast output from all COI ASVs (minimum identity = 96%)
calvert2021_combined.CO1.tar.gz -> Full, raw, bioinformatics output from Evan
