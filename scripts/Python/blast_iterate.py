#!/usr/bin/env python3
import os
from Bio.Blast import NCBIWWW
NCBIWWW.email = "benjamin.millard-martin@mail.mcgill.ca"
from Bio import Entrez
Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"
Entrez.email = "benmillardmartin@gmail.com"

print(os.getcwd())

blastn_results = []
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import SearchIO

#option one
for record in SeqIO.parse("asv_table_sub.fasta", "fasta"):
    print(record.id)


#option two
record_iterator = SeqIO.parse("asv_table_sub.fasta", "fasta")
print("1")
print(record_iterator)
print("2")
for record in record_iterator:
    entry = str(">" + i.description + "\n" + i.seq)
    f1 = open("test.txt", "w")
    f1.write(entry)
    f1.close()
    f2 = open("test.txt", "r")
    blastn_cline = NCBIWWW.qblast("blastn", 
                               "nt",
                               "test.txt", 
                               hitlist_size = 20, 
                               format_type = 'XML')    
    res = blastn_cline()
    blastn_results.append(res)
    f2.close()
SearchIO.write(blastn_results, "results.tab", "blast-tab")  # write to tabular file
(3, 4, 239, 277, 277)
