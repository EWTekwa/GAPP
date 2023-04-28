#!/usr/bin/env python3
import os
from Bio.Blast import NCBIWWW
NCBIWWW.email = "benjamin.millard-martin@mail.mcgill.ca"
from Bio import Entrez
Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"
Entrez.email = "benmillardmartin@gmail.com"

print(os.getcwd())

sequence_data = open("12S_ASV_sequences.length_var_subset.fasta").read() 
print(sequence_data) 
print(type(sequence_data))


print("data in")
result_handle = NCBIWWW.qblast("blastn", 
                               "nt",
                               sequence_data, 
                               hitlist_size = 20, 
                               format_type = 'XML',
                               entrez_query = "txid7742 [ORGN]")

# plants (taxid:3193) Arthropoda (taxid:6656) bacteria (taxid:2) vertebrates (taxid:7742)

print(result_handle)
with open('blast.xml', 'w') as save_file: 
  blast_results = result_handle.read() 
  save_file.write(blast_results)
print("done")

#writing tab delimited blast results
from Bio import SearchIO
qresults = SearchIO.parse("blast.xml", "blast-xml")  # read XML file
SearchIO.write(qresults, "results.tab", "blast-tab")  # write to tabular file
(3, 4, 239, 277, 277)
