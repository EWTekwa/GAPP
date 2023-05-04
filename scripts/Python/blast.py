#!/usr/bin/env python3
import os
from Bio.Blast import NCBIWWW
NCBIWWW.email = "benjamin.millard-martin@mail.mcgill.ca"
from Bio import Entrez
Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"
Entrez.email = "benmillardmartin@gmail.com"

print(os.getcwd())

fasta = open("top500_consensus.fasta").read()


print("data in")
result_handle = NCBIWWW.qblast("blastn", 
                               "nt",
                               fasta, 
                               hitlist_size = 20, 
                               format_type = 'XML',
                               entrez_query = "txid3193[ORGN]")

# plants (taxid:3193) Arthropoda (taxid:6656) bacteria (taxid:2) vertebrates (taxid:7742)

print(result_handle)
with open('blast_a.xml', 'w') as save_file: 
  blast_results = result_handle.read() 
  save_file.write(blast_results)
print("done")

#writing tab delimited blast results
from Bio import SearchIO
qresults = SearchIO.parse("blast_a.xml", "blast-xml")  # read XML file
SearchIO.write(qresults, "results_a.tab", "blast-tab")  # write to tabular file

