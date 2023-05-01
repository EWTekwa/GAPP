#!/usr/bin/env python3
import os
from Bio.Blast import NCBIWWW
NCBIWWW.email = "benjamin.millard-martin@mail.mcgill.ca"
from Bio import Entrez
Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"
Entrez.email = "benmillardmartin@gmail.com"

print(os.getcwd())

#blastn_all = []
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import SearchIO

fasta = SeqIO.parse("top500_consensus.fasta", "fasta")


record_iterator = SeqIO.parse("asv_table_sub.fasta", "fasta")
print("1")
#print(record_iterator)
print("2")
for record in record_iterator:
    entry = str(">" + record.description + "\n" + record.seq)
    print("entry")
    blastn_cline = NCBIWWW.qblast("blastn", 
                               "nt",
                               entry, 
                               hitlist_size = 2, 
                               format_type = 'XML',
                               entrez_query= "txid7898 [ORGN]")    
#print(blastn_cline)
    with open('blast1.xml', 'w') as save_file: 
        blast_results = blastn_cline.read() #append
        save_file.write(blast_results)
    print("3")
    qresults = SearchIO.parse("blast1.xml", "blast-xml")  # read XML file and parse
    SearchIO.write(qresults, "results4.tab", "blast-tab") # save as .tab
    with open("results4.tab", "r") as f2:
        data = f2.read()
    print(data)
    with open("blastn_all.tab", "a") as f:
        f.write(data)
    print("done")


