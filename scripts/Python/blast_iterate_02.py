#!/usr/bin/env python3

#       This code takes 10 consensus sequences from the 10 orders that were
#   present in the top500 blast search, and which are found in the region (see
#   R script 05b),and blasts them against NCBI nucleotide databases filtered to 
#   each species (R script 05a) found in the region (OBIS). A returned hit with a reasonable
#   percent identity, coverage, and e-value indicate that the primer targeted region for that
#   taxa exists in the NCBI database. It does not indicate how well the primer
#   works. It can be read into R using script 05c. Also, requires that a species has an
#   NCBI taxonomy ID. Those without are discarded. For Mifish it is ~1

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

fasta = open("top500_consensus_sub.fasta").read()
print(fasta)
with open("NCBI_IDs_sub.txt", "r") as f2:
    IDstring = f2.read()
print(IDstring)

IDs = IDstring.split(sep= ", ")
print(IDs)

#with open('ncbi_IDs.txt') as csvfile:
#    IDs = f2.read()


for i in IDs:
    blastn_cline = NCBIWWW.qblast("blastn", 
                               "nt",
                               fasta, 
                               hitlist_size = 1, 
                               format_type = 'XML',
                               entrez_query= i)
    with open('blast1.xml', 'w') as save_file: 
        blast_results = blastn_cline.read() #append
        save_file.write(blast_results)
    print("3")
    qresults = SearchIO.parse("blast1.xml", "blast-xml")  # read XML file and parse
    print("4")
    SearchIO.write(qresults, "results.tab", "blast-tab") # save as .tab

    infile = open("results.tab")
    outfile = open("results1.tab", "w")
    for line in infile:
        col = i + "\t" + line
        columns = col.split("\t")
        print(columns)
        print(col)
        outfile.write("\t".join(columns)+"\n")

    with open("results1.tab", "r") as f2:
        data = f2.read()

    print(i)
    with open("blastn_all.tab", "a") as f:
        f.write(data)
    print("done")
print("done")



#record_iterator = SeqIO.parse("asv_table_sub.fasta", "fasta")
#print("1")
##print(record_iterator)
#print("2")
#for record in record_iterator:
#    entry = str(">" + record.description + "\n" + record.seq)
#    print("entry")
#    blastn_cline = NCBIWWW.qblast("blastn", 
#                               "nt",
#                               entry, 
#                               hitlist_size = 2, 
#                               format_type = 'XML',
#                               entrez_query= "txid7898 [ORGN]")    
##print(blastn_cline)
#    with open('blast1.xml', 'w') as save_file: 
#        blast_results = blastn_cline.read() #append
#        save_file.write(blast_results)
##    print("3")
#    qresults = SearchIO.parse("blast1.xml", "blast-xml")  # read XML file and parse
#    SearchIO.write(qresults, "results4.tab", "blast-tab") # save as .tab
#    with open("results4.tab", "r") as f2:
#        data = f2.read()
#    print(data)
#    with open("blastn_all.tab", "a") as f:
#        f.write(data)
#    print("done")


