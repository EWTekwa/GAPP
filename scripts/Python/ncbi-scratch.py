#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
this is a basic template for a python script bc I'm lazy
I also do the same thing a lot. Anyway, notes go here

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
import logging
from pyprojroot import here  # like here in R
import os

# common libraries and login
from Bio.Blast import NCBIWWW

from Bio import Entrez
from Bio import SeqIO
Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"


# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development

LOG_FILENAME = 'blast-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions

def get_search_terms(termsfile):
    with open(termsfile) as f:
        # we just need the first line, technically
        # make it into a list
        first = f.readline()
        #data_list = first.split(",")
        #logging.info(f"{data_list}")
        #return(data_list)
        logging.info(f"{first}")
        return(str(first))
        yield classname(name, tlen)

# basic-body


if __name__ == "__main__":
    # set emails
    # NCBIWWW.email = "kate.sheridan@mail.mcgill.ca"
    Entrez.email = "kate.sheridan@mail.mcgill.ca"
    # testing queries
    handle = Entrez.egquery(term="Gadus")
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"] == "nuccore":
            print(row["Count"])

    # downloading records
    handle = Entrez.esearch(db="nuccore", term="Gadus AND 12S",
                            retmax=12, usehistory="y")
    record = Entrez.read(handle)
    # note i may need to change this to accession number
    gi_list = record["IdList"]
    gi_str = ",".join(gi_list)
    handle = Entrez.efetch(db="nuccore", id=gi_str,
                           rettype="gb", retmode="text", usehistory="y")
    records = SeqIO.parse(handle, "gb")
    records
    for record in records:
        print(f"{record.name}")

    # downloading sequences
    search_handle = Entrez.esearch(db="nucleotide", term="Gadus and 12S",
                                   retmax=12, idtype="acc", usehistory="y"
                                   )
    search_results = Entrez.read(search_handle)
    search_handle.close()
    acc_list = search_results["IdList"]
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    batch_size = 3
    out_handle = open("testing/gadus.fasta", "w")
for start in range(0, count, batch_size):
    end = min(count, start + batch_size)
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
    )
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()

print("hello world")


# write code here
# don't forget logging.debug('message' or object) periodically


# end
logging.info('Fin!')
