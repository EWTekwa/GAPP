#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script queries NCBI using Biopython.
It requires a list of UIDS formatted: txid###[Orgn],txid###[Orgn]
The list is queried then 3 sequences matching the query are written to a FASTA file.
Refining the query process is still in progress.

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
#import os

# common libraries
from Bio import Entrez

from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'blast-log2.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions

# extract our search from the one we set up
def get_search_terms(termsfile):
    with open(termsfile) as f:
        # we just need the first line, technically
        # make it into a list
        first = f.readline()
        data_list = first.split(",")
        logging.info(f"{data_list}")
        return(data_list)

# now we want to use the search terms interatively


def fasta_by_acc(query_list, outfile_path):
    with open(here(outfile_path+"test.fasta"), "w") as out:
        for q in query_list:
            logging.info(q + "is searching")
            handle = Entrez.esearch(
                db='nucleotide', term=q + "AND 12S NOT UNVERIFIED [word] NOT PREDICTED [word]", usehistory="y", retmax="3", idtype="acc")
            search_results = Entrez.read(handle)
            handle.close()
            webenv = search_results["WebEnv"]
            querykey = search_results["QueryKey"]
            acc_list = search_results["IdList"]
            logging.info(acc_list)
            # if we didn't get any ids, don't write
            if not acc_list:
                pass
            else:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="fasta",
                    retmode="text",
                    webenv=webenv,
                    query_key=querykey,
                    idtype="acc",
                    )
                data = fetch_handle.read()
                fetch_handle.close()
                out.write(data)


# basic-body

if __name__ == "__main__":

    # set emails
    Entrez.email = "kate.sheridan@mail.mcgill.ca"
    Entrez.api_key = "330c561fca772258d2c869e4caefe8c98b09"
    # Path to fasta output,
    # name of fasta file currently defined in function (change)
    fasta_out_path = ("./processeddata/species_lists/ncbi/")
    # path to txt for search terms
    search_terms = (
        "processeddata/species_lists/ncbi/20230501_region-fish-ncbisearch.txt")
    # search_terms = ("processeddata/species_lists/region/ncbi_IDs.txt")
    # get search terms
    query = get_search_terms(search_terms)
    fasta_by_acc(query, fasta_out_path)


# end
logging.info('Fin!')
