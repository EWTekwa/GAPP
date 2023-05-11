#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script extracts mitofish information from CRABS output

This assumes you've run CRABS with mitofish and assigned taxonomy to the original database.

By extracting all taxonomy from this output, it is possible to eliminate species that CRABS will extract from mitofish, since it is written in such a way that the whole mitofish database is always pulled without filter. Mitofish guarantees a 12S sequence of a length not too long for cutadapt, so these sequences should be preferred.

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
import os

# common libraries
from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'mitofish-extract-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions


def extractMitofish(crabs_output):
    with open(crabs_output, 'r') as f:
        for line in f:
            # no headers
            # tab separated
            # extract column 0, 2, 6, 8
            # accession, ?uid?, family, species

# basic-body

if __name__ == "__main__":

    print("hello world")
    # read in fasta from crabs
    # adjust fasta headers
    crabs_mitofish = ("processeddata/species_lists/ncbi/20230511_mitofishdb_tax.tsv")
    extractMitofish(crabs_mitofish)



# end
logging.info('Fin!')
