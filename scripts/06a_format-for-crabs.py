#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This updates our output from 05 to read into crabs

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
import os
import re

# common libraries
from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'script-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')

# set up functions
def extractMitofish(crabs_output, ready_out_path):
    with open(crabs_output, 'r') as f:
        with open(ready_out_path + '20230511_test_crabs-ready.fasta', 'w') as out:
            for line in f:
                if line.startswith('>'):
                    # extract ID
                    accn = line[1:-1]
                    # modify ID
                    accn_mod = re.sub('\.[1-9]', '', line)
                    #logging.info(accn_mod)
                    # writes accnmod newline then sequence then newline
                    out.write(accn_mod)
                else:
                    out.write(line)

# basic-body

if __name__ == "__main__":

    # read in fasta from crabs: pga results
    # adjust fasta headers
    ncbi_database = ("processeddata/crabs/20230504_test.fasta")
    crabsready_out_path = ("./processeddata/crabs/")
    extractMitofish(ncbi_database, crabsready_out_path)


# end
logging.info('Fin!')
