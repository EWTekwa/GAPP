#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script formats CRABS output to our specifications


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

LOG_FILENAME = 'crabs-refine-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions


def extractMitofish(crabs_output):
    with open(crabs_output, 'r') as f:
        for line in f:
            if line.startswith('>'):
                #modify line
                logging.info('line')
            else:
                pass)

# basic-body

if __name__ == "__main__":

    print("hello world")
    # read in fasta from crabs
    # adjust fasta headers
    crabs_mitofish = "path/to/fasta"
    extractMitofish(crabs_mitofish)



# end
logging.info('Fin!')
