#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This updates our output from 06b to match taxonomy
Currently input from custom databases gets 'CRABS:' appended.
This is potentially a bug in CRABS to report in github issues.

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

LOG_FILENAME = 'format-taxonomy-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions


def updateTaxonomy(crabs_output):
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
    crabs_fasta = "path/to/fasta"
    updateTaxonomy(crabs_fasta)



# end
logging.info('Fin!')
