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

# common libraries
import csv
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


# set up classes


class classname(object):

    def __init__(self, foo1, foo2):
        self.foo1 = foo1
        self.foo2 = foo2

# set up functions


def namefunction(foo):
    with open(foo) as f:
        for line in f:
            if line.startswith('@'):
                pass
            else:
                fields = line.strip().split()
                name = fields[0]
                tlen = fields[8]
                yield classname(name, tlen)

# basic-body


if __name__ == "__main__":

    print("hello world")
    # read in fasta from ncbi
    # adjust fasta headers
    


# write code here
# don't forget logging.debug('message' or object) periodically


# end
logging.info('Fin!')
