#!/bin/bash

# starts in project root
# move to folder for data download/processing
cd processeddata/crabs

# taxonomy
# This file needs to be downloaded if we use
# NCBI etc data, even a custom database
crabs db_download --source taxonomy
