#!/bin/bash

crabs -h

crabs db_import --input ./processeddata/species_lists/ncbi/20230504_test.fasta --output output.fasta --seq_header species --delim ' '

crabs insilico_pcr --input output.fasta --output pcrout.fasta --fwd GTCGGTAAAA --rev CATAGTGGGGTA --error 5

crabs pga --input pcrout.fasta --output pgaout.fasta --database output.fasta --fwd GTCGGTAAAA --rev CATAGTGGGGTA --speed medium --percid 0.95 --coverage 0.95 --filter_method relaxed
