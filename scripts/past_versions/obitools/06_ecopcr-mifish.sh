#!/bin/bash






# ecoPCR and obitools testing
obi test

# obitools is working; obicount shows this as 7847 sequences
# obicount processeddata/species_lists/ncbi/20230504_test.fasta



#obiconvert -d ./processeddata/species_lists/ncbi/ --fasta --ecopcrdb-output=sequencesdb 20230504_test.fasta --skip-on-error
