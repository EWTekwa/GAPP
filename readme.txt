readme

Scripts
Pipeline

- 00-species-lists.Rmd
	Generates species lists from results, no taxonomic correction. 
	
- 01a_OBIS-species-list.Rmd
	This script pulls species data from OBIS to generate a species list for a defined region. By ecoregion section uses shapefiles from the MEOW nearshore ecoregions and PPOW pelagic provinces. Then uses fishbase to validate fish from OBIS checklist
	
	Output saved in processeddata>species_lists>region>[date]_obis-[taxa]_nep, if validated by fishbase: [date]_obis-fish_fishbase_nep:
`accepted_name` -- The most up to date name from WoRMS
`verbatim_name` -- column coalesced from subspecies and species, includes unaccepted names.
`taxonomic_status` -- status of `verbatim_name`
`taxon_rank` -- level of identification
taxonomy: species, subspecies, genus, family, order, class, phylum (invertebrates only)
`is_marine`, `is_brackish`, `is_freshwater`: true/false
`taxon_id` -- should be WoRMS AphiaID for lowest found taxon
`genus_id` -- (invertebrates only) to help look for NCBI IDs
`ncbi_id` -- ID for taxon in NCBI
`region` -- region found in, from MEOW or PPOW
`fishbase_name` -- name validated in fishbase, if applicable

- 01b_GBIF-species-list.Rmd
	Currently removed from the pipeline.

- 02_top500-taxize-worrms.Rmd
	Takes top 500 and extracts, cleans species names, and populates with higher taxonomy for downstream analysis. 
	For fish, first only fish groups are selected, then taxize get_ids_() is run for 'gbif' and 'ncbi'. This removes some of the dirty names. Then the resulting list is run through taxize::get_worms_id_(sci_com, marine_only = FALSE, accepted = FALSE) to get WoRMS IDs to fetch records from worrms. Any unmatched taxa is paired with higher taxonomy for downstream analysis by extracting the genus from the verbatim_name.
	
	Output saved in processeddata>species_lists>from_data>[date]_[survey]_[taxa]-top500-worrms:
`found_taxa` - single column of found species names coalesced in this order worms > gbif > ncbi > verbatim_name. To be used downstream.
`verbatim_name` - verbatim from top500
`worms_name` - sciname from worrms
`gbif_name` - canonicalname from gbif
`ncbi_name` - scientificname from ncbi
`ncbi_commonname` - common name from ncbi
`worms_aphiaid` - unique ID for WoRMS
`gbif_usagekey` - unique ID for gbif
`ncbi_uid` - unique ID for ncbi
`gbif_status` - accepted/synonym etc from gbif
`worms_status` - accepted/synonym etc from WoRMS
taxonomy - kingdom, phylum, class, order, family, genus -- from WoRMS
`freshwater_only` - 0/1, if 1, 'freshwater' was the only field true out of 'marine', 'brackish', 'freshwater' on WoRMS

	
- 03_species-in-region.Rmd
	This script uses eulerr and arsenal to generate both proportional venn diagrams (euler plots) and a comparison of an input and whether its in the region.
	Output of comparisons saved in:
		processeddata>species_lists>comparisons>[date]_[survey]_[taxa]-[compared list]_compared-region for added binary column:
			For top 500:
			Same output as 02 with one additional column:
			`in_region` -- yes - intersects with our region list, no - does not interesect with region list, NA - unable to be compared
		processeddata>species_lists>comparisons>[date]_[survey]_[taxa]-[compared list]_out-of-region:
			`source` -- origin of taxa
			`found_taxa` -- from found_taxa column in upstream files
	Output of eulerr saved as .png in:
		output>region_comparison>[date]_[survey]_[taxa]-[compared list]_region-venn
		

- 04a_ncbi-search.Rmd
	Generates region species list in Genus_species txt format
		Output saved in processeddata>specieslists>region>[date]_[group]sp.txt
	generates NCBI search and saves as .txt file
		Output saved in processeddata>specieslists>ncbi>[date]_[parameters]-ncbisearch.txt
	generates NCBI search filtered by presence in mitofish
		output saved in processeddata>specieslists>ncbi>[date]_[parameters]-ncbisearch-nomitofish.txt
		
		Note 04_extract-mitofish-crabs.py in scripts>python for the origin of the mitofish list
	
	
 - 05 blast iterate
 
 
 
 - 06 crabs; 
 To get CRABS working: clone from https://github.com/gjeunen/reference_database_creator/tree/main, add path to crabs command, also requires vsearch, cutadapt, muscle, and the python packages argparse, biopython, tqdm, numpy, matplotlib, pandas. May need to `brew install wget` on mac.
 Citation: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13741
 
 To work through 06:
 - 06a_format-for-crabs.py 
 	Takes our custom NCBI output and reformats the header to be read in properly by CRABS
 
 - 06b_crabs_mitofish-custom-ncbi_pcr.sh
 	Prepares taxonomy
 	Loads in mitofish
 	Loads in our custom NCBI results
 	Merges the databases
 	Runs in silico silico pcr on the merged results
 	PGA/pairwise global alignment to catch trimmed sequences
 
 - 06c_format-for-taxonomy.py
 	reads in in-silico-pcr results
 	updates headers from custom db to match crabs taxonomy file
 		note this is due to a 'bug' in crabs to report as a github issue later
 	
 - 06d_crabs_dereplicate-filter.sh




 
 Overview of CRABS (remove later when our workflow is updated):
 
 step 1; databases
 	crabs can search multiple databases
 	At this step you need to download a taxonomy file too
 
 step 2: in house import
 	you can import your own files
 
 step 3: merge dbs
 	if download + custom or multiple downloads
 
 step 4: in silico pcr
 	can also do pga; pairwise global alignment
 	
 step 5: reassign taxonomy
 	Uses taxonomy files downloaded from step 1 plus --missing for custom lineages
 	
 step 6: dereplicate
	removes duplicate sequences/species
		strict: only unique sequences
		single_species: one sequence per species
		uniq_species: all unique sequences for each species
		
 step 7: filter/clean
	removes based on parameters:
	minimum length: '--minlen'
	maximum length: '--maxlen'
	number of ambiguous bases: '--maxns'
	environmental sequences: '--enviro' discard yes/no
	unspecified species name: '--species' discard yes/no
	missing taxonomic information: '--nans'	 # of unspecified taxonomic levels
		
 step 8: visualize
 	diversity = barplot of #sp and # sequences by taxonomic rank in database
 	amplicon length = frequency plot of length
 	db completeness = provide a txt with species names (include underscore), it gives you:
 		name of the species of interest
 		if species is present in the reference database (indicated by a 1 or 0)
 		number of species in the reference database that share the same genus
 		number of species in the genus according to the NCBI taxonomy
 		percentage of species in the genus present in the reference database
 		number of species in the reference database that share the same family
 		number of species in the family according to the NCBI taxonomy
 		percentage of species in the family present in the reference database
 		list of species sharing the same genus in the reference database
 		list of species sharing the same family in the reference database
 	primer_efficiency = barplot showing % basepair occurrences on each site of binder and fasta of sequences that contributed to plot
 	
 	
 step 9: export
 	This step is if you want ot use the DB for eDNA blasting
 
 
 
 
 