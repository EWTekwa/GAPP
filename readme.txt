readme

Scripts
Pipeline

- 00-species-lists.Rmd
	Generates species lists from results, no taxonomic correction. 
	
- 01a_OBIS-species-list.Rmd
	This script pulls species data from OBIS to generate a species list for a defined region. By ecoregion section uses shapefiles from the MEOW nearshore ecoregions and PPOW pelagic provinces. Then uses fishbase to validate fish from OBIS checklist
	Output saved in processeddata>species_lists>region:
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


- 02_top500-taxize-worrms.Rmd
	Takes top 500 and extracts and cleans species names
	
	
- 03_species-in-region.Rmd
	