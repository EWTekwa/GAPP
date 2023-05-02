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
	
	