# Initial Date: 12.20.22

# Project: Macroinvertebrate Inferred Temperature

# summary: Calculate MTTI for new data

#				


# packages used
library(tidyverse)
library(rioja)






#####
#####
#			bring in taxonomy lookup table 
#####
#####

taxa <- read.csv('_TAXA_TRANSLATOR_MASTER_20221201.csv')



# modify to make easier to join/link with bug data files
taxa.otu <- taxa %>%
	select(Taxon_orig, OTU_MTTI) 



#####
#####

#		bring in bug data

#####
#####

bugs_7977 <- read.csv('_Bio_ALL_7977.csv') # non-wadeable sites removed

# limit the bug file to three necessary columns: sample ID, taxa, counts (relative abundance)
bugs.all_lim <- bugs_7977 %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(sample.id = UniqueID_v2) %>%
	rename(Taxon_v2 = TaxaID_v2)
	

#####
#####

#		join bugs and OTUs
#####
#####

bugs.all_otu <- bugs.all_lim %>%
	left_join(taxa.otu, by = c('Taxon_v2' = 'Taxon_orig')) 


bugs.all_otu <- bugs.all_otu %>%
	select(sample.id, OTU_MTTI, RA) %>%
	filter(OTU_MTTI != 'DNI') 
# sum RA's across all OTUs--should see a reduction in rows
bugs.all_otu<-plyr::ddply(.data = bugs.all_otu, c('sample.id', 'OTU_MTTI'), 
										plyr::summarize, RA=sum(RA))



# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

bugs.all_wide <- bugs.all_otu %>% 
	pivot_wider(id_cols = sample.id, names_from = OTU_MTTI, values_from = RA,
					values_fn = sum) 
bugs.all_wide[is.na(bugs.all_wide)] <- 0
bugs.all_wide <-	column_to_rownames(bugs.all_wide, 'sample.id')






# bring in the saved MTTI

load('MTTI_dec2022.Rdata')


@@@@@ carefully think about what we want to title this model for distribution
		currently titled "wa_high.res"


mtti.all.samps <- predict(wa_high.res, newdata=bugs.all_wide, sse=TRUE, nboot=100,
									match.data=TRUE, verbose=TRUE)

mtti.all.samps_scores <- as.data.frame(mtti.all.samps$fit)


mtti.all.samps_scores$sample.id <- row.names(mtti.all.samps_scores) 
mtti.all.samps_scores <- mtti.all.samps_scores %>%
	select(sample.id, everything()) 
row.names(mtti.all.samps_scores) <- NULL



