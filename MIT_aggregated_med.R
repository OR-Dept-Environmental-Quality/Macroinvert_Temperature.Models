# Author: Shannon Hubler
# Initial Date: 11.10.21

# Project: Macroinvertebrate Inferred Temperature

# summary: Building on thermal traits work, as part of the Western OR/WA BCG workgroup.
#			  Develop Weighted Averaging models to infer seasonal 7-d maximum temperatures,
#			  derived from NorWeSTs MWMT metric.

#				AGGREGATED taxonomic resolution models


# packages used
library(tidyverse)
library(rioja)

 
											#####
											#####
											#			bring in taxonomy lookup table -- with OTUs for different models
											#####
											#####
											
											taxa <- read.csv('Taxa_translator_FINAL.csv')
											
											# OTU = taxa translator to bug data files  -- use this for aggregation
											# OTU_high.res = highest taxonomic resolution (genus/species as much as possible)
											# OTU_med.res = medium taxa resolution--compromises for groups where majority of info is at the next level up from 'high.res'
											# OTU_low.res = lowest taxa resolution--many groups moved to tribe/subfamily/family/order (mostly midges, oligochaetes, etc.)
											
											# '666' = used to filter out and drop taxa not included in model,
											
											# modify to make easier to join/link with bug data files
											taxa.otu <- taxa %>%
												select(TaxaID_v2, OTU) 
											
											@@@@ dont need this because Sean included it in the agg data files?
												
												
#####
#####

#		bring in bug data

#####
#####
# calibration AGG data set, used to build the models (n = 3875)
cal <- read.csv('FINAL_AGG_DATA/CAL_3875_AGG_W_META.csv') 
cal <- cal %>%
	rename(sample.id = UniqueID_v2)

# independent validation AGG data set, used to get better error  estimates of model performance (n = 629)
val <- read.csv('FINAL_AGG_DATA/VAL_629_AGG_W_META.csv')  
val <- val %>%
	rename(sample.id = UniqueID_v2)

#####
#
#		CALIBRATION
#
#####
## calculate total abundance from raw data

tot.abund<-aggregate(cal$Count,  list(sample.id=cal$sample.id), sum)  
colnames(tot.abund)[colnames(tot.abund)=="x"] <- "tot.abund"


# join back into CAL
cal.agg <- cal %>%
	left_join(tot.abund, by = c('sample.id'))

# calculate Relative Abundance
rel.abund<-cal.agg[, c("sample.id", "OTU_med.res", "Count")]
rel.abund <- rel.abund %>% 
	left_join(tot.abund, by = c('sample.id')) 


# calculate RA
rel.abund<-plyr::ddply(.data = rel.abund, c('sample.id', 'OTU_med.res'), 
							  plyr::summarize, rel.abund=Count/tot.abund)

# sum RA's across all OTUs--should see a reduction in rows
rel.abund<-plyr::ddply(.data = rel.abund, c('sample.id', 'OTU_med.res'), 
							  plyr::summarize, rel.abund=sum(rel.abund))

# check for duplicates, across sample.id and OTU
rel.abund[duplicated(rel.abund[,1:2]),]


# filter for MED RES OTUs and remove unwanted taxa


cal.agg_RA <- rel.abund %>%
	select(sample.id, OTU_med.res, rel.abund) %>%
	filter(OTU_med.res != 666) 

# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

# high res
bug.cal.agg_med_wide <- cal.agg_RA %>% 
	pivot_wider(id_cols = sample.id, names_from = OTU_med.res, values_from = rel.abund,
					values_fn = sum) 
bug.cal.agg_med_wide[is.na(bug.cal.agg_med_wide)] <- 0
bug.cal.agg_med_wide <-	column_to_rownames(bug.cal.agg_med_wide, 'sample.id')




####
# repeat for VAL
####



@@@@@@
	@@@@@@@@@
	@@@@@@@@@@@@
	@@@@@@@@@@@@@@@@@@
	
	
	
	#####
#####

#     create ENV data frame (MWMT)

#####
#####

env <- read.csv('CAL_3875_JS.csv')

env.cal <- env %>%
	select(UniqueID_v2, MWMT_final) %>%
	rename(sample.id = UniqueID_v2) %>%
	distinct_all()


env.cal <-	column_to_rownames(env.cal, 'sample.id')


####
# repeat for VAL
####










#####
#####

#		Create inference models with CAL data

#####
#####


# high res
spec <- bug.cal.agg_med_wide
env <- env.cal

wa_med.res <- WA(y=spec, x=env, mono=TRUE, tolDW = TRUE, use.N2=TRUE, tol.cut=.01, 
						check.data=TRUE, lean=FALSE)


wa_med.res # 280 taxa, RMSE (inv/cla) = 2.3/2.7, r2 = 0.682 (both), max bias 6.2/3.2

crossval(wa_med.res, cv.method="loo", verbose=TRUE, ngroups=10,
			nboot=100, h.cutoff=0, h.dist=NULL)
# crossval results almost exactly the same as original

