# Author: Shannon Hubler
# Initial Date: 11.10.21

# Project: Macroinvertebrate Inferred Temperature

# summary: Building on thermal traits work, as part of the Western OR/WA BCG workgroup.
#			  Develop Weighted Averaging models to infer seasonal 7-d maximum temperatures,
#			  derived from NorWeSTs MWMT metric.

#				LOW taxonomic resolution models

 
# packages used
library(tidyverse)
library(rioja)


#####
#####
#			bring in taxonomy lookup table -- with OTUs for different models
#####
#####

taxa <- read.csv('OTU_final_894.csv')

# OTU = taxa translator to bug data files
# OTU_high.res = highest taxonomic resolution (genus/species as much as possible)
# OTU_med.res = medium taxa resolution--compromises for groups where majority of info is at the next level up from 'high.res'
# OTU_low.res = lowest taxa resolution--many groups moved to tribe/subfamily/family/order (mostly midges, oligochaetes, etc.)

# '666' = used to filter out and drop taxa not included in model,

# modify to make easier to join/link with bug data files
taxa.otu <- taxa %>%
	select(OTU, OTU_high.res, OTU_med.res, OTU_low.res) %>%
	rename(TaxaID_v2 = OTU) # rename OTU to match taxa ID column in bug data files




#####
#####

#		bring in bug data

#####
#####

cal <- read.csv('CAL_3875.csv') # calibration data set, used to build the models (n = 3875)
val <- read.csv('VAL_629.csv')  # independent validation data set, used to get better error 
# estimates of model performance (n = 629)


# 'TaxaID_v2' = taxa ID, link this to OTU level for the appropriate model


#####
#####

#		join bugs and OTUs, create multiple bug files for each OTU resolution level

#####
#####

cal.lim <- cal %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(site.id = UniqueID_v2)


bug.cal_low.res <- cal.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) %>%
	select(site.id, OTU_low.res, RA) %>%
	filter(OTU_low.res != 666) 

# need to crosstab the bug data (turn into a wide format) so that OTUs are columns


# low res
bug.cal_low_wide <- bug.cal_low.res %>% 
	pivot_wider(id_cols = site.id, names_from = OTU_low.res, values_from = RA,
					values_fn = sum) 
bug.cal_low_wide[is.na(bug.cal_low_wide)] <- 0

####
# repeat for VAL
####
val.lim <- val %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(site.id = UniqueID_v2)


bug.val_low.res <- val.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) %>%
	select(site.id, OTU_low.res, RA) %>%
	filter(OTU_low.res != 666) 

# need to crosstab the bug data (turn into a wide format) so that OTUs are columns


# low res
bug.val_low_wide <- bug.val_low.res %>% 
	pivot_wider(id_cols = site.id, names_from = OTU_low.res, values_from = RA,
					values_fn = sum) 
bug.val_low_wide[is.na(bug.val_low_wide)] <- 0



#####
#####

#     create ENV data frame (MWMT)

#####
#####

env.cal <- cal %>%
	select(UniqueID_v2, MWMT_final) %>%
	rename(site.id = UniqueID_v2) %>%
	distinct_all()




#####
#####

#		Create inference models with CAL data

#####
#####





####
# low res
####
spec <- bug.cal_low_wide[,-1]
env <- env.cal[,-1]

wa_low.res <- WA(y=spec, x=env, mono=FALSE, tolDW = FALSE, use.N2=TRUE, tol.cut=.01, 
					  check.data=TRUE, lean=FALSE)

wa_low.res # 181 taxa , RMSE 2.4/3.0, r2 = 0.638, max bias = 7.6/4.9

crossval.low <- crossval(wa_low.res, cv.method="loo", verbose=TRUE, ngroups=10,
								 nboot=100, h.cutoff=0, h.dist=NULL)
# crossval results almost exactly the same as original

performance(wa_low.res)

names(wa_low.res)
wa_low.res$fitted.values

predict(wa_low.res, newdata=bug.val_low_wide, sse=FALSE, nboot=100,
		  match.data=TRUE, verbose=TRUE)

# plot inferred vs MWMT
plot(wa_low.res, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="classical",
	  xlab="", ylab="", ylim=c(0,35), xlim=c(0,35), add.ref=TRUE,
	  add.smooth=FALSE)

# plot residuals -- INVERSE shows BIAS, CLASSICAL doesn't
plot(wa_low.res, resid=TRUE, xval=FALSE, tolDW=FALSE, deshrink="classical",
	  xlab="", ylab="", ylim=c(-15,15), xlim=c(0,35), add.ref=TRUE,
	  add.smooth=FALSE)








