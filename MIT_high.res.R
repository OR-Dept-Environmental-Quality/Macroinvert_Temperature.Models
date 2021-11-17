# Author: Shannon Hubler
# Initial Date: 11.10.21

# Project: Macroinvertebrate Inferred Temperature

# summary: Building on thermal traits work, as part of the Western OR/WA BCG workgroup.
#			  Develop Weighted Averaging models to infer seasonal 7-d maximum temperatures,
#			  derived from NorWeSTs MWMT metric.

#				HIGH taxonomic resolution models


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

cal <- read.csv('CAL_3875_JS.csv') # calibration data set, used to build the models (n = 3875)
val <- read.csv('VAL_629_JS.csv')  # independent validation data set, used to get better error 
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

hist(cal.lim$RA)
qqnorm(cal.lim$RA ,main="QQ plot of RA data",pch=19)
qqline(cal.lim$RA)

hist(asin(sqrt(cal.lim$RA)))
qqnorm(asin(sqrt(cal.lim$RA)) ,main="QQ plot of RA data",pch=19)
qqline(asin(sqrt(cal.lim$RA)))



bug.cal_high.res <- cal.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) %>%
	select(site.id, OTU_high.res, RA) %>%
	filter(OTU_high.res != 666) 


bug.cal_high.res$RA.trans <- asin(sqrt(bug.cal_high.res$RA))

# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

# high res
bug.cal_high_wide <- bug.cal_high.res %>% 
	pivot_wider(id_cols = site.id, names_from = OTU_high.res, values_from = RA,
					values_fn = sum) 
bug.cal_high_wide[is.na(bug.cal_high_wide)] <- 0
bug.cal_high_wide <-	column_to_rownames(bug.cal_high_wide, 'site.id')


	# create a wide dataset for RA.trans
	bug.cal_high_wide_RA.trans <- bug.cal_high.res %>% 
		pivot_wider(id_cols = site.id, names_from = OTU_high.res, values_from = RA.trans,
						values_fn = sum) 
	bug.cal_high_wide_RA.trans[is.na(bug.cal_high_wide_RA.trans)] <- 0
	bug.cal_high_wide_RA.trans <-	column_to_rownames(bug.cal_high_wide_RA.trans, 'site.id')




####
# repeat for VAL
####
val.lim <- val %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(site.id = UniqueID_v2)

bug.val_high.res <- val.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) %>%
	select(site.id, OTU_high.res, RA) %>%
	filter(OTU_high.res != 666) 

bug.val_high.res$RA.trans <- asin(sqrt(bug.val_high.res$RA))


# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

# high res
bug.val_high_wide <- bug.val_high.res %>% 
	pivot_wider(id_cols = site.id, names_from = OTU_high.res, values_from = RA,
					values_fn = sum) 
bug.val_high_wide[is.na(bug.val_high_wide)] <- 0
bug.val_high_wide <-	column_to_rownames(bug.val_high_wide, 'site.id')

			# create a data set for transformed RA
			bug.val_high_wide_RA.trans <- bug.val_high.res %>% 
				pivot_wider(id_cols = site.id, names_from = OTU_high.res, values_from = RA.trans,
								values_fn = sum) 
			bug.val_high_wide_RA.trans[is.na(bug.val_high_wide_RA.trans)] <- 0
			bug.val_high_wide_RA.trans <-	column_to_rownames(bug.val_high_wide_RA.trans, 'site.id')


#####
#####

#     create ENV data frame (MWMT)

#####
#####

env.cal <- cal %>%
	select(UniqueID_v2, MWMT_final) %>%
	rename(site.id = UniqueID_v2) %>%
	distinct_all()


env.cal <-	column_to_rownames(env.cal, 'site.id')


hist(env.cal$MWMT_final)
qqnorm(env.cal$MWMT_final ,main="QQ plot of MWMT data",pch=19)
qqline(env.cal$MWMT_final)



#####
#####

#		Create inference models with CAL data

#####
#####


# high res
spec <- bug.cal_high_wide
env <- env.cal

wa_high.res <- WA(y=spec, x=env, mono=TRUE, tolDW = TRUE, use.N2=TRUE, tol.cut=.01, 
							check.data=TRUE, lean=FALSE)

		# look at WA-PLS
		wa.pls_high.res <- WAPLS(y=spec, x=env, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE,
		check.data=TRUE)
		# minimal to no improvements
		
		

wa_high.res # 342 taxa, RMSE (inv/cla) = 2.2/2.7, r2 = 0.687 (both), max bias 6.4/3.7

crossval(wa_high.res, cv.method="loo", verbose=TRUE, ngroups=10,
			nboot=100, h.cutoff=0, h.dist=NULL)
		# crossval results almost exactly the same as original



performance(wa_high.res)

names(wa_high.res)
wa_high.res$fitted.values
wa_high.res$coefficients

WA.resid_high <- residuals(wa_high.res)	
hist(WA.resid_high[,2])
WA.resid_high <- (WA.resid_high[,2])
WA.resid.high <- as.data.frame(WA.resid_high)
WA.resid.high <- WA.resid.high %>%
	mutate(sample.id = rownames(WA.resid.high))

# plot inferred vs MWMT
plot(wa_high.res, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="classical",
	  xlab="", ylab="", ylim=c(0,40), xlim=c(0,40), add.ref=TRUE,
	  add.smooth=TRUE)

# plot residuals -- INVERSE shows BIAS, CLASSICAL doesn't
plot(wa_high.res, resid=TRUE, xval=FALSE, tolDW=FALSE, deshrink="classical",
	  xlab="", ylab="", ylim=c(-15,15), xlim=c(0,35), add.ref=TRUE,
	  add.smooth=TRUE)

				# high res---RA.trans
				spec <- bug.cal_high_wide_RA.trans
				env <- env.cal
				
				wa_high.res_RA.trans <- WA(y=spec, x=env, mono=FALSE, tolDW = FALSE, use.N2=TRUE, tol.cut=.01, 
										check.data=TRUE, lean=FALSE)
				
				wa_high.res_RA.trans # 342 taxa, RMSE (inv/cla) = 2.2/2.7, r2 = 0.687 (both), max bias 6.4/3.7
				
				crossval(wa_high.res_RA.trans, cv.method="loo", verbose=TRUE, ngroups=10,
							nboot=100, h.cutoff=0, h.dist=NULL)
				# crossval results almost exactly the same as original
				
				
				
				performance(wa_high.res_RA.trans) # RMSE 2.1/2.5, r2 0.73/0.73, max bias 6.4/4.0
				
				names(wa_high.res_RA.trans)
				wa_high.res$fitted.values_RA.trans
				wa_high.res$coefficients_RA.trans
				wa_high.res$cv.summary_RA.trans
				
				
				predict(wa_high.res_RA.trans, newdata=bug.val_high_wide_RA.trans, sse=FALSE, nboot=100,
						  match.data=TRUE, verbose=TRUE)
				
				# plot inferred vs MWMT
				plot(wa_high.res_RA.trans, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="classical",
					  xlab="MWMT", ylab="Inferred temp (MWMT)", ylim=c(0,40), xlim=c(0,40), add.ref=TRUE,
					  add.smooth=TRUE)
				
				# plot residuals -- INVERSE shows BIAS, CLASSICAL doesn't
				plot(wa_high.res_RA.trans, resid=TRUE, xval=FALSE, tolDW=FALSE, deshrink="classical",
					  xlab="MWMT", ylab="residuals", ylim=c(-15,15), xlim=c(0,35), add.ref=TRUE,
					  add.smooth=TRUE)

# make predictions for VAL dataset				
				
wa_high.res_VAL <- predict(wa_high.res, newdata=bug.val_high_wide, sse=TRUE, nboot=100,
									match.data=TRUE, verbose=TRUE)

# get RMSE for VAL data
wa_high.res_VAL$v2.boot # Inv = 2.285071, Cla = 2.717451 (1000 boots)
								# [1] 2.282864 2.715980 (10,000 boots)
wa_high.res_VAL$fit
			@@@@ whats the difference between fit and fit.boot?
	

				
				
# join MWMT and Inferred values, plus site specific data

site.data <- read.csv('site.data_PrimVal_JS.csv')
site.data <- site.data %>%
				rename('sample.id' = 'UniqueID_v2')


site.data_residuals <- site.data %>%
							 inner_join(WA.resid.high, by = c('sample.id'))


site.data_residuals$eco3 <- plyr::revalue(site.data_residuals$eco3, c("Blue Mountains" = "B.Mts",                                   
			"Cascades" = "Casc", "Coast Range" = "Co.Ra", "Columbia Plateau" = "Col.Pl",
			"Eastern Cascades Slopes and Foothills" = "E.Casc",
			"Klamath Mountains/California High North Coast Range" = "Kl.Mts", 
			"North Cascades" = "N.Casc", "Northern Basin and Range" = "N.Ba.Ra",
			"Northern Rockies" = "N.Rock", "Puget Lowland" = "Pu.Low",
			"Snake River Plain" = "Sna.Riv.Pl", "Willamette Valley" = "Wil.Val"))     
	
# residuals ~ longitude
p <- ggplot(data=site.data_residuals, aes(x=long, y=WA.resid_high))
p + geom_point()+ geom_smooth(method="lm")
								 
								 
# residuals ~ latitude
p <- ggplot(data=site.data_residuals, aes(x=lat, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")						 
								 
# residuals ~ slope
p <- ggplot(data=site.data_residuals, aes(x=slope_nhd, y=WA.resid_high))
p + geom_point()	+ geom_smooth(method="lm")										 
								 
								 
# residuals ~ watershed area
p <- ggplot(data=site.data_residuals, aes(x=wsarea_km2, y=WA.resid_high))
p + geom_point()	+ geom_smooth(method="lm")	+ xlim(c(0, 10000))									 
			
	@@@@ some very large watrersheds....????
	

# residuals ~ elevation
p <- ggplot(data=site.data_residuals, aes(x=elev_m, y=WA.resid_high))
p + geom_point()	+ geom_smooth(method="lm")				
								 
								 
# residuals ~ IWI
p <- ggplot(data=site.data_residuals, aes(x=IWI, y=WA.resid_high))
p + geom_point()					



# residuals ~ eocregions
p <- ggplot(data=site.data_residuals, aes(x=eco3, y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 0, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 
								 
								 
							 
								 
								 
