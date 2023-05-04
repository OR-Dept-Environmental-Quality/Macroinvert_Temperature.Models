# Author: Shannon Hubler
# Initial Date: 4.30.2023
# updates:
					#  changes from BSTI v 1.0:
					#	1) no longer using C2 software
					#	2) added in WA ECY fines data, using same sites as BCG Thermal work
					

# Project: BSTI_v 2.0 (Biological Sediment Tolerance Index)

# summary: Building on thermal traits work, as part of the Western OR/WA BCG workgroup.
#			  Develop Weighted Averaging models to represent fine sediment niche for macro assemblage
#			  

#				

 
# packages used
library(tidyverse)
library(rioja)
library(RODBC) 








#####
#####
#				SITE level data, for each SAMPLE
#####
#####

site.data <- read.csv('input files/_SiteData_7977_20220213.csv')%>%
	rename(sample.id = UniqueID_v2)

site.data <- site.data[order(site.data$sample.id),]	


# bring in CAL/VAL designations, based on Mark's stratified random sampling

calval <- read.csv('FINES/FN_Cal.Val.csv')%>%
	rename(sample.id = UniqueID_v2) %>%
	select(sample.id, Cal.Val, PCT_FN)

# join to site data
site.data_calval <- calval %>%
	left_join(site.data, by = "sample.id") %>%
	select(-old_calval)
	

# FINES transformations: 

hist(log(site.data_calval$PCT_FN/(1-site.data_calval$PCT_FN)))

fn.log <- (log(site.data_calval$PCT_FN+1))
fn <- ((site.data_calval$PCT_FN))
fn.sqrt <- sqrt(site.data_calval$PCT_FN)

library(ggplot2)
ggplot(site.data_calval, aes(sample=sqrt(PCT_FN)))+stat_qq()+theme_bw()
ggplot(site.data_calval, aes(sample=log10(PCT_FN+1)))+stat_qq()+theme_bw()
ggplot(site.data_calval, aes(sample=asin(sqrt(PCT_FN/100))))+stat_qq()+theme_bw()
ggplot(log(site.data_calval, aes(sample=PCT_FN/(1-PCT_FN)))) + stat_qq()+theme_bw()


site.data_calval$pct.fn_log <- log10(site.data_calval$PCT_FN+1)


#####
#####
#			bring in taxonomy lookup table 
#####
#####

# install "BiomonTools" package
					
					#install if needed
					if(!require(remotes)){install.packages("remotes")}  
					install_github("leppott/BioMonTools", force=TRUE)

				# Download first then read.
				## URL BioMonTools
				url_bmt_base <- 'https://github.com/leppott/BioMonTools_SupportFiles/raw/main/data'
				url_taxa_official_pick <- file.path(url_bmt_base, "taxa_official", "ORWA_TaxaTranslator_20230228.csv")
				taxonomy <- read.csv(url_taxa_official_pick)
							# old csv copy locally: taxa <- read.csv('ORWA_TaxaTranslator_20230112_SLH.updated.csv') \
				

				# 'DNI' = used to filter out and drop taxa not included in model,
				
				
taxonomy <- read.csv('FINES/OTU_bsti_v2.csv')				
				
									

# limit to necessary fields to void messy joins
taxa.otu <- taxonomy %>%
	select(Taxon_orig, OTU_BSTI.v2) #%>%
	#rename(TaxaID_v2 = Bio_ALL_TaxaID_v2) %>% # rename OTU to match taxa ID column in bug data files
	#rename(OTU_high.res = MTTI_OTU_highres_v2) # renames MTTI_OTU.... to simplify



#####
#####

#		bring in bug data

#####
#####

bugs_7977 <- read.csv('input files/_Bio_ALL_7977.csv')  # non-wadeable sites removed

bugs_7977 <- bugs_7977 %>%
	rename(Taxon_orig = TaxaID_v2) %>%
	rename(sample.id = UniqueID_v2) %>%
	select(-old_calval)

# limit to cal/val data
bugs_cal.val <- calval %>%
	left_join(bugs_7977, by = 'sample.id', multiple = "all")


#####
#####

#		join bugs and OTUs, filter out 'DNI' taxa, sum across OTUs within a sample

#####
#####

# join
bugs_otu <- bugs_cal.val %>%
	left_join(taxa.otu, by='Taxon_orig') %>% # join dataframes
	filter(OTU_BSTI.v2 != 'DNI')						# filter out DNI taxa


# sum RA's across all OTUs--should see a reduction in rows.  
# Also limits to the following: dataset (CAL/VAL/not), sample, OTU, (summed) RA


bugs_otu_sum.RA<-plyr::ddply(.data = bugs_otu, c('Cal.Val', 'sample.id', 'OTU_BSTI.v2'), 
										plyr::summarize, RA=sum(RA), Count=sum(Count))



# create a dataframe of total counts and relative abundances
# important! = after applying OTUs 
abunds <- bugs_otu_sum.RA %>% 			
	group_by(sample.id, Cal.Val) %>% 
	summarize(tot.abund = sum(Count),
				 tot.RA = sum(RA))	
							

#####
#####

#	Prepare data sets for modeling
#	need to crosstab the bug data (turn into a wide format) so that OTUs are columns
# then split into separate CAl and VAl datasets 

bugs_cross <- bugs_otu_sum.RA %>% 
	pivot_wider(id_cols = c(Cal.Val, sample.id), names_from = OTU_BSTI.v2, values_from = RA,
					values_fn = sum) 

bugs_cross[is.na(bugs_cross)] <- 0 # change all NA values to zeros (no presence in the sample)
								# # Exploratory
								# 	# limit samples by total abundances
								# 	bugs_cross.300 <- bugs_cross %>%
								# 		filter(sample.id %in% abunds.300$sample.id)
								# 
								# 	bugs_cross.150 <- bugs_cross %>%
								# 		filter(sample.id %in% abunds.150$sample.id)
								# 
								# 	bugs_cross.500 <- bugs_cross %>%
								# 		filter(sample.id %in% abunds.500$sample.id)
			


####
####

#		Split into CAL and VAL datasets
bug.cal <- bugs_cross %>%
	filter(Cal.Val == 'CAL') %>% # calibration data set, used to build the models (starting 'n' = 3658)
	select(!(Cal.Val))



bug.val <- bugs_cross %>%
	filter(Cal.Val == 'VAL') %>% # calibration data set, used to build the models (starting 'n' = 3658)
	select(!(Cal.Val))




			

#####
#####

#     create ENV data frame (MWMT)

#####
#####

# cal  


env.cal <- site.data_calval %>%
	filter(Cal.Val == 'CAL') %>%
	select(sample.id, pct.fn_log) %>%		
	distinct_all()

env.cal <- env.cal %>%
	filter(sample.id %in% bug.cal$sample.id)


# move sample ID column to rownames--necessary for modeling
bug.cal <-	column_to_rownames(bug.cal, 'sample.id') 
env.cal <-	column_to_rownames(env.cal, 'sample.id')		

					

# val --------> based on ANNUAL MWMT
env.val <- site.data_calval %>%
	filter(Cal.Val == 'VAL') %>%
	select(sample.id, pct.fn_log) %>%    
	distinct_all()

env.val <- env.val %>%
	filter(sample.id %in% bug.val$sample.id)

bug.val <-	column_to_rownames(bug.val, 'sample.id') 
env.val <-	column_to_rownames(env.val, 'sample.id')


							

#####
#####

#		Create inference models with CAL data

#####
#####


spec <- bug.cal

env <- env.cal		
					

wa_BSTI.v2_apr23 <- WA(y=spec, x=env, mono=TRUE, tolDW = TRUE, use.N2=TRUE, tol.cut=.01, 
							check.data=TRUE, lean=FALSE)

		# look at WA-PLS
		#wa.pls_high.res <- WAPLS(y=spec, x=env, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE,
		#check.data=TRUE)
		# minimal to no improvements
		
		
# WA model results
wa_BSTI.v2_apr23 
				

crossval(wa_BSTI.v2_apr23, cv.method="bootstrap", verbose=TRUE, ngroups=10,
			nboot=1000, h.cutoff=0, h.dist=NULL)
		# crossval results almost exactly the same as original





							# names(wa_high.res)
							# wa_high.res$fitted.values
							# wa_high.res$coefficients

WA.resid <- residuals(wa_BSTI.v2_apr23)	
							# hist(WA.resid_high[,2])
WA.resid_Cla.tol <- (WA.resid[,5])
WA.resid_Cla.tol <- as.data.frame(WA.resid_Cla.tol)
				



##
	### Save the WA model to be used later for getting MTTI for new samples
	save(wa_MTTI.mar23, file=('wa_MTTI.mar23.Rdata'))
##

				

# ######
				
# make predictions for VAL dataset				
				
######
				
	wa_BSTI.v2_apr23_VAL <- predict(wa_BSTI.v2_apr23, newdata=bug.val, sse=TRUE, nboot=100,
									match.data=TRUE, verbose=TRUE)

#### get RMSE and r2 for VAL data


# extract VAL predictions as a standalone dataframe, so can join with env.val and calculate RMSEP
VAL.predicts <- as.data.frame(wa_BSTI.v2_apr23_VAL$fit)
VAL.predicts <- cbind(VAL.predicts, env.val)




	
VAL.rmsep <- VAL.predicts %>%
	mutate(VAL_WA.inv_rmsep = sqrt(mean((VAL.predicts$WA.inv - VAL.predicts$pct.fn_log)^2))) %>%
	mutate(VAL_WA.cla_rmsep = sqrt(mean((VAL.predicts$WA.cla - VAL.predicts$pct.fn_log)^2)))	%>%	 
	mutate(VAL_WA.m_rmsep = sqrt(mean((VAL.predicts$WA.m - VAL.predicts$pct.fn_log)^2)))		 %>%
	mutate(VAL_WA.inv.tol_rmsep = sqrt(mean((VAL.predicts$WA.inv.tol - VAL.predicts$pct.fn_log)^2)))	%>%	 
	mutate(VAL_WA.cla.tol_rmsep = sqrt(mean((VAL.predicts$WA.cla.tol - VAL.predicts$pct.fn_log)^2)))	%>%	 
	mutate(VAL_WA.m.tol_rmsep = sqrt(mean((VAL.predicts$WA.m.tol - VAL.predicts$pct.fn_log)^2)))	%>%
	mutate(VAL.WA.inv_rsq = cor(VAL.predicts$WA.inv, VAL.predicts$pct.fn_log)^2) 	%>%
	mutate(VAL.WA.cla_rsq = cor(VAL.predicts$WA.cla, VAL.predicts$pct.fn_log)^2)		%>%
	mutate(VAL.WA.m_rsq = cor(VAL.predicts$WA.m, VAL.predicts$pct.fn_log)^2)			%>%
	mutate(VAL.WA.inv.tol_rsq = cor(VAL.predicts$WA.inv.tol, VAL.predicts$pct.fn_log)^2)	%>%
	mutate(VAL.WA.cla.tol_rsq = cor(VAL.predicts$WA.cla.tol, VAL.predicts$pct.fn_log)^2)	%>%
	mutate(VAL.WA.m.tol_rsq = cor(VAL.predicts$WA.m.tol, VAL.predicts$pct.fn_log)^2)

VAL.rmsep <- VAL.rmsep[,8:19]
VAL.rmsep <- distinct(VAL.rmsep)
row.names(VAL.rmsep) <- NULL

# calculate Max Bias for VAL
				# From ter Braak and Juggins 1993: For estimation of the maximum bias, the sampling interval (0, 100) was subdivided 
				# into 10 equal intervals, the bias per interval calculated and the (signed) maximum of the 10 values calculated. 

VAL.maxbias <- VAL.predicts %>%
	select(WA.cla.tol, pct.fn_log) %>%
	mutate(val.resid = WA.cla.tol - pct.fn_log)

val.interval = (max(VAL.maxbias$MWMT_annual) - min(VAL.maxbias$MWMT_annual))/10 # range is 5.9 - 29.6....interval is every 2.4 degrees

#average bias per interval

	
# create a new column to assign interval grouping
VAL.maxbias$interval <- ifelse(VAL.maxbias$pct.fn_log <= 8.3, 'int.1',   
		ifelse(VAL.maxbias$pct.fn_log > 8.3 & VAL.maxbias$pct.fn_log <= 10.7, 'int.2',  
			ifelse(VAL.maxbias$pct.fn_log > 10.7 & VAL.maxbias$pct.fn_log <= 13.1, 'int.3',
				ifelse(VAL.maxbias$pct.fn_log > 13.1 & VAL.maxbias$pct.fn_log <= 15.5, 'int.4',
					ifelse(VAL.maxbias$pct.fn_log > 15.5 & VAL.maxbias$pct.fn_log <= 17.9, 'int.5',
						ifelse(VAL.maxbias$pct.fn_log > 17.9 & VAL.maxbias$pct.fn_log <= 20.3, 'int.6',
							ifelse(VAL.maxbias$pct.fn_log > 20.3 & VAL.maxbias$pct.fn_log <= 22.7, 'int.7',
								ifelse(VAL.maxbias$pct.fn_log > 22.7 & VAL.maxbias$pct.fn_log <= 25.1, 'int.8', 
									 ifelse(VAL.maxbias$pct.fn_log > 25.1 & VAL.maxbias$pct.fn_log <= 27.5, 'int.9',
										'int.10'))))))))) 
									 
VAL.ave.bias <- (VAL.maxbias %>%  # as.data.table
	group_by(interval) %>%
	summarise(VAL.ave.bias = mean(val.resid)))
###

#				Get predicted MTTI for CAL and VAL

###

cal_wa.cla.tol <- as.data.frame(wa_MTTI.mar23$fitted.values)
cal_wa.cla.tol <- as.data.frame(wa_MTTI.mar23$fitted.values[,5])
cal_wa.cla.tol <- cal_wa.cla.tol %>% rename(MTTI = `wa_MTTI.mar23$fitted.values[, 5]`)
cal_wa.cla.tol$Model <- "CAL"

val_wa.cla.tol <- VAL.predicts[,-c(1:4,6)]
val_wa.cla.tol <- VAL.predicts %>%
	select(WA.cla.tol) %>%
	rename(MTTI = WA.cla.tol)
val_wa.cla.tol$Model <- "VAL"

MTTI_cal.val <- rbind(cal_wa.cla.tol, val_wa.cla.tol)

MTTI_cal.val$sample.id <- row.names(MTTI_cal.val)

MTTI_cal.val <- MTTI_cal.val %>%
	left_join(site.data_calval, by=c('sample.id'))



					# check for duplicates
					# n_occur <- data.frame(table(site.data_calval$UniqueID_v2))
					
					# n_occur[n_occur$Freq > 1,]
					
					# MTTI_cal.val[MTTI_cal.val$sample.id %in% n_occur$Var1[n_occur$Freq > 1],]




write.csv(MTTI_cal.val, 'MTTI_cal.val_mar23.csv')

###########
##########

# Plot CAL and VAL together


#########
##########

# plot inferred vs MWMT
par(mar=c(4.3,4.3,1,1))
plot(wa_MTTI.mar23, resid=FALSE, xval=FALSE, tolDW=TRUE, deshrink="classical",
	  xlab="MWMT", ylab="MTTI", main="", 
	  ylim=c(-5,40), xlim=c(0,40), add.ref=TRUE)#, add.smooth=TRUE

par(new=TRUE)
plot(x=VAL.predicts$MWMT_annual, y=VAL.predicts$WA.cla.tol,
	  xlab="", ylab="", main="", pch=2, col='red',
	  ylim=c(-5,40), xlim=c(0,40), axes=FALSE)
abline(a=0, b=1)

# plot residuals -- INVERSE shows BIAS, CLASSICAL doesn't
par(mar=c(4.3,4.3,1,1))
plot(wa_MTTI.mar23, resid=TRUE, xval=FALSE, tolDW=TRUE, deshrink="classical",
	  xlab="MWMT", ylab="residuals (MTTI)", ylim=c(-15,15), xlim=c(0,35), add.ref=TRUE,
	   main='')#add.smooth=TRUE,


par(new=TRUE)	
plot(x=VAL.predicts$MWMT_annual, y=(VAL.predicts$MWMT_annual - VAL.predicts$WA.cla.tol), 
	  pch=2,  xlab="", col='red',ylab="", main="", ylim=c(-15,15), xlim=c(0,35),
	  axes=FALSE)
abline(h=0)


# histogram MWMT val and cal



library(ggplot2)


ggplot(MTTI_cal.val, aes(x=MWMT_TolAnal, color=dataset, fill=dataset)) +
	geom_histogram(position="identity", alpha=0.8)+
	scale_fill_manual(values=c("light blue", "black")) +
	theme(legend.position="none") +
	labs(,x="MWMT", y = "Count")


######


#			Look for potential BIAS in predictions - natural gradients, disturbance, etc.

#####

# join MWMT and Inferred values, plus site specific data

eco.cnt <- site.data_calval %>%
	group_by(dataset, eco3) %>%
	summarize(count = n()) %>%
	pivot_wider(names_from = dataset, values_from = count)

write.csv(eco.cnt, 'ecoregion_count.csv')

# numerical summaries -- need to remove categorical vars
site.data_num <- site.data_calval %>%
	# rename(sample.id = UniqueID_v2) %>%
	select(sample.id, siteid, dataset, lat, long, daynum, Year, MWMT_TolAnal,  slope_nhd, wsarea_km2, 
			 elev_m, IWI, PctUrbWs, PctAgWs, PctForWs, RdDensWs, Tmean8110Ws, Precip8110Ws)

site.data_num <- site.data_num %>%
	rename(MWMT_Annual = MWMT_TolAnal) #%>%
	# rename(slope = slope_nhd) %>%
	# rename(area = wsarea_km2) %>%
	# rename(elevation = elev_m) %>%
	# rename(temperature = Tmean8110Ws) %>%
	# rename(precipitation = Precip8110Ws)

site.data_num.long <- site.data_num %>%
	pivot_longer(!c(sample.id, siteid, dataset)  , names_to = "variable", values_to = "value")


sum.stats <- site.data_num.long %>%                               # Summary by group using dplyr
					group_by(dataset, variable) %>% 
					summarize(min = min(value, na.rm = TRUE),
								 q1 = quantile(value, 0.25, na.rm=TRUE),
								 median = median(value, na.rm=TRUE),
								 mean = mean(value, na.rm=TRUE),
								 q3 = quantile(value, 0.75, na.rm=TRUE),
								 max = max(value, na.rm=TRUE))


	
write.csv(sum.stats, 'sum.stats.csv')
	
# boxplots
	site.data_num.long %>% 
	ggplot(aes(x= dataset, y=value, fill=dataset)) +
	geom_boxplot() +
	#geom_jitter(width=0.1,alpha=0.2) +
	xlab("")+ 
	facet_wrap(~variable, scales = "free") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


	
	
# PCA

library(factoextra)

# limit dframe to numeric variables = use 'site.data_num'

								# remove columns with "zero variance"
								#which(apply(site.data_num, 2, var)==0)  
# remove incomplete cases
site.data_pca <- site.data_num[complete.cases(site.data_num), ]

# remove extreme WSarea outlier
site.data_pca <- site.data_pca %>%
	filter(sample.id != 'SampEPA_5080_1')


#@@@@@@@@@@@@@@@@@@@@@@
# data transformations
#@@@@@@@@@@@@@@@@@@@@@@								#colnames(site.data_pca)

								source("transform.variables[2_asinsqrt].r")
								transform.view2(site.data_pca)
									
								
								
								#         lat  							none            
								#         long          				none
								#         daynum        				none
								#         Year          				none
								#		    MWMT_Annual   				none
								#	       slope_nhd              	log10 + 1
								#         wsarea_km2             	log10 
								#         elev_m              		sqrt
								#         IWI           				none
								#         PctUrbWs            		log10 + 1 
								#         PctAgWs             		log10 + 1
								#         PctForWs                  asin(sqrt)
								#         RdDensWs            		log10 + 1
								#         Tmean8110Ws         		none
								#         Precip8110Ws        		sqrt

site.data_pca.trans <- site.data_pca %>%
	mutate(slope.log = log10(slope_nhd + 1)) %>%
	mutate(wsarea.log = log10(wsarea_km2))   %>%
	mutate(elev.sqrt = sqrt(elev_m)) %>%
	mutate(pct.urban_log = log10(PctUrbWs + 1)) %>%
	mutate(pct.ag_log = log10(PctAgWs + 1))   %>%
	mutate(pct.forest_asin = asin(sqrt(PctForWs/100)))  %>%
	mutate(road.den_log = log10(RdDensWs + 1)) %>%
	mutate(Precip_sqrt = sqrt(Precip8110Ws))
	
	
site.data_pca.trans <- site.data_pca.trans %>%
	select(sample.id, siteid, dataset, lat, long, Year, daynum, MWMT_Annual, slope.log, 
			 wsarea.log, elev.sqrt, IWI, pct.urban_log, pct.ag_log, pct.forest_asin, 
			 road.den_log, Tmean8110Ws, Precip_sqrt)

# remove ag, roads, urban, forest ---> captured by IWI, makes plot less messy, remove lat/long = not true env gradients			
site.data_pca.trans_LIM <- site.data_pca.trans %>%
	select(sample.id, siteid, dataset, Year, daynum, MWMT_Annual, slope.log, 
			  wsarea.log, elev.sqrt, IWI, Tmean8110Ws, Precip_sqrt) %>%
	rename(MWMT = MWMT_Annual, slope = slope.log,
			 elevation = elev.sqrt, precipitation = Precip_sqrt, year = Year,
			 area = wsarea.log, "day of year" = daynum, "air temp" = Tmean8110Ws)
	

site.data_pca.trans_LIM.long <- site.data_pca.trans_LIM %>%
	pivot_longer(!c(sample.id, siteid, dataset)  , names_to = "variable", values_to = "value")

												# try to reorder the levels of 'variable'
												site.data_pca.trans_LIM.long$variable <- as.factor(site.data_pca.trans_LIM.long$variable)
												site.data_pca.trans_LIM.long$variable <- relevel(site.data_pca.trans_LIM.long$variable, 
													c( "MWMT","air temp","sqrt(precipitation)", "sqrt(elevation)", "IWI", "Year", "day of year"))
												
												
												site.data_pca.trans_LIM.long <- site.data_pca.trans_LIM.long %>%
													mutate(variable = reorder("MWMT","air temp","sqrt(precipitation)", 
																					  "sqrt(elevation)", "IWI", "Year", "day of year"))
												
											

site.data_pca.trans_LIM.long %>% 
	ggplot(aes(x= dataset, y=value, fill=dataset)) +
	geom_boxplot() +
	#geom_jitter(width=0.1,alpha=0.2) +
	xlab("")+ ylab("")+
	facet_wrap(~variable, scales = "free") +
	theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
	theme(legend.position = 'none')

# histograms for CAL and VAL + MWMT	
p <- ggplot(site.data_pca.trans_LIM, aes(x=MWMT, fill = dataset, color = dataset)) +
	geom_histogram( alpha = 0.7)	+ theme(legend.position='none') 

p +	scale_color_manual(values=c("blue", "red"))+
	scale_fill_manual(values=c("light blue", "orange"))



	
								#rownames(site.data_complete) <- site.data_complete[,1]
								#site.data_pca <- site.data_complete[,-c(1:3)]


# run pca 
library(FactoMineR)  
library("factoextra")
library("corrplot")

# center and scale, before PCA
site.data_pca.trans_LIM_c.s=scale(site.data_pca.trans_LIM[,c(4:12)],center=T,scale=T)

pca.allsites <- PCA(site.data_pca.trans_LIM_c.s, graph = FALSE, scale.unit = TRUE)



summary(pca.allsites)


				# str(pca.allsites)
				# print(pca.allsites)
eig.val <- get_eigenvalue(pca.allsites) # Extract the eigenvalues/variances of principal components
fviz_eig(pca.allsites, addlabels = TRUE, ylim = c(0, 50) ) # Visualize the eigenvalues

# Extract the results for individuals and variables, respectively.
ind <- get_pca_ind(pca.allsites); 
var <- get_pca_var(pca.allsites) 

# Visualize the results individuals and variables, respectively.
fviz_pca_ind(pca.allsites)
fviz_pca_var(pca.allsites, col.var = "blue") 

fviz_pca_biplot(pca.allsites) # Make a biplot of individuals and variables.


corrplot(var$cos2, is.corr=FALSE)  # ???????

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(pca.allsites, choice = "var", axes = 1:2)


# Color by cos2 values: quality on the factor map
fviz_pca_var(pca.allsites, col.var = "cos2",
				 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
				 repel = TRUE # Avoid text overlapping
)


fviz_pca_var(pca.allsites, alpha.var = "cos2") # Change the transparency by cos2 values

corrplot(var$contrib, is.corr=FALSE)  

fviz_contrib(pca.allsites, choice = "var", axes = 1, top = 15)
fviz_contrib(pca.allsites, choice = "var", axes = 2, top = 15)
fviz_contrib(pca.allsites, choice = "var", axes = 3, top = 15)
fviz_contrib(pca.allsites, choice = "var", axes = 1:3, top = 15)


# The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
fviz_pca_var(pca.allsites, col.var = "contrib",
				 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)


fviz_pca_ind(pca.allsites, col.ind = "cos2", #axes = 2:3,
				 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
				 repel = TRUE # Avoid text overlapping (slow if many points)
)


grp <- as.factor(site.data_pca.trans_LIM$dataset)
# Color variables by groups


fviz_pca_ind(pca.allsites,
				 col.ind = grp, # color by groups
				 palette = c("#00AFBB",  "#FC4E07"),
				 addEllipses = FALSE, # Concentration ellipses
				 ellipse.type = "confidence",
				 legend.title = "Groups",
				 repel = TRUE
)



fviz_pca_biplot(pca.allsites, axes=1:2,
				 col.ind = grp, # color by groups
				 palette = c("steelblue1",  "red"),#c("#00AFBB",  "#FC4E07")
				 # habillage = site.data_pca$calval,
				 # addEllipses = TRUE, # Concentration ellipses
				 # ellipse.type = "confidence",
				 legend.title = "Dataset",
				 repel = TRUE, #label='var', 
				 geom.ind = c("point"), geom.var = c("arrow", "text"),
				 col.var = 'black',
				 alpha.ind = 1, # transparency of individuals 
				 title = ""
)


fviz_pca_biplot(pca.allsites, axes=c(2,3),
					 col.ind = grp, # color by groups
					 palette = c("steelblue1",  "red"),#c("#00AFBB",  "#FC4E07")
					 # habillage = site.data_pca$calval,
					 # addEllipses = TRUE, # Concentration ellipses
					 # ellipse.type = "confidence",
					 legend.title = "Dataset",
					 repel = TRUE, #label='var', 
					 geom.ind = c("point"), geom.var = c("arrow", "text"),
					 col.var = 'black',
					 alpha.ind = 1, # transparency of individuals 
					 title=''
)

fviz_pca_biplot(pca.allsites, axes=c(3,4),
					 col.ind = grp, # color by groups
					 palette = c("steelblue1",  "red"),#c("#00AFBB",  "#FC4E07")
					 # habillage = site.data_pca$calval,
					 # addEllipses = TRUE, # Concentration ellipses
					 # ellipse.type = "confidence",
					 legend.title = "Dataset",
					 repel = TRUE, #label='var', 
					 geom.ind = c("point"), geom.var = c("arrow", "text"),
					 col.var = 'black',
					 alpha.ind = 1, # transparency of individuals 
					 title=''
)



######

# BIAS = look for patterns in CAL predictions--do we see bias associated with other variables?
		
######

WA.resid_Cla.tol_2 <- WA.resid_Cla.tol %>%
 	mutate(sample.id = rownames(WA.resid_Cla.tol)) # get site ID as a column for join




WA.resid_Cla.tol_2$sample.id <- row.names(WA.resid_Cla.tol_2)

site.data_residuals <- site.data %>%
							 rename(sample.id = sample.id)  %>%
							 inner_join(WA.resid_Cla.tol_2, by = c('sample.id'))



site.data_residuals <- site.data_residuals %>%
	left_join(abunds, by='sample.id') %>%
	select(-dataset)


				#site.data_residuals <- site.data_residuals %>%
				#								filter(tot.abund > 499)
										# site.data_residuals$eco3 <- plyr::revalue(site.data_residuals$eco3, c("Blue Mountains" = "B.Mts",                                   
										# 			"Cascades" = "Casc", "Coast Range" = "Co.Ra", "Columbia Plateau" = "Col.Pl",
										# 			"Eastern Cascades Slopes and Foothills" = "E.Casc",
										# 			"Klamath Mountains/California High North Coast Range" = "Kl.Mts", 
										# 			"North Cascades" = "N.Casc", "Northern Basin and Range" = "N.Ba.Ra",
										# 			"Northern Rockies" = "N.Rock", "Puget Lowland" = "Pu.Low",
										# 			"Snake River Plain" = "Sna.Riv.Pl", "Willamette Valley" = "Wil.Val"))     
	
# residuals ~ longitude
p <- ggplot(data=site.data_residuals, aes(x=long, y=WA.resid_Cla.tol))
p + geom_point()+ geom_smooth(method="lm")
								 
								 
# residuals ~ latitude
p <- ggplot(data=site.data_residuals, aes(x=lat, y=WA.resid_Cla.tol))
p + geom_point()		+ geom_smooth(method="lm")						 
								 
# residuals ~ slope
p <- ggplot(data=site.data_residuals, aes(x=sqrt(slope_nhd), y=WA.resid_Cla.tol))
p + geom_point()	+ geom_smooth(method="lm")										 
								 
								 
# residuals ~ watershed area
p <- ggplot(data=site.data_residuals, aes(x=log10(wsarea_km2), y=WA.resid_Cla.tol))
p + geom_point()	+ geom_smooth(method="lm")	+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 


# residuals ~ elevation
p <- ggplot(data=site.data_residuals, aes(x=elev_m, y=WA.resid_Cla.tol))
p + geom_point()	+ geom_smooth(method="lm")	+
	geom_hline(yintercept=0, linetype="dashed", color = "red")			
								 
								 
# residuals ~ IWI
p <- ggplot(data=site.data_residuals, aes(x=IWI, y=WA.resid_Cla.tol))
p + geom_point()		+ geom_smooth(method="lm")	+
	geom_hline(yintercept=0, linetype="dashed", color = "red")	


@@@@@@@@@ MTTI higher than observed MWMT at disturbed sites--makes sense?  
@@@@@@@@@ More tolerant taxa = higher WAopt 

# residuals ~ daynum
p <- ggplot(data=site.data_residuals, aes(x=daynum, y=WA.resid_Cla.tol))
p + geom_point()		+ geom_smooth(method="lm")		+
	geom_hline(yintercept=0, linetype="dashed", color = "red")			


# residuals ~ year
p <- ggplot(data=site.data_residuals, aes(x=Year, y=WA.resid_Cla.tol))
p + geom_point()		+ geom_smooth(method="lm")		+
	geom_hline(yintercept=0, linetype="dashed", color = "red")			



# residuals ~ ecocregions
p <- ggplot(data=site.data_residuals, aes(x=eco3, y=WA.resid_Cla.tol))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 0, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 
								 

# residuals ~ Collection Method
p <- ggplot(data=site.data_residuals, aes(x=as.factor(CollMeth), y=WA.resid_Cla.tol))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 0, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 


# residuals ~ Source
p <- ggplot(data=site.data_residuals, aes(x=as.factor(SourceEntity), y=WA.resid_Cla.tol))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 90, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 




# residuals ~ Lab
p <- ggplot(data=site.data_residuals, aes(x=as.factor(lab), y=WA.resid_Cla.tol))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 90, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 

# resids vs total abundance
p <- ggplot(data=site.data_residuals, aes(x=log10(tot.abund), y=WA.resid_Cla.tol))								 
p + geom_point()		+ geom_smooth(method="lm")		+
	geom_hline(yintercept=0, linetype="dashed", color = "red")			


# resids vs Relative abundance
p <- ggplot(data=site.data_residuals, aes(x=tot.RA, y=WA.resid_Cla.tol))								 
p + geom_point()		+ geom_smooth(method="lm")		+
	geom_hline(yintercept=0, linetype="dashed", color = "red")			


# RA vs tot.abund
p <- ggplot(data=site.data_residuals, aes(x=log10(tot.abund), y=tot.RA))								 
p + geom_point()





# REFERENCE/DISTURBANCE: ODEQ data only

# bring in reference designations and disturbance metrics
one_rule_all <- read.csv("//deqlab1/GIS_WA/Project_Working_Folders/Reference/2020/_Final outputs/one.table_rule.all.csv")
translator <- read.csv('C://Users/shubler/Desktop/ORWA BCG_thermal/Inference modeling/translator_ODEQ.csv')
translator <- translator %>% select(UniqueID_v2, MLocID, Date) %>% rename(sample.id = UniqueID_v2) 

or.ref <- translator %>%
	left_join(one_rule_all, by = 'MLocID') 
	

or.ref_residuals <- or.ref %>%
	inner_join(WA.resid.high_2, by = c('sample.id'))

# residuals ~ road density
p <- ggplot(data=or.ref_residuals, aes(x=rdden_km_km2, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	+ 
	geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)


# residuals ~ road crossings
p <- ggplot(data=or.ref_residuals, aes(x=xings_km2, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	+ 
	geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)

# residuals ~ % Ag
p <- ggplot(data=or.ref_residuals, aes(x=P_AgLand, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	+ 
	geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)


# residuals ~ Code 21
p <- ggplot(data=or.ref_residuals, aes(x=P_Urban21Land, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	+ 
	geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)


# residuals ~ Disturbance Score (Google Earth screens--only REF CANDIDATE sites get this)
p <- ggplot(data=or.ref_residuals, aes(x=Disturb.score, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	+ 
	geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)


# residuals ~ OR Reference designation
p <- ggplot(data=or.ref_residuals, aes(x=as.factor(Ref2020_FINAL), y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 90, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 

@@@@@
	@@@@@@@  Can we use RF to describe leading factors in residuals?  Natural vs Disturbance more important?
@@@@@
	
library(randomForest)
rf <- randomForest(WA.resid_high~ Lat_DD+Long_DD+Eco3+rdden_km_km2+xings_km2+P_AgLand+P_Urban21Land+mines+grvl_mn_km2+
			P_canal+WorE+Ref2020_FINAL, or.ref_residuals, proximity=TRUE) 


@@@@@ needs to be tied back to natural gradients

print(rf)
varImpPlot(rf,
			  sort = T,
			  n.var = 10,
			  main = "Top 10 - Variable Importance")
importance(rf)

#####

# Segment Wise RMSEP: https://quantpalaeo.wordpress.com/2014/06/10/uneven-sampling-of-the-gradient-and-segment-wise-rmsep/

			# k sets the column of the model output to use. For example, with WA, the first column 
			# is for the inverse deshrinking model, and the second for the classical deshrinking model. 
			# Run performance(mod) to see which model is in which column. 

####

#Now a function to calculate and plot segment-wise RMSEP.


segmentwise.rmse <- function(mod, ng = 10, k = 5, plot = TRUE, ...){ 
	if(is.null(mod$residuals.cv)){
		if(class(mod) == "MAT"){
			r <- mod$fitted.values[,k]-mod$x
			perf<-performance(mod)$object
		}else{
			stop("Need cross-validated model to calculate RMSEP")
		}
	}else{
		r <- mod$residuals.cv[,k]
		perf<-performance(mod)$crossval
	}  
	breaks <- seq(min(mod$x),max(mod$x), length=ng)
	envcut <- cut(mod$x,breaks=breaks, include.lowest=TRUE)
	segRMSEP <- tapply(r,envcut,function(x)sqrt(mean(x^2)))
	allsegRMSEP <- sqrt(mean(segRMSEP^2))
	
	if(plot){
		hist(mod$x, breaks=breaks, col="grey70", border=NA, ...)
		par(new=T)
		mid<-((c(breaks,NA)+c(NA,breaks))/2)
		mid<-mid[!is.na(mid)]
		plot(mid,segRMSEP, type="n", xlim=par()$usr[1:2],xaxs="i", yaxt="n", ylab="", xlab="", col=2)
		lines(breaks,c(segRMSEP[1],segRMSEP), type="S", col=3)
		axis(4)
		mtext("RMSEP", side=4, line=1.5)
		abline(h = perf[k,"RMSE"], col=2, lty=2)
		abline(h = allsegRMSEP, col=4, lty=2) 
	}
	list(breaks = breaks, segRMSEP = segRMSEP, allsegRMSEP = allsegRMSEP)
}

# Now we can use this function.
#modWA<-crossval(WA(spec, env))
modWA <- crossval(wa_high.res, cv.method="bootstrap", verbose=TRUE, ngroups=10,
			nboot=1000, h.cutoff=0, h.dist=NULL)


x11(4.5,4.5); par(mar=c(3,3,1,3), mgp=c(1.5, 0.5, 0))
segmentwise.rmse(modWA, k=5, main="WA.cla.tol", xlab="MWMT")   # k = model type, 5 = classical, tol down-weighted



palaeoSig











####

# Alternative statistical models for transfer functions

####

# need to convert env to vector, not dataframe

env.vector <- env.cal$MWMT_annual

IKFA(y=spec, x=env.vector, nFact = 5, IsPoly = FALSE, IsRot = TRUE,
	  ccoef = 1:nFact, check.data=TRUE, lean=FALSE)
		# results: 4 factors to get RMSE similar to WA.  R2 poor.  Max bias high.

MAT(y=spec, x=env.vector, dist.method="bray", k=5, lean=TRUE)
		# results: RMSE similar to WA.  R2 equivalent at NO4.  Max bias equal to INV deshrinking.

mrlc <- MLRC(y=spec, x=env.vector, check.data=TRUE, lean=FALSE, n.cut=5, verbose=TRUE)
		# results: RMSE similar to WA.  R2 similar.  Max bias lower than WA-cla (2.8)

	summary(mrlc)
	plot(mrlc, add.ref = TRUE, add.smooth = TRUE) # no bias at extremes
	crossval(mrlc, cv.method="lgo", verbose=TRUE, ngroups=10,
				nboot=100, h.cutoff=0, h.dist=NULL)

mr <- MR(y=spec, x=env.vector, check.data=TRUE, lean=FALSE)
		# RMSE a fair bit lower.  R2 a bit higher.  But max bias = WA.inv
	summary(mr)
	plot(mr, add.ref = TRUE, add.smooth = TRUE) # bias at extremes


wapls <-	WAPLS(y=spec, x=env.cal, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE,
			check.data=TRUE)
		# results: RMSE improves with 3 or more, R2 improves a bit.  Max bias = WA.inv
		summary(wapls)
		plot(wapls, add.ref = TRUE, add.smooth = TRUE) # bias at extremes
		crossval(wapls, cv.method="lgo", verbose=TRUE, ngroups=10,
					nboot=100, h.cutoff=0, h.dist=NULL)



		
		
		
		
##################################################
		###############################
		##############################
		
		Compare MTTI and thermal metrics
		
		##############################
		###############################
#####################################################
		

mets <- read.csv('Metric_BioALL_wSampleInfo_20220927.csv')
		
mets <- mets %>%
	rename('sample.id' = 'UniqueID_v2')

mets.cal <- mets %>%
	filter(dataset == 'CAL')
		


##########	
	# Explore relationships between MTTI and metrics
##########

##
	# Number of taxa
##

p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_stenocold))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_cold))
p + geom_point()+ geom_smooth(method="lm")		



		
p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_cool))
p + geom_point()+ geom_smooth(method="lm")		


p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_warm))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_stenowarm))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_eury))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_cowa))
p + geom_point()+ geom_smooth(method="lm")		

              
p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_na))
p + geom_point()+ geom_smooth(method="lm")		


p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="lm")		


# steno.cold to cool
p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="lm")	


# cool.warm to steno.warm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="lm")	

# warm to steno.warm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=nt_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="lm")	



## 
	# Percent individuals
##

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_stenocold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_cold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_cool))
p + geom_point()+ geom_smooth(method="loess")	



p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_warm))
p + geom_point()+ geom_smooth(method="loess")	


# steno.warm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_eury))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_cowa))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_na))
p + geom_point()+ geom_smooth(method="loess")	

# steno.cold + cold
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="loess")	

# steno.cold, cold, cool
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="loess")	

# cool.warm, warm, steno.warm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

# warm, steno.warm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pi_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


# what if we make a ratio?
p <- ggplot(data=mets.cal, aes(x=MTTI, y=c(pi_ti_stenocold_cold/pi_ti_warm_stenowarm)/100))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=c(pi_ti_warm_stenowarm/pi_ti_stenocold_cold)/100))
p + geom_point()+ geom_smooth(method="loess")	


# +/- ?
p <- ggplot(data=mets.cal, aes(x=MTTI, y=c(pi_ti_stenocold_cold-pi_ti_warm_stenowarm)))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=c(pi_ti_warm_stenowarm-pi_ti_stenocold_cold)))
p + geom_point()+ geom_smooth(method="loess")	


##
	# Percent taxa
##

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_stenocold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_cold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_cool))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_warm))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_eury))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_cowa))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_na))
p + geom_point()+ geom_smooth(method="loess")	

# steno.col, cold
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="loess")	


# steno.col, cold, cool
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="loess")	


# cool.warm, warm, stenowarm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


# warm, stenowarm
p <- ggplot(data=mets.cal, aes(x=MTTI, y=pt_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	



##################################	
##################################
##	
# Explore relationships between MWMT and metrics
##
#################################
#################################



###
	# Number of taxa
###

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_stenocold))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_cold))
p + geom_point()+ geom_smooth(method="lm")		




p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_cool))
p + geom_point()+ geom_smooth(method="lm")		


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_warm))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_stenowarm))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_eury))
p + geom_point()+ geom_smooth(method="lm")		

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_cowa))
p + geom_point()+ geom_smooth(method="lm")		


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_na))
p + geom_point()+ geom_smooth(method="lm")		


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="lm")		


# steno.cold to cool
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="loess")	


# cool.warm to steno.warm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

# warm to steno.warm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=nt_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

# what if we make a ratio?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(nt_ti_stenocold_cold/nt_ti_warm_stenowarm)/100))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(nt_ti_warm_stenowarm/nt_ti_stenocold_cold)/100))
p + geom_point()+ geom_smooth(method="loess")	

# +/- ?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(nt_ti_stenocold_cold-nt_ti_warm_stenowarm)))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(nt_ti_warm_stenowarm-nt_ti_stenocold_cold)))
p + geom_point()+ geom_smooth(method="loess")




## 
# Percent individuals
##

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_stenocold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_cold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_cool))
p + geom_point()+ geom_smooth(method="loess")	



p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_warm))
p + geom_point()+ geom_smooth(method="loess")	


# steno.warm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_eury))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_cowa))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_na))
p + geom_point()+ geom_smooth(method="loess")	

# steno.cold + cold
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="loess")	

# steno.cold, cold, cool
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="loess")	

# cool.warm, warm, steno.warm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

# warm, steno.warm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pi_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


# what if we make a ratio?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pi_ti_stenocold_cold/pi_ti_warm_stenowarm)/100))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pi_ti_warm_stenowarm/pi_ti_stenocold_cold)/100))
p + geom_point()+ geom_smooth(method="loess")	

# +/- ?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pi_ti_stenocold_cold-pi_ti_warm_stenowarm)))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pi_ti_warm_stenowarm-pi_ti_stenocold_cold)))
p + geom_point()+ geom_smooth(method="loess")	


##
# Percent taxa
##

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_stenocold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_cold))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_cool))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_warm))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_eury))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_cowa))
p + geom_point()+ geom_smooth(method="loess")	


p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_na))
p + geom_point()+ geom_smooth(method="loess")	

# steno.col, cold
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_stenocold_cold))
p + geom_point()+ geom_smooth(method="loess")	


# steno.col, cold, cool
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_stenocold_cold_cool))
p + geom_point()+ geom_smooth(method="loess")	


# cool.warm, warm, stenowarm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_cowa_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	


# warm, stenowarm
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=pt_ti_warm_stenowarm))
p + geom_point()+ geom_smooth(method="loess")	

# what if we make a ratio?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pt_ti_stenocold_cold/pt_ti_warm_stenowarm)/100))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pt_ti_warm_stenowarm/pt_ti_stenocold_cold)/100))
p + geom_point()+ geom_smooth(method="loess")	

# +/- ?
p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pt_ti_stenocold_cold-pt_ti_warm_stenowarm)))
p + geom_point()+ geom_smooth(method="loess")	

p <- ggplot(data=mets.cal, aes(x=MWMT_TolAnal, y=c(pt_ti_warm_stenowarm-pt_ti_stenocold_cold)))
p + geom_point()+ geom_smooth(method="loess")	


######################################
######################################
#####
			# regression model for MWMT vs metrics
#####
######################################
######################################


mod <- lm(mets.cal$MWMT_TolAnal ~ mets.cal$pi_ti_stenocold_cold_cool + mets.cal$pi_ti_warm_stenowarm)
summary(mod)
			# r2 = 0.57
			# correlation between these two is r = -0.65
		
mod <- lm(mets.cal$MWMT_TolAnal ~ mets.cal$pi_ti_stenocold + mets.cal$pi_ti_cold + mets.cal$pi_ti_cool + mets.cal$pi_ti_warm + mets.cal$pi_ti_stenowarm)
summary(mod)
# r2 = 0.61

@@@@ problem is these are all likely highly correlated.  Would need to go through formal process for regression model


######################################
######################################
#####
#			RF/CART model for MWMT vs metrics - can we develop a model, using metrics only,
#			that classifies a site into a bio-thermal regime
#####
######################################
######################################

library(randomForest)

# pare down mets dataframe to predictors and response (MWMT)

met.cal_RF <- mets.cal %>%
	select(MWMT_TolAnal, nt_total:pt_ti_na) # L3_eco_name, slope_nhd, wsarea_km2, elev_m, ICI, IWI --- these metrics have mising values and RF crashes
	# currently including lat, long, daynum, CollDate, daynum--as a way of looking to see how much these relate to MWMT compared to bio.temp metrics

RF_mwmt.mets <- randomForest(MWMT_TolAnal ~ ., data=met.cal_RF, importance=TRUE,  keep.forest=TRUE)
RF_mwmt.mets # 74.6% variance explained 

imp <- importance(RF_mwmt.mets) #, conditional = TRUE)  
RF_mwmt.mets$importance

plot(RF_mwmt.mets)
varImpPlot(RF_mwmt.mets)   

getTree(RF_mwmt.mets)


@@@ partial dependence plots
partialPlot(RF_mwmt.mets, met.cal_RF, MWMT_TolAnal, x.var = pt_ti_stenowarm)

impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op <- par(mfrow=c(2,2))
for (i in seq_along(impvar)) {
	partialPlot(RF_mwmt.mets, met.cal_RF, impvar[i], xlab=impvar[i],
					main=paste("Partial Dependence on", impvar[i])) # , ylim=c(0,30)
}
par(op)
