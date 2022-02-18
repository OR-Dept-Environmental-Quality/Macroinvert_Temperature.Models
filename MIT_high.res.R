# Author: Shannon Hubler
# Initial Date: 11.10.21

# Project: Macroinvertebrate Inferred Temperature

# summary: Building on thermal traits work, as part of the Western OR/WA BCG workgroup.
#			  Develop Weighted Averaging models to infer seasonal 7-d maximum temperatures,
#			  derived from NorWeSTs MWMT metric.

#				HIGH taxonomic resolution models
#				Check model performance
#				Look for bias

 
# packages used
library(tidyverse)
library(rioja)
library(RODBC) 








#####
#####
#				SITE level data, for each SAMPLE
#####
#####

site.data <- read.csv('_SiteData_forShannon_20220213.csv')

site.data <- site.data[order(site.data$UniqueID_v2),]	


#####
#####
#			bring in taxonomy lookup table 
#####
#####

	
	
taxa <- read.csv('_OTU_worksheet_20220216.csv')

	# OTU_high.res = highest taxonomic resolution (genus/species as much as possible)
	
	# '666' = used to filter out and drop taxa not included in model,

# modify to make easier to join/link with bug data files
taxa.otu <- taxa %>%
	select(Bio_ALL_TaxaID_v2, MTTI_OTU_highres_v2) %>%
	rename(TaxaID_v2 = Bio_ALL_TaxaID_v2) %>% # rename OTU to match taxa ID column in bug data files
	rename(OTU_high.res = MTTI_OTU_highres_v2) # renames MTTI_OTU.... to simplify



#####
#####

#		bring in bug data

#####
#####

bugs_7977 <- read.csv('_Bio_ALL_7977.csv') # non-wadeable sites removed

#Split into CAL and VAL datasets
cal <- bugs_7977 %>%
	filter(dataset == 'CAL') # calibration data set, used to build the models (starting 'n' = 3658)



val <- bugs_7977 %>%
	filter(dataset == 'VAL') # independent validation data set, used to get better error 
									 # estimates of model performance (starting 'n' = 603)


	


#####
#####

#		join bugs and OTUs
#####
#####

cal.lim <- cal %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(sample.id = UniqueID_v2)

				# hist(cal.lim$RA)
				# qqnorm(cal.lim$RA ,main="QQ plot of RA data",pch=19)
				# qqline(cal.lim$RA)

			


bug.cal_high.res <- cal.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) 


bug.cal_high.res <- bug.cal_high.res %>%
	select(sample.id, OTU_high.res, RA) %>%
	filter(OTU_high.res != 666) 
# sum RA's across all OTUs--should see a reduction in rows
bug.cal_high.res<-plyr::ddply(.data = bug.cal_high.res, c('sample.id', 'OTU_high.res'), 
							  plyr::summarize, RA=sum(RA))


bug.cal_high.res_totRA <- bug.cal_high.res %>%
				group_by(sample.id) %>%
				summarise(tot.RA = sum(RA))				

hist(bug.cal_high.res_totRA$tot.RA)

f <- bug.cal_high.res_totRA %>%
	filter (tot.RA < 0.4)


# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

bug.cal_high_wide <- bug.cal_high.res %>% 
	pivot_wider(id_cols = sample.id, names_from = OTU_high.res, values_from = RA,
					values_fn = sum) 
bug.cal_high_wide[is.na(bug.cal_high_wide)] <- 0
bug.cal_high_wide <-	column_to_rownames(bug.cal_high_wide, 'sample.id')


								



####
# repeat for VAL
####
val.lim <- val %>%
	select(UniqueID_v2, TaxaID_v2, RA) %>%
	rename(sample.id = UniqueID_v2)

bug.val_high.res <- val.lim %>%
	left_join(taxa.otu, by = c('TaxaID_v2')) %>%
	select(sample.id, OTU_high.res, RA) %>%
	filter(OTU_high.res != 666) 

# sum RA's across all OTUs--should see a reduction in rows
bug.val_high.res<-plyr::ddply(.data = bug.val_high.res, c('sample.id', 'OTU_high.res'), 
										plyr::summarize, RA=sum(RA))



# need to crosstab the bug data (turn into a wide format) so that OTUs are columns

bug.val_high_wide <- bug.val_high.res %>% 
	pivot_wider(id_cols = sample.id, names_from = OTU_high.res, values_from = RA,
					values_fn = sum) 
bug.val_high_wide[is.na(bug.val_high_wide)] <- 0
bug.val_high_wide <-	column_to_rownames(bug.val_high_wide, 'sample.id')

							

#####
#####

#     create ENV data frame (MWMT)

#####
#####

# cal  -----> based on ANNUAL MWMT
env.cal <- site.data %>%
	filter(calval == 'CAL') %>%
	select(UniqueID_v2, MWMT_TolAnal) %>%		#@@@@@@@@@@@@@@@@@@@@@@@@@@@
	rename(sample.id = UniqueID_v2) %>%
	rename(MWMT_annual = 'MWMT_TolAnal') %>%   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	distinct_all()


env.cal <-	column_to_rownames(env.cal, 'sample.id')			

							# # is there model improvement if we drop low and high temp sites--low 'n' and creates high model errors?
							# env.cal_10.25 <- env.cal %>%
							# 	filter(MWMT_annual > 9.99 & MWMT_annual < 25.01)
							# bug.cal_high_wide_10.25 <- left_join(rownames_to_column(env.cal_10.25), rownames_to_column(bug.cal_high_wide), by=c("rowname"))
							# rownames(bug.cal_high_wide_10.25) <- bug.cal_high_wide_10.25$rowname
							# bug.cal_high_wide_10.25 <- bug.cal_high_wide_10.25 %>%
							# 	select(-c(rowname, MWMT_annual))
							



							# hist(env.cal$MWMT_annual)
							# qqnorm(env.cal$MWMT_annual ,main="QQ plot of MWMT data",pch=19)
							# qqline(env.cal$MWMT_annual)

# val --------> based on ANNUAL MWMT
env.val <- site.data %>%
	filter(calval == 'VAL') %>%
	select(UniqueID_v2, MWMT_TolAnal) %>%    #@@@@@@@@@@@@@@@@@@@@@@@@@
	rename(sample.id = UniqueID_v2) %>%
	rename(MWMT_annual = 'MWMT_TolAnal') %>%       #@@@@@@@@@@@@@@@@@@@@@@@@@@@
	distinct_all()


env.val <-	column_to_rownames(env.val, 'sample.id')


							

#####
#####

#		Create inference models with CAL data

#####
#####


spec <- bug.cal_high_wide
env <- env.cal		




							# spec_10.25 <- bug.cal_high_wide_10.25
							# env_10.25 <- env.cal_10.25

wa_high.res <- WA(y=spec, x=env, mono=TRUE, tolDW = TRUE, use.N2=TRUE, tol.cut=.01, 
							check.data=TRUE, lean=FALSE)

		# look at WA-PLS
		# wa.pls_high.res <- WAPLS(y=spec, x=env, npls=5, iswapls=TRUE, standx=FALSE, lean=FALSE,
		# check.data=TRUE)
		# minimal to no improvements
		
		

wa_high.res 
				

crossval(wa_high.res, cv.method="bootstrap", verbose=TRUE, ngroups=10,
			nboot=1000, h.cutoff=0, h.dist=NULL)
		# crossval results almost exactly the same as original



performance(wa_high.res)

							# names(wa_high.res)
							# wa_high.res$fitted.values
							# wa_high.res$coefficients

WA.resid_high <- residuals(wa_high.res)	
							# hist(WA.resid_high[,2])
WA.resid_high <- (WA.resid_high[,2])
WA.resid.high <- as.data.frame(WA.resid_high)
				# WA.resid.high <- WA.resid.high %>%
				# 	mutate(sample.id = rownames(WA.resid.high))

# plot inferred vs MWMT
plot(wa_high.res, resid=FALSE, xval=FALSE, tolDW=TRUE, deshrink="classical",
	  xlab="MWMT", ylab="MTSI (inferred MWMT)", main="Classical deshrinking: tolerance downweighted", 
	  ylim=c(0,40), xlim=c(0,40), add.ref=TRUE, add.smooth=TRUE)

# plot residuals -- INVERSE shows BIAS, CLASSICAL doesn't
plot(wa_high.res, resid=TRUE, xval=FALSE, tolDW=TRUE, deshrink="classical",
	  xlab="MWMT", ylab="residuals (MTSI)", ylim=c(-15,15), xlim=c(0,35), add.ref=TRUE,
	  add.smooth=TRUE, main='Classical deshrinking: tolerance downweighted')

				

# ######
				
# make predictions for VAL dataset				
				
######
				
wa_high.res_VAL <- predict(wa_high.res, newdata=bug.val_high_wide, sse=TRUE, nboot=100,
									match.data=TRUE, verbose=TRUE)

#### get RMSE and r2 for VAL data


# extract VAL predictions as a standalone dataframe, so can join with env.val and calculate RMSEP
VAL.predicts <- as.data.frame(wa_high.res_VAL$fit)
VAL.predicts <- cbind(VAL.predicts, env.val)
	


VAL.rmsep <- VAL.predicts %>%
	mutate(VAL_WA.inv_rmsep = sqrt(mean((VAL.predicts$WA.inv - VAL.predicts$MWMT_annual)^2))) %>%
	mutate(VAL_WA.cla_rmsep = sqrt(mean((VAL.predicts$WA.cla - VAL.predicts$MWMT_annual)^2)))	%>%	 
	mutate(VAL_WA.m_rmsep = sqrt(mean((VAL.predicts$WA.m - VAL.predicts$MWMT_annual)^2)))		 %>%
	mutate(VAL_WA.inv.tol_rmsep = sqrt(mean((VAL.predicts$WA.inv.tol - VAL.predicts$MWMT_annual)^2)))	%>%	 
	mutate(VAL_WA.cla.tol_rmsep = sqrt(mean((VAL.predicts$WA.cla.tol - VAL.predicts$MWMT_annual)^2)))	%>%	 
	mutate(VAL_WA.m.tol_rmsep = sqrt(mean((VAL.predicts$WA.m.tol - VAL.predicts$MWMT_annual)^2)))	%>%
	mutate(VAL.WA.inv_rsq = cor(VAL.predicts$WA.inv, VAL.predicts$MWMT_annual)^2) 	%>%
	mutate(VAL.WA.cla_rsq = cor(VAL.predicts$WA.cla, VAL.predicts$MWMT_annual)^2)		%>%
	mutate(VAL.WA.m_rsq = cor(VAL.predicts$WA.m, VAL.predicts$MWMT_annual)^2)			%>%
	mutate(VAL.WA.inv.tol_rsq = cor(VAL.predicts$WA.inv.tol, VAL.predicts$MWMT_annual)^2)	%>%
	mutate(VAL.WA.cla.tol_rsq = cor(VAL.predicts$WA.cla.tol, VAL.predicts$MWMT_annual)^2)	%>%
	mutate(VAL.WA.m.tol_rsq = cor(VAL.predicts$WA.m.tol, VAL.predicts$MWMT_annual)^2)

VAL.rmsep <- VAL.rmsep[,8:19]
VAL.rmsep <- distinct(VAL.rmsep)
row.names(VAL.rmsep) <- NULL


######

#			Look for potential BIAS in predictions - natural gradients, disturbance, etc.

#####
					
# join MWMT and Inferred values, plus site specific data
		
site.data_calval <- site.data %>%
	filter(calval != 'NOT_CALVAL')

eco.cnt <- site.data_calval %>%
	group_by(calval, eco3) %>%
	summarize(count = n()) %>%
	pivot_wider(names_from = calval, values_from = count)

# numerical summaries -- need to remove categorical vars
site.data_num <- site.data_calval %>%
	rename(sample.id = UniqueID_v2) %>%
	select(sample.id, siteid, calval, lat, long, daynum, Year, MWMT_TolAnal,  slope_nhd, wsarea_km2, elev_m, IWI, 
			 PctUrbWs, PctAgWs, PctForWs, RdDensWs, Tmean8110Ws, Precip8110Ws)

site.data_num <- site.data_num %>%
	rename(MWMT_Annual = MWMT_TolAnal)


site.data_num.long <- site.data_num %>%
	pivot_longer(!c(sample.id, siteid, calval)  , names_to = "variable", values_to = "value")


sum.stats <- site.data_num.long %>%                               # Summary by group using dplyr
					group_by(calval, variable) %>% 
					summarize(min = min(value, na.rm = TRUE),
								 q1 = quantile(value, 0.25, na.rm=TRUE),
								 median = median(value, na.rm=TRUE),
								 mean = mean(value, na.rm=TRUE),
								 q3 = quantile(value, 0.75, na.rm=TRUE),
								 max = max(value, na.rm=TRUE))


	
write.csv(sum.stats, 'sum.stats.csv')
	
# boxplots
	site.data_num.long %>% 
	ggplot(aes(x= calval, y=value, fill=calval)) +
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


@@@@@@@@@@@@@@@@@@@@@@
# data transformations
@@@@@@@@@@@@@@@@@@@@@@								#colnames(site.data_pca)

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
	select(sample.id, siteid, calval, lat, long, Year, daynum, MWMT_Annual, slope.log, 
			 wsarea.log, elev.sqrt, IWI, pct.urban_log, pct.ag_log, pct.forest_asin, 
			 road.den_log, Tmean8110Ws, Precip_sqrt)

				

								#rownames(site.data_complete) <- site.data_complete[,1]
								#site.data_pca <- site.data_complete[,-c(1:3)]


# run pca 
library(FactoMineR)  
library("factoextra")
library("corrplot")


pca.allsites <- PCA(site.data_pca.trans[ , c(4:18)], graph = FALSE, scale.unit = TRUE)

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
fviz_contrib(pca.allsites, choice = "var", axes = 1:2, top = 15)


# The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
fviz_pca_var(pca.allsites, col.var = "contrib",
				 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)


fviz_pca_ind(pca.allsites, col.ind = "cos2", 
				 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
				 repel = TRUE # Avoid text overlapping (slow if many points)
)


grp <- as.factor(site.data_pca$calval)
# Color variables by groups
fviz_pca_var(pca.allsites, col.var = grp, 
				 palette = c("#0073C2FF", "#EFC000FF"),
				 legend.title = "Cluster")


fviz_pca_ind(pca.allsites,
				 col.ind = grp, # color by groups
				 palette = c("#00AFBB",  "#FC4E07"),
				 addEllipses = FALSE, # Concentration ellipses
				 ellipse.type = "confidence",
				 legend.title = "Groups",
				 repel = TRUE
)

@@@ how to unlabel all points (individuals)?
	
@@@ how to add variable biplot as well?  It owuld be nice to see what is driving the sites in upper-left.











######

# BIAS = look for patterns in CAL predictions--do we see bias associated with other variables?
			@@@@ do we need to look at VAL as well?
######

WA.resid.high_2 <- WA.resid.high %>%
 	mutate(sample.id = rownames(WA.resid.high)) # get site ID as a column for join



site.data_residuals <- site.data %>%
							 rename(sample.id = UniqueID_v2)  %>%
							 inner_join(WA.resid.high_2, by = c('sample.id'))


										# site.data_residuals$eco3 <- plyr::revalue(site.data_residuals$eco3, c("Blue Mountains" = "B.Mts",                                   
										# 			"Cascades" = "Casc", "Coast Range" = "Co.Ra", "Columbia Plateau" = "Col.Pl",
										# 			"Eastern Cascades Slopes and Foothills" = "E.Casc",
										# 			"Klamath Mountains/California High North Coast Range" = "Kl.Mts", 
										# 			"North Cascades" = "N.Casc", "Northern Basin and Range" = "N.Ba.Ra",
										# 			"Northern Rockies" = "N.Rock", "Puget Lowland" = "Pu.Low",
										# 			"Snake River Plain" = "Sna.Riv.Pl", "Willamette Valley" = "Wil.Val"))     
	
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
p + geom_point()	+ geom_smooth(method="lm")	+ xlim(c(0, 25000))	

@@@@@ problematic?
	

# residuals ~ elevation
p <- ggplot(data=site.data_residuals, aes(x=elev_m, y=WA.resid_high))
p + geom_point()	+ geom_smooth(method="lm")				
								 
								 
# residuals ~ IWI
p <- ggplot(data=site.data_residuals, aes(x=IWI, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")	


@@@@@@@@@ MTTI higher than observed MWMT at disturbed sites--makes sense?  
@@@@@@@@@ More tolerant taxa = higher WAopt 

# residuals ~ daynum
p <- ggplot(data=site.data_residuals, aes(x=daynum, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")				


# residuals ~ year
p <- ggplot(data=site.data_residuals, aes(x=Year, y=WA.resid_high))
p + geom_point()		+ geom_smooth(method="lm")				



# residuals ~ ecocregions
p <- ggplot(data=site.data_residuals, aes(x=eco3, y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 0, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 
								 
# residuals ~ Year
p <- ggplot(data=site.data_residuals, aes(x=as.factor(Year), y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 90, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 
		
		p + geom_point()		+ geom_smooth(method="lm")				

# residuals ~ Collection Method
p <- ggplot(data=site.data_residuals, aes(x=as.factor(CollMeth), y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 0, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 


# residuals ~ Source
p <- ggplot(data=site.data_residuals, aes(x=as.factor(SourceEntity), y=WA.resid_high))								 
p + geom_boxplot() + theme(axis.text.x=element_text(size=8, angle = 90, vjust=.5))+
	geom_hline(yintercept=0, linetype="dashed", color = "red") 



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



