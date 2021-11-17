
## Load libraries
# library(rioja)
library(analogue)
# Note: rioja and analogue mask each other's functions,
# so it's advisable to not load both at the same time.

# Clean workspace
rm(list=ls()); gc()

# Load RA matrix
RA=read.csv('RA.csv')
RA=as.matrix(RA)
dim(RA) # 3875  809
RA[1:5,1:5]

# Load observed environmental data - Sample-specific
obsX=read.csv('ObsEnvData.csv')
obsX=obsX$MWMT_93_11
length(obsX) # 3875
head(obsX)

# WA - Taxa-specific
fit=analogue::wa(RA,env=obsX,deshrink='classical') # deshrink option not important in manual method.
wa=fit$wa.optima
head(wa)
length(wa) # 809

# Inferred model performance
perf=analogue::performance(fit)
perf
#  RMSE       R2 Avg.Bias Max.Bias 
# 3.859    0.078    0.000  -11.485

# Cross-validation inferred model performance
# Options include bootstrapping or k-fold
# Requires "deshrink" argument to be set appropriately above.
set.seed(123)
cv=analogue::crossval(fit,method='kfold',nfold=10,folds=10) # "folds" = how many times to repeat the k-folds analysis.
cv$performance
#           R2       avgBias   maxBias    RMSEP   RMSEP2       s1       s2
# 1 0.02546348 -0.0001541593 -11.71365 3.994871 3.999148 0.184911 3.994871

						###
						###   SLH: see what happens when use the 'rioja' package
						###
						
						library(rioja)
						
						
						spec <- RA
						env <- obsX
						wa.norwest_ALL <- WA(y=spec, x=env, mono=FALSE, tolDW = FALSE, use.N2=TRUE, tol.cut=.01, 
											  check.data=TRUE, lean=FALSE)
						
						crossval(wa.norwest_ALL, cv.method="loo", verbose=TRUE, ngroups=10,
									nboot=100, h.cutoff=0, h.dist=NULL)
						
						
						plot(wa.norwest_ALL, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="inverse",
							  xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE,
							  add.smooth=FALSE)
						
						
# SLH 9/28/21 -- bring in different FULL OR/WA data set, check model performance
orwa <- read.csv('MacroInferTemp_core_primary_20210817_v2.csv')

# reshape the data
library(tidyverse)
orwa.env <- orwa %>%
	select(UniqueID_v2,  MWMT_final)
orwa.env <- orwa %>%
						spread(key=UniqueID_v2, value=MWMT_final)




#### Practitioner would start here.
# Practitioner provides his RA.
# Practitioner uses our WAopt (wa) and our deshrinking coefs


# Raw Inferred X - Sample-specific
infx=((RA %*% wa) / rowSums(RA))[, 1] # Matrix multiplication
length(infx) # 3875
head(infx)
plot(obsX~infx); abline(0,1,col=2) # Raw inferred values

# Deshrinking model
infx_Deshrink=analogue::deshrink(obsX,wa.env=infx,type='inverse')
coefs=infx_Deshrink$coefficients
coefs
#                wa.env 
# -58.340382   4.267188

fit=lm(obsX~infx) # Manual deshrinking model (inverse)
coefs=coef(fit)
coefs
# (Intercept)        infx 
#  -58.340382    4.267188

# Deshrink Inferred X - Sample-specific
infx_Deshrink=infx_Deshrink$env
infx_Deshrink=coefs[1] + coefs[2]*infx # Manual method
length(infx_Deshrink) # 3875
head(infx_Deshrink)

# Plot
plot(obsX~infx); abline(0,1,col=2)          # Raw inferred values
plot(obsX~infx_Deshrink); abline(0,1,col=2) # Deshrunk inferred values
