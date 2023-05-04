
###################################*
### MWMT - Split Cal/Val ####
###################################*


###################################*
#### Contact Info ####
###################################*

# Mark Fernandez (mark.fernandez@tetratech.com)



###################################*
### Citations ####
###################################*

citation()
# R Core Team. 2023. R: A language and environment for statistical computing. R Foundation for Statistical Computing,Vienna,Austria. URL https://www.R-project.org/.




###################################.
### Load Packages ####
###################################.

# Load libraries
# options(digits=7L)
par(mar=c(4.5,4.5,1,1),bty='l')
library(openxlsx)    # read.xlsx()
library(dplyr)       # slice_sample()
library(data.table)  # []
library(ggplot2)     # ggplot()
library(ggmap)       # ggmap()
theme_set(theme_classic())
theme_update(legend.position='bottom')
# detach('package:aaa',unload=T)




# Clean workspace
rm(list=ls()); gc()

# Set na.strings
na.strings=c(NA,'NA','N/A','#N/A','na','-88','',' ','None','none','<Null>')
seed=27709L
pointsize=9L
txtSize=9L

# labels for plot ticks
labelsPretty=function(x){
  formatC(x,digits=NULL,big.mark=',',format='fg')
}


# Set arguments
formals(melt.data.table)$na.rm=T
formals(melt.data.table)$variable.factor=F
formals(dcast.data.table)$fill=NA
formals(dcast.data.table)$value.var='value'
formals(fread)$sep=','
formals(fread)$skip=0
formals(fread)$header=T
formals(fread)$check.names=F
formals(fread)$na.strings=na.strings
formals(png)$units='in'
formals(png)$res=300
formals(png)$type='cairo'
formals(png)$pointsize=pointsize
formals(expand.grid)$KEEP.OUT.ATTRS=F
formals(expand.grid)$stringsAsFactors=F
formals(rbindlist)$use.names=T
formals(rbindlist)$fill=T
formals(read.xlsx)$startRow=1L
formals(read.xlsx)$check.names=F
formals(read.xlsx)$na.strings=na.strings





###################################*
### Load Data ####
###################################*

# SiteDate
rm(SiteDate)
SiteDate=read.xlsx('Data/VALdesignations_forMark_20230308.xlsx',sheet='CAL_samples')
setDT(SiteDate)
dim(SiteDate)                 # 3658 xx
uniqueN(SiteDate$UniqueID_v2) # 3658  - Same :)

# Do NOT rename fields :)

# ID
setorder(SiteDate,UniqueID_v2)
SiteDate[,ID:=.I]

# Check
SiteDate[,.N,by=dataset]
SiteDate[,.N,by=NORWEST_Unit]
SiteDate[,.N,by=eco3]
plot(lat~long,data=SiteDate)

# ranges
tmp=SiteDate[,.(Min=min(MWMT_TolAnal),
                Max=max(MWMT_TolAnal)),
             by=.(NORWEST_Unit)]
tmp
#    NORWEST_Unit      Min      Max
# 1:       MidCol  5.85000 30.58000
# 2:     MidSnake  9.10000 30.77000
# 3:      ORCoast  3.59000 28.70000
# 4:      WACoast  5.42000 26.79000
# 5:         SCOR  9.45000 29.09000
# 6:     UpColYak  3.60000 27.42000
# 7:      SpoKoot 10.88832 25.32726

# MWMT_Bin
SiteDate$MWMT_Bin=NULL
SiteDate[,MWMT_BinChar:=cut(MWMT_TolAnal,breaks=5,ordered_result=T),
         by=.(NORWEST_Unit)]
SiteDate[,MWMT_Bin:=cut(MWMT_TolAnal,breaks=5,ordered_result=T,labels=F),
         by=.(NORWEST_Unit)]

# counts - to ensure a more accurate Val count per strata
rm(counts)
counts=SiteDate[,.(nGroup=.N),
                by=.(NORWEST_Unit,MWMT_Bin)]
dim(counts) # 35 xx
counts[,nGroupVal:=round(nGroup*0.10)]
counts[,nGroupVal:=as.integer(nGroupVal)]

# Try to make each group have at least 1 Val sample
counts[nGroupVal==0,.N] # 5 groups
counts[nGroupVal==0,nGroupVal:=1]
counts[nGroup==1,nGroupVal:=0]
counts[,sum(nGroupVal)/sum(nGroup)] # 0.1008748  ~ 10%

# Stratified random sample
if(F){ # Doesn't guarantee each group has a sample. So, do a manual loop.
  rm(tmp)
  tmp=SiteDate %>%
    group_by(NORWEST_Unit,MWMT_Bin) %>%
    slice_sample(n=nGroupVal,replace=F)
  setDT(tmp)
  SiteDate[,Val:=0L]
  SiteDate[UniqueID_v2 %in% tmp$UniqueID_v2,Val:=1]
} # END if(F)

# Stratified random sample
SiteDate[,Val:=0L]
i=1L
# i=26
set.seed(seed)
for(i in 1:nrow(counts)){ # Loop across groups
  
  # Subset
  rm(tmp)
  tmp=SiteDate[NORWEST_Unit==counts[i,NORWEST_Unit] & MWMT_Bin==counts[i,MWMT_Bin],ID]
  # length(tmp) # 194
  
  # Check
  if(counts[i,nGroupVal]==0) next
  
  # Sample
  SiteDate[ID %in% sample(tmp,size=counts[i,nGroupVal],replace=F),Val:=1]
  
} # END loop across groups
warnings()

# Check
SiteDate[,.N,by=Val]
#    Val    N
# 1:   0 3289
# 2:   1  369
SiteDate[,sum(Val)/.N]              # 0.1008748
counts[,sum(nGroupVal)/sum(nGroup)] # 0.1008748

# dataset
SiteDate[,dataset:='CAL']
SiteDate[Val==1,dataset:='VAL']
SiteDate[,.N,by=dataset]

#### Tables
# table
rm(df)
df=SiteDate[,.N,by=.(NORWEST_Unit,MWMT_Bin,dataset,Val)]
dim(df) # 69 xx

# All
rm(dt)
dt=dcast(df,NORWEST_Unit~MWMT_Bin,fun.aggregate=sum,value.var='N',fill=0L)
fwrite(dt,'Tables/Counts_All.csv')

# Cal
rm(dt)
dt=dcast(df[Val==0,],NORWEST_Unit~MWMT_Bin,fun.aggregate=sum,value.var='N',fill=0L)
fwrite(dt,'Tables/Counts_Cal.csv')

# Val
rm(dt)
dt=dcast(df[Val==1,],NORWEST_Unit~MWMT_Bin,fun.aggregate=sum,value.var='N',fill=0L)
fwrite(dt,'Tables/Counts_Val.csv')

# Export
setorder(SiteDate,NORWEST_Unit,UniqueID_v2)
fwrite(SiteDate,'Tables/CalValSplit_ExtraFields.csv')

# Export
setorder(SiteDate,UniqueID_v2)
fwrite(SiteDate[,.(UniqueID_v2,dataset)],'Tables/CalValSplit.csv')





###################################*
### Density Plots ####
###################################*

# Density Plots
# png('Plots/DensityPlots.png',width=20,height=10)
ggplot(SiteDate,aes(x=MWMT_TolAnal))+
  geom_density(adjust=0.9,fill='grey')+
  geom_rug()+
  labs(y='Density',x='MWMT (C)',title=NULL)+
  facet_grid(dataset~NORWEST_Unit,scales='fixed')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(text=element_text(size=txtSize))+
  theme(legend.position='none')
# dev.off()




###################################*
### Boxlots ####
###################################*

# Boxlots
png('Plots/Boxplots.png',width=6.5,height=7)
set.seed(seed)
ggplot(SiteDate,aes(y=MWMT_TolAnal,x=dataset,color=dataset))+
  geom_jitter(position=position_jitter(w=0.1,h=0),alpha=0.8)+
  geom_boxplot(outlier.shape=NA,fill=NA)+
  
  labs(y='Density',x='MWMT (C)',title=NULL)+
  facet_wrap(~NORWEST_Unit,nrow=3,scales='fixed')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(text=element_text(size=txtSize))+
  theme(legend.position='none')
dev.off()




###################################*
### Basemap ####
###################################*

# Run once
if(F){
  # Set API key - Once per R session
  # register_google(key='AIzaSyBtViuE2LzvfOKVg2XQob8T8K3eZ3Dqx1s',write=T)
  has_google_key() # TRUE
  
  # Download basemap
  tmp=c(-120.0,46.0)
  tmp1=6
  basemap=list()
  # basemap$terrainStamen=get_map(location=tmp,zoom=tmp1,source='stamen',maptype='terrain')
  # basemap$terrain=get_map(location=tmp,zoom=tmp1,source='google',maptype='terrain')
  # basemap$hybrid=get_map(location=tmp,zoom=tmp1,source='google',maptype='hybrid')
  basemap$satellite=get_map(location=tmp,zoom=tmp1,source='google',maptype='satellite')
  # basemap$roadmap=get_map(location=tmp,zoom=tmp1,source='google',maptype='roadmap')
  
  # Save map
  save(basemap,file='basemap.RData')
} # END if(F)





###################################*
### Maps ####
###################################*

# Map
setorder(SiteDate,dataset)
png('Plots/Map.png',width=6.5,height=8.2,units='in',res=300,type='cairo',pointsize=pointsize)
ggmap(basemap$satellite)+ # terrain terrainStamen satellite roadmap hybrid watercolor
  geom_point(data=SiteDate,mapping=aes(x=long,y=lat,fill=dataset),shape=24,size=2,alpha=0.8)+
  coord_map(xlim=c(-125.9,-115.5),ylim=c(41.1,49.7))+
  labs(x='Longitude',y='Latitude',title=NULL,fill=NULL)+
  theme(legend.position='bottom')
dev.off()

# Map - By MWMT - Dataset by shape isn't distinguishable :(
png('Plots/Map_MWMT.png',width=6.5,height=8.2,units='in',res=300,type='cairo',pointsize=pointsize)
ggmap(basemap$satellite)+ # terrain terrainStamen satellite roadmap hybrid watercolor
  geom_point(data=SiteDate,mapping=aes(x=long,y=lat,color=MWMT_TolAnal),size=1.5,alpha=0.8)+
  coord_map(xlim=c(-125.9,-115.5),ylim=c(41.1,49.7))+
  labs(x='Longitude',y='Latitude',title=NULL,color='MWMT (C)',shape=NULL)+
  scale_color_distiller(palette='YlOrBr',direction=1)+ # OrRd PuOr Oranges YlOrBr
  theme(legend.position='bottom')
dev.off()


























