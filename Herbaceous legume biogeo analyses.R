###HERBACEOUS LEGUME BIOGEOGRAPHY###
#Author: Ben Taylor
#bentonneiltaylor@gmail.com

###These analyses aim to understand the broad biogeographic (primarily latitudinal)
###patterns of herbaceous N-fixing legumes. To do this, we analyze a large
###network of broadly distributed grassland experiments from the NutNet and CoRRE networks

###########################################################################
### Read in the data and load R packages ##################################
setwd("C:/Users/Owner/Desktop/Work Files/Research/Herbaceous Legume Biogeography/Data Analysis & Figure Files/herbaceous-N-fixer-biogoegraphy")
nnplt<-read.csv("NutNet_plot_data_31-Aug-2018.csv")#reads in NutNet plot data
nnab<-read.csv("NutNet_raw_abundance_31-Aug-2018.csv")#reads in NutNet abundance data
cplt<-read.csv("CORRE_site_summary.csv")#reads in CoRRE plot data
cab<-read.csv("CORRE_raw_abundance.csv")#reads in CoRRE abundance data
ctrt<-read.csv("CORRE_treatment_summary.csv")#reads in CoRRE treatment data
gex<-read.csv("GEx_All data_Genera Species cleaned_BNT_2019.csv")
gexmeta<-read.csv("GEx_MetaData_forAnalysis_v3_July2016.csv")
nfsp<-read.csv("N fixer species list_Duncan.csv")#reads in a species list of confirmed N-fixers

#loads necessary libraries
library(bbmle)
library(vegan)
library(mapdata)
library(lme4)
library(tidyverse)
library(maps)
library(ggplot2)
source("C:\\RCode\\Homegrown Functions.R")

###########################################################################
### Getting CoRRE Data into shape #########################################
###########################################################################
cab$site.plt.yr<-with(cab, paste(site_code,"_",project_name,"_",block,"_",plot_id,"_",calendar_year))
cab$genus_species<-as.factor(paste(cab$genus,"_",cab$species))
cab<-merge(cab,cplt,by="site_code")
cab<-merge(cab,ctrt[,-1],by=c("site_code","project_name","calendar_year","treatment_year","treatment","community_type"),
           all.x=T,all.y=F)

#### Identifying N-fixing Species ##################################################################
nfsp$genus_species<-as.factor(paste(tolower(nfsp$Genus),"_",tolower(nfsp$Species)))
nfsp$Genus<-tolower(nfsp$Genus)
fixlst<-nfsp[nfsp$Nod.==1,c("Genus","Nod.")]
colnames(fixlst)<-c("genus","N_fixer")
fixlst<-unique(fixlst)
fixlst$genus<-as.factor(trimws(fixlst$genus))


cab<-merge(cab,fixlst,by="genus",all.x=T,all.y=F)
cab<-cab[!duplicated(cab$X),]
cab$N_fixer<-ifelse(is.na(cab$N_fixer),yes=0,no=cab$N_fixer)
####################################################################################################

#### Calculating total abundance and relative abundance of N fixers in each plot ###################
cplts<-unique(cab$site.plt.yr)
cpab<-NULL
for(i in 1:length(cplts)){
  print(i)
  temp<-cab[cab$site.plt.yr==cplts[i],]
  temp<-temp[!is.na(temp$site.plt.yr),]
  temp<-droplevels(temp)
  tempdf<-data.frame("site.plt.yr"=unique(temp$site.plt.yr),
                     "treatment_year"=unique(temp$treatment_year),
                     "totab"=sum(temp$abundance, na.rm=T),
                     "fixab"=sum(temp[temp$N_fixer==1,]$abundance, na.rm=T))
  tempdf$fixra<-with(tempdf, (fixab/totab)*100)
  cpab<-rbind(cpab,tempdf)
}

cabplots<-cab[,c("site.plt.yr","site_code")]
cabplots<-cabplots[!duplicated(cabplots$site.plt.yr)&!is.na(cabplots$site.plt.yr),]
cr<-merge(cabplots,cpab,by="site.plt.yr",all.x=T,all.y=T)

cr<-merge(cr,cplt,by="site_code")
#####################################################################################################

#### Calculating species richness for fixers and non-fixers in each plot ############################
cplts<-unique(cab$site.plt.yr)
cabrich<-NULL
for(p in 1:length(cplts)){
  print(p)
  temp<-cab[cab$site.plt.yr==cplts[p],]
  tempdf<-data.frame("site.plt.yr"=temp$site.plt.yr,
                     #"treatment_year"=unique(temp$treatment_year),
                     "tot.rich"=length(unique(temp$genus_species)),
                     "fix.rich"=length(unique(temp[temp$N_fixer==1,]$genus_species)))
  cabrich<-rbind(cabrich,tempdf)
}
cabrich$fixrr<-with(cabrich, (fix.rich/tot.rich)*100)
cabrich<-cabrich[!duplicated(cabrich),]

#####################################################################################################

### Calculating N limitation for each site ##########################################################
cabtrt<-cab[cab$treatment_year==1,]
cabsites<-unique(cabtrt$site_code)
cabsitelim<-NULL
for(i in 1:length(cabsites)){
  print(i)
  temp1<-cabtrt[cabtrt$site_code==cabsites[i],]#makes a temporary dataframe with an individual site and just the 1st year of trt data
  ctempcon<-temp1[temp1$n==0&temp1$p==0&temp1$k==0&temp1$CO2==0&temp1$precip==0&temp1$temp==0&temp1$mow_clip==0&temp1$burn==0&temp1$herb_removal==0,]
  ctempconoth<-temp1[temp1$n==0,]
  ctempN<-temp1[temp1$n>0&temp1$p==0&temp1$k==0&temp1$CO2==0&temp1$precip==0&temp1$temp==0&temp1$mow_clip==0&temp1$burn==0&temp1$herb_removal==0,]
  ctempNoth<-temp1[temp1$n>1,]
  limdat<-as.data.frame(unique(temp1$site_code))
  colnames(limdat)<-"site_code"
  limdat$ctrl.cov<-(sum(ctempcon$abundance)/length(unique(ctempcon$plot_id)))
  limdat$ctrloth.cov<-(sum(ctempconoth$abundance)/length(unique(ctempconoth$plot_id)))
  limdat$N.cov<-(sum(ctempN$abundance)/length(unique(ctempN$plot_id)))
  limdat$Noth.cov<-(sum(ctempNoth$abundance)/length(unique(ctempNoth$plot_id)))
  cabsitelim<-rbind(cabsitelim,limdat)
}
cabsitelim$NRR<-with(cabsitelim, ((N.cov-ctrl.cov)/ctrl.cov)*100)
cabsitelim$NothRR<-with(cabsitelim, ((Noth.cov-ctrloth.cov)/ctrloth.cov)*100)
cabsitelim<-cabsitelim[!is.na(cabsitelim$NRR),]

csl.mrg<-cabsitelim[,c(1,6)]
#####################################################################################################

### Calculating P limitation for each site ##########################################################
cabtrt<-cab[cab$treatment_year==1,]
cabsites<-unique(cabtrt$site_code)
cabsitePlim<-NULL
for(i in 1:length(cabsites)){
  print(i)
  temp1<-cabtrt[cabtrt$site_code==cabsites[i],]#makes a temporary dataframe with an individual site and just the 1st year of trt data
  ctempcon<-temp1[temp1$n==0&temp1$p==0&temp1$k==0&temp1$CO2==0&temp1$precip==0&temp1$temp==0&temp1$mow_clip==0&temp1$burn==0&temp1$herb_removal==0,]
  ctempconoth<-temp1[temp1$p==0,]
  ctempP<-temp1[temp1$p>0&temp1$n==0&temp1$k==0&temp1$CO2==0&temp1$precip==0&temp1$temp==0&temp1$mow_clip==0&temp1$burn==0&temp1$herb_removal==0,]
  ctempPoth<-temp1[temp1$p>1,]
  Plimdat<-as.data.frame(unique(temp1$site_code))
  colnames(Plimdat)<-"site_code"
  Plimdat$ctrl.cov<-(sum(ctempcon$abundance)/length(unique(ctempcon$plot_id)))
  Plimdat$ctrloth.cov<-(sum(ctempconoth$abundance)/length(unique(ctempconoth$plot_id)))
  Plimdat$P.cov<-(sum(ctempP$abundance)/length(unique(ctempP$plot_id)))
  Plimdat$Poth.cov<-(sum(ctempPoth$abundance)/length(unique(ctempPoth$plot_id)))
  cabsitePlim<-rbind(cabsitePlim,Plimdat)
}
cabsitePlim$PRR<-with(cabsitePlim, ((P.cov-ctrl.cov)/ctrl.cov)*100)
cabsitePlim$PothRR<-with(cabsitePlim, ((Poth.cov-ctrloth.cov)/ctrloth.cov)*100)
cabsitePlim<-cabsitePlim[!is.na(cabsitePlim$PRR),]

csPl.mrg<-cabsitePlim[,c(1,6)]
#####################################################################################################

#### Merging in richness and site nutrient limitation columns #######################################
cr<-merge(cr,csl.mrg,by="site_code",all.x=T,all.y=F)
cr<-merge(cr,csPl.mrg,by="site_code",all.x=T,all.y=F)
cr<-merge(cr,cabrich, by="site.plt.yr",all.x=T,all.y=F)
cr$dataset<-"CoRRE"

colnames(cr)[c(10,11,15,16)]<-c("LAT","LON","Nlim","Plim")
cr<-cr[,c(1,2,3,4,5,6,10,11,13,14,15,16,17,18,19,20)]
cr$abs.LAT<-abs(cr$LAT)
#####################################################################################################

#### Subsetting for just the control plots or the pre-treatmnet data ################################
cr.pt<-cr[cr$treatment_year==0,]
#cabctl<-dat[dat$treatment%in%c("T0F0","control","UnwarmedControl","N0P0S0",
#                               "ambient","0N0P","Control","N0F0","AcAt",
#                               "N0M0","AMBIENT","N0","amb","XXX","0_CONTROL",
#                               "N0P0","Cont","ambient_control","CONTROL","N0B0",
#                               "Reference"),]
#cabctl2<-dat[dat$nutrients==0&dat$light==0&dat$carbon==0&dat$water==0&dat$other_manipulation==0&dat$other_trt==0,]
#cabctl<-rbind(cabpt,cabctl,cabctl2)
#cabctl<-unique(cabctl)
#cabctl<-cabctl[!is.na(cabctl$site.plt.yr),]
######################################################################################################

###########################################################################
### Getting GEx Data into shape ###########################################
###########################################################################
gex$site.plt.yr<-with(gex, paste(site,"_",block,"_",plot,"_",year))
gex$site.plt<-with(gex, paste(site,"_",block,"_",plot))
gexmeta.mrg<-gexmeta[,c("site","Final.Lat","Final.Long","precip","grazing.pressure")]
gex<-merge(gex,gexmeta.mrg, by="site",all.x=T)

gex<-merge(gex,fixlst,by="genus",all.x=T,all.y=F)
gex$N_fixer<-ifelse(is.na(gex$N_fixer),yes=0,no=gex$N_fixer)

#Exporting File for Sally
#gex_sally<-gex[,c(3,2,4:10,20,1,11:19)]
#write.csv(gex_sally, file="GEx_Names Cleaned and N fix_BNT_2019.csv", row.names=F)

gex.g<-gex[gex$trt=="G",]

#### Calculating total abundance and relative abundance of N fixers in each plot (only uses the grazed, i.e. unfenced, data)
gplts<-unique(gex.g$site.plt.yr)
gpab<-NULL
for(i in 1:length(gplts)){
  print(i)
  temp<-gex.g[gex.g$site.plt.yr==gplts[i],]
  temp<-temp[!is.na(temp$site.plt.yr),]
  temp<-droplevels(temp)
  tempdf<-data.frame("site.plt.yr"=unique(temp$site.plt.yr),
                     "site.plt"=unique(temp$site.plt),
                     "year"=unique(temp$year),
                     "totab"=sum(temp$relcov, na.rm=T),
                     "fixab"=sum(temp[temp$N_fixer==1,]$relcov, na.rm=T))
  tempdf$fixra<-with(tempdf, (fixab/totab)*100)
  gpab<-rbind(gpab,tempdf)
}


#### Calculating Grazing response ratio for each site from 1st year of exclosure ##########
# This does not currently work because the plant data are all in "relative cover" 
#gplts<-unique(gex$site.plt)
#gsiteGlim<-NULL
#for(i in 1:length(gplts)){
#  print(i)
#  temp<-gex[gex$site.plt==gplts[i],]
#  temp2<-temp[temp$exage==min(temp$exage),]
#  tempg<-temp2[temp2$trt=="G",]
#  tempu<-temp2[temp2$trt=="U",]
#  Glimdat<-as.data.frame(unique(temp2$site.plt))
#  colnames(Glimdat)<-"site.plt"
#  Glimdat$grazed<-sum(tempg$relcov)
#  Glimdat$ungrazed<-sum(tempu$relcov)
#  gsiteGlim<-rbind(gsiteGlim,Glimdat)
#}
############################################################################################

gplots1<-gex.g[,c("site.plt.yr","site")]
gplots1<-gplots1[!duplicated(gplots1$site.plt.yr)&!is.na(gplots1$site.plt.yr),]
gplots<-merge(gplots1,gpab, by="site.plt.yr",all.x=T,all.y=T)

gplt.info<-gex.g[,c("site","Final.Lat","Final.Long","precip")]
gplt.info<-unique(gplt.info)
g<-merge(gplots,gplt.info,by="site",all.x=T,all.y=T)
g$MAT<-NA
g$Nlim<-NA
g$Plim<-NA

##########################################################################################

#calculating species richness for fixers and non-fixers in just the grazed section of each plot
gplts<-unique(gex.g$site.plt.yr)
grich<-NULL
for(p in 1:length(gplts)){
  print(p)
  temp<-gex.g[gex.g$site.plt.yr==gplts[p],]
  tempdf<-data.frame("site.plt.yr"=unique(temp$site.plt.yr),
                     "tot.rich"=length(unique(temp$genus_species)),
                     "fix.rich"=length(unique(temp[temp$N_fixer==1,]$genus_species)))
  grich<-rbind(grich,tempdf)
}
grich$fixrr<-with(grich, (fix.rich/tot.rich)*100)
g<-merge(g,grich,by="site.plt.yr",all.x=T,all.y=T)

#### Extracting just the first year of data for each plot #################################
g.y1<-NULL
siteplots<-unique(g$site.plt)
for(i in 1:length(siteplots)){
  print(i)
  temp<-g[g$site.plt==siteplots[i],]
  temp2<-temp[temp$year==min(temp$year),]
  g.y1<-rbind(g.y1,temp2)
}
###########################################################################################

g<-g[,c(1,2,4:9,11,10,12:16)]
g.y1<-g.y1[,c(1,2,4:9,11,10,12:16)]
g$dataset<-"GEx"
g.y1$dataset<-"GEx"
g$abs.LAT<-abs(g$Final.Lat)
g.y1$abs.LAT<-abs(g.y1$Final.Lat)

clnms<-names(cr)
colnames(g)<-clnms
colnames(g.y1)<-clnms

#### Adding GEx and CoRRE data together to create the "dat" object ############################
dat<-rbind(cr,g)

#### Adding first year of CoRRE and GEx data together into the "dat.pt" object #################
dat.pt<-rbind(cr.pt,g.y1)#this strings together the pre-treatment data from CoRRE and the control data from year 1 of GEx


###########################################################################
### Getting NutNet Data into shape ########################################
###########################################################################
nnab$site.plt<-with(nnab, paste(site_code,"_",block,"_",plot))
nnab$site.yr<-with(nnab, paste(site_code,"_",year))
nnab$site.plt.yr<-with(nnab, paste(site.plt,"_",year))
nnabc<-nnab[nnab$trt=="Control",]#pulls out just the control (unfertilized) plots
nnabpt<-nnab[nnab$year_trt==0,]#pulls out all plots in the pre-treatment (before fertilization) year
plots<-data.frame("site.plt.yr"=sort(unique(nnabpt$site.plt.yr)))
fixab<-data.frame("site.plt.yr"=with(nnabpt[nnabpt$N_fixer==1,],sort(unique(site.plt.yr))),
  "fixab"<-with(nnabpt[nnabpt$N_fixer==1,], tapply(max_cover,site.plt.yr,sum)))

nn<-merge(plots,fixab, by="site.plt.yr",all.x=T,all.y=T)
colnames(nn)<-c("site.plt.yr","fixab")
nn$fixab<-ifelse(is.na(nn$fixab),yes=0,no=nn$fixab)

#calculating species richness for fixers and non-fixers in each plot
nnplots<-unique(nnabpt$site.plt.yr)
nnrich<-NULL
for(p in 1:length(nnplots)){
  print(p)
  temp<-nnabpt[nnabpt$site.plt.yr==nnplots[p],]
  tempdf<-data.frame("site.plt.yr"=temp$site.plt.yr,
                     "tot.rich"=length(unique(temp$Taxon)),
                     "fix.rich"=length(unique(temp[temp$N_fixer==1,]$Taxon)))
  nnrich<-rbind(nnrich,tempdf)
}
nnrich$fix.rr<-with(nnrich, (fix.rich/tot.rich)*100)
nnrich<-nnrich[!duplicated(nnrich),]

nnplt$site.plt.yr<-with(nnplt, paste(site_code,"_",block,"_",plot,"_",year))
nnmrg<-nnplt[,c(91,1,3:8,45,13:17,38,29,32,18,21,39,52,56,63,64,67:90)]

nn<-merge(nn,nnmrg,by="site.plt.yr",all.x=T,all.y=F)
nn<-merge(nn, nnrich, by="site.plt.yr", all.x=T, all.y=F)
nn$fixrpc<-with(nn, (fixab/total_cover)*100)#calculates the percent cover of N fixers relativized to the total percent cover (often >100%) of the plot
nn$fixrpc<-ifelse(nn$fixab==0, yes=0, no=nn$fixrpc)#several plots are NA for fixer relative % cover because total % cover is NA, 
#so here I'm assigning 0's for the plots that have no fixers (we know it's 0% fixer relative % abundance regardless of what total % abundance is)
nn<-nn[!duplicated(nn$site.plt.yr),]
nn<-nn[!is.na(nn$site_code),]

## Calculating N limitation for each site ###
nnabtrt<-nnab[nnab$year_trt==1,]
sites<-unique(nnabtrt$site_code)
sitelim<-NULL
for(s in 1:length(sites)){
  print(s)
  temp1<-nnab[nnab$site_code==sites[s]&nnab$year_trt==1&nnab$live==1,]#makes a temporary dataframe with an individual site and just the 1st year of trt data
  tempcon<-temp1[temp1$trt=="Control",]
  tempconK<-temp1[temp1$trt%in%c("Control","K"),]
  tempN<-temp1[temp1$trt=="N",]
  tempNK<-temp1[temp1$trt%in%c("N","NK"),]
  tempP<-temp1[temp1$trt=="P",]
  tempPK<-temp1[temp1$trt%in%c("P","PK"),]
  limdat<-as.data.frame(unique(temp1$site_code))
  colnames(limdat)<-"site_code"
  limdat$ctrl.cov<-(sum(tempcon$max_cover)/length(unique(tempcon$plot)))
  limdat$ctrlK.cov<-(sum(tempconK$max_cover)/length(unique(tempconK$plot)))
  limdat$N.cov<-(sum(tempN$max_cover)/length(unique(tempN$plot)))
  limdat$NK.cov<-(sum(tempNK$max_cover)/length(unique(tempNK$plot)))
  limdat$P.cov<-(sum(tempP$max_cover)/length(unique(tempP$plot)))
  limdat$PK.cov<-(sum(tempPK$max_cover)/length(unique(tempPK$plot)))
  sitelim<-rbind(sitelim,limdat)
}
sitelim$NRR<-with(sitelim, ((N.cov-ctrl.cov)/ctrl.cov)*100)
sitelim$NKRR<-with(sitelim, ((NK.cov-ctrlK.cov)/ctrlK.cov)*100)
sitelim$PRR<-with(sitelim, ((P.cov-ctrl.cov)/ctrl.cov)*100)
sitelim$PKRR<-with(sitelim, ((PK.cov-ctrlK.cov)/ctrlK.cov)*100)
sitelim<-sitelim[!is.na(sitelim$NRR),]
sitelim<-sitelim[!is.na(sitelim$PRR),]
#sitelim$msRR<-with(sitelim, ((N.mass-ctrl.mass)/ctrl.mass)*100)

sitelim.mrg<-sitelim[,c(1,8,10)]
colnames(sitelim.mrg)<-c("site_code","Nlim","Plim")
nn<-merge(nn,sitelim.mrg,by="site_code",all.x=T,all.y=F)


###########################################################################
### Combining NutNet & CoRRE data #########################################
###########################################################################
nncomb<-nn[,c(2,1,5,22,3,53,13,14,19,17,54,55,50,51,52)]
nncomb$dataset<-"NutNet"
nncomb$abs.LAT<-abs(nncomb$latitude)
colnames(nncomb)<-clnms


dat2<-rbind(dat,nncomb)
dat2.pt<-rbind(dat.pt,nncomb)

dat.NAM<-dat[dat$LON>-180&dat$LON< -50&dat$LAT>0,]
dat2.NAM<-dat2[dat2$LON>-180&dat2$LON< -50&dat2$LAT>0,]

datpt.NAM<-dat.pt[dat.pt$LON>-180&dat.pt$LON< -50&dat.pt$LAT>0,]
dat2pt.NAM<-dat2.pt[dat2.pt$LON>-180&dat2.pt$LON< -50&dat2.pt$LAT>0,]

###########################################################################
### Making a Site-Level data set #########################################
###########################################################################
stfxra<-with(dat.pt, tapply(fixra, site_code, mean, na.rm=T))
dat.pt<-droplevels(dat.pt)
sitedat<-data.frame("site_code"=sort(unique(dat.pt$site_code)),
                    "fixra"=with(dat.pt, tapply(fixra, site_code, mean, na.rm=T)),
                    "fixrr"=with(dat.pt, tapply(fixrr, site_code, mean, na.rm=T)),
                    "MAT"=with(dat.pt, tapply(MAT, site_code, mean, na.rm=T)),
                    "abs.LAT"=with(dat.pt, tapply(abs.LAT, site_code, mean, na.rm=T)),
                    "Nlim"=with(dat.pt, tapply(Nlim, site_code, mean, na.rm=T))
                    )

stfxra2<-with(dat2.pt, tapply(fixra, site_code, mean, na.rm=T))
dat2.pt<-droplevels(dat2.pt)
sitedat2<-data.frame("site_code"=sort(unique(dat2.pt$site_code)),
                    "fixra"=with(dat2.pt, tapply(fixra, site_code, mean, na.rm=T)),
                    "fixrr"=with(dat2.pt, tapply(fixrr, site_code, mean, na.rm=T)),
                    "MAT"=with(dat2.pt, tapply(MAT, site_code, mean, na.rm=T)),
                    "abs.LAT"=with(dat2.pt, tapply(abs.LAT, site_code, mean, na.rm=T)),
                    "Nlim"=with(dat2.pt, tapply(Nlim, site_code, mean, na.rm=T))
)

###########################################################################
### What is the Data Distribution of N-fixer Relative Abundance? ##########
###########################################################################
hist(dat2.pt$fixra)#looks 0-inflated lognormal
hist(dat2.pt[dat2.pt$fixra>0,]$fixra)#Definitely 0-inflated
hist(log(dat2.pt[dat2.pt$fixra>0,]$fixra))#And the non-zero data are lognormal. So yes, 0-inflated lognormal

###########################################################################
### Writing data frames to files to use in analyses and figures ###########
###########################################################################

### Combined Data Files ############################################
write.csv(dat2, file="Processed Grassland Data_3 Datasets_All Years.csv", row.names=F)
write.csv(dat2.pt, file="Processed Grassland Data_3 Datasets_Pretreatment Year.csv", row.names=F)
write.csv(dat, file="Processed Grassland Data_CoRRE and GEx_All Years.csv", row.names=F)
write.csv(dat.pt, file="Processed Grassland Data_CoRRE and GEx_Pretreatment Year.csv", row.names=F)

### Individual Dataset Files ########################################
write.csv(cr, file="Processed CoRRE Data_All Years.csv", row.names=F)
write.csv(cr.pt, file="Processed CoRRE Data_Pretreatment Year.csv", row.names=F)
write.csv(g, file="Processed GEx Data_All Years.csv", row.names=F)
write.csv(g.y1, file="Processed NutNet Data_First Year Only.csv", row.names=F)
write.csv(nncomb, file="Processed NutNet Data_Pretreatment Year.csv", row.names=F)

### North American Data Files #######################################
write.csv(dat.NAM, file="Processed North American Data_CoRRE and GEX_All Years.csv", row.names=F)
write.csv(dat2.NAM, file="Processed North American Data_3 Datasets_All Years.csv", row.names=F)
write.csv(datpt.NAM, file="Processed North American Data_CoRRE and GEX_All Years.csv", row.names=F)
write.csv(dat2pt.NAM, file="Processed North American Data_CoRRE and GEX_Pretreatment Year.csv", row.names=F)


###########################################################################
### What is the latitudinal pattern of N-fixer Relative Abundance? ########
###########################################################################
with(dat, summary(lm(fixra~abs.LAT)))
with(dat, summary(lm(fixra~abs.LAT+site_code)))
with(dat, summary(lmer(fixra~abs.LAT+(1|site_code))))

with(sitedat, summary(lm(log(fixra+1)~abs.LAT)))

### What is the latitudinal pattern for North American N-fixers? ##########
with(datNAM, summary(lm(fixra~abs.LAT)))
with(datNAM, summary(lm(fixra~abs.LAT+site_code)))
with(datNAM, summary(lmer(fixra~abs.LAT+(1|site_code))))

###########################################################################
### What is the latitudinal pattern of N-fixer Relative Abundance? ########
### INCLUDING NUTNET DATA #################################################
###########################################################################
with(dat2, summary(lm(fixra~abs.LAT)))
with(dat2, summary(lm(fixra~abs.LAT+site_code)))
with(dat2, summary(lmer(fixra~abs.LAT+(1|site_code))))

with(sitedat2, summary(lm(log(fixra+1)~abs.LAT)))

### What is the latitudinal pattern for North American N-fixers? ##########
with(dat2NAM, summary(lm(fixra~abs.LAT)))
with(dat2NAM, summary(lm(fixra~abs.LAT+site_code)))
with(dat2NAM, summary(lmer(fixra~abs.LAT+(1|site_code))))

###########################################################################
### What is the latitudinal pattern of N-fixer Relative Richness? #########
###########################################################################
with(dat, summary(lm(fixrr~abs.LAT)))
with(dat, summary(lm(fixrr~abs.LAT+site_code)))
with(dat, summary(lmer(fixrr~abs.LAT+(1|site_code))))

### What is the latitudinal pattern for North American N-fixers? ##########
with(datNAM, summary(lm(fixrr~abs.LAT)))
with(datNAM, summary(lm(fixrr~abs.LAT+site_code)))
with(datNAM, summary(lmer(fixrr~abs.LAT+(1|site_code))))

###########################################################################
### What is the latitudinal pattern of N-fixer Relative Richness? #########
### INCLUDING NUTNET DATA #################################################
###########################################################################
with(dat2, summary(lm(fixrr~abs.LAT)))
with(dat2, summary(lm(fixrr~abs.LAT+site_code)))
with(dat2, summary(lmer(fixrr~abs.LAT+(1|site_code))))

### What is the latitudinal pattern for North American N-fixers? ##########
with(dat2NAM, summary(lm(fixrr~abs.LAT)))
with(dat2NAM, summary(lm(fixrr~abs.LAT+site_code)))
with(dat2NAM, summary(lmer(fixrr~abs.LAT+(1|site_code))))

###########################################################################
### Are N-fixer Relative Abundance and Relative Richness Correlated? ######
###########################################################################
with(dat, summary(lm(fixra~fixrr))) #Yes, but R2 is only .33
with(datNAM, summary(lm(fixra~fixrr+site_code)))

###########################################################################
### Does N-fixer Relative Abundance and Richness Increase with MAT? #######
###########################################################################
#relative abundance
with(dat, summary(lm(fixra~MAT+site_code)))
with(datNAM, summary(lm(fixra~MAT*LAT)))

#relative richness
with(dat, summary(lm(fixrr~MAT+site_code)))
with(datNAM, summary(lm(fixrr~MAT+site_code)))

###########################################################################
### Does N-fixer Relative Abundance increase with site N limitation? ######
###########################################################################
#relative abundance
with(dat[!is.na(dat$Nlim),], summary(lm(fixra~Nlim)))
with(dat[!is.na(dat$Nlim),], summary(lm(fixra~Nlim+site_code)))
with(dat[!is.na(dat$Nlim),], summary(lmer(fixra~Nlim+(1|site_code))))

with(datNAM[!is.na(datNAM$Nlim),], summary(lm(fixra~Nlim+site_code)))

#relative richness
with(dat[!is.na(dat$Nlim),], summary(lm(fixrr~Nlim+site_code)))
with(datNAM[!is.na(datNAM$Nlim),], summary(lm(fixrr~Nlim+site_code)))

          
###########################################################################################
###########################################################################################
### FIGURES ###############################################################################
###########################################################################################
###########################################################################################

###########################################################################
### Mapping the data distribution #########################################
###########################################################################
#All the datapoints together
#png(filename = "Map of Data Distribution_all.png", width=8, height=4, units="in", res=300)
map('world')
with(dat2, points(x=LON, y=LAT, bg='forestgreen', pch=21, cex=1.2))
#dev.off()

#Now with GEx and CoRRE split out
png(filename = "Map of Data Distribution_3 datasets.png", width=14, height=7, units="in", res=300)
map('world')
with(dat2[dat2$dataset=="NutNet",], points(x=LON, y=LAT, bg='forestgreen', pch=21, cex=1.2))
with(dat2[dat2$dataset=="GEx",], points(x=LON, y=LAT, bg='blue', pch=21, cex=1.2))
with(dat2[dat2$dataset=="CoRRE",], points(x=LON, y=LAT, bg='red', pch=21, cex=1.2))
legend(x=43, y=-28,legend=c("NutNet","CoRRE", "GEx"), pt.bg=c("forestgreen","Red","Blue"), pch=21, bty="n", cex=1.2)
dev.off()

###########################################################################
### Fixer Abundance by Latitude ###########################################
###########################################################################
dat$latbins<-cut(dat$abs.LAT,c(seq(30,75,1)),labels=c(seq(30.5,75,1)))
latbndat<-data.frame("LAT"=seq(30.5,75,1),
                     "fixra"=with(dat, tapply(fixra,latbins,mean, na.rm=T)))

#png(filename = "CoRRE N-Fixer Relative Abundance.png", width=8, height=6, units="in", res=300)
latfull<-ggplot(dat, aes(x=abs.LAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="red")+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=30,xend=30,y=0,yend=70),colour="black")+
  geom_segment(aes(x=30,xend=75,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))
latfull
#dev.off()
#Mean N-fixer relative % cover ranges from ~1% in the tropics up to ~13% at 49 degrees lat.

###########################################################################
### Fixer Abundance by Latitude in North America ##########################
###########################################################################
datNAM$latbins<-cut(datNAM$abs.LAT,c(seq(27,65,1)),labels=c(seq(27.5,64.5,1)))
NAMlatbndat<-data.frame("LAT"=seq(27.5,64.5,1),
                     "fixra"=with(datNAM, tapply(fixra,latbins,median, na.rm=T)))

#png(filename = "North American N-Fixer Relative Abundance.png", width=6, height=6, units="in", res=300)
latNA<-ggplot(datNAM, aes(x=abs.LAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=NAMlatbndat, mapping=aes(x=LAT,y=fixra),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=25,xend=25,y=0,yend=70),colour="black")+
  geom_segment(aes(x=25,xend=65,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("North America")+
  theme(plot.title=element_text(hjust=.5))
latNA
#dev.off()

png(filename = "North American N-Fixer Relative Abundance.png", width=6, height=10, units="in", res=300)
multiplot(latfull, latNA, cols = 1)
dev.off()
###########################################################################
### NutNet Fixer Abundance by Latitude ####################################
###########################################################################
nn$latbins<-cut(abs(nn$latitude),c(seq(0,80,1)),labels=c(seq(1,80,1)))
latbndat<-data.frame("latitude"=seq(1,80,1),
                     "fixrpc"=with(nn, tapply(fixrpc,latbins,median, na.rm=T)))

#png(filename = "NutNet Fixer Relative Abundance.png", width=8, height=6, units="in", res=300)
ggplot(nn, aes(x=abs(latitude),y=fixrpc))+
  geom_point(size=2,colour="grey")+
  geom_point(data=latbndat, mapping=aes(x=latitude,y=fixrpc),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Percent Cover (%)")+
  ggtitle("NutNet N-fixer Relative Percent Cover by Latitude")
#dev.off()
#Mean N-fixer relative % cover ranges from ~1% in the tropics up to ~13% at 49 degrees lat.

###########################################################################
### CoRRE Fixer Abundance by Latitude #####################################
###########################################################################
cr$latbins<-cut(abs(cr$latitude),c(seq(34,54,1)),labels=c(seq(35,54,1)))
latbndat<-data.frame("latitude"=seq(35,54,1),
                     "fixra"=with(cr, tapply(fixra,latbins,median, na.rm=T)))

#png(filename = "CoRRE Fixer Abundance.png", width=8, height=6, units="in", res=300)
ggplot(cr, aes(x=abs(latitude),y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=latbndat, mapping=aes(x=latitude,y=fixra),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=34,xend=34,y=0,yend=70),colour="black")+
  geom_segment(aes(x=34,xend=55,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("CoRRE N-fixer Abundance by Latitude")
#dev.off()

###########################################################################
### Fixer Richness by Fixer Abundance #####################################
###########################################################################
rva<-ggplot(dat, aes(x=fixrr, y=fixra))+
  geom_point(size=2, colour="grey")+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black")+
  geom_segment(aes(x=0,xend=70,y=0,yend=0),colour="black")+
  xlab(expression("N-Fixer Relative Richness (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("N-fixer Richness vs. Abundance")+
  theme(plot.title=element_text(hjust=.5))
###########################################################################
### Fixer Richness by Latitude ###########################################
###########################################################################
dat$latbins<-cut(dat$abs.LAT,c(seq(0,80,1)),labels=c(seq(1,80,1)))
latbndat<-data.frame("LAT"=seq(1,80,1),
                     "fixrr"=with(dat, tapply(fixrr,latbins,median, na.rm=T)))

#png(filename = "N-Fixer Relative Richness.png", width=8, height=6, units="in", res=300)
richfull<-ggplot(dat, aes(x=abs.LAT,y=fixrr))+
  geom_point(size=2,colour="grey")+
  geom_point(data=latbndat, mapping=aes(x=LAT,y=fixrr),size=4)+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Richness (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))
#dev.off()
#Mean N-fixer relative % cover ranges from ~1% in the tropics up to ~13% at 49 degrees lat.

###########################################################################
### Fixer Richness by Latitude in North America ##########################
###########################################################################
datNAM$latbins<-cut(datNAM$abs.LAT,c(seq(27,49,1)),labels=c(seq(28,49,1)))
NAMlatbndat<-data.frame("LAT"=seq(28,49,1),
                        "fixrr"=with(datNAM, tapply(fixrr,latbins,median, na.rm=T)))

#png(filename = "North American N-Fixer Relative Richness.png", width=6, height=6, units="in", res=300)
richNA<-ggplot(datNAM, aes(x=abs.LAT,y=fixrr))+
  geom_point(size=2,colour="grey")+
  geom_point(data=NAMlatbndat, mapping=aes(x=LAT,y=fixrr),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=25,xend=25,y=0,yend=70),colour="black")+
  geom_segment(aes(x=25,xend=50,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Richness (%)")+
  ggtitle("North America")+
  theme(plot.title=element_text(hjust=.5))
#dev.off()

png(filename = "Richness vs. Abundance.png", width=6, height=12, units="in", res=300)
multiplot(rva,richfull,richNA, cols=1)
dev.off()

###########################################################################
### Fixer Abundance by Temperature ###########################################
###########################################################################
dat$tempbins<-cut(dat$MAT,c(seq(-8,28,1)),labels=c(seq(-7,28,1)))
tempbndat<-data.frame("MAT"=seq(-7,28,1),
                     "fixra"=with(dat, tapply(fixra,tempbins,median, na.rm=T)))

tempfull<-ggplot(dat, aes(x=MAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=tempbndat, mapping=aes(x=MAT,y=fixra),size=4)+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-8,xend=-8,y=0,yend=80),colour="black")+
  geom_segment(aes(x=-8,xend=30,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))

###########################################################################
### Fixer Richness by Latitude in North America ##########################
###########################################################################
datNAM$tempbins<-cut(datNAM$MAT,c(seq(-1,23,1)),labels=c(seq(0,23,1)))
tempbndat<-data.frame("MAT"=seq(0,23,1),
                      "fixra"=with(datNAM, tapply(fixra,tempbins,median, na.rm=T)))

tempNA<-ggplot(datNAM, aes(x=MAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=tempbndat, mapping=aes(x=MAT,y=fixra),size=4)+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-1,xend=-1,y=0,yend=80),colour="black")+
  geom_segment(aes(x=-1,xend=25,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("North America")+
  theme(plot.title=element_text(hjust=.5))

png(filename = "Abundance vs. Temperature.png", width=6, height=10, units="in", res=300)
multiplot(tempfull,tempNA, cols=1)
dev.off()

###########################################################################
### Fixer Abundance by N limitation ###########################################
###########################################################################
nlimdat<-dat[!is.na(dat$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-78,160,1)),labels=c(seq(-77,160,1)))
limbndat<-data.frame("Nlim"=seq(-77,160,1),
                     "fixra"=with(nlimdat, tapply(fixra,bins,median, na.rm=T)))

limfull<-ggplot(nlimdat, aes(x=Nlim,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=limbndat, mapping=aes(x=Nlim,y=fixra),size=4)+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-80,xend=-80,y=0,yend=60),colour="black")+
  geom_segment(aes(x=-80,xend=160,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))

###########################################################################
### Fixer Abundance by N limitation in North America ##########################
###########################################################################
nlimdat<-datNAM[!is.na(datNAM$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-38,160,1)),labels=c(seq(-37,160,1)))
limbndat<-data.frame("Nlim"=seq(-37,160,1),
                     "fixra"=with(nlimdat, tapply(fixra,bins,median, na.rm=T)))

limNA<-ggplot(nlimdat, aes(x=Nlim,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=limbndat, mapping=aes(x=Nlim,y=fixra),size=4)+
  geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-40,xend=-40,y=0,yend=60),colour="black")+
  geom_segment(aes(x=-40,xend=160,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("North America")+
  theme(plot.title=element_text(hjust=.5))

png(filename = "Abundance vs. N limitation.png", width=6, height=10, units="in", res=300)
multiplot(limfull,limNA, cols=1)
dev.off()


###########################################################################
### Fixer Richness by N limitation ###########################################
###########################################################################
nlimdat<-dat[!is.na(dat$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-78,160,1)),labels=c(seq(-77,160,1)))
limbndat<-data.frame("Nlim"=seq(-77,160,1),
                     "fixrr"=with(nlimdat, tapply(fixrr,bins,median, na.rm=T)))

#png(filename = "N-Fixer Relative Richness by N Limitation.png", width=8, height=6, units="in", res=300)
ggplot(nlimdat, aes(x=Nlim,y=fixrr))+
  geom_point(size=2,colour="grey")+
  geom_point(data=limbndat, mapping=aes(x=Nlim,y=fixrr),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-80,xend=-80,y=0,yend=70),colour="black")+
  geom_segment(aes(x=-80,xend=160,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Richness (%)")+
  ggtitle("N-fixer Richness by N Limitation")
#dev.off()

###########################################################################
### Fixer Richness by N limitation in North America ##########################
###########################################################################
nlimdat<-datNAM[!is.na(datNAM$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-38,160,1)),labels=c(seq(-37,160,1)))
limbndat<-data.frame("Nlim"=seq(-37,160,1),
                     "fixrr"=with(nlimdat, tapply(fixrr,bins,median, na.rm=T)))

png(filename = "North American N-Fixer Relative Richness by N Limitation.png", width=8, height=6, units="in", res=300)
ggplot(nlimdat, aes(x=Nlim,y=fixrr))+
  geom_point(size=2,colour="grey")+
  geom_point(data=limbndat, mapping=aes(x=Nlim,y=fixrr),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-40,xend=-40,y=0,yend=70),colour="black")+
  geom_segment(aes(x=-40,xend=160,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=70),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Richness (%)")+
  ggtitle("North American N-fixer Abundance by N Limitation")
#dev.off()



###########################################################################
### NutNet Fixer Abundance by Soil %N ####################################
###########################################################################
nn$nbins<-cut(nn$pct_N,c(seq(0,1.55,.1)),labels=c(seq(.1,1.55,.1)))
nbndat<-data.frame("pct_N"=seq(.1,1.55,.1),
                     "fixrpc"=with(nn, tapply(fixrpc,nbins,mean, na.rm=T)))

#png(filename = "NutNet Fixer Abundance by Soil N.png", width=8, height=6, units="in", res=300)
ggplot(nn, aes(x=pct_N,y=fixrpc))+
  geom_point(size=2,colour="grey")+
  geom_point(data=nbndat, mapping=aes(x=pct_N,y=fixrpc),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black")+
  geom_segment(aes(x=0,xend=1.6,y=0,yend=0),colour="black")+
  xlab(expression("Soil N (%)"))+
  ylab("N-Fixer Percent Cover (%)")+
  ggtitle("NutNet N-fixer Percent Cover by Soil %N")
#dev.off()

###########################################################################
### NutNet Fixer Abundance by N Deposition ####################################
###########################################################################
nn$depbins<-cut(nn$N_Dep,c(seq(0,20,1)),labels=c(seq(1,20,1)))
depbndat<-data.frame("N_Dep"=seq(1,20,1),
                     "fixrpc"=with(nn, tapply(fixrpc,depbins,mean, na.rm=T)))

#png(filename = "NutNet Fixer Abundance by N Deposition.png", width=8, height=6, units="in", res=300)
ggplot(nn[nn$N_Dep<30,], aes(x=N_Dep,y=fixrpc))+
  geom_point(size=2,colour="grey")+
  geom_point(data=depbndat, mapping=aes(x=N_Dep,y=fixrpc),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=75),colour="black")+
  geom_segment(aes(x=0,xend=20,y=0,yend=0),colour="black")+
  xlab(expression("N Deposition"~(kg~N~ha^{-1}~yr^{-1})))+
  ylab("N-Fixer Percent Cover (%)")+
  ggtitle("NutNet N-fixer Percent Cover by N Deposition")
#dev.off()

###########################################################################
### NutNet Fixer Abundance by Site N Limitation ###########################
###########################################################################
nn$limbins<-cut(nn$RR,c(seq(-31,134,1)),labels=c(seq(-30,134,1)))
nlimdat<-data.frame("N_Lim"=seq(-30,134,1),
                     "fixrpc"=with(nn, tapply(fixrpc,limbins,mean, na.rm=T)))

#png(filename = "NutNet Fixer Abundance by N Limitation.png", width=8, height=6, units="in", res=300)
ggplot(nn, aes(x=RR,y=fixrpc))+
  geom_point(size=2,colour="grey")+
  geom_point(data=nlimdat, mapping=aes(x=N_Lim,y=fixrpc),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-30,xend=-30,y=0,yend=60),colour="black")+
  geom_segment(aes(x=-30,xend=130,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black",linetype="dashed")+
  xlab(expression("Site N Limitation (% Response to N Fert.)"))+
  ylab("N-Fixer Percent Cover (%)")+
  ggtitle("NutNet N-fixer Percent Cover by N Limitation")
#dev.off()

###########################################################################
### NutNet Fixer Abundance by Light ####################################
###########################################################################
nn$parbins<-cut(nn$Ambient_PAR,c(seq(0,3000,100)),labels=c(seq(1,3000,100)))
parbndat<-data.frame("PAR"=seq(1,3000,100),
                     "fixab"=with(nn, tapply(fixab,parbins,mean, na.rm=T)))

#png(filename = "NutNet Fixer Abundance by PAR.png", width=8, height=6, units="in", res=300)
ggplot(nn, aes(x=Ambient_PAR,y=fixab))+
  geom_point(size=2,colour="grey")+
  geom_smooth(method="lm",se=T)+
  geom_point(data=parbndat, mapping=aes(x=PAR,y=fixab),size=4)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=100),colour="black")+
  geom_segment(aes(x=0,xend=3000,y=0,yend=0),colour="black")+
  xlab(expression("PAR (mmol)"))+
  ylab("N-Fixer Percent Cover (%)")+
  ggtitle("NutNet N-fixer Percent Cover by Light")
#dev.off()