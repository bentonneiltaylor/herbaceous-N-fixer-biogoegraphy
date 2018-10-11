###HERBACEOUS LEGUME BIOGEOGRAPHY###
#Author: Ben Taylor
#bentonneiltaylor@gmail.com

###These analyses aim to understand the broad biogeographic (primarily latitudinal)
###patterns of herbaceous N-fixing legumes. To do this, we analyze a large
###network of broadly distributed grassland experiments from the NutNet and CoRRE networks

###########################################################################
### Read in the data and load R packages ##################################
setwd("C:/Users/Benton/Desktop/Work Files/Research/Herbaceous Legume Biogeography/Data Analysis & Figure Files")
nnplt<-read.csv("NutNet_plot_data_31-Aug-2018.csv")#reads in NutNet plot data
nnab<-read.csv("NutNet_raw_abundance_31-Aug-2018.csv")#reads in NutNet abundance data
cplt<-read.csv("CORRE_site_summary.csv")#reads in CoRRE plot data
cab<-read.csv("CORRE_raw_abundance.csv")#reads in CoRRE abundance data
ctrt<-read.csv("CORRE_treatment_summary.csv")#reads in CoRRE treatment data
nfsp<-read.csv("N fixer species list_Duncan.csv")#reads in a species list of confirmed N-fixers

library(ggplot2);library(bbmle);library(vegan);library(maps);library(mapdata)#loads necessary libraries
source("C:\\RCode\\Homegrown Functions.R")
###########################################################################
### Getting NutNet Data into shape ########################################
###########################################################################
nnab$site.plt<-with(nnab, paste(site_code,"_",plot))
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

nnplt$site.plt.yr<-with(nnplt, paste(site_code,"_",plot,"_",year))
nnmrg<-nnplt[,c(91,1,3:8,45,13:17,38,29,32,18,21,39,52,56,63,64,67:90)]

nn<-merge(nn,nnmrg,by="site.plt.yr",all.x=T,all.y=F)
nn$fixrpc<-with(nn, (fixab/total_cover)*100)#calculates the percent cover of N fixers relativized to the total percent cover (often >100%) of the plot
nn$fixrpc<-ifelse(nn$fixab==0, yes=0, no=nn$fixrpc)#several plots are NA for fixer relative % cover because total % cover is NA, 
#so here I'm assigning 0's for the plots that have no fixers (we know it's 0% fixer relative % abundance regardless of what total % abundance is)
nn<-nn[!duplicated(nn$site.plt.yr),]
nn<-nn[!is.na(nn$site_code),]

### Calculating N limitation for each site ###
finyr<-as.data.frame(aggregate(nnab$year, by=list(nnab$site_code),max))
colnames(finyr)<-c("site_code","year")
sites<-finyr$site_code
sitelim<-NULL
for(s in 1:length(sites)){
  print(s)
  temp1<-nnplt[nnplt$site_code==sites[s],]
  temp2<-temp1[temp1$year==max(temp1$year),]
  limdat<-as.data.frame(unique(temp2$site_code))
  limdat$ctrl.cov<-mean(temp2[temp2$trt=="Control",]$total_cover)
  limdat$ctrl.se<-sefun(temp2[temp2$trt=="Control",]$total_cover)
  limdat$N.cov<-mean(temp2[temp2$trt=="N",]$total_cover)
  limdat$N.se<-sefun(temp2[temp2$trt=="N",]$total_cover)
  limdat$ctrl.mass<-mean(temp2[temp2$trt=="Control",]$total_mass)
  limdat$ctrl.ms.se<-sefun(temp2[temp2$trt=="Control",]$total_mass)
  limdat$N.mass<-mean(temp2[temp2$trt=="N",]$total_mass)
  limdat$N.ms.se<-sefun(temp2[temp2$trt=="N",]$total_mass)
  sitelim<-rbind(sitelim,limdat)
}
sitelim$covRR<-with(sitelim, ((N.cov-ctrl.cov)/ctrl.cov)*100)
sitelim$msRR<-with(sitelim, ((N.mass-ctrl.mass)/ctrl.mass)*100)
colnames(sitelim)[1]<-"site_code"
sitelim.mrg<-sitelim[,c(1,10,11)]
nn<-merge(nn,sitelim.mrg,by="site_code",all.x=T,all.y=F)

###########################################################################
### Getting CoRRE Data into shape #########################################
###########################################################################
cab$site.plt.yr<-with(cab, paste(site_code,"_",plot_id,"_",calendar_year))
cab$genus_species<-as.factor(paste(cab$genus,"_",cab$species))
cab<-merge(cab,cplt,by="site_code")
cab<-merge(cab,ctrt[,-1],by=c("site_code","project_name","calendar_year","treatment_year","treatment","community_type"),
           all.x=T,all.y=F)

nfsp$genus_species<-as.factor(paste(tolower(nfsp$Genus),"_",tolower(nfsp$Species)))
fixlst<-nfsp[nfsp$Nod.==1,c("genus_species","Nod.")]
colnames(fixlst)<-c("genus_species","N_fixer")
fixlst$genus_species<-as.factor(trimws(fixlst$genus_species))

cab<-merge(cab,fixlst,by="genus_species",all.x=T,all.y=F)
cab<-cab[!duplicated(cab$X),]
cab$N_fixer<-ifelse(is.na(cab$N_fixer),yes=0,no=cab$N_fixer)

cab<-cab[cab$treatment_year==0,]

cabplots<-cab[,c("site.plt.yr","site_code")]
cabplots<-cabplots[!duplicated(cabplots$site.plt.yr)&!is.na(cabplots$site.plt.yr),]
cabplots$totab<-with(cab, tapply(abundance,site.plt.yr,sum))
cabfixab<-data.frame("site.plt.yr"=with(cab[cab$N_fixer==1,],sort(unique(site.plt.yr))),
                     "fixab"<-with(cab[cab$N_fixer==1,], tapply(abundance,site.plt.yr,sum)))
cr<-merge(cabplots,cabfixab,by="site.plt.yr",all.x=T,all.y=T)
colnames(cr)<-c("site.plt.yr","site_code","totab","fixab")
cr$fixab<-ifelse(is.na(cr$fixab),yes=0,no=cr$fixab)
cr$fixrab<-with(cr, (fixab/totab)*100)
cr<-merge(cr,cplt,by="site_code")

###########################################################################
### Mapping the data distribution #########################################
###########################################################################
#png(filename = "Map of Data Distribution.png", width=8, height=4, units="in", res=300)
map('worldHires')
with(nn, points(x=longitude, y=latitude, col='blue', pch=16, cex=.8))
with(cr, points(x=longitude, y=latitude, col='red', pch=16, cex=.8))
legend(x=43, y=-28,legend=c("CoRRE", "NutNet"), col=c("Red","Blue"), pch=16, bty="n", cex=1)
#dev.off()

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
png(filename = "CoRRE Fixer Abundance.png", width=8, height=6, units="in", res=300)
ggplot(cr, aes(x=abs(latitude),y=fixab))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=latbndat, mapping=aes(x=latitude,y=fixab),size=4)+
  geom_smooth(method="lm",se=T)+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=20,xend=20,y=0,yend=700),colour="black")+
  geom_segment(aes(x=20,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Abundance (Units?)")+
  ggtitle("CoRRE N-fixer Abundance by Latitude")
dev.off()

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