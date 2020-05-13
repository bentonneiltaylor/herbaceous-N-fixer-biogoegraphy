### FIGURE CODE FOR GRASSLAND N-FIXER BIOGEOGRAPHY #####

#### Load packages and data files ####
library(ggplot2)
library(grid)
library(ggpubr)
library(maps)
library(mapdata)

datfull.pt<-read.csv("Processed Grassland Data_3 Datasets_Pretreatment Year.csv")
dun.trees<-read.csv("Menge 2017_tree data_BNT.csv")
sc<-read.csv("Gridcell-Level Control and Pretreatment Data_CoRRE and GEx.csv")

########################################################
#### HYPOTHETICAL FIGURE OF LAT PATTERNS ###############
########################################################
latseq<-seq(0,90,.1)
nlim.hyp.seq<-(0+.5*latseq)
ra.hyp.seq<-(0+.5*latseq)

hyp1<-ggplot()+
  geom_line(aes(x=nlim.hyp.seq,y=ra.hyp.seq), colour="black", size=2)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=50,y=0,yend=0),colour="black")+
  xlab(expression("Nitrogen Limitation"))+
  ylab("N-fixer Relative Abundance")

hyp2<-ggplot()+
  geom_line(aes(x=latseq,y=nlim.hyp.seq), colour="black", size=2)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=90,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("Nitrogen Limitation")

hyp3<-ggplot()+
  geom_line(aes(x=latseq,y=ra.hyp.seq), colour="black", size=2)+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=90,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-fixer Relative Abundance")

#png(filename="Hyp Relationships between Nlim Lat and Fixer RA.png",height=11,width=5,units="in",res=300)
ggarrange(hyp1,hyp2,hyp3,nrow=3,ncol=1)
#dev.off()

###########################################################################
### Mapping the data distribution #########################################
###########################################################################
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
world <- ne_countries(scale = "medium", returnclass = "sf")

#All the datapoints together
#png(filename = "Map of Data Distribution_all.png", width=8, height=4, units="in", res=300)
map('world')
with(sc, points(x=LON, y=LAT, bg='forestgreen', pch=21, cex=1.2))
#dev.off()

#Now with GEx and CoRRE split out
png(filename = "Map of Control Data Distribution_3 datasets.png", width=14, height=7, units="in", res=300)
map('world')
#with(sc[sc$dataset=="NutNet",], points(x=LON, y=LAT, bg='forestgreen', pch=21, cex=1.2))
with(sc[sc$dataset=="GEx",], points(x=LON, y=LAT, bg='blue', pch=21, cex=1.2))
with(sc[sc$dataset=="CoRRE",], points(x=LON, y=LAT, bg='red', pch=21, cex=1.2))
legend(x=43, y=-28,legend=c("CoRRE", "GEx"), pt.bg=c("Red","Blue"), pch=21, bty="n", cex=1.2)
#legend(x=43, y=-28,legend=c("NutNet","CoRRE", "GEx"), pt.bg=c("forestgreen","Red","Blue"), pch=21, bty="n", cex=1.2)
dev.off()

png(filename = "Map of Control Data Distribution_3 datasets_V2.png", width=14, height=7, units="in", res=300)
ggplot(data=world)+
  theme(panel.background=element_rect(fill="aliceblue", color="aliceblue"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  geom_sf(color="black", fill="antiquewhite")+
  #geom_point(data=sc[sc$dataset=="NutNet",], mapping=aes(x=LON,y=LAT),size=3,shape=21,fill='forestgreen')+
  #geom_point(data=sc[sc$dataset=="GEx",], mapping=aes(x=LON,y=LAT),size=3,shape=21,fill='blue')+
  geom_point(data=sc, mapping=aes(x=LON,y=LAT),size=3,shape=21,fill='forestgreen')+
  ylab(expression("Latitude "*degree*""))+
  xlab(expression("Longitude "*degree*""))
dev.off()

########################################################
#### GLOBAL LATITUDINAL PATTERN FOR ALL  DATASETS #####
########################################################
sc$latbins<-cut(sc$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
ziln.median.fun<-function(x){(1-(length(x[x==0])/length(x)))*exp(mean(log(x[x>0]),na.rm=T))}
latbndat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(sc, tapply(fixra,latbins,mean, na.rm=T)))
latbndat2<-data.frame("LAT"=seq(0.5,80,1),
                      "fixra"=with(sc, tapply(fixra,latbins,ziln.median.fun)))
#l.xseq<-seq(0,80,.1)
#l.yseq<-(1-(coef(abs.LAT_bestmod)[3]+(coef(abs.LAT_bestmod)[4]*l.xseq)))*exp(coef(abs.LAT_bestmod)[1])
#l.predcrv<-data.frame("x"=l.xseq,"y"=l.yseq)

latfull<-ggplot(sc, aes(x=abs.LAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="red")+
  geom_point(data=latbndat2, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="green")+
  #geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=35),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))
  #geom_line(aes(x=x,y=y), data=l.predcrv, colour="black")
png(filename = "3 Dataset Global N-Fixer Relative Abundance by Latitude.png", width=8, height=6, units="in", res=300)
latfull
dev.off()
#Mean N-fixer relative % cover ranges from ~1% in the tropics up to ~13% at 49 degrees lat.

########################################################
#### COMPARING TREE AND GRASSLAND PATTERNS #############
########################################################
dun.trees$latbins<-cut(dun.trees$LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
d.tr.dat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(dun.trees, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="Trees")
latbndat$data<-"Grasslands"
comp.dat<-rbind(d.tr.dat,latbndat)

comp.fig<-ggplot(comp.dat, aes(x=LAT,y=fixra))+
  geom_point(aes(shape=data, fill=data, size=5))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.key=element_blank(),
        legend.position=c(.8,1),legend.text=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  geom_segment(aes(x=0,xend=0,y=0,yend=30),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  scale_shape_manual(name="",labels=c("Trees","Grasslands"),
                     values=c(21,24))+
  scale_fill_manual(name="", labels=c("Trees","Grasslands"),
                    values=c("dodgerblue","forestgreen"))+
  geom_line(aes(x=x,y=y), data=l.predcrv, colour="forestgreen",size=1.3,)+
  geom_smooth(data=comp.dat[comp.dat$data=="Trees",],se=F,colour="dodgerblue",span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  coord_cartesian(xlim=c(-1,80),ylim=c(-1,30))

png(filename = "Trees vs. Grasslands Comparison Fig.png", width=8, height=6, units="in", res=300)
comp.fig
dev.off()

  
########################################################
#### N-FIXER ABUNDANCE VS. TEMPERATURE #################
########################################################
sc$tempbins<-cut(sc$MAT,c(seq(-12,28,1)),labels=c(seq(-11,28,1)))
tempbn.mn<-data.frame("MAT"=seq(-11,28,1),
                      "fixra"=with(sc, tapply(fixra,tempbins,mean, na.rm=T)))
tempbn.mdn<-data.frame("MAT"=seq(-11,28,1),
                      "fixra"=with(sc, tapply(fixra,tempbins,ziln.median.fun)))

t.xseq<-seq(-12,28,.1)
t.yseq<-(1-alogitfn(coef(MAT_bestmod)[3]+(coef(MAT_bestmod)[4]*t.xseq)))*exp(coef(MAT_bestmod)[1])
t.predcrv<-data.frame("x"=t.xseq,"y"=t.yseq)

tempfull<-ggplot(sc, aes(x=MAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=tempbn.mn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=tempbn.mdn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill="green")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-12,xend=-12,y=0,yend=30),colour="black")+
  geom_segment(aes(x=-12,xend=30,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=t.predcrv, colour="black")+
  coord_cartesian(xlim=c(-12,30),ylim=c(0,30),expand=F)

png(filename = "3-dataset N-fixer Abundance vs. Temperature.png", width=8, height=6, units="in", res=300)
tempfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. N LIMITATION ################
########################################################
nlimdat<-sc[!is.na(sc$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-25,85,1)),labels=c(seq(-24.5,85,1)))
limbn.mn<-data.frame("Nlim"=seq(-24.5,85,1),
                     "fixra"=with(nlimdat, tapply(fixra,bins,mean, na.rm=T)))
limbn.mdn<-data.frame("Nlim"=seq(-24.5,85,1),
                     "fixra"=with(nlimdat, tapply(fixra,bins,ziln.median.fun)))

nlim.xseq<-seq(-25,85,.1)
nlim.yseq<-(1-(coef(N.lim_bestmod)[4])*exp(coef(N.lim_bestmod)[2]+(coef(N.lim_bestmod)[3]*nlim.xseq)))
nlim.predcrv<-data.frame("x"=nlim.xseq,"y"=nlim.yseq)
#predmedian<-(1-coef(N.lim_bestmod)[3])*exp(coef(N.lim_bestmod)[2])

limfull<-ggplot(nlimdat, aes(x=Nlim,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=limbn.mdn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="green")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-25,xend=-25,y=0,yend=20),colour="black")+
  geom_segment(aes(x=-25,xend=85,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=20),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  #geom_abline(slope=0,intercept=predmedian,colour="black")
  geom_line(aes(x=x,y=y), data=nlim.predcrv, colour="black")
#png(filename = "3-dataset N-fixer Abundance vs. N Limitation.png", width=8, height=6, units="in", res=300)
limfull
#dev.off()

limfull+
  coord_cartesian(xlim=c(-40,60),ylim=c(0,15))

########################################################
#### N-FIXER ABUNDANCE VS. P LIMITATION ################
########################################################
plimdat<-datfull.pt[!is.na(datfull.pt$Plim),]
plimdat$bins<-cut(plimdat$Plim,c(seq(-97,350,1)),labels=c(seq(-96.5,350,1)))
plimbn.mn<-data.frame("Plim"=seq(-96.5,350,1),
                     "fixra"=with(plimdat, tapply(fixra,bins,mean, na.rm=T)))
plimbn.mdn<-data.frame("Plim"=seq(-96.5,350,1),
                      "fixra"=with(plimdat, tapply(fixra,bins,ziln.median.fun)))

plim.xseq<-seq(-97,350,.1)
plim.yseq<-(1-(coef(Plim_bestmod)[5]+(coef(Plim_bestmod)[6]*plim.xseq)))*exp(coef(Plim_bestmod)[2]+(coef(Plim_bestmod)[3]*plim.xseq)+(coef(Plim_bestmod)[4]*plim.xseq^2))
plim.predcrv<-data.frame("x"=plim.xseq,"y"=plim.yseq)

Plimfull<-ggplot(plimdat, aes(x=Plim,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=plimbn.mn, mapping=aes(x=Plim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=plimbn.mdn, mapping=aes(x=Plim,y=fixra),size=4,shape=21,fill="green")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-100,xend=-100,y=0,yend=60),colour="black")+
  geom_segment(aes(x=-100,xend=110,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black",linetype="dashed")+
  xlab(expression("P Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=plim.predcrv, colour="black")+
  coord_cartesian(xlim=c(-100,110))
png(filename = "3-dataset N-fixer Abundance vs. P Limitation.png", width=8, height=6, units="in", res=300)
Plimfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. N-FIXER RICHNESS ############
########################################################

### !!! Currently Incorrect - fitted curve is artifically divided by 5 !!!

sc$rrbins<-cut(sc$fixrr,c(seq(0,100,1)),labels=c(seq(.5,100,1)))
rrbn.mn<-data.frame("fixrr"=seq(.5,100,1),
                     "fixra"=with(sc, tapply(fixra,rrbins,mean, na.rm=T)))
rrbn.mdn<-data.frame("fixrr"=seq(.5,100,1),
                      "fixra"=with(sc, tapply(fixra,rrbins,ziln.median.fun)))

rr.xseq<-seq(0,60,.1)
rr.yseq<-(1-alogitfn(coef(fixrr_bestmod)[5]+(coef(fixrr_bestmod)[6]*rr.xseq)))*exp(coef(fixrr_bestmod)[2]+(coef(fixrr_bestmod)[3]*rr.xseq)+(coef(fixrr_bestmod)[4]*rr.xseq^2))
rr.predcrv<-data.frame("x"=rr.xseq,"y"=rr.yseq)

rva<-ggplot(sc, aes(x=fixrr, y=fixra))+
  geom_point(size=2, colour="grey")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_point(data=rrbn.mn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=rrbn.mdn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill="green")+
  geom_segment(aes(x=0,xend=0,y=0,yend=75),colour="black")+
  geom_segment(aes(x=0,xend=70,y=0,yend=0),colour="black")+
  xlab(expression("N-Fixer Relative Richness (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("N-fixer Richness vs. Abundance")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=(y)), data=rr.predcrv, colour="black")+
  coord_cartesian(xlim=c(0,70),ylim=c(0,75))

png(filename = "3-dataset N-fixer Abundance vs. N-Fixer Relative Richness.png", width=8, height=6, units="in", res=300)
rva
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. PRECIPITATION ###############
########################################################
sc$precbins<-cut(sc$MAP,c(seq(45,2566,10)),labels=c(seq(47.5,2566,10)))
precbn.mn<-data.frame("MAP"=seq(47.5,2566,10),
                      "fixra"=with(sc, tapply(fixra,precbins,mean, na.rm=T)))
precbn.mdn<-data.frame("MAP"=seq(47.5,2566,10),
                       "fixra"=with(sc, tapply(fixra,precbins,ziln.median.fun)))

prec.xseq<-seq(45,2566,.1)
prec.yseq<-(1-(coef(MAP_bestmod)[5]+(coef(MAP_bestmod)[6]*prec.xseq)))*exp(coef(MAP_bestmod)[2]+(coef(MAP_bestmod)[3]*prec.xseq)+(coef(MAP_bestmod)[4]*prec.xseq^2))
prec.predcrv<-data.frame("x"=prec.xseq,"y"=prec.yseq)
predmedian<-(1-coef(MAP_bestmod)[3])*exp(coef(MAP_bestmod)[2])

precfull<-ggplot(sc, aes(x=MAP,y=fixra))+
  geom_point(size=2,colour="grey")+
  geom_point(data=precbn.mn, mapping=aes(x=MAP,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=precbn.mdn, mapping=aes(x=MAP,y=fixra),size=4,shape=21, fill="green")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=45,xend=45,y=0,yend=75),colour="black")+
  geom_segment(aes(x=45,xend=2566,y=0,yend=0),colour="black")+
  xlab(expression("MAP (mm)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  #geom_line(aes(x=x,y=(y)), data=prec.predcrv, colour="black")+
  geom_abline(slope=0, intercept=predmedian, colour="black")+
  coord_cartesian(xlim=c(45,2566),ylim=c(0,75),expand=F)

png(filename = "3-dataset N-fixer Abundance vs. Precipitation.png", width=8, height=6, units="in", res=300)
precfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. HERBIVORY ###################
########################################################

herb.mdn.plt<-ggplot(data=g.herbdf[g.herbdf$fix.ra>0,],aes(x=graz.trt,y=fix.ra,fill=graz.trt))+
  geom_boxplot(notch=F,fill=c("grey40","grey70"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"))+
  scale_x_discrete(label=c("Grazed","Ungrazed"))+
  annotate("text",x=1.5,y=25,label=c("n.s."),size=5)+
  xlab(expression("Grazer Presence"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_segment(aes(x=0.5,xend=0.5,y=0,yend=25),colour="black")+
  geom_segment(aes(x=0.5,xend=2.5,y=0,yend=0),colour="black")

png(filename = "Herbivory Effects on N-fixer Abundance.png", width=6, height=6, units="in", res=300)
herb.mdn.plt
dev.off()

########################################################
#### MULTI-PANEL FIG OF ECOLOGICAL DRIVERS #############
########################################################

lat.mdn.plt<-ggplot()+
  #geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="red")+
  geom_point(data=latbndat2, mapping=aes(x=LAT,y=fixra),size=5, shape=21, fill="grey")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=100),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=l.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(0,80),ylim=c(-1,18))

rich.mdn.plt<-ggplot()+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  #geom_point(data=rrbn.mn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=rrbn.mdn, mapping=aes(x=fixrr,y=fixra),size=5,shape=21,fill="grey")+
  geom_segment(aes(x=0,xend=0,y=0,yend=50),colour="black")+
  geom_segment(aes(x=0,xend=75,y=0,yend=0),colour="black")+
  xlab(expression("N-Fixer Relative Richness (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=(y)), data=rr.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-1,75),ylim=c(-1,50))

Nlim.mdn.plt<-ggplot()+
  #geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=limbn.mdn, mapping=aes(x=Nlim,y=fixra),size=5,shape=21,fill="grey")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-80,xend=-80,y=0,yend=60),colour="black")+
  geom_segment(aes(x=-80,xend=160,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=60),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=nlim.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-80,80),ylim=c(-1,11))

Plim.mdn.plt<-ggplot()+
  #geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=plimbn.mdn, mapping=aes(x=Plim,y=fixra),size=5,shape=21,fill="grey")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-100,xend=-100,y=0,yend=20),colour="black")+
  geom_segment(aes(x=-100,xend=110,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=20),colour="black",linetype="dashed")+
  xlab(expression("P Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=plim.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-100,110),ylim=c(0,20))

t.mdn.plt<-ggplot()+
  #geom_point(data=tempbn.mn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=tempbn.mdn, mapping=aes(x=MAT,y=fixra),size=5,shape=21, fill="grey")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="none",
        plot.margin = margin(1.5,.8,.8,.8,"cm"))+
  geom_segment(aes(x=-5,xend=-5,y=0,yend=15),colour="black")+
  geom_segment(aes(x=-5,xend=30,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=t.predcrv, colour="black", size=1.2)+
  coord_cartesian(xlim=c(-5,32),ylim=c(0,15),expand=F)

prec.mdn.plt<-ggplot()+
  geom_point(data=precbn.mdn, mapping=aes(x=MAP,y=fixra),size=5,shape=21, fill="grey")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="none",
        plot.margin = margin(1.5,.8,.8,.8,"cm"))+
  geom_segment(aes(x=45,xend=45,y=0,yend=30),colour="black")+
  geom_segment(aes(x=45,xend=2566,y=0,yend=0),colour="black")+
  xlab(expression("MAP (mm)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=prec.predcrv, colour="black", size=1.2)+
  coord_cartesian(xlim=c(45,2566),ylim=c(0,30),expand=F)

herb.mdn.plt<-ggplot(data=g.herbdf[g.herbdf$fix.ra>0,],aes(x=graz.trt,y=fix.ra,fill=graz.trt))+
  geom_boxplot(notch=F,fill=c("grey40","grey70"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"))+
  scale_x_discrete(label=c("Grazed","Ungrazed"))+
  annotate("text",x=1.5,y=25,label=c("n.s."),size=5)+
  xlab(expression("Grazer Presence"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_segment(aes(x=0.5,xend=0.5,y=0,yend=25),colour="black")+
  geom_segment(aes(x=0.5,xend=2.5,y=0,yend=0),colour="black")

png(filename = "Ecological Drivers of Grassland N fixers.png", width=10, height=12, units="in", res=300)
ggarrange(t.mdn.plt,prec.mdn.plt,Nlim.mdn.plt,Plim.mdn.plt,rich.mdn.plt,herb.mdn.plt,ncol=2,nrow=3)
dev.off()

########################################################
#### LATITUDINAL PATTERN OF N LIMITATION ###############
########################################################
Nlimlat.lm<-with(datfull.pt, lm(Nlim~abs.LAT))
Nlimlat<-ggplot(data=datfull.pt, aes(x=abs.LAT, y=Nlim))+
  geom_point(size=5, shape=21, fill="grey30")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=-80,yend=165),colour="black")+
  geom_segment(aes(x=0,xend=77,y=-80,yend=-80),colour="black")+
  geom_segment(aes(x=0,xend=77,y=0,yend=0),colour="black", linetype="dashed")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab(expression("N Limitation RR (%)"))+
  geom_smooth(method="lm",color="black")+
  coord_cartesian(xlim=c(0,77),ylim=c(-80,165),expand=F)

png(filename = "N-limitation Severity by Latitude.png", width=8, height=6, units="in", res=300)
Nlimlat
dev.off()