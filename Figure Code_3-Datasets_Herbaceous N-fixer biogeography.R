### FIGURE CODE FOR GRASSLAND N-FIXER BIOGEOGRAPHY #####
### All 3 Datasets #####################################

#### Load packages and data files ####
library(ggplot2)
library(grid)
library(ggpubr)
library(maps)
library(mapdata)
library(RColorBrewer)
library(wesanderson)

sc<-read.csv("Gridcell-Level Control and Pretreatment Data_3 Datasets.csv")
dun.trees<-read.csv("Menge 2017_tree data_BNT.csv")
crre<-read.csv("Processed CoRRE Data_Pretreatment Year.csv")
gx<-read.csv("Processed GEx Data_First Year Only.csv")
nut<-read.csv("Processed NutNet Data_Pretreatment Year.csv")
sdat3<-read.csv("Site-Level Control and Pretreatment Data_3 Datasets.csv")

sc.nam<-read.csv("Gridcell-Level Control and Pretreatment Data_3 Datasets NORTH AMERICA.csv")

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

##Now with GEx and CoRRE split out
#png(filename = "Map of Control Data Distribution_3 datasets.png", width=14, height=7, units="in", res=300)
#map('world')
#with(sc[sc$dataset=="NutNet",], points(x=LON, y=LAT, bg='forestgreen', pch=21, cex=1.2))
#with(sc[sc$dataset=="GEx",], points(x=LON, y=LAT, bg='blue', pch=21, cex=1.2))
#with(sc[sc$dataset=="CoRRE",], points(x=LON, y=LAT, bg='red', pch=21, cex=1.2))
#legend(x=43, y=-28,legend=c("CoRRE", "GEx"), pt.bg=c("Red","Blue"), pch=21, bty="n", cex=1.2)
#legend(x=43, y=-28,legend=c("NutNet","CoRRE", "GEx"), pt.bg=c("forestgreen","Red","Blue"), pch=21, bty="n", cex=1.2)
#dev.off()

z.pal<-wes_palette("Zissou1", 100, type = "continuous")
png(filename = "Map of Control Data Distribution_Gridcell level 3 Datasets.png", width=14, height=7, units="in", res=300)
ggplot(data=world)+
  theme(panel.background=element_rect(fill="aliceblue", color="aliceblue"))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  geom_sf(color="black", fill="antiquewhite")+
  #geom_point(data=sc[sc$dataset=="CoRRE",], mapping=aes(x=LON,y=LAT),size=3,shape=21,fill='forestgreen')+
  #geom_point(data=sc[sc$dataset=="GEx",], mapping=aes(x=LON,y=LAT),size=3,shape=21,fill='blue')+
  geom_point(data=sc, mapping=aes(x=LON,y=LAT,fill=fixra),size=4.5,shape=21)+
  scale_fill_gradientn(colours = z.pal) +
  labs(fill="N-fixer Relative Abundance (%)")+
  theme(legend.position = "top")+
  ylab(expression("Latitude "*degree*""))+
  xlab(expression("Longitude "*degree*""))
dev.off()

########################################################
#### LATITUDINAL PATTERN FOR INDIVIDUAL  DATASETS ######
########################################################
clrs<-brewer.pal(n = 5,name = "Set1")
ggplot(sdat3[sdat3$dataset=="GEx",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
                axis.text.y=element_text(size=20, colour="black"),
                legend.position="none")+
  geom_point(size=2,colour=clrs[4])+
  ggtitle(label="GEx")+
  geom_text(aes(label=site_code),hjust=0, vjust=0)+
  coord_cartesian(xlim=c(0,80))

ggplot(sdat3[sdat3$dataset=="CoRRE",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        legend.position="none")+
  geom_point(size=2,colour=clrs[3])+
  ggtitle(label="CoRRE")+
  geom_text(aes(label=site_code),hjust=0, vjust=0)+
  coord_cartesian(xlim=c(0,80))

ggplot(sdat3[sdat3$dataset=="NutNet",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        legend.position="none")+
  geom_point(size=2,colour=clrs[5])+
  ggtitle(label="NutNet")+
  geom_text(aes(label=site_code),hjust=0, vjust=0)+
  coord_cartesian(xlim=c(0,80))

g.lat.lbl<-ggplot(sdat3[sdat3$dataset=="GEx",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        legend.position="none")+
  geom_point(size=4,colour=clrs[4])+
  ggtitle(label="GEx")+
  coord_cartesian(xlim=c(0,80))
c.lat.lbl<-ggplot(sdat3[sdat3$dataset=="CoRRE",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        legend.position="none")+
  geom_point(size=4,colour=clrs[3])+
  ggtitle(label="CoRRE")+
  coord_cartesian(xlim=c(0,80))
n.lat.lbl<-ggplot(sdat3[sdat3$dataset=="NutNet",], aes(x=abs.LAT,y=fixra))+
  theme_classic()+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        legend.position="none")+
  geom_point(size=4,colour=clrs[5])+
  ggtitle(label="NutNet")+
  coord_cartesian(xlim=c(0,80))


png(filename = "Latitudinal Pattern for Each Dataset with Labelled Sites.png", width=8, height=12, units="in", res=300)
ggarrange(g.lat.lbl,c.lat.lbl,n.lat.lbl,ncol=1,nrow=3)
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
l.xseq<-seq(0,80,.1)
l.yseq<-(1-alogitfn(coef(abs.LAT_bestmod2)[3]))*exp(coef(abs.LAT_bestmod2)[2])
l.predcrv<-data.frame("x"=l.xseq,"y"=l.yseq)

latfull<-ggplot(sc, aes(x=abs.LAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="black")+
  geom_point(data=latbndat2, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="black")+
  #geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=35),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("Global - 3 Datasets")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=l.predcrv, colour="black", size=1.6)
png(filename = "3 Dataset Global N-Fixer Relative Abundance by Latitude.png", width=8, height=6, units="in", res=300)
latfull
dev.off()
#Mean N-fixer relative % cover ranges from ~1% in the tropics up to ~13% at 49 degrees lat.

###############################################################
#### NORTH AMERICAN LATITUDINAL PATTERN FOR ALL  DATASETS #####
###############################################################
sc.nam$latbins<-cut(sc.nam$abs.LAT,c(seq(20,70,1)),labels=c(seq(20.5,70,1)))
ziln.median.fun<-function(x){(1-(length(x[x==0])/length(x)))*exp(mean(log(x[x>0]),na.rm=T))}
latbndat<-data.frame("LAT"=seq(20.5,70,1),
                     "fixra"=with(sc.nam, tapply(fixra,latbins,mean, na.rm=T)))
latbndat2<-data.frame("LAT"=seq(20.5,70,1),
                      "fixra"=with(sc.nam, tapply(fixra,latbins,ziln.median.fun)))
l.xseq<-seq(20,70,.1)
l.yseq<-(1-alogitfn(coef(abs.LAT_bestmod2.NAM)[3]+coef(abs.LAT_bestmod2.NAM)[4]*l.xseq))*exp(coef(abs.LAT_bestmod2.NAM)[1])
l.predcrv<-data.frame("x"=l.xseq,"y"=l.yseq)

latfull.NAM<-ggplot(sc.nam, aes(x=abs.LAT,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="black")+
  geom_point(data=latbndat2, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="black")+
  #geom_smooth(method="lm",se=T, linetype="dashed")+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=20,xend=20,y=0,yend=35),colour="black")+
  geom_segment(aes(x=20,xend=70,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  ggtitle("North America - 3 Datasets")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=l.predcrv, colour="black", size=1.6)
png(filename = "3 Dataset North American N-Fixer Relative Abundance by Latitude.png", width=8, height=6, units="in", res=300)
latfull.NAM
dev.off()



########################################################
#### COMPARING TREE AND GRASSLAND PATTERNS #############
########################################################
clrs<-brewer.pal(n = 5,name = "Set1")
crre$latbins<-cut(crre$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
c.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(crre, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="CoRRE")
gx$latbins<-cut(gx$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
g.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(gx, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="GEx")
nut$latbins<-cut(nut$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
n.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(nut, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="NutNet")
dun.trees$latbins<-cut(dun.trees$LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
d.tr.dat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(dun.trees, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="Trees")
latbndat$data<-"Grasslands"
comp.dat<-rbind(d.tr.dat,latbndat,c.latbndat,g.latbndat,n.latbndat)

comp.fig<-ggplot(comp.dat[comp.dat$data%in%c("Trees","Grasslands"),], aes(x=LAT,y=fixra))+
  geom_point(aes(shape=data, fill=data, size=5))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.key=element_blank(),
        legend.position=c(.9,1),legend.text=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  geom_segment(aes(x=0,xend=0,y=0,yend=30),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  scale_shape_manual(name="",labels=c("Trees","Grasslands"),
                     values=rep(21,5))+
  scale_fill_manual(name="", labels=c("Trees","Grasslands"),
                    values=c(clrs))+
  geom_line(aes(x=x,y=y), data=l.predcrv, colour=clrs[2],size=1.3,)+
  geom_smooth(data=comp.dat[comp.dat$data=="Trees",],se=F,colour=clrs[1],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  #geom_smooth(data=comp.dat[comp.dat$data=="Grasslands",],se=F,colour=clrs[2],span=1,size=1.3,
  #            method="glm",method.args=list(family=poisson))+
  coord_cartesian(xlim=c(-1,80),ylim=c(-1,30))

png(filename = "Trees vs. Grasslands Comparison Fig.png", width=8, height=6, units="in", res=300)
comp.fig
dev.off()

########################################################
#### COMPARING TREE AND 3 GRASSLAND DATASETS ###########
########################################################
clrs<-brewer.pal(n = 5,name = "Set1")
crre$latbins<-cut(crre$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
c.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                       "fixra"=with(crre, tapply(fixra,latbins,mean, na.rm=T)),
                       "data"="CoRRE")
gx$latbins<-cut(gx$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
g.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                       "fixra"=with(gx, tapply(fixra,latbins,mean, na.rm=T)),
                       "data"="GEx")
nut$latbins<-cut(nut$abs.LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
n.latbndat<-data.frame("LAT"=seq(0.5,80,1),
                       "fixra"=with(nut, tapply(fixra,latbins,mean, na.rm=T)),
                       "data"="NutNet")
dun.trees$latbins<-cut(dun.trees$LAT,c(seq(0,80,1)),labels=c(seq(0.5,80,1)))
d.tr.dat<-data.frame("LAT"=seq(0.5,80,1),
                     "fixra"=with(dun.trees, tapply(fixra,latbins,mean, na.rm=T)),
                     "data"="Trees")
latbndat$data<-"Grasslands"
comp.dat<-rbind(d.tr.dat,latbndat,c.latbndat,g.latbndat,n.latbndat)

l.xseq<-seq(0,80,.1)
l.yseq<-(1-alogitfn(coef(abs.LAT_bestmod2)[3]))*exp(coef(abs.LAT_bestmod2)[2])
l.yseq.nn<-(1-alogitfn(coef(abs.LAT_bestmod.nn)[3]+(coef(abs.LAT_bestmod.nn)[4]*l.xseq)))*exp(coef(abs.LAT_bestmod.nn)[1])
l.yseq.gx<-(1-alogitfn(coef(abs.LAT_bestmod.gex)[3]+(coef(abs.LAT_bestmod.gex)[4]*l.xseq)))*exp(coef(abs.LAT_bestmod.gex)[1])
l.yseq.cr<-(1-alogitfn(coef(abs.LAT_bestmod.cr)[3]))*exp(coef(abs.LAT_bestmod.cr)[2])
l.predcrv<-data.frame("x"=l.xseq,"all"=l.yseq,"nn"=l.yseq.nn,"gx"=l.yseq.gx,"cr"=l.yseq.cr)

comp.3d.fig_smooths<-ggplot(comp.dat, aes(x=LAT,y=fixra))+
  geom_point(aes(shape=data, fill=data, size=5))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.key=element_blank(),
        legend.position=c(.9,1),legend.text=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  geom_segment(aes(x=0,xend=0,y=0,yend=30),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  scale_shape_manual(name="",labels=c("Trees","Grasslands","CoRRE","GEx","NutNet"),
                     values=rep(21,5))+
  scale_fill_manual(name="", labels=c("Trees","Grasslands","CoRRE","GEx","NutNet"),
                    values=c(clrs))+
  #geom_line(aes(x=x,y=all), data=l.predcrv, colour=clrs[2],size=1.3,linetype="dashed")+
  #geom_line(aes(x=x,y=cr), data=l.predcrv, colour=clrs[3],size=1.3,linetype="dashed")+
  #geom_line(aes(x=x,y=gx), data=l.predcrv, colour=clrs[4],size=1.3,linetype="dashed")+
  #geom_line(aes(x=x,y=nn), data=l.predcrv, colour=clrs[5],size=1.3,linetype="dashed")+
  geom_smooth(data=comp.dat[comp.dat$data=="Trees",],se=F,colour=clrs[1],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  geom_smooth(data=comp.dat[comp.dat$data=="Grasslands",],se=F,colour=clrs[2],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  geom_smooth(data=comp.dat[comp.dat$data=="CoRRE",],se=F,colour=clrs[3],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  geom_smooth(data=comp.dat[comp.dat$data=="GEx",],se=F,colour=clrs[4],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  geom_smooth(data=comp.dat[comp.dat$data=="NutNet",],se=F,colour=clrs[5],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  coord_cartesian(xlim=c(-1,80),ylim=c(-1,30))

png(filename = "Trees vs. 3-Dataset Grasslands Comparison Fig_Smooths.png", width=8, height=6, units="in", res=300)
comp.3d.fig_smooths
dev.off()

comp.3d.fig_mle<-ggplot(comp.dat, aes(x=LAT,y=fixra))+
  geom_point(aes(shape=data, fill=data, size=5))+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.key=element_blank(),
        legend.position=c(.9,1),legend.text=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  geom_segment(aes(x=0,xend=0,y=0,yend=30),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  scale_shape_manual(name="",labels=c("Trees","Grasslands","CoRRE","GEx","NutNet"),
                     values=rep(21,5))+
  scale_fill_manual(name="", labels=c("Trees","Grasslands","CoRRE","GEx","NutNet"),
                    values=c(clrs))+
  geom_line(aes(x=x,y=all), data=l.predcrv, colour=clrs[2],size=1.3)+
  geom_line(aes(x=x,y=cr), data=l.predcrv, colour=clrs[3],size=1.3)+
  geom_line(aes(x=x,y=gx), data=l.predcrv, colour=clrs[4],size=1.3)+
  geom_line(aes(x=x,y=nn), data=l.predcrv, colour=clrs[5],size=1.3)+
  geom_smooth(data=comp.dat[comp.dat$data=="Trees",],se=F,colour=clrs[1],span=1,size=1.3,
              method="glm",method.args=list(family=poisson))+
  coord_cartesian(xlim=c(-1,80),ylim=c(-1,30))

png(filename = "Trees vs. 3-Dataset Grasslands Comparison Fig_MLE fits.png", width=8, height=6, units="in", res=300)
comp.3d.fig_mle
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. TEMPERATURE #################
########################################################
sc$tempbins<-cut((sc$MAT),c(seq(-12,29,1)),labels=c(seq(-11.5,29,1)))
tempbn.mn<-data.frame("MAT"=seq(-11.5,29,1),
                      "fixra"=with(sc, tapply(fixra,tempbins,mean, na.rm=T)))
tempbn.mdn<-data.frame("MAT"=seq(-11.5,29,1),
                       "fixra"=with(sc, tapply(fixra,tempbins,ziln.median.fun)))

t.xseq<-seq(-12,29,.1)
t.yseq<-(1-alogitfn(coef(MAT_bestmod2)[3]))*exp(coef(MAT_bestmod2)[2])
t.predcrv<-data.frame("x"=t.xseq,"y"=t.yseq)

tempfull<-ggplot(sc, aes(x=(MAT),y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=tempbn.mn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=tempbn.mdn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill=clrs[2])+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-12,xend=-12,y=0,yend=33),colour="black")+
  geom_segment(aes(x=-12,xend=30,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=t.predcrv, colour="black",size=2)+
  #geom_abline(slope=0,intercept=temp.predmedian,colour="black", size=2)+
  coord_cartesian(xlim=c(-12,30),ylim=c(0,33),expand=F)

png(filename = "N-fixer Abundance vs. Temperature_3-Datasets.png", width=8, height=6, units="in", res=300)
tempfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. PRECIPITATION ###############
########################################################
sc$precbins<-cut(sc$MAP,c(seq(80,2900,10)),labels=c(seq(85,2900,10)))
precbn.mn<-data.frame("MAP"=seq(85,2900,10),
                      "fixra"=with(sc, tapply(fixra,precbins,mean, na.rm=T)))
precbn.mdn<-data.frame("MAP"=seq(85,2900,10),
                       "fixra"=with(sc, tapply(fixra,precbins,ziln.median.fun)))

prec.xseq<-seq(80,2900,.1)
prec.yseq<-(1-alogitfn(coef(MAP_bestmod2)[3]))*exp(coef(MAP_bestmod2)[2])
prec.predcrv<-data.frame("x"=prec.xseq,"y"=prec.yseq)

precfull<-ggplot(sc, aes(x=MAP,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=precbn.mn, mapping=aes(x=MAP,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=precbn.mdn, mapping=aes(x=MAP,y=fixra),size=4,shape=21, fill=clrs[2])+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=80,xend=80,y=0,yend=35),colour="black")+
  geom_segment(aes(x=80,xend=1900,y=0,yend=0),colour="black")+
  xlab(expression("MAP (mm)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=(y)), data=prec.predcrv, colour="black",size=2)+
  #geom_abline(slope=0, intercept=predmedian, colour="black")+
  coord_cartesian(xlim=c(80,1900),ylim=c(0,35),expand=F)

png(filename = "N-fixer Abundance vs. Precipitation_3 Datasets.png", width=8, height=6, units="in", res=300)
precfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. N LIMITATION ################
########################################################
nlimdat<-sc[!is.na(sc$Nlim),]
nlimdat$bins<-cut(nlimdat$Nlim,c(seq(-26,85,1)),labels=c(seq(-25.5,85,1)))
limbn.mn<-data.frame("Nlim"=seq(-25.5,85,1),
                     "fixra"=with(nlimdat, tapply(fixra,bins,mean, na.rm=T)))
limbn.mdn<-data.frame("Nlim"=seq(-25.5,85,1),
                      "fixra"=with(nlimdat, tapply(fixra,bins,ziln.median.fun)))

nlim.xseq<-seq(-26,85,.1)
nlim.yseq<-(1-alogitfn(coef(N.lim_bestmod2)[4]))*exp(coef(N.lim_bestmod2)[2]+(coef(N.lim_bestmod2)[3]*nlim.xseq))
nlim.predcrv<-data.frame("x"=nlim.xseq,"y"=nlim.yseq)

Nlimfull<-ggplot(nlimdat, aes(x=Nlim,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=limbn.mdn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill=clrs[2])+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-26,xend=-26,y=0,yend=20),colour="black")+
  geom_segment(aes(x=-26,xend=90,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=20),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=nlim.predcrv, colour="black", size=2)+
  xlim(-26,85)+
  coord_cartesian(xlim=c(-26,85),ylim=c(0,20))

png(filename = "N-fixer Abundance vs. N Limitation_3-Datasets.png", width=8, height=6, units="in", res=300)
Nlimfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. P LIMITATION ################
########################################################
plimdat<-sc[!is.na(sc$Plim),]
plimdat$bins<-cut(plimdat$Plim,c(seq(-97,152,1)),labels=c(seq(-96.5,152,1)))
plimbn.mn<-data.frame("Plim"=seq(-96.5,152,1),
                      "fixra"=with(plimdat, tapply(fixra,bins,mean, na.rm=T)))
plimbn.mdn<-data.frame("Plim"=seq(-96.5,152,1),
                       "fixra"=with(plimdat, tapply(fixra,bins,ziln.median.fun)))

plim.xseq<-seq(-97,152,.1)
plim.yseq<-(1-alogitfn(coef(Plim_bestmod2)[3]))*exp(coef(Plim_bestmod2)[2])
plim.predcrv<-data.frame("x"=plim.xseq,"y"=plim.yseq)

Plimfull<-ggplot(plimdat, aes(x=Plim,y=fixra))+
  geom_point(size=2,colour="grey")+
  #geom_point(data=plimbn.mn, mapping=aes(x=Plim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=plimbn.mdn, mapping=aes(x=Plim,y=fixra),size=4,shape=21,fill=clrs[2])+
  theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-100,xend=-100,y=0,yend=20),colour="black")+
  geom_segment(aes(x=-100,xend=155,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=20),colour="black",linetype="dashed")+
  xlab(expression("P Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("Global")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=y), data=plim.predcrv, colour="black", size=2)+
  coord_cartesian(xlim=c(-100,155))

png(filename = "N-fixer Abundance vs. P Limitation_3-Datasets.png", width=8, height=6, units="in", res=300)
Plimfull
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. N-FIXER RICHNESS ############
########################################################
sc$rrbins<-cut(sc$fixrr,c(seq(0,60,1)),labels=c(seq(.5,60,1)))
rrbn.mn<-data.frame("fixrr"=seq(.5,60,1),
                    "fixra"=with(sc, tapply(fixra,rrbins,mean, na.rm=T)))
rrbn.mdn<-data.frame("fixrr"=seq(.5,60,1),
                     "fixra"=with(sc, tapply(fixra,rrbins,ziln.median.fun)))

rr.xseq<-seq(0,60,.1)
rr.yseq<-(1-alogitfn(coef(fixrr_bestmod)[5]+(coef(fixrr_bestmod)[6]*rr.xseq)))*exp(coef(fixrr_bestmod)[2]+(coef(fixrr_bestmod)[3]*rr.xseq)+(coef(fixrr_bestmod)[4]*rr.xseq^2))
rr.predcrv<-data.frame("x"=rr.xseq,"y"=rr.yseq)

rva<-ggplot(sc, aes(x=fixrr, y=fixra))+
  geom_point(size=2, colour="grey")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  #geom_point(data=rrbn.mn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=rrbn.mdn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill=clrs[2])+
  geom_segment(aes(x=0,xend=0,y=0,yend=33),colour="black")+
  geom_segment(aes(x=0,xend=60,y=0,yend=0),colour="black")+
  xlab(expression("N-Fixer Relative Richness (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("N-fixer Richness vs. Abundance")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=(y)), data=rr.predcrv, colour="black",size=2)+
  coord_cartesian(xlim=c(0,60),ylim=c(0,33))

png(filename = "N-fixer Abundance vs. N-Fixer Relative Richness_3-Datasets.png", width=8, height=6, units="in", res=300)
rva
dev.off()

########################################################
#### N-FIXER ABUNDANCE VS. HERBIVORY ###################
########################################################
glimdat<-sc[!is.na(sc$Glim),]
glimdat$glimbins<-cut(glimdat$Glim,c(seq(-3.7,3.1,.1)),labels=c(seq(-3.65,3.1,.1)))
glimbn.mn<-data.frame("Glim"=seq(-3.65,3.1,.1),
                      "fixra"=with(glimdat, tapply(fixra,glimbins,mean, na.rm=T)))
glimbn.mdn<-data.frame("Glim"=seq(-3.65,3.1,.1),
                       "fixra"=with(glimdat, tapply(fixra,glimbins,ziln.median.fun)))

glim.xseq<-seq(-3.7,3.1,.001)
glim.yseq<-(1-alogitfn(coef(Glim_bestmod2)[5]))*exp(coef(Glim_bestmod2)[2]+(coef(Glim_bestmod2)[3]*glim.xseq)+(coef(Glim_bestmod2)[4]*glim.xseq^2))
glim.predcrv<-data.frame("x"=glim.xseq,"y"=glim.yseq)

Glimfull<-ggplot(glimdat, aes(x=Glim, y=fixra))+
  geom_point(size=2, colour="grey")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  #geom_point(data=glimbn.mn, mapping=aes(x=Glim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=glimbn.mdn, mapping=aes(x=Glim,y=fixra),size=4,shape=21,fill=clrs[2])+
  geom_segment(aes(x=0,xend=0,y=0,yend=33),colour="black",linetype="dashed")+
  geom_segment(aes(x=-3.7,xend=-3.7,y=0,yend=33),colour="black")+
  geom_segment(aes(x=-3.7,xend=3.1,y=0,yend=0),colour="black")+
  xlab(expression("Grazer Preference for N fixers"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("N-fixer Richness vs. Abundance")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=(y)), data=glim.predcrv, colour="black",size=2)+
  coord_cartesian(xlim=c(-3.7,3.1),ylim=c(0,33))

png(filename = "Glim Effects on N-fixer Abundance_3-Datasets.png", width=6, height=6, units="in", res=300)
Glimfull
dev.off()

#### Plotting GEx data only with designation for Herbivore Type
g.anml$herb.type.num<-ifelse(g.anml$herbivore.type=="grazer", yes=1, no=3)
g.anml$herb.type.num<-ifelse(g.anml$herbivore.type=="grazer_browser", yes=2, no=g.anml$herb.type.num)
glimdat2<-g.anml[!is.na(g.anml$Glim),]
glimdat2$glimbins<-cut(glimdat2$Glim,c(seq(-3.7,3.1,.1)),labels=c(seq(-3.65,3.1,.1)))
glimbn.mn2<-data.frame("Glim"=seq(-3.65,3.1,.1),
                      "fixra"=with(glimdat2, tapply(fixra,glimbins,mean, na.rm=T)),
                      "herb.type"=with(glimdat2, tapply(herb.type.num,glimbins,mean, na.rm=T)))
glimbn.mdn2<-data.frame("Glim"=seq(-3.65,3.1,.1),
                       "fixra"=with(glimdat2, tapply(fixra,glimbins,ziln.median.fun)),
                       "herb.type"=with(glimdat2, tapply(herb.type.num,glimbins,mean, na.rm=T)))

#glim.xseq<-seq(-3.7,3.1,.001)
#glim.yseq<-(1-alogitfn(coef(Glim_bestmod2)[5]))*exp(coef(Glim_bestmod2)[2]+(coef(Glim_bestmod2)[3]*glim.xseq)+(coef(Glim_bestmod2)[4]*glim.xseq^2))
#glim.predcrv<-data.frame("x"=glim.xseq,"y"=glim.yseq)

Glim_herb.type<-ggplot(glimdat2, aes(x=Glim, y=fixra))+
  #geom_point(size=2, colour="grey")+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="top")+
  #geom_point(data=glimbn.mn, mapping=aes(x=Glim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=glimbn.mdn2, mapping=aes(x=Glim,y=fixra, fill=herb.type),size=4,shape=21)+
  scale_fill_gradientn(colours=z.pal,breaks=c(1,3),labels=c("Grazer","Browser"))+
  labs(fill="Herbivore Type")+
  geom_segment(aes(x=0,xend=0,y=0,yend=37),colour="black",linetype="dashed")+
  geom_segment(aes(x=-3.7,xend=-3.7,y=0,yend=37),colour="black")+
  geom_segment(aes(x=-3.7,xend=3.1,y=0,yend=0),colour="black")+
  xlab(expression("Herbivore Preference for N fixers"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #ggtitle("N-fixer Richness vs. Abundance")+
  theme(plot.title=element_text(hjust=.5))+
  geom_line(aes(x=x,y=(y)), data=glim.predcrv, colour="black",size=2)+
  coord_cartesian(xlim=c(-3.7,3.1),ylim=c(0,37))
png(filename = "Glim plot with herbivore type colored.png", width=8, height=6, units="in", res=300)
Glim_herb.type
dev.off()

#### Interaction plot of herbivore type and herbivory pressure on N-fixer Abundance #####################
#glimdat3<-glimdat2
#glimdat3$
v1<-ggplot(glimdat2, aes(x=as.factor(herb.type.num), y=fixra,fill=as.factor(grazing.pressure)))+
  geom_boxplot()+
  labs(fill="Grazing Pressure", x="Herbivore Type", y="N-fixer Relative Abundance (%)")+
  ylim(c(0,25))+
  theme_classic()+
  scale_x_discrete(breaks=c("1","2","3"),labels=c("Grazer", "Mixed", "Browser"))+
  theme(text=element_text(size=18, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), legend.position="top")
amend<-ggplot_build(v1)
amend$data[[1]]$x[7:8]<-c(2.75,3.25)
amend$data[[1]]$xmax[7:8]<-c(2.6375,3.1375)
amend$data[[1]]$xmin[7:8]<-c(2.8625,3.3625)
library(grid)
png(filename = "Interaction of Herbivore Type and Grazing Pressure.png", width=8, height=6, units="in", res=300)
grid.draw(ggplot_gtable(amend))
dev.off()

#herb.mdn.plt<-ggplot(data=g.herbdf[g.herbdf$fix.ra>0,],aes(x=graz.trt,y=fix.ra,fill=graz.trt))+
#  geom_boxplot(notch=F,fill=c("grey40","grey70"))+
#  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
#  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=15, colour="black"),
#        axis.text.y=element_text(size=15, colour="black"))+
#  scale_x_discrete(label=c("Grazed","Ungrazed"))+
#  annotate("text",x=1.5,y=25,label=c("n.s."),size=5)+
#  xlab(expression("Grazer Presence"))+
#  ylab("N-Fixer Relative Abundance (%)")+
#  geom_segment(aes(x=0.5,xend=0.5,y=0,yend=25),colour="black")+
#  geom_segment(aes(x=0.5,xend=2.5,y=0,yend=0),colour="black")

########################################################
#### MULTI-PANEL FIG OF ECOLOGICAL DRIVERS #############
########################################################

lat.mdn.plt<-ggplot()+
  #geom_point(data=latbndat, mapping=aes(x=LAT,y=fixra),size=4, shape=21, fill="red")+
  geom_point(data=latbndat2, mapping=aes(x=LAT,y=fixra),size=5, shape=21, fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=0,xend=0,y=0,yend=30),colour="black")+
  geom_segment(aes(x=0,xend=80,y=0,yend=0),colour="black")+
  xlab(expression("abs(latitude "*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=all), data=l.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(0,80),ylim=c(-1,30))

rich.mdn.plt<-ggplot()+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  #geom_point(data=rrbn.mn, mapping=aes(x=fixrr,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=rrbn.mdn, mapping=aes(x=fixrr,y=fixra),size=5,shape=21,fill=clrs[2])+
  geom_segment(aes(x=0,xend=0,y=0,yend=22),colour="black")+
  geom_segment(aes(x=0,xend=60,y=0,yend=0),colour="black")+
  xlab(expression("N-Fixer Relative Richness (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=(y)), data=rr.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-1,60),ylim=c(-1,22))

Nlim.mdn.plt<-ggplot()+
  #geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=limbn.mdn, mapping=aes(x=Nlim,y=fixra),size=5,shape=21,fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-33,xend=-33,y=0,yend=15),colour="black")+
  geom_segment(aes(x=-33,xend=85,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=15),colour="black",linetype="dashed")+
  xlab(expression("N Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=nlim.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-33,85),ylim=c(-1,15))

Plim.mdn.plt<-ggplot()+
  #geom_point(data=limbn.mn, mapping=aes(x=Nlim,y=fixra),size=4,shape=21,fill="red")+
  geom_point(data=plimbn.mdn, mapping=aes(x=Plim,y=fixra),size=5,shape=21,fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), legend.position="none")+
  geom_segment(aes(x=-100,xend=-100,y=0,yend=15),colour="black")+
  geom_segment(aes(x=-100,xend=151,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=15),colour="black",linetype="dashed")+
  xlab(expression("P Limitation RR (%)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_segment(aes(x=-100,xend=151,y=plim.predmedian,yend=plim.predmedian),colour="black",size=1.2)+
  #geom_line(aes(x=x,y=y), data=plim.predcrv, colour="black",size=1.2)+
  coord_cartesian(xlim=c(-100,151),ylim=c(0,15))

t.mdn.plt<-ggplot()+
  #geom_point(data=tempbn.mn, mapping=aes(x=MAT,y=fixra),size=4,shape=21, fill="red")+
  geom_point(data=tempbn.mdn, mapping=aes(x=MAT,y=fixra),size=5,shape=21, fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="none",
        plot.margin = margin(1.5,.8,.8,.8,"cm"))+
  geom_segment(aes(x=-12,xend=-12,y=0,yend=33),colour="black")+
  geom_segment(aes(x=-12,xend=30,y=0,yend=0),colour="black")+
  xlab(expression("MAT ("*degree*")"))+
  ylab("N-Fixer Relative Abundance (%)")+
  #geom_segment(aes(x=-6,xend=28,y=temp.predmedian,yend=temp.predmedian),colour="black",size=1.2)+
  geom_line(aes(x=x,y=y), data=t.predcrv, colour="black", size=1.2)+
  coord_cartesian(xlim=c(-12,30),ylim=c(-1,34),expand=F)

prec.mdn.plt<-ggplot()+
  geom_point(data=precbn.mdn, mapping=aes(x=MAP,y=fixra),size=5,shape=21, fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="none",
        plot.margin = margin(1.5,.8,.8,.8,"cm"))+
  geom_segment(aes(x=85,xend=85,y=0,yend=36),colour="black")+
  geom_segment(aes(x=85,xend=2900,y=0,yend=0),colour="black")+
  xlab(expression("MAP (mm)"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=prec.predcrv, colour="black", size=1.2)+
  coord_cartesian(xlim=c(85,2900),ylim=c(0,36.5),expand=F)

glim.mdn.plt<-ggplot()+
  geom_point(data=glimbn.mdn, mapping=aes(x=Glim,y=fixra),size=5,shape=21, fill=clrs[2])+
  theme(text=element_text(size=15, colour="black"),axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"), 
        legend.position="none",
        plot.margin = margin(1.5,.8,.8,.8,"cm"))+
  geom_segment(aes(x=-3.7,xend=-3.7,y=0,yend=33),colour="black")+
  geom_segment(aes(x=-3.7,xend=3.1,y=0,yend=0),colour="black")+
  geom_segment(aes(x=0,xend=0,y=0,yend=33),colour="black",linetype="dashed")+
  xlab(expression("Grazer Preference for N fixers"))+
  ylab("N-Fixer Relative Abundance (%)")+
  geom_line(aes(x=x,y=y), data=glim.predcrv, colour="black", size=1.2)+
  coord_cartesian(xlim=c(-3.7,3.1),ylim=c(0,33.5),expand=F)

png(filename = "Ecological Drivers of Grassland N fixers_3-Datasets.png", width=10, height=12, units="in", res=300)
ggarrange(t.mdn.plt,prec.mdn.plt,Nlim.mdn.plt,Plim.mdn.plt,rich.mdn.plt,glim.mdn.plt,ncol=2,nrow=3)
dev.off()
