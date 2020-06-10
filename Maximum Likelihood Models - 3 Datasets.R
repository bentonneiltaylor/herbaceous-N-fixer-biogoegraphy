##############################################################################################
#### MAXIMUM LIKELIHOOD MODELS FOR N-FIXER RELATIVE ABUNDANCE ACROSS ECOLOGICAL GRADIENTS ####
#### Uses "pre-treatment" data for all 3 datasets (NutNet, CoRRE, and GEx) ####

library(bbmle)

dat.mle<-read.csv("Gridcell-Level Control and Pretreatment Data_3 Datasets.csv")
dat.mle$logfixra<-log(dat.mle$fixra)
dat.mle$MAP<-as.numeric(as.character(dat.mle$MAP))

sdat3<-read.csv("Site-Level Control and Pretreatment Data_3 Datasets.csv")
nn.mle<-sdat3[sdat3$dataset=="NutNet",]
nn.mle$logfixra<-log(nn.mle$fixra)
cr.mle<-sdat3[sdat3$dataset=="CoRRE",]
cr.mle$logfixra<-log(cr.mle$fixra)
g.mle<-sdat3[sdat3$dataset=="GEx",]
g.mle$logfixra<-log(g.mle$fixra)

###########################################################################
###########################################################################
### What is the latitudinal pattern of N-fixer Relative Abundance? ########
###########################################################################
###########################################################################

###########################################################################
### First, All Grassland Data Put Together ################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$abs.LAT),]
rv<-sub$fixra
ind<-sub$abs.LAT  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_rv_null)
abs.LAT_bestmod2<-fit_rv_null

###########################################################################
### Now For Each Dataset Individually - NutNet, GEx, and CoRRE ############
###########################################################################

### NUTNET FIRST ################################################
### Subsetting the data and assigning variables for this question  
sub<-nn.mle[!is.na(nn.mle$abs.LAT),]
rv<-sub$fixra
ind<-sub$abs.LAT  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_p_lin)
abs.LAT_bestmod.nn<-fit_p_lin

#################################################################
### NOW GEx #####################################################
### Subsetting the data and assigning variables for this question  
sub<-g.mle[!is.na(g.mle$abs.LAT),]
rv<-sub$fixra
ind<-sub$abs.LAT  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_p_lin)
abs.LAT_bestmod.gex<-fit_p_lin

### NOW CoRRE ################################################
### Subsetting the data and assigning variables for this question  
sub<-cr.mle[!is.na(cr.mle$abs.LAT),]
rv<-sub$fixra
ind<-sub$abs.LAT  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_rv_null)
abs.LAT_bestmod.cr<-fit_rv_null





###########################################################################
###########################################################################
### Now Assessing Ecological Drivers of N-fixer Abundance #################
###########################################################################
###########################################################################

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with Temperature? ################
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$MAT),]
rv<-sub$fixra
ind<-sub$MAT  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<-.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=.00005,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_p_quad,sub),AICc(fit_u_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_rv_null)
MAT_bestmod2<-fit_rv_null

with(sub, summary(lm(fixra~MAT)))#Basic regression agrees with no relationship

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with Precipitation? ##############
###########################################################################
###########################################################################

# !!!!! NOT ALL MODELS ARE WORKING CURRENTLY !!!! 

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$MAP),]
rv<-sub$fixra
ind<-sub$MAP  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<--.003
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=.003),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=.5,d2=.5),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub)))
### Summary for the best-fit model ###
summary(fit_rv_null)
MAP_bestmod2<-fit_rv_null

with(sub, summary(lm(fixra~MAP)))#Basic regression agrees - no significant relationship

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with N Limitation? ###############
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$Nlim),]
rv<-sub$fixra
ind<-sub$Nlim  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_u_lin)
N.lim_bestmod2<-fit_u_lin

with(sub, summary(lm(fixra~Nlim)))#Basic regression shows non-significant decline with N lim

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with P Limitation? ###############
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$Plim),]
rv<-sub$fixra
ind<-sub$Plim  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_rv_null)
Plim_bestmod2<-fit_rv_null

with(sub, summary(lm(fixra~Plim)))#Basic regression agrees with no significant relationshp

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with N-fixer Relative Richness? ##
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$fixrr),]
rv<-sub$fixra
ind<-sub$fixrr  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_p.u_quad)
fixrr_bestmod2<-fit_p.u_quad

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with Grazing Limitation? #########
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$Glim),]
rv<-sub$fixra
ind<-sub$Glim  

#Variables and starting values
pstart<-.5
ustart<-log(3)
sstart<-log(1.7)

cstart<-log(2)
bstart<-log(1.1)
b2start<-.01
estart<-.5
dstart<-.1
d2start<--.005

#Flat mean ################################################################
ziln_null <- function(sdlogrv,ulogrv,p_no_rv,rv_dat){
  J <- sum(rv==0)
  K <- sum(rv>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}
fit_rv_null <- mle2(ziln_null,
                    start=list(sdlogrv=sstart,ulogrv=ustart,p_no_rv=pstart),
                    data=list(rv_dat=rv))
print(summary(fit_rv_null))
predmedian<-(1-coef(fit_rv_null)[3])*exp(coef(fit_rv_null)[2])#for plotting

#Linear variation in just u ################################################
ziln_u_lin <- function(sdlogrv,c,b,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_lin <- mle2(ziln_u_lin,
                  start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_lin))

#Linear variation in just p ################################################
ziln_p_lin <- function(ulogrv,sdlogrv,e,d,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_lin <- mle2(ziln_p_lin,
                  start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart),
                  data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_lin))

#Linear variation in both p and u ##########################################
ziln_p.u_lin <- function(sdlogrv,c,b,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_lin <- mle2(ziln_p.u_lin,
                    start=list(sdlogrv=sstart,c=cstart,b=bstart,e=estart,d=dstart),
                    data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_lin))

#Quadratic variation in just u ################################################
ziln_u_quad<- function(sdlogrv,c,b,b2,p_no_rv,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  J <- sum(rv_dat==0)
  K <- sum(rv_dat>0)
  -J*log(p_no_rv) - K*log(1-p_no_rv) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u_quad <- mle2(ziln_u_quad,
                   start=list(sdlogrv=sstart,p_no_rv=pstart,c=cstart,b=bstart,b2=b2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u_quad))

#Quadratic variation in just p ################################################
ziln_p_quad<- function(ulogrv,sdlogrv,e,d,d2,rv_dat,ind_dat){
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv,sd=exp(sdlogrv),log=TRUE))
}

fit_p_quad <- mle2(ziln_p_quad,
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=dstart,d2=d2start),
                   data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p_quad))

#Quadratic variation in u,linear variation in p ################################
ziln_p.u_quad <- function(sdlogrv,c,b,b2,e,d,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_quad <- mle2(ziln_p.u_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_quad))

#Linear variation in u,quadratic variation in p ################################
ziln_u.p_quad <- function(sdlogrv,c,b,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_u.p_quad <- mle2(ziln_u.p_quad,
                     start=list(sdlogrv=sstart,c=cstart,b=bstart,d=dstart,d2=d2start,e=estart),
                     data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_u.p_quad))

#Quadratic variation in p and u ################################
ziln_p.u_dualquad <- function(sdlogrv,c,b,b2,e,d,d2,rv_dat,ind_dat){
  ulogrv<-c+b*ind+b2*ind^2
  p_no_rv<-alogitfn(e+d*ind+d2*ind^2)
  -sum(log(p_no_rv[rv_dat==0])) - sum(log(1 - p_no_rv[rv_dat>0])) - 
    sum(dlnorm(rv_dat[rv_dat>0],mean=ulogrv[rv_dat>0],sd=exp(sdlogrv),log=TRUE))
}

fit_p.u_dualquad <- mle2(ziln_p.u_dualquad,
                         start=list(sdlogrv=sstart,c=cstart,b=bstart,b2=b2start,d=dstart,d2=d2start,e=estart),
                         data=list(rv_dat=rv,ind_dat=ind))
print(summary(fit_p.u_dualquad))

# Comparing the models #########################################################
deltaAICs(c(AICc(fit_rv_null,sub),AICc(fit_u_lin,sub),AICc(fit_p_lin,sub),
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p_quad,sub),
            AICc(fit_p.u_quad,sub),AICc(fit_u.p_quad,sub),AICc(fit_p.u_dualquad,sub)))
### Summary for the best-fit model ###
summary(fit_u_quad)
Glim_bestmod2<-fit_u_quad

###########################################################################
###########################################################################
### END ###################################################################
###########################################################################
###########################################################################