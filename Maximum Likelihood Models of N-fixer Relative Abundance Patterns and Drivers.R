#### MAXIMUM LIKELIHOOD MODELS FOR N-FIXER RELATIVE ABUNDANCE ACROSS ECOLOGICAL GRADIENTS ####
#### Uses "pre-treatment" data for all 3 datasets (NutNet, CoRRE, and GEx) ####

library(bbmle)

dat.mle<-read.csv("Gridcell-Level Control and Pretreatment Data_CoRRE and GEx.csv")
dat.mle$logfixra<-log(dat.mle$fixra)
dat.mle$MAP<-as.numeric(as.character(dat.mle$MAP))
###########################################################################
###########################################################################
### What is the latitudinal pattern of N-fixer Relative Abundance? ########
###########################################################################
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
summary(fit_u_lin)
abs.LAT_bestmod<-fit_u_lin

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance Vary with Temperature? ################
###########################################################################
###########################################################################

### Subsetting the data and assigning variables for this question  
sub<-dat.mle[!is.na(dat.mle$MAT),]
rv<-sub$fixra
ind<-sub$MAT/10  

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
MAT_bestmod<-fit_rv_null

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
                   start=list(sdlogrv=sstart,ulogrv=ustart,e=estart,d=.003,d2=.000005),
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
            AICc(fit_p.u_lin,sub),AICc(fit_u_quad,sub),AICc(fit_p.u_quad,sub)))
### Summary for the best-fit model ###
summary(fit_u_lin)
MAP_bestmod<-fit_u_lin

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
summary(fit_rv_null)
N.lim_bestmod<-fit_rv_null

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
Plim_bestmod<-fit_rv_null

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
                     start=list(sdlogrv=sstart,c=cstart,b=.001,b2=.001,d=dstart,e=.1),
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
fixrr_bestmod<-fit_p.u_quad

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
Glim_bestmod<-fit_u_quad

###########################################################################
###########################################################################
### Does N-fixer Relative Abundance change with herbivory in GEx? #########
###########################################################################
###########################################################################

#calculating N-fixer relative abundance and relative richness in grazed and ungrazed plots in GEx
gplts<-unique(gex$site.plt.yr)
gherb<-NULL
for(p in 1:length(gplts)){
  print(p)
  temp<-gex[gex$site.plt.yr==gplts[p],]
  u<-temp[temp$trt=="U",]
  gr<-temp[temp$trt=="G",]
  tempdf<-data.frame("site.plt.yr"=unique(temp$site.plt.yr),
                     "grazed.tot.relcov"=sum(gr$relcov, na.rm=T),
                     "ungrazed.tot.relcov"=sum(u$relcov, na.rm=T),
                     "grazed.tot.rich"=length(unique(gr$genus_species, na.rm=T)),
                     "ungrazed.tot.rich"=length(unique(u$genus_species, na.rm=T)),
                     "grazed.fix.ra"=(sum(gr[gr$N_fixer==1,]$relcov, na.rm=T)/sum(gr$relcov, na.rm=T))*100,
                     "ungrazed.fix.ra"=(sum(u[u$N_fixer==1,]$relcov, na.rm=T)/sum(u$relcov, na.rm=T))*100,
                     "grazed.fix.rr"=(length(unique(gr[gr$N_fixer==1,]$genus_species, na.rm=T))/length(unique(gr$genus_species, na.rm=T)))*100,
                     "ungrazed.fix.rr"=(length(unique(u[u$N_fixer==1,]$genus_species, na.rm=T))/length(unique(u$genus_species, na.rm=T)))*100)
  gherb<-rbind(gherb,tempdf)
}

g2<-merge(g,gherb,by="site.plt.yr",all.x=T,all.y=T)


#This creates a dataframe where grazed and ungrazed data are split out by rows to run anovas on grazing treatment
g.herbdf<-data.frame("site.plt.yr"=c(g2$site.plt.yr,g2$site.plt.yr),
                     "graz.trt"=c(rep("U",776),rep("G",776)),
                     "tot.relcov"=c(g2$ungrazed.tot.relcov,g2$grazed.tot.relcov),
                     "tot.rich"=c(g2$ungrazed.tot.rich,g2$grazed.tot.rich),
                     "fix.ra"=c(g2$ungrazed.fix.ra,g2$grazed.fix.ra),
                     "fix.rr"=c(g2$ungrazed.fix.rr,g2$grazed.fix.rr))

##########################################################################
# Running MLE models to test for effect of Herbivory on fixer relative abundance 
# assuming 0-inflated lognormal data

# Adding columns of 0s and 1s for each treatment
g.herbdf$U <- g.herbdf$G <- 0
g.herbdf[g.herbdf$graz.trt=="U",]$U <- 1
g.herbdf[g.herbdf$graz.trt=="G",]$G <- 1

fixra<-g.herbdf$fix.ra
p <- length(fixra[fixra==0])/length(fixra) #gets the probability of NO FIXATION
u<-mean(log(fixra[fixra>0])) # log space mean of SNF rate if there fixation > 0
s <- sd(log(fixra[fixra>0])) # log space sd of SNF rate if there fixation > 0
knownmedian <- (1-p)*exp(u)

# Starting values
pstart <- 0.25
ustart <- log(5)
sstart <- log(3)

#First we'll build a null model that fits the same mean to both grazed and ungrazed plots
fixra_null <- function(sdlogfix,ulogfix,p_no_fix,fix_dat){
  J <- sum(fixra==0)
  K <- sum(fixra>0)
  -J*log(p_no_fix) - K*log(1-p_no_fix) - sum(dlnorm(fix_dat[fix_dat>0],mean=ulogfix,sd=exp(sdlogfix),log=TRUE))
}
fit_fixra_null <- mle2(fixra_null,
                       start=list(sdlogfix=sstart,ulogfix=ustart,p_no_fix=pstart),
                       data=list(fix_dat=fixra))
print(summary(fit_fixra_null))

#Now a model that only fits a different MEAN for the different GRAZING treatments
fixra_u_G <- function(u_u, u_g, U, G, sdlogfix, p_no_fix,fix_dat){
  u_vec <- u_u*U + u_g*G
  J <- sum(fix_dat==0)
  K <- sum(fix_dat>0)
  -J*log(p_no_fix) - K*log(1-p_no_fix) - sum(dlnorm(fix_dat[fix_dat>0],mean=u_vec,sd=exp(sdlogfix),log=TRUE))
}
fit_fixra_u_G <- mle2(fixra_u_G,
                      start=list(sdlogfix=sstart,u_u=ustart,u_g=ustart,p_no_fix=pstart),
                      data=list(fix_dat=g.herbdf$fix.ra, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixra_u_G))

#Now a model where PROBABILITY OF FIXING is different between GRAZING treatments
fixra_p_G <- function(p_u, p_g, U, G, sdlogfix, ulogfix,fix_dat){
  p_vec <- p_u*U + p_g*G
  -sum(log(p_vec[fixra==0])) - sum(log(1-p_vec[fixra>0])) - sum(dlnorm(fix_dat[fix_dat>0],mean=ulogfix,sd=exp(sdlogfix),log=TRUE))
}
fit_fixra_p_G <- mle2(fixra_p_G,
                      start=list(sdlogfix=sstart,p_u=pstart,p_g=pstart,ulogfix=ustart),
                      data=list(fix_dat=g.herbdf$fix.ra, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixra_p_G))

#Now a model where BOTH PROBABILITY AND MEAN are different between GRAZING treatments
fixra_up_G <- function(u_u, u_g, p_u, p_g, U, G, sdlogfix, p_no_fix,fix_dat){
  u_vec <- u_u*U + u_g*G
  p_vec <- p_u*U + p_g*G
  -sum(log(p_vec[fixra==0])) - sum(log(1-p_vec[fixra>0]))- sum(dlnorm(fix_dat[fix_dat>0],mean=u_vec,sd=exp(sdlogfix),log=TRUE))
}
fit_fixra_up_G <- mle2(fixra_up_G,
                       start=list(sdlogfix=sstart,u_u=ustart,u_g=ustart,p_u=pstart, p_g=pstart),
                       data=list(fix_dat=g.herbdf$fix.ra, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixra_up_G))

deltaAICs(c(AICc(fit_fixra_null,fixra),AICc(fit_fixra_u_G,fixra),AICc(fit_fixra_p_G,fixra),AICc(fit_fixra_up_G,fixra)))
# It looks like the null model is the best - so no effect of grazers on either the probability of N-fixers being present
# or the mean relative abundance of N fixers when they are present

##########################################################################
# Running MLE models to test for effect of Herbivory on fixer relative species richness assuming 0-inflated lognormal data

fixrr<-g.herbdf$fix.rr
p <- length(fixrr[fixrr==0])/length(fixrr) #gets the probability of NO FIXERS
u<-mean(log(fixrr[fixrr>0])) # log space mean of relative richness if fixrr > 0
s <- sd(log(fixrr[fixrr>0])) # log space sd of relative richness if there fixrr > 0
knownmedian <- (1-p)*exp(u)

# Starting values
pstart <- 0.25
ustart <- log(6)
sstart <- log(2)

#First we'll build a null model that fits the same mean to both grazed and ungrazed plots
fixrr_null <- function(sdlogfix,ulogfix,p_no_fix,fix_dat){
  J <- sum(fixrr==0)
  K <- sum(fixrr>0)
  -J*log(p_no_fix) - K*log(1-p_no_fix) - sum(dlnorm(fix_dat[fix_dat>0],mean=ulogfix,sd=exp(sdlogfix),log=TRUE))
}
fit_fixrr_null <- mle2(fixrr_null,
                       start=list(sdlogfix=sstart,ulogfix=ustart,p_no_fix=pstart),
                       data=list(fix_dat=fixrr))
print(summary(fit_fixrr_null))

#Now a model that only fits a different MEAN for the different GRAZING treatments
fixrr_u_G <- function(u_u, u_g, U, G, sdlogfix, p_no_fix,fix_dat){
  u_vec <- u_u*U + u_g*G
  J <- sum(fix_dat==0)
  K <- sum(fix_dat>0)
  -J*log(p_no_fix) - K*log(1-p_no_fix) - sum(dlnorm(fix_dat[fix_dat>0],mean=u_vec,sd=exp(sdlogfix),log=TRUE))
}
fit_fixrr_u_G <- mle2(fixrr_u_G,
                      start=list(sdlogfix=sstart,u_u=ustart,u_g=ustart,p_no_fix=pstart),
                      data=list(fix_dat=g.herbdf$fix.rr, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixrr_u_G))

#Now a model where PROBABILITY OF FIXING is different between GRAZING treatments
fixrr_p_G <- function(p_u, p_g, U, G, sdlogfix, ulogfix,fix_dat){
  p_vec <- p_u*U + p_g*G
  -sum(log(p_vec[fixra==0])) - sum(log(1-p_vec[fixra>0])) - sum(dlnorm(fix_dat[fix_dat>0],mean=ulogfix,sd=exp(sdlogfix),log=TRUE))
}
fit_fixrr_p_G <- mle2(fixrr_p_G,
                      start=list(sdlogfix=sstart,p_u=pstart,p_g=pstart,ulogfix=ustart),
                      data=list(fix_dat=g.herbdf$fix.rr, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixrr_p_G))

#Now a model where BOTH PROBABILITY AND MEAN are different between GRAZING treatments
fixrr_up_G <- function(u_u, u_g, p_u, p_g, U, G, sdlogfix, p_no_fix,fix_dat){
  u_vec <- u_u*U + u_g*G
  p_vec <- p_u*U + p_g*G
  -sum(log(p_vec[fixra==0])) - sum(log(1-p_vec[fixra>0]))- sum(dlnorm(fix_dat[fix_dat>0],mean=u_vec,sd=exp(sdlogfix),log=TRUE))
}
fit_fixrr_up_G <- mle2(fixrr_up_G,
                       start=list(sdlogfix=sstart,u_u=ustart,u_g=ustart,p_u=pstart, p_g=pstart),
                       data=list(fix_dat=g.herbdf$fix.rr, U=g.herbdf$U, G=g.herbdf$G))
print(summary(fit_fixrr_up_G))

deltaAICs(c(AICc(fit_fixrr_null,fixrr),AICc(fit_fixrr_u_G,fixrr),AICc(fit_fixrr_p_G,fixrr),AICc(fit_fixrr_up_G,fixrr)))
# It looks like the null model is the best - so no effect of grazers on either the probability of N-fixers being present
# or the mean relative richness of N fixers when they are present

