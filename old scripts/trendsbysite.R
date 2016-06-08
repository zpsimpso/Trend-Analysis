# install.packages("bisoreg")
# library("bisoreg")

install.packages("bisoreg", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("bisoreg", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
library("bisoreg", lib.loc="~/Documents/R_PACKAGES")

#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("F:/TRENDS/UWRB TRENDS/UWRB/Kings/")
getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
NO3temp<-data.frame(concentrations$Date,concentrations$NO3)
names(NO3temp)[names(NO3temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for NO3 is 0.004 mg/L
NO3temp$concentrations.NO3 <- ifelse(NO3temp$concentrations.NO3 == 0, 0.004, NO3temp$concentrations.NO3)
NO3temp<-na.omit(NO3temp)



NO3<-merge(NO3temp,siteflow, by="Date",all.NO3temp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
NO3$Flow<-ifelse(NO3$Flow == 0, 0.001, NO3$Flow)
#note that the log function is actually natural log
lnC<-log(NO3$concentrations.NO3)
lnQ<-log(NO3$Flow)
lnNO3<-data.frame(NO3$Date,NO3$Flow,NO3$concentrations.NO3,lnC,lnQ)
#make a subset that omits the NA's
lnNO3subs<-na.omit(lnNO3)

write.csv(NO3,"NO3.csv")
write.csv(lnNO3subs,"lnNO3subs.csv")

#start trend analysis process

lnq.lims=range(lnNO3subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnNO3subs$lnQ,lnNO3subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnNO3subs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$NO3.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))


######### repeat but use 'special' span value
lnq.lims=range(lnNO3subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnNO3subs$lnQ,lnNO3subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnNO3subs$lnQ,lnNO3subs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))




####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

NO3results.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnNO3subs$NO3.Flow,Conc=lnNO3subs$NO3.concentrations.NO3,lnQ=lnNO3subs$lnQ,lnC=lnNO3subs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(NO3results.frame,"NO3results.csv")




########################################################################
######  REPEAT FOR TN  #################################

#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
TNtemp<-data.frame(concentrations$Date,concentrations$TN)
names(TNtemp)[names(TNtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TN is 0.012 mg/L
TNtemp$concentrations.TN <- ifelse(TNtemp$concentrations.TN == 0, 0.012, TNtemp$concentrations.TN)
TNtemp<-na.omit(TNtemp)

TN<-merge(TNtemp,siteflow, by="Date",all.TNtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TN$Flow<-ifelse(TN$Flow == 0, 0.001, TN$Flow)

#note that the log function is actually natural log
lnC<-log(TN$concentrations.TN)
lnQ<-log(TN$Flow)
lnTN<-data.frame(TN$Date,TN$Flow,TN$concentrations.TN,lnC,lnQ)
#make a subset that omits the NA's
lnTNsubs<-na.omit(lnTN)

write.csv(TN,"TN.csv")
write.csv(lnTNsubs,"lnTNsubs.csv")

#start trend analysis process

lnq.lims=range(lnTNsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTNsubs$lnQ,lnTNsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnTNsubs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$TN.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))


######### repeat but use 'special' span value
lnq.lims=range(lnTNsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTNsubs$lnQ,lnTNsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnTNsubs$lnQ,lnTNsubs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))




####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
x[1:len]
}
makePaddedDataFrame <- function(l,...){
maxlen <- max(sapply(l,length))
data.frame(lapply(l,na.pad,len=maxlen),...)
}



TNresults.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnTNsubs$TN.Flow,Conc=lnTNsubs$TN.concentrations.TN,lnQ=lnTNsubs$lnQ,lnC=lnTNsubs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(TNresults.frame,"TNresults.csv")


#################################################################################################
#####################  REPEAT FOR SRP ###########################################################

#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
SRPtemp<-data.frame(concentrations$Date,concentrations$SRP)
names(SRPtemp)[names(SRPtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SRP is 0.002 mg/L
SRPtemp$concentrations.SRP <- ifelse(SRPtemp$concentrations.SRP == 0, 0.002, SRPtemp$concentrations.SRP)
SRPtemp<-na.omit(SRPtemp)

SRP<-merge(SRPtemp,siteflow, by="Date",all.SRPtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
SRP$Flow<-ifelse(SRP$Flow == 0, 0.001, SRP$Flow)


#note that the log function is actually natural log
lnC<-log(SRP$concentrations.SRP)
lnQ<-log(SRP$Flow)
lnSRP<-data.frame(SRP$Date,SRP$Flow,SRP$concentrations.SRP,lnC,lnQ)
#make a subset that omits the NA's
lnSRPsubs<-na.omit(lnSRP)

write.csv(SRP,"SRP.csv")
write.csv(lnSRPsubs,"lnSRPsubs.csv")

#start trend analysis process

lnq.lims=range(lnSRPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSRPsubs$lnQ,lnSRPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnSRPsubs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$SRP.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))

######### repeat but use 'special' span value
lnq.lims=range(lnSRPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSRPsubs$lnQ,lnSRPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnSRPsubs$lnQ,lnSRPsubs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))

####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}



SRPresults.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnSRPsubs$SRP.Flow,Conc=lnSRPsubs$SRP.concentrations.SRP,lnQ=lnSRPsubs$lnQ,lnC=lnSRPsubs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(SRPresults.frame,"SRPresults.csv")

##############################################################################
###############  REPEAT FOR TP  #######################################


#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
TPtemp<-data.frame(concentrations$Date,concentrations$TP)
names(TPtemp)[names(TPtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TP is 0.003 mg/L
TPtemp$concentrations.TP <- ifelse(TPtemp$concentrations.TP == 0, 0.003, TPtemp$concentrations.TP)
TPtemp<-na.omit(TPtemp)


TP<-merge(TPtemp,siteflow, by="Date",all.TPtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TP$Flow<-ifelse(TP$Flow == 0, 0.001, TP$Flow)

#note that the log function is actually natural log
lnC<-log(TP$concentrations.TP)
lnQ<-log(TP$Flow)
lnTP<-data.frame(TP$Date,TP$Flow,TP$concentrations.TP,lnC,lnQ)
#make a subset that omits the NA's
lnTPsubs<-na.omit(lnTP)

write.csv(TP,"TP.csv")
write.csv(lnTPsubs,"lnTPsubs.csv")

#start trend analysis process

lnq.lims=range(lnTPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTPsubs$lnQ,lnTPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnTPsubs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$TP.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))



######### repeat but use 'special' span value
lnq.lims=range(lnTPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTPsubs$lnQ,lnTPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnTPsubs$lnQ,lnTPsubs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))


####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}



TPresults.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnTPsubs$TP.Flow,Conc=lnTPsubs$TP.concentrations.TP,lnQ=lnTPsubs$lnQ,lnC=lnTPsubs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(TPresults.frame,"TPresults.csv")


##################################################################################################
###########################  REPEAT FOR TSS ######################################################




#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
TSStemp<-data.frame(concentrations$Date,concentrations$TSS)
names(TSStemp)[names(TSStemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TSS is 2.880 mg/L
TSStemp$concentrations.TSS <- ifelse(TSStemp$concentrations.TSS == 0, 2.880, TSStemp$concentrations.TSS)
TSStemp<-na.omit(TSStemp)



TSS<-merge(TSStemp,siteflow, by="Date",all.TSStemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TSS$Flow<-ifelse(TSS$Flow == 0, 0.001, TSS$Flow)
#note that the log function is actually natural log
lnC<-log(TSS$concentrations.TSS)
lnQ<-log(TSS$Flow)
lnTSS<-data.frame(TSS$Date,TSS$Flow,TSS$concentrations.TSS,lnC,lnQ)
#make a subset that omits the NA's
lnTSSsubs<-na.omit(lnTSS)

write.csv(TSS,"TSS.csv")
write.csv(lnTSSsubs,"lnTSSsubs.csv")

#start trend analysis process

lnq.lims=range(lnTSSsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTSSsubs$lnQ,lnTSSsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnTSSsubs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$TSS.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))


######### repeat but use 'special' span value
lnq.lims=range(lnTSSsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTSSsubs$lnQ,lnTSSsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnTSSsubs$lnQ,lnTSSsubs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))


####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}



TSSresults.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnTSSsubs$TSS.Flow,Conc=lnTSSsubs$TSS.concentrations.TSS,lnQ=lnTSSsubs$lnQ,lnC=lnTSSsubs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(TSSresults.frame,"TSSresults.csv")


#############################################################################################
##################   REPEAT FOR Cl  #########################################################




#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
Cltemp<-data.frame(concentrations$Date,concentrations$Cl)
names(Cltemp)[names(Cltemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for Cl is 0.034 mg/L
Cltemp$concentrations.Cl <- ifelse(Cltemp$concentrations.Cl == 0, 0.034,Cltemp$concentrations.Cl)
Cltemp<-na.omit(Cltemp)

Cl<-merge(Cltemp,siteflow, by="Date",all.Cltemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
Cl$Flow<-ifelse(Cl$Flow == 0, 0.001, Cl$Flow)
#note that the log function is actually natural log
lnC<-log(Cl$concentrations.Cl)
lnQ<-log(Cl$Flow)
lnCl<-data.frame(Cl$Date,Cl$Flow,Cl$concentrations.Cl,lnC,lnQ)
#make a subset that omits the NA's
lnClsubs<-na.omit(lnCl)

write.csv(Cl,"Cl.csv")
write.csv(lnClsubs,"lnClsubs.csv")

#start trend analysis process

lnq.lims=range(lnClsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnClsubs$lnQ,lnClsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnClsubs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$Cl.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))


######### repeat but use 'special' span value
lnq.lims=range(lnClsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnClsubs$lnQ,lnClsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnClsubs$lnQ,lnClsubs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))



####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}



Clresults.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnClsubs$Cl.Flow,Conc=lnClsubs$Cl.concentrations.Cl,lnQ=lnClsubs$lnQ,lnC=lnClsubs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(Clresults.frame,"Clresults.csv")



########################################################################################################
######################## REPEAT FOR SO4  ###############################################################



#clear environment
rm(list=ls())


getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make separate files for each constituent + flow
SO4temp<-data.frame(concentrations$Date,concentrations$SO4)
names(SO4temp)[names(SO4temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SO4 is 0.038 mg/L
SO4temp$concentrations.SO4 <- ifelse(SO4temp$concentrations.SO4 == 0, 0.038,SO4temp$concentrations.SO4)
SO4temp<-na.omit(SO4temp)

SO4<-merge(SO4temp,siteflow, by="Date",all.SO4temp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
SO4$Flow<-ifelse(SO4$Flow == 0, 0.001, SO4$Flow)
#note that the log function is actually natural log
lnC<-log(SO4$concentrations.SO4)
lnQ<-log(SO4$Flow)
lnSO4<-data.frame(SO4$Date,SO4$Flow,SO4$concentrations.SO4,lnC,lnQ)
#make a subset that omits the NA's
lnSO4subs<-na.omit(lnSO4)

write.csv(SO4,"SO4.csv")
write.csv(lnSO4subs,"lnSO4subs.csv")

#start trend analysis process

lnq.lims=range(lnSO4subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSO4subs$lnQ,lnSO4subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Local Regression")

#make sure that there are no 0 values for concentrations
#those lead to inf/-inf values for lnC

dat<-read.csv("lnSO4subs.csv",header=TRUE)
attach(dat)

fit=loess(lnC~lnQ,span=0.5,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)

#get residuals
loess.residuals <- fit$residuals

#transform dates
newdates <- as.POSIXct(dat$SO4.Date)
trend.frame <- data.frame(newdates,loess.residuals)
#plot the trend
plot(newdates,loess.residuals)
FACreg <-  lm(loess.residuals~newdates)
abline(FACreg)
#check these sweet ass statistics on the trend
anova.trends <- anova(FACreg)
anova.trends

#get % change per year
#the slope coefficient on FACreg is in log units per second

m <- FACreg$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))



######### repeat but use 'special' span value
lnq.lims=range(lnSO4subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSO4subs$lnQ,lnSO4subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

#for the cross validation procedure, use set.seed to make the RNG consistent
#hail Satan

set.seed(666)

fitspecial.wrp<- loess.wrapper(lnSO4subs$lnQ,lnSO4subs$lnC,span.vals=seq(0.2,0.80,by=0.05),folds=10)
#plot loess line
spanloessspecial <- fitspecial.wrp$pars$span

fitspecial=loess(lnC~lnQ,span=spanloessspecial,data=dat)
lines(lnq.grid,predict(fit,data.frame(lnQ=lnq.grid)), col="red",lwd=2)
spanloessspecial

specialloess.residuals <- fitspecial$residuals
specialtrend.frame <- data.frame(newdates, specialloess.residuals)
plot(newdates, specialloess.residuals)
specialFACreg <- lm(specialloess.residuals~newdates)
abline(specialFACreg)
#check statistics
specialanova.trends <- anova(specialFACreg)
specialanova.trends

m2 <- specialFACreg$coefficients[2]
specialpercentslope <- (((exp(m2)-1)*60*60*24*365*100))


####write results into file?


####create functions to make data frame with columns of unequal lengths
na.pad <- function(x,len){
  x[1:len]
}
makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}



SO4results.frame<-data.frame(makePaddedDataFrame(list(Dates=newdates,Flow=lnSO4subs$SO4.Flow,Conc=lnSO4subs$SO4.concentrations.SO4,lnQ=lnSO4subs$lnQ,lnC=lnSO4subs$lnC,Reg_Residuals=loess.residuals,Special_Residuals=specialloess.residuals,Reg_Pvalue=anova.trends$"Pr(>F)",Special_Span=spanloessspecial,Special_Pvalue=specialanova.trends$"Pr(>F)",Reg_Perc_Change_yr=percentslope,Special_Perc_Change_yr=specialpercentslope)))



write.csv(SO4results.frame,"SO4results.csv")



