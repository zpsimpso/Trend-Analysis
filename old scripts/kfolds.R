# install.packages("bisoreg")
# library("bisoreg")

install.packages("bisoreg", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("bisoreg", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
library("bisoreg", lib.loc="~/Documents/R_PACKAGES")

#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/UWRB/Kings")
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

######### do trend analysis with custom span value
lnq.lims=range(lnNO3subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnNO3subs$lnQ,lnNO3subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

customfit <- function(kfolds, iterations) {
  y <-vector(mode="numeric", length=iterations)
  for (i in 1:iterations)
  
  fitspecial.wrp <- loess.wrapper(lnNO3subs$lnQ,lnNO3subs$lnC,span.vals=seq(0.2,0.8,by=0.05),folds=kfolds)
  customspan <- fitspecial.wrp$pars$span
  y[i] <- customspan
  
}

try5 <- customfit(5, 20)


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
