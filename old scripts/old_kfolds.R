# install.packages("bisoreg")
# library("bisoreg")

install.packages("bisoreg", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("bisoreg", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
library("bisoreg", lib.loc="~/Documents/R_PACKAGES")

#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("E:/TRENDS/UWRB TRENDS/UWRB/Kings/2009-2015, no 2010 data/")
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

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC


###################### FUNCTIONS
#create function to find custom span values at k folds for specified iterations

customfit <- function(kfolds, iterations) {
  y <-vector(mode="numeric", length=iterations)
  for (i in 1:iterations) {
    fitspecial.wrp <- loess.wrapper(FUNlnQ,FUNlnC,span.vals=seq(0.2,0.9,by=0.05),folds=kfolds)
    customspan <- fitspecial.wrp$pars$span
    y[i] <- customspan
    print(i)
  }
  return(y)
}

#SET NUMBER OF ITERATIONS HERE
iternum = 50


foldsfunction <- function(lowfolds, highfolds){
  
  
  ncol.cal =  (highfolds - lowfolds) +1 
 
  
  df <- data.frame(matrix(NA, nrow = iternum, ncol = ncol.cal))
  colnames(df)<-paste("folds", c(lowfolds:highfolds))
  
  
  for (k in lowfolds:highfolds) {
    data = as.data.frame(customfit(k, iternum))
    df[k+(1-lowfolds)] <- data
     print(k)
  }
  return(df)
}

#function for finding the mode of the data
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
##############################



NO3test5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(NO3test5to20)

NO3test21to30 <- foldsfunction(21,30)
# 
# par(mfrow=c(2,1))
# boxplot(test5to20)
# boxplot(test21to30)


#make sure Mode is capitalized
apply(NO3test5to20, 2, Mode)
apply(NO3test21to30, 2, Mode)

NO3test5to30 <- data.frame(NO3test5to20,NO3test21to30)
boxplot(NO3test5to30)
NO3test5to30modes <- apply(NO3test5to30, 2, Mode)
NO3test5to30stdevs <- apply(NO3test5to30, 2, sd)

NO3test5to30stats <- data.frame((5:30), NO3test5to30modes, NO3test5to30stdevs)

