# make sure you have bisoreg installed
# install.packages("bisoreg")
# library("bisoreg")


#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("F:/TRENDS/IRW TRENDS/IRW/Watts/")
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


NO3test5to30 <- data.frame(NO3test5to20,NO3test21to30)
boxplot(NO3test5to30)
NO3test5to30modes <- apply(NO3test5to30, 2, Mode)
NO3test5to30stdevs <- apply(NO3test5to30, 2, sd)

NO3test5to30stats <- data.frame((5:30), NO3test5to30modes, NO3test5to30stdevs)

write.csv(NO3test5to30, "NO35to30.csv")



################################################################################
####REPEAT FOR TN####


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


######### do trend analysis with custom span value
lnq.lims=range(lnTNsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTNsubs$lnQ,lnTNsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnTNsubs$lnQ
FUNlnC <- lnTNsubs$lnC



TNtest5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(TNtest5to20)

TNtest21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


TNtest5to30 <- data.frame(TNtest5to20,TNtest21to30)
boxplot(TNtest5to30)
TNtest5to30modes <- apply(TNtest5to30, 2, Mode)
TNtest5to30stdevs <- apply(TNtest5to30, 2, sd)

TNtest5to30stats <- data.frame((5:30), TNtest5to30modes, TNtest5to30stdevs)

write.csv(TNtest5to30, "TN5to30.csv")




################################################################################
####################REPEAT FOR SRP #############################################




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


######### do trend analysis with custom span value
lnq.lims=range(lnSRPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSRPsubs$lnQ,lnSRPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnSRPsubs$lnQ
FUNlnC <- lnSRPsubs$lnC



SRPtest5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(SRPtest5to20)

SRPtest21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


SRPtest5to30 <- data.frame(SRPtest5to20,SRPtest21to30)
boxplot(SRPtest5to30)
SRPtest5to30modes <- apply(SRPtest5to30, 2, Mode)
SRPtest5to30stdevs <- apply(SRPtest5to30, 2, sd)

SRPtest5to30stats <- data.frame((5:30), SRPtest5to30modes, SRPtest5to30stdevs)

write.csv(SRPtest5to30, "SRP5to30.csv")




###############################################################################
#################### REPEAT FOR TP ############################################




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


######### do trend analysis with custom span value
lnq.lims=range(lnTPsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTPsubs$lnQ,lnTPsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnTPsubs$lnQ
FUNlnC <- lnTPsubs$lnC



TPtest5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(TPtest5to20)

TPtest21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


TPtest5to30 <- data.frame(TPtest5to20,TPtest21to30)
boxplot(TPtest5to30)
TPtest5to30modes <- apply(TPtest5to30, 2, Mode)
TPtest5to30stdevs <- apply(TPtest5to30, 2, sd)

TPtest5to30stats <- data.frame((5:30), TPtest5to30modes, TPtest5to30stdevs)

write.csv(TPtest5to30, "TP5to30.csv")





######################################################################################
################# REPEAT FOR TSS #####################################################



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


######### do trend analysis with custom span value
lnq.lims=range(lnTSSsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnTSSsubs$lnQ,lnTSSsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnTSSsubs$lnQ
FUNlnC <- lnTSSsubs$lnC



TSStest5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(TSStest5to20)

TSStest21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


TSStest5to30 <- data.frame(TSStest5to20,TSStest21to30)
boxplot(TSStest5to30)
TSStest5to30modes <- apply(TSStest5to30, 2, Mode)
TSStest5to30stdevs <- apply(TSStest5to30, 2, sd)

TSStest5to30stats <- data.frame((5:30), TSStest5to30modes, TSStest5to30stdevs)

write.csv(TSStest5to30, "TSS5to30.csv")



##########################################################################################
################## REPEAT FOR Cl #########################################################





#make separate files for each constituent + flow
Cltemp<-data.frame(concentrations$Date,concentrations$Cl)
names(Cltemp)[names(Cltemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for Cl is 0.034 mg/L
Cltemp$concentrations.Cl <- ifelse(Cltemp$concentrations.Cl == 0, 0.034, Cltemp$concentrations.Cl)
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


######### do trend analysis with custom span value
lnq.lims=range(lnClsubs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnClsubs$lnQ,lnClsubs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnClsubs$lnQ
FUNlnC <- lnClsubs$lnC



Cltest5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(Cltest5to20)

Cltest21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


Cltest5to30 <- data.frame(Cltest5to20,Cltest21to30)
boxplot(Cltest5to30)
Cltest5to30modes <- apply(Cltest5to30, 2, Mode)
Cltest5to30stdevs <- apply(Cltest5to30, 2, sd)

Cltest5to30stats <- data.frame((5:30), Cltest5to30modes, Cltest5to30stdevs)

write.csv(Cltest5to30, "Cl5to30.csv")




######################################################################################
################# ONE MO TIME FOR SO4 ################################################





#make separate files for each constituent + flow
SO4temp<-data.frame(concentrations$Date,concentrations$SO4)
names(SO4temp)[names(SO4temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SO4 is 0.038 mg/L
SO4temp$concentrations.SO4 <- ifelse(SO4temp$concentrations.SO4 == 0, 0.038, SO4temp$concentrations.SO4)
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


######### do trend analysis with custom span value
lnq.lims=range(lnSO4subs$lnQ)
lnq.grid=seq(from=lnq.lims[1],to=lnq.lims[2])
plot(lnSO4subs$lnQ,lnSO4subs$lnC,xlim=lnq.lims,cex=0.5,col="darkgrey")
title("Special LOESS")

FUNlnQ <- lnSO4subs$lnQ
FUNlnC <- lnSO4subs$lnC



SO4test5to20 <- foldsfunction(5,20)
#look at boxplot
boxplot(SO4test5to20)

SO4test21to30 <- foldsfunction(21,30)


#make sure Mode is capitalized


SO4test5to30 <- data.frame(SO4test5to20,SO4test21to30)
boxplot(SO4test5to30)
SO4test5to30modes <- apply(SO4test5to30, 2, Mode)
SO4test5to30stdevs <- apply(SO4test5to30, 2, sd)

SO4test5to30stats <- data.frame((5:30), SO4test5to30modes, SO4test5to30stdevs)

write.csv(SO4test5to30, "SO45to30.csv")

############################################

fiveto30foldsSTATS <- data.frame(NO3test5to30stats, TNtest5to30stats, SRPtest5to30stats, TPtest5to30stats, TSStest5to30stats, Cltest5to30stats, SO4test5to30stats)
write.csv(fiveto30foldsSTATS, "fiveto30foldsSTATS.csv")

############################

#we gonna try to make a 'long' dataset of all the constituents and all the folds (with all 50 iterations)
#hwatever


#reshape2 package
TESTNO3test5to30 <- NO3test5to30
TESTNO3test5to30$id <- 1:nrow(TESTNO3test5to30)
meltNO3 <- melt(TESTNO3test5to30, id.vars="id", variable.name="folds", value.name="span")
meltNO3$const <- "NO3"

TESTTNtest5to30 <- TNtest5to30
TESTTNtest5to30$id <- 1:nrow(TESTTNtest5to30)
meltTN <- melt(TESTTNtest5to30, id.vars="id", variable.name="folds", value.name="span")
meltTN$const <- "TN"

TESTSRPtest5to30 <- SRPtest5to30
TESTSRPtest5to30$id <- 1:nrow(TESTSRPtest5to30)
meltSRP <- melt(TESTSRPtest5to30, id.vars="id", variable.name="folds", value.name="span")
meltSRP$const <- "SRP"

TESTTPtest5to30 <- TPtest5to30
TESTTPtest5to30$id <- 1:nrow(TESTTPtest5to30)
meltTP <- melt(TESTTPtest5to30, id.vars="id", variable.name="folds", value.name="span")
meltTP$const <- "TP"

TESTTSStest5to30 <- TSStest5to30
TESTTSStest5to30$id <- 1:nrow(TESTTSStest5to30)
meltTSS <- melt(TESTTSStest5to30, id.vars="id", variable.name="folds", value.name="span")
meltTSS$const <- "TSS"

TESTCltest5to30 <- Cltest5to30
TESTCltest5to30$id <- 1:nrow(TESTCltest5to30)
meltCl <- melt(TESTCltest5to30, id.vars="id", variable.name="folds", value.name="span")
meltCl$const <- "Cl"

TESTSO4test5to30 <- SO4test5to30
TESTSO4test5to30$id <- 1:nrow(TESTSO4test5to30)
meltSO4 <- melt(TESTSO4test5to30, id.vars="id", variable.name="folds", value.name="span")
meltSO4$const <- "SO4"

THEWHOLEDAMNTHING <- data.frame(meltNO3, meltTN, meltSRP, meltTP, meltTSS, meltCl, meltSO4)
THEWHOLEDAMNTHING$REALID <- 1:nrow(THEWHOLEDAMNTHING)

write.csv(THEWHOLEDAMNTHING, "THEWHOLEDAMNTHING.csv")


MELT_5to30FOLDS <- melt(THEWHOLEDAMNTHING, id.vars="REALID", variable.name="const", value.name="span")


