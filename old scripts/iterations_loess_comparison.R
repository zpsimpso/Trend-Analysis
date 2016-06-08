#compare 10 folds of 10 iterations to 10 folds of 1000 iterations
#UWRB EDITION

#NEW FUNCTION - loess.wrapperMSE
#exactly like original, but evaluates fit based on MSE, not MAE
#can return the fit, span selected, and MSE of that span

loess.wrapperMSE <- function(x, y, span.vals = seq(0.1, 1, by = 0.05), folds, iteration){
  mse <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, span = span)
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals) {
    ii <- ii + 1
    #make the kfolds procedure consistent across all spans
    set.seed(iteration) #can also change to a constant if preferred
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    mse[ii] <- mean((y[fltr] - y.cv[fltr])^2)
  }
  span <- span.vals[which.min(mse)]
  out <- loess(y ~ x, span = span)
  return(list(fit = out, span = span, MSE = (min(mse)), MSE.list = mse))
}


customfitSPAN <- function(kfolds, iterations) {
  y<-numeric(length(iterations))
  for (i in 1:iterations) {
    fitspecial.wrp <- loess.wrapperMSE(FUNlnQ,FUNlnC,folds=kfolds, iteration = i)
    iter <- fitspecial.wrp$span
    y[i]<-iter
    print(i)
  }
  median_f<-median(y)
  return(median_f)
}


#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers
setwd("F:/TRENDS/UWRB TRENDS/UWRB/Kings/2009-2015, no 2010 data/")
getwd()
SITE<-'Kings'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")

###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/RC45/")
getwd()
SITE<-'RC45'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/TB62/")
getwd()
SITE<-'TB62'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/TBT/")
getwd()
SITE<-'TBT'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WEC/")
getwd()
SITE<-'WEC'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WFWR/")
getwd()
SITE<-'WFWR'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")

###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WR45/")
getwd()
SITE<-'WR45'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
###################################################################
##################################################################
###################################################################
setwd("F:/TRENDS/UWRB TRENDS/UWRB/Wyman/")
getwd()
SITE<-'Wyman'

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

FUNlnQ<-lnNO3subs$lnQ
FUNlnC<-lnNO3subs$lnC

NO3_10_span <- customfitSPAN(10, 10)
NO3_1000_span <- customfitSPAN(35, 10)

############################ TN
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

FUNlnQ<-lnTNsubs$lnQ
FUNlnC<-lnTNsubs$lnC

TN_10_span <- customfitSPAN(10, 10)
TN_1000_span <- customfitSPAN(35, 10)

############################################## SRP
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

FUNlnQ<-lnSRPsubs$lnQ
FUNlnC<-lnSRPsubs$lnC

SRP_10_span <- customfitSPAN(10, 10)
SRP_1000_span <- customfitSPAN(35, 10)


################################### TP
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

FUNlnQ<-lnTPsubs$lnQ
FUNlnC<-lnTPsubs$lnC

TP_10_span <- customfitSPAN(10, 10)
TP_1000_span <- customfitSPAN(35, 10)

########################### TSS
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

FUNlnQ<-lnTSSsubs$lnQ
FUNlnC<-lnTSSsubs$lnC

TSS_10_span <- customfitSPAN(10, 10)
TSS_1000_span <- customfitSPAN(35, 10)

################################# Cl
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

FUNlnQ<-lnClsubs$lnQ
FUNlnC<-lnClsubs$lnC

Cl_10_span <- customfitSPAN(10, 10)
Cl_1000_span <- customfitSPAN(35, 10)

############################ SO4
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

FUNlnQ<-lnSO4subs$lnQ
FUNlnC<-lnSO4subs$lnC

SO4_10_span <- customfitSPAN(10, 10)
SO4_1000_span <- customfitSPAN(35, 10)

#############################################
#make site data frame with both spans across all constituents
SITE

df <- data.frame(matrix(NA, ncol = 4))
colnames(df) <- c('Site', 'Const', 'I_10', 'I_1000')
df[1:7,1]<-SITE
df[1:7,2] <- c('NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', 'SO4')
df[1,3]<-NO3_10_span
df[1,4]<-NO3_1000_span
df[2,3]<-TN_10_span
df[2,4]<-TN_1000_span
df[3,3]<-SRP_10_span
df[3,4]<-SRP_1000_span
df[4,3]<-TP_10_span
df[4,4]<-TP_1000_span
df[5,3]<-TSS_10_span
df[5,4]<-TSS_1000_span
df[6,3]<-Cl_10_span
df[6,4]<-Cl_1000_span
df[7,3]<-SO4_10_span
df[7,4]<-SO4_1000_span

write.csv(df, "TENBYFOLD_COMPARISON.csv")
