
#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("G:/UWRB/TBT/")
getwd()

#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

#make estimation file
estfile <- data.frame(siteflow$Date,1200,siteflow$Flow)
colnames(estfile)[1] <- "Date"
colnames(estfile)[2] <- "Time"
colnames(estfile)[3] <- "FullFlow"

#replace 0 flows
estfile$FullFlow<-ifelse(estfile$FullFlow == 0, 0.001, estfile$FullFlow)

#reformat dates
estfile$Date <- as.character(as.Date(estfile$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(estfile, "est_file.txt", row.names = FALSE, quote = FALSE, sep="\t")




######################## NITRATE  ########################




#make separate files for each constituent + flow
NO3temp<-data.frame(concentrations$Date,concentrations$NO3)
names(NO3temp)[names(NO3temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for NO3 is 0.004 mg/L
#let LOADEST use <values for MLE
NO3temp$concentrations.NO3 <- ifelse(NO3temp$concentrations.NO3 == 0, "<0.004", NO3temp$concentrations.NO3)
NO3temp<-na.omit(NO3temp)

NO3<-merge(NO3temp,siteflow, by="Date",all.NO3temp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
NO3$Flow<-ifelse(NO3$Flow == 0, 0.001, NO3$Flow)

#make calib file
calib_NO3_temp <- data.frame(NO3$Date,NO3$Flow,NO3$concentrations.NO3)

#add 'Time' column
calib_NO3 <- data.frame(calib_NO3_temp[,1],1200,calib_NO3_temp[,2:3])
colnames(calib_NO3)[2] <- "Time"
colnames(calib_NO3)[1] <- "Date"
#omit NA's
calib_NO3 <- na.omit(calib_NO3)
#reformat dates
calib_NO3$Date <- as.character(as.Date(calib_NO3$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_NO3, "NO3_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")


################# TOTAL NITROGEN #########################


#make separate files for each constituent + flow
TNtemp<-data.frame(concentrations$Date,concentrations$TN)
names(TNtemp)[names(TNtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TN is 0.012 mg/L
#let LOADEST use <values for MLE
TNtemp$concentrations.TN <- ifelse(TNtemp$concentrations.TN == 0, "<0.012", TNtemp$concentrations.TN)
TNtemp<-na.omit(TNtemp)

TN<-merge(TNtemp,siteflow, by="Date",all.TNtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TN$Flow<-ifelse(TN$Flow == 0, 0.001, TN$Flow)

#make calib file
calib_TN_temp <- data.frame(TN$Date,TN$Flow,TN$concentrations.TN)

#add 'Time' column
calib_TN <- data.frame(calib_TN_temp[,1],1200,calib_TN_temp[,2:3])
colnames(calib_TN)[2] <- "Time"
colnames(calib_TN)[1] <- "Date"
#omit NA's
calib_TN <- na.omit(calib_TN)
#reformat dates
calib_TN$Date <- as.character(as.Date(calib_TN$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_TN, "TN_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")




################# SOLUBLE REACTIVE PHOSPHORUS #########################


#make separate files for each constituent + flow
SRPtemp<-data.frame(concentrations$Date,concentrations$SRP)
names(SRPtemp)[names(SRPtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SRP is 0.002 mg/L
#let LOADEST use <values for MLE
SRPtemp$concentrations.SRP <- ifelse(SRPtemp$concentrations.SRP == 0, "<0.002", SRPtemp$concentrations.SRP)
SRPtemp<-na.omit(SRPtemp)

SRP<-merge(SRPtemp,siteflow, by="Date",all.SRPtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
SRP$Flow<-ifelse(SRP$Flow == 0, 0.001, SRP$Flow)

#make calib file
calib_SRP_temp <- data.frame(SRP$Date,SRP$Flow,SRP$concentrations.SRP)

#add 'Time' column
calib_SRP <- data.frame(calib_SRP_temp[,1],1200,calib_SRP_temp[,2:3])
colnames(calib_SRP)[2] <- "Time"
colnames(calib_SRP)[1] <- "Date"
#omit NA's
calib_SRP <- na.omit(calib_SRP)
#reformat dates
calib_SRP$Date <- as.character(as.Date(calib_SRP$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_SRP, "SRP_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")




################# TOTAL PHOSPHORUS #########################


#make separate files for each constituent + flow
TPtemp<-data.frame(concentrations$Date,concentrations$TP)
names(TPtemp)[names(TPtemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TP is 0.003 mg/L
#let LOADEST use <values for MLE
TPtemp$concentrations.TP <- ifelse(TPtemp$concentrations.TP == 0, "<0.003", TPtemp$concentrations.TP)
TPtemp<-na.omit(TPtemp)

TP<-merge(TPtemp,siteflow, by="Date",all.TPtemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TP$Flow<-ifelse(TP$Flow == 0, 0.001, TP$Flow)

#make calib file
calib_TP_temp <- data.frame(TP$Date,TP$Flow,TP$concentrations.TP)

#add 'Time' column
calib_TP <- data.frame(calib_TP_temp[,1],1200,calib_TP_temp[,2:3])
colnames(calib_TP)[2] <- "Time"
colnames(calib_TP)[1] <- "Date"
#omit NA's
calib_TP <- na.omit(calib_TP)
#reformat dates
calib_TP$Date <- as.character(as.Date(calib_TP$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_TP, "TP_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")



################# TOTAL SUSPENDED SOLIDS #########################


#make separate files for each constituent + flow
TSStemp<-data.frame(concentrations$Date,concentrations$TSS)
names(TSStemp)[names(TSStemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for TSS is 2.880 mg/L
#let LOADEST use <values for MLE
TSStemp$concentrations.TSS <- ifelse(TSStemp$concentrations.TSS == 0, "<2.880", TSStemp$concentrations.TSS)
TSStemp<-na.omit(TSStemp)

TSS<-merge(TSStemp,siteflow, by="Date",all.TSStemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
TSS$Flow<-ifelse(TSS$Flow == 0, 0.001, TSS$Flow)

#make calib file
calib_TSS_temp <- data.frame(TSS$Date,TSS$Flow,TSS$concentrations.TSS)

#add 'Time' column
calib_TSS <- data.frame(calib_TSS_temp[,1],1200,calib_TSS_temp[,2:3])
colnames(calib_TSS)[2] <- "Time"
colnames(calib_TSS)[1] <- "Date"
#omit NA's
calib_TSS <- na.omit(calib_TSS)
#reformat dates
calib_TSS$Date <- as.character(as.Date(calib_TSS$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_TSS, "TSS_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")



################# CHLORIDE #########################


#make separate files for each constituent + flow
Cltemp<-data.frame(concentrations$Date,concentrations$Cl)
names(Cltemp)[names(Cltemp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for Cl is 0.034 mg/L
#let LOADEST use <values for MLE
Cltemp$concentrations.Cl <- ifelse(Cltemp$concentrations.Cl == 0, "<0.034", Cltemp$concentrations.Cl)
Cltemp<-na.omit(Cltemp)

Cl<-merge(Cltemp,siteflow, by="Date",all.Cltemp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
Cl$Flow<-ifelse(Cl$Flow == 0, 0.001, Cl$Flow)

#make calib file
calib_Cl_temp <- data.frame(Cl$Date,Cl$Flow,Cl$concentrations.Cl)

#add 'Time' column
calib_Cl <- data.frame(calib_Cl_temp[,1],1200,calib_Cl_temp[,2:3])
colnames(calib_Cl)[2] <- "Time"
colnames(calib_Cl)[1] <- "Date"
#omit NA's
calib_Cl <- na.omit(calib_Cl)
#reformat dates
calib_Cl$Date <- as.character(as.Date(calib_Cl$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_Cl, "Cl_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")



################# SULFATE #########################


#make separate files for each constituent + flow
SO4temp<-data.frame(concentrations$Date,concentrations$SO4)
names(SO4temp)[names(SO4temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SO4 is 0.038 mg/L
#let LOADEST use <values for MLE
SO4temp$concentrations.SO4 <- ifelse(SO4temp$concentrations.SO4 == 0, "<0.038", SO4temp$concentrations.SO4)
SO4temp<-na.omit(SO4temp)

SO4<-merge(SO4temp,siteflow, by="Date",all.SO4temp=TRUE,all.siteflow=FALSE,sort=TRUE)
#in case of sites that can have 0 flow, replace 0 values with 0.001 cfs
SO4$Flow<-ifelse(SO4$Flow == 0, 0.001, SO4$Flow)

#make calib file
calib_SO4_temp <- data.frame(SO4$Date,SO4$Flow,SO4$concentrations.SO4)

#add 'Time' column
calib_SO4 <- data.frame(calib_SO4_temp[,1],1200,calib_SO4_temp[,2:3])
colnames(calib_SO4)[2] <- "Time"
colnames(calib_SO4)[1] <- "Date"
#omit NA's
calib_SO4 <- na.omit(calib_SO4)
#reformat dates
calib_SO4$Date <- as.character(as.Date(calib_SO4$Date, "%Y/%m/%d", tz= "GMT"),"%Y%m%d")

#write out the calib file
write.table(calib_SO4, "SO4_calib.txt", row.names = FALSE, quote = FALSE, sep="\t")
















