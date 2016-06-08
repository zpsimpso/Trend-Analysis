

#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("F:/TRENDS/UWRB TRENDS/UWRB/Wyman/")
getwd()

site <- "WR at Wyman"
#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)

ReadInWQData <- function(siteflow, concentrations){
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
  
  ######################################
  
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
  #######################################
  
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
  #########################################
  
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
  
  ############################################
  
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
  
  ############################################
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
  
  ###############################################
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
  ############################################
  return(list(lnNO3subs, lnTNsubs, lnSRPsubs, lnTPsubs, lnTSSsubs, lnClsubs, lnSO4subs))
  
}

lnNO3subs <- ReadInWQData(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData(siteflow, concentrations)[[7]]
########################################################

#ggfilename <- paste(site,const,".pdf",sep="_")
#ggsave(ggfilename, TSS_FACplot, device="pdf", dpi=500)

#############################################################
const <- "NO3"
FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnNO3subs)
lnNO3subs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnNO3subs[7] <- fitmse$fit$residuals

NO3plot <- ggplot(data = lnNO3subs, aes(x=NO3.Date, y=NO3.concentrations.NO3))
NO3plot <- NO3plot + geom_point() + labs(x="Year", y=expression(NO3~(mg~L^{-1})))
# +annotate("text", x=lnNO3subs[30,1], y=750, size=9, label="A")
NO3plot #I had to use a data point from the dataframe as the x for the label

FA_NO3plot <- ggplot(data=lnNO3subs, aes(x=lnQ, y=lnC)) + geom_point()
FA_NO3plot <- FA_NO3plot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(NO3)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_NO3plot

Trendline <- lm(V7~NO3.Date, data=lnNO3subs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))


if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
#text=element_text(size=rel(0.9))

FACplot <- ggplot(data=lnNO3subs, aes(x=NO3.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y=expression(paste(NO[3]," FACs")))
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnNO3subs[30,1], xend=lnNO3subs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnNO3subs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnNO3subs[30,1], xend=lnNO3subs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnNO3subs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnNO3subs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnNO3subs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnNO3subs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnNO3subs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
NO3_FACplot <- FACplot
#############################################################
const <- "TN"
FUNlnQ <- lnTNsubs$lnQ
FUNlnC <- lnTNsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnTNsubs)
lnTNsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTNsubs[7] <- fitmse$fit$residuals

TNplot <- ggplot(data = lnTNsubs, aes(x=TN.Date, y=TN.concentrations.TN))
TNplot <- TNplot + geom_point() + labs(x="Year", y=expression(TN~(mg~L^{-1})))
# +annotate("text", x=lnTNsubs[30,1], y=750, size=9, label="A")
TNplot #I had to use a data point from the dataframe as the x for the label

FA_TNplot <- ggplot(data=lnTNsubs, aes(x=lnQ, y=lnC)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FA_TNplot <- FA_TNplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(TN)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_TNplot

Trendline <- lm(V7~TN.Date, data=lnTNsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnTNsubs, aes(x=TN.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y="TN FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTNsubs[30,1], xend=lnTNsubs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnTNsubs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnTNsubs[30,1], xend=lnTNsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnTNsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnTNsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnTNsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnTNsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnTNsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
TN_FACplot <- FACplot
#############################################################
const <- "SRP"
FUNlnQ <- lnSRPsubs$lnQ
FUNlnC <- lnSRPsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnSRPsubs)
lnSRPsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnSRPsubs[7] <- fitmse$fit$residuals

SRPplot <- ggplot(data = lnSRPsubs, aes(x=SRP.Date, y=SRP.concentrations.SRP))
SRPplot <- SRPplot + geom_point() + labs(x="Year", y=expression(SRP~(mg~L^{-1})))
# +annotate("text", x=lnSRPsubs[30,1], y=750, size=9, label="A")
SRPplot #I had to use a data point from the dataframe as the x for the label

FA_SRPplot <- ggplot(data=lnSRPsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_SRPplot <- FA_SRPplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(SRP)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_SRPplot

Trendline <- lm(V7~SRP.Date, data=lnSRPsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnSRPsubs, aes(x=SRP.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y="SRP FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnSRPsubs[30,1], xend=lnSRPsubs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnSRPsubs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnSRPsubs[30,1], xend=lnSRPsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnSRPsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnSRPsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnSRPsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnSRPsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnSRPsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
SRP_FACplot <- FACplot
#############################################################
const <- "TP"
FUNlnQ <- lnTPsubs$lnQ
FUNlnC <- lnTPsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnTPsubs)
lnTPsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTPsubs[7] <- fitmse$fit$residuals

TPplot <- ggplot(data = lnTPsubs, aes(x=TP.Date, y=TP.concentrations.TP))
TPplot <- TPplot + geom_point() + labs(x="Year", y=expression(TP~(mg~L^{-1})))
# +annotate("text", x=lnTPsubs[30,1], y=750, size=9, label="A")
TPplot #I had to use a data point from the dataframe as the x for the label

FA_TPplot <- ggplot(data=lnTPsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_TPplot <- FA_TPplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(TP)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_TPplot

Trendline <- lm(V7~TP.Date, data=lnTPsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnTPsubs, aes(x=TP.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y="TP FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTPsubs[30,1], xend=lnTPsubs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnTPsubs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnTPsubs[30,1], xend=lnTPsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnTPsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnTPsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnTPsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnTPsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnTPsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
TP_FACplot <- FACplot
#############################################################
const <- "TSS"
FUNlnQ <- lnTSSsubs$lnQ
FUNlnC <- lnTSSsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnTSSsubs)
lnTSSsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTSSsubs[7] <- fitmse$fit$residuals

TSSplot <- ggplot(data = lnTSSsubs, aes(x=TSS.Date, y=TSS.concentrations.TSS))
TSSplot <- TSSplot + geom_point() + labs(x="Year", y=expression(TSS~(mg~L^{-1})))
# +annotate("text", x=lnTSSsubs[30,1], y=750, size=9, label="A")
TSSplot #I had to use a data point from the dataframe as the x for the label

FA_TSSplot <- ggplot(data=lnTSSsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_TSSplot <- FA_TSSplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(TSS)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_TSSplot

Trendline <- lm(V7~TSS.Date, data=lnTSSsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnTSSsubs, aes(x=TSS.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y="TSS FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTSSsubs[30,1], xend=lnTSSsubs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnTSSsubs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnTSSsubs[30,1], xend=lnTSSsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnTSSsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnTSSsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
TSS_FACplot <- FACplot
#############################################################
const <- "Cl"
FUNlnQ <- lnClsubs$lnQ
FUNlnC <- lnClsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnClsubs)
lnClsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnClsubs[7] <- fitmse$fit$residuals

Clplot <- ggplot(data = lnClsubs, aes(x=Cl.Date, y=Cl.concentrations.Cl))
Clplot <- Clplot + geom_point() + labs(x="Year", y=expression(Cl~(mg~L^{-1})))
# +annotate("text", x=lnClsubs[30,1], y=750, size=9, label="A")
Clplot #I had to use a data point from the dataframe as the x for the label

FA_Clplot <- ggplot(data=lnClsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_Clplot <- FA_Clplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(Cl)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_Clplot

Trendline <- lm(V7~Cl.Date, data=lnClsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnClsubs, aes(x=Cl.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y="Cl FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnClsubs[30,1], xend=lnClsubs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnClsubs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnClsubs[30,1], xend=lnClsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnClsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnClsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnClsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnClsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnClsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
Cl_FACplot <- FACplot
#############################################################
const <- "SO4"
FUNlnQ <- lnSO4subs$lnQ
FUNlnC <- lnSO4subs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnSO4subs)
lnSO4subs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnSO4subs[7] <- fitmse$fit$residuals

SO4plot <- ggplot(data = lnSO4subs, aes(x=SO4.Date, y=SO4.concentrations.SO4))
SO4plot <- SO4plot + geom_point() + labs(x="Year", y=expression(SO4~(mg~L^{-1})))
# +annotate("text", x=lnSO4subs[30,1], y=750, size=9, label="A")
SO4plot #I had to use a data point from the dataframe as the x for the label

FA_SO4plot <- ggplot(data=lnSO4subs, aes(x=lnQ, y=lnC)) + geom_point()
FA_SO4plot <- FA_SO4plot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(SO4)")+ annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
#+annotate("text", x=2.5, y=6, size=9, label="B") 
FA_SO4plot

Trendline <- lm(V7~SO4.Date, data=lnSO4subs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]
anova.trends <- anova(Trendline)
m <- summary.lm$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))

if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
FACplot <- ggplot(data=lnSO4subs, aes(x=SO4.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)))
FACplot <- FACplot + labs(x="Year", y=expression(paste(SO[4]," FACs")))
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnSO4subs[30,1], xend=lnSO4subs[50,1], y=-2.7, yend=-2.7, colour="blue")+
  annotate("text", x=lnSO4subs[60,1], y=-2.7, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot + annotate("segment", x=lnSO4subs[30,1], xend=lnSO4subs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnSO4subs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot + ggtitle(site) + stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnSO4subs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
# FACplot + geom_text(x=lnSO4subs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnSO4subs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) + scale_y_continuous(limits=c(-3.3,3.3))
#FACplot + annotate("text", x=lnSO4subs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot
SO4_FACplot <- FACplot

###############################################
#pdffilename <- paste(site,"_","FACplots",".pdf",sep="")
#pdf(file = pdffilename)
# pngfilename <- paste(site, "_", "FACplots%02d", ".png", sep="")
# png(file = pngfilename, pointsize = 40) #pointsize doesn't seem to change anything
# NO3_FACplot
# TN_FACplot
# SRP_FACplot
# TP_FACplot
# #TSS_FACplot
# #Cl_FACplot
# #SO4_FACplot
# 
# dev.off()

tifffilename <- paste(site, "_", "FACplots%02d", ".tiff", sep="")
tiff(file = tifffilename, res=144, width = 960, height = 960, compression = "none") 
NO3_FACplot #oh hell yeah son look at dat quality
TN_FACplot
SRP_FACplot
TP_FACplot
#TSS_FACplot
#Cl_FACplot
#SO4_FACplot

dev.off()

