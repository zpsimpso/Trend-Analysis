#WYman is BULLSHIT I tell ya, B U L L S H I T

###########
#Wyman

setwd("F:/TRENDS/UWRB TRENDS/UWRB/Wyman/")
siteflow <- read.csv(file = "flow.csv")
Wyman_flow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")




#can't separate bflow, too many gaps
#use 'cut' files for everything up to 12-17-2014
siteflow <- read.csv(file = "flow_cut.csv")
Wyman_flow <- read.csv(file = "flow_cut.csv")
concentrations <- read.csv(file = "concentrations_cut.csv")

lnNO3subs <- ReadInWQData2(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData2(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData2(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData2(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData2(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData2(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData2(siteflow, concentrations)[[7]]

Wyman_SRP <- lnSRPsubs
Wyman_TP <- lnTPsubs
Wyman_NO3 <- lnNO3subs
Wyman_TN  <- lnTNsubs

#separate bflow
Wyman_flow$bf <- Eckhardt(siteflow$Flow, filter_parameter = 0.96, passes=1)$bt
Wyman_flow$ro <- Eckhardt(siteflow$Flow, filter_parameter = 0.96, passes=1)$qft
Wyman_flow$BFF <- Wyman_flow$bf/Wyman_flow$Flow
Wyman_flow$Date <- as.POSIXct(Wyman_flow$Date)


Wyman_SRP <- merge(lnSRPsubs, Wyman_flow, by.x = "SRP.Date", by.y = "Date")
Wyman_TP <- merge(lnTPsubs, Wyman_flow, by.x = "TP.Date", by.y = "Date")
Wyman_NO3 <- merge(lnNO3subs, Wyman_flow, by.x = "NO3.Date", by.y = "Date")
Wyman_TN <- merge(lnTNsubs, Wyman_flow, by.x = "TN.Date", by.y = "Date")

ggplot(data=Wyman_SRP, aes(x=BFF, y=lnC)) + geom_point()
ggplot(data=Wyman_TP, aes(x=BFF, y=lnC)) + geom_point()

Wyman_SRP_storm_cut <- data.frame(Wyman_SRP[which(Wyman_SRP$BFF <= 0.6), ])
Wyman_SRP_base_cut <- data.frame(Wyman_SRP[which(Wyman_SRP$BFF > 0.6), ])
Wyman_TP_storm_cut <- data.frame(Wyman_TP[which(Wyman_TP$BFF <= 0.6), ])
Wyman_TP_base_cut <- data.frame(Wyman_TP[which(Wyman_TP$BFF > 0.6), ])
Wyman_NO3_storm_cut <- data.frame(Wyman_NO3[which(Wyman_NO3$BFF <= 0.6), ])
Wyman_NO3_base_cut <- data.frame(Wyman_NO3[which(Wyman_NO3$BFF > 0.6), ])
Wyman_TN_storm_cut <- data.frame(Wyman_TN[which(Wyman_TN$BFF <= 0.6), ])
Wyman_TN_base_cut <- data.frame(Wyman_TN[which(Wyman_TN$BFF > 0.6), ])
#delete last 3 columns (7-9)
Wyman_TN_storm_cut$bf<-NULL
Wyman_TN_storm_cut$ro<-NULL
Wyman_TN_storm_cut$BFF<-NULL
Wyman_TN_base_cut$bf<-NULL
Wyman_TN_base_cut$ro<-NULL
Wyman_TN_base_cut$BFF<-NULL



#read in stuff after 12-17-14
siteflow <- read.csv(file = "flow_after.csv")
Wyman_flow <- read.csv(file = "flow_after.csv")
concentrations <- read.csv(file = "concentrations_after.csv")

lnNO3subs <- ReadInWQData2(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData2(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData2(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData2(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData2(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData2(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData2(siteflow, concentrations)[[7]]

Wyman_flow$Date <- as.POSIXct(Wyman_flow$Date)
Wyman_SRP <- merge(lnSRPsubs, Wyman_flow, by.x = "SRP.Date", by.y = "Date")
Wyman_TP <- merge(lnTPsubs, Wyman_flow, by.x = "TP.Date", by.y = "Date")
Wyman_NO3 <- merge(lnNO3subs, Wyman_flow, by.x = "NO3.Date", by.y = "Date")
Wyman_TN <- merge(lnTNsubs, Wyman_flow, by.x = "TN.Date", by.y = "Date")
#make identity vector for storm vs base
storm_id <- concentrations$STORM

Wyman_SRP_storm_after <- data.frame(Wyman_SRP[which(storm_id == "YES"), ])
Wyman_SRP_base_after <- data.frame(Wyman_SRP[which(storm_id == "NO"), ])
Wyman_TP_storm_after <- data.frame(Wyman_TP[which(storm_id == "YES"), ])
Wyman_TP_base_after <- data.frame(Wyman_TP[which(storm_id == "NO"), ])
Wyman_NO3_storm_after <- data.frame(Wyman_NO3[which(storm_id == "YES"), ])
Wyman_NO3_base_after <- data.frame(Wyman_NO3[which(storm_id == "NO"), ])
Wyman_TN_storm_after <- data.frame(Wyman_TN[which(storm_id == "YES"), ])
Wyman_TN_base_after <- data.frame(Wyman_TN[which(storm_id == "NO"), ])

#now merge back together
Wyman_SRP_storm <- rbind(Wyman_SRP_storm_cut, Wyman_SRP_storm_after)
Wyman_SRP_base <- rbind(Wyman_SRP_base_cut, Wyman_SRP_base_after)
Wyman_TP_storm <- rbind(Wyman_TP_storm_cut, Wyman_TP_storm_after)
Wyman_TP_base <- rbind(Wyman_TP_base_cut, Wyman_TP_base_after)
Wyman_NO3_storm <- rbind(Wyman_NO3_storm_cut, Wyman_NO3_storm_after)
Wyman_NO3_base <- rbind(Wyman_NO3_base_cut, Wyman_NO3_base_after)
Wyman_TN_storm <- rbind(Wyman_TN_storm_cut, Wyman_TN_storm_after)
Wyman_TN_base <- rbind(Wyman_TN_base_cut, Wyman_TN_base_after)
#FUCKING FINALLY
write.csv(Wyman_SRP_storm, "Wyman_SRP_storm.csv") #for safe keeping
write.csv(Wyman_SRP_base, "Wyman_SRP_base.csv")
write.csv(Wyman_TP_storm, "Wyman_TP_storm.csv")
write.csv(Wyman_TP_base, "Wyman_TP_base.csv")
write.csv(Wyman_NO3_storm, "Wyman_NO3_storm.csv")
write.csv(Wyman_NO3_base, "Wyman_NO3_base.csv")
write.csv(Wyman_TN_storm, "Wyman_TN_storm.csv")
write.csv(Wyman_TN_base, "Wyman_TN_base.csv")



avgTP <- mean(Wyman_TP$TP.concentrations.TP) #0.092
avgSRP <- mean(Wyman_SRP$SRP.concentrations.SRP) #0.010 
gm_mean(Wyman_SRP$SRP.concentrations.SRP) #0.007
gm_mean(Wyman_TP$TP.concentrations.TP) #0.049


ggplot(Wyman_SRP, aes(x=SRP.Date, y=lnC))+geom_point()
Wyman_SRP_fit <- loess.wrapperMSE(Wyman_SRP$lnQ, Wyman_SRP$lnC, folds=10, iteration=10)
Wyman_SRP_fit$span #o.45
Wyman_SRP$FACs <- Wyman_SRP_fit$fit$residuals
ggplot(Wyman_SRP, aes(SRP.Date, FACs))+geom_point()
Kendall(Wyman_SRP$SRP.Date, Wyman_SRP$FACs) #0.15
SRP_lm <- lm(FACs~SRP.Date, data=Wyman_SRP)
summary(SRP_lm) #0.35

ggplot(Wyman_TP, aes(x=TP.Date, y=lnC))+geom_point()
Wyman_TP_fit <- loess.wrapperMSE(Wyman_TP$lnQ, Wyman_TP$lnC, folds=10, iteration=10)
Wyman_TP_fit$span #0.40
Wyman_TP$FACs <- Wyman_TP_fit$fit$residuals
ggplot(Wyman_TP, aes(TP.Date, FACs))+geom_point()
Kendall(Wyman_TP$TP.Date, Wyman_TP$FACs) #0.26
TP_lm <- lm(FACs~TP.Date, data=Wyman_TP)
summary(TP_lm) #0.12

ggplot(Wyman_NO3, aes(x=NO3.Date, y=lnC))+geom_point()
Wyman_NO3_fit <- loess.wrapperMSE(Wyman_NO3$lnQ, Wyman_NO3$lnC, folds=10, iteration=10)
Wyman_NO3_fit$span #0.7
Wyman_NO3$FACs <- Wyman_NO3_fit$fit$residuals
ggplot(Wyman_NO3, aes(NO3.Date, FACs))+geom_point()
Kendall(Wyman_NO3$NO3.Date, Wyman_NO3$FACs) #0.29
NO3_lm <- lm(FACs~NO3.Date, data=Wyman_NO3)
summary(NO3_lm) #0.80

ggplot(Wyman_TN, aes(x=TN.Date, y=lnC))+geom_point()
Wyman_TN_fit <- loess.wrapperMSE(Wyman_TN$lnQ, Wyman_TN$lnC, folds=10, iteration=10)
Wyman_TN_fit$span #0.7
Wyman_TN$FACs <- Wyman_TN_fit$fit$residuals
ggplot(Wyman_TN, aes(TN.Date, FACs))+geom_point()
Kendall(Wyman_TN$TN.Date, Wyman_TN$FACs) #0.003
TN_lm <- lm(FACs~TN.Date, data=Wyman_TN)
summary(TN_lm) #0.0023

Wyman_TNplot <- ggplot(Wyman_TN, aes(x=TN.Date, y=FACs))+geom_point(size=3)+zachs_theme + x_scale
Wyman_TNplot <- Wyman_TNplot + labs(x="Year", y="Flow-adjusted TN", title="All Streamflow Data at White River") + scale_y_continuous(limits=c(-2.2,2.2))
pval <- summary(TN_lm)$coefficients[2,4] #0.29 but Kendall's is 0.003
m <- summary(TN_lm)$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
perc_slope_val <- bquote(.(round(percentslope,2)))
perc_slope_char <- as.character(as.expression(substitute("%"~year^-1 ~ "=" ~ s, list(s = perc_slope_val))))
Wyman_TNplot <- Wyman_TNplot + stat_smooth(method=lm, se=FALSE, colour="darkturquoise") + annotate("text", x=as.POSIXct("2009-09-30"), y=-1.5, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(5.5))
Wyman_TNplot


