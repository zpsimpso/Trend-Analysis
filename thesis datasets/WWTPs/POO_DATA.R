#THIS IS THE SCRIPT FOR LOOKING AT WWTP SHIT DISCHARGE
setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/WWTPs/")
noland <- read.csv("Noland_WWTP.csv")
#as.Date(noland$DATE, "%m-%d-%Y", tz) #FUCKING DAMN IT
noland$DATE <- as.POSIXct(as.Date(noland$DATE, "%m-%d-%Y", tz="US"))
oneday <- 24*3600
noland$DATE <- noland$DATE + oneday #SHIFTED BY ONE DAY, F U C K
noland_shitplot <- ggplot(noland, aes(x=DATE, y=EffTP))+geom_point()
noland_shitplot #INTERESTING
noland$loads <- (noland$EffQ)*(noland$EffTP)*(0.02832*3600*24*0.001)#UNIT CONVERSION BULLSHIT
noland_shitplot <- ggplot(noland, aes(x=DATE, y=log(loads)))+geom_point()
noland_shitplot #LOOKIN' GOOD, TIM


###################### NACA
#library(zoo)
NACA <- read.csv("NACA_WWTP.csv")
NACA$Date <- as.Date(as.yearmon(NACA$Date))
NACAplot <- ggplot(data=NACA, aes(x=Date, y=EffTP))+geom_point()
NACAplot
NACA$AvgTPload <- (NACA$EffTP*NACA$EffQ)*(0.02832*3600*24*0.001)*30
NACA_tpload <- ggplot(data=NACA, aes(x=Date, y=AvgTPload)) + geom_point()
NACA_tpload

############## Rogers
Rogers <- read.csv("Rogers_WWTP.csv")
Rogers$DATE <- as.POSIXct(Rogers$DATE)
Rogers_TP <- ggplot(data=Rogers, aes(x=DATE, y=EffTP))+geom_point()
Rogers_TP

#find starting month/year of data
start<-as.yearmon(paste(strftime(min(Rogers$DATE), "%Y"),"-", strftime(min(Rogers$DATE), "%m"), sep = ""))
end <- as.yearmon(paste(strftime(max(Rogers$DATE), "%Y"),"-", strftime(max(Rogers$DATE), "%m"), sep = ""))

monthvector <- seq.Date(from=as.Date(start), to=as.Date(end), by="month")
monthframe <- as.data.frame(matrix(data = NA, nrow = length(monthvector), ncol = 3))
monthframe[,1] <- monthvector
colnames(monthframe)[1:3]<-c("Date", "EffTP", "EffQ")

for (i in 1:length(monthvector)){
  a <- data.frame(Rogers[which(as.yearmon(strftime(Rogers$DATE, "%Y-%m")) == as.yearmon(monthvector[i])),])
  meanTP<-mean(a$EffTP, na.rm = TRUE) #there are some NA's
  meanQ <- mean(a$EffQ, na.rm = TRUE)
  monthframe[i,2:3]<-c(meanTP, meanQ)
}
start_year <- strftime(start, "%Y")
end_year <- strftime(end, "%Y")
start_Jan <- as.yearmon(paste(start_year,"-", "1", sep = ""))
end_Jan <- as.yearmon(paste("2017","-", "1", sep = "")) #to make it cover all later points
#plot
monthframe$Date <- as.yearmon(monthframe$Date)
ggplot(monthframe, aes(x=as.Date(Date), y=EffTP))+geom_point() +
  scale_x_date(breaks=seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="12 month"),labels=date_format("%b %Y"), minor_breaks = seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="6 month"))
write.csv(monthframe,"Rogers_late.csv")

#Merge with Thad's data
Rogers_full <- read.csv("Rogers_1996-2016.csv")
Rogers_full$Date <- as.POSIXct(Rogers_full$Date)

start<-as.yearmon(paste(strftime(min(Rogers_full$Date), "%Y"),"-", strftime(min(Rogers_full$Date), "%m"), sep = ""))
end <- as.yearmon(paste(strftime(max(Rogers_full$Date), "%Y"),"-", strftime(max(Rogers_full$Date), "%m"), sep = ""))
start_year <- strftime(start, "%Y")
end_year <- strftime(end, "%Y")
start_Jan <- as.yearmon(paste(start_year,"-", "1", sep = ""))
end_Jan <- as.yearmon(paste("2017","-", "1", sep = "")) #to make it cover all later points
#plot
Rogers_full$Date <- as.yearmon(Rogers_full$Date)
ggplot(Rogers_full, aes(x=as.Date(Date), y=TP_mgL))+geom_point() +
  scale_x_date(breaks=seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="12 month"),labels=date_format("%Y"), minor_breaks = seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="6 month"))



Rogers <- na.omit(data.frame(IRW_WWTP[which(IRW_WWTP$Site == "Rogers"), ]))
Rogers$Date <- as.Date(Rogers$Date)
Rogers$monthindex <- seq(1,241) #tree don't like date data i guess

#Rog_tree <- rpart("TP ~ Date", data=Rogers, method = "anova")
#I set minsize and mincut to 40 and 20; This seems appropriate since defaults were chasing all sorts of variance
Rog_tree <- tree("TP ~ monthindex", data=Rogers, control = tree.control(nobs = length(Rogers$TP), mincut = 20, minsize = 40))
plot(Rog_tree)
text(Rog_tree, cex=0.75)


Rog_tree_cv <- cv.tree(Rog_tree, best=4)
plot(Rog_tree_cv)
text(Rog_tree_cv, cex=0.75)

plot(Rogers$monthindex,Rogers$TP)
partition.tree(Rog_tree_cv, ordvars = "monthindex", add=TRUE, cex=0.3)
# 
# summary(Rog_tree_cv) #the numerator in the res mean dev expression is SSE
# Rog_tree_cv$frame$splits #splits with "" are the end nodes/leaves/whatever
# Rog_tree_cv$frame$dev #deviance of each split
# sum(Rog_tree_cv$frame$dev[3:7])-Rog_tree_cv$frame$dev[5] #use only the ones identified as the actual split, not dev of parent

#WOW ZACH OR YOU COULD USE THIS:
Tree_SSE <- deviance(Rog_tree_cv) #FUCK I'M DUMB
Total_SS <- sum((Rogers$TP - mean(Rogers$TP))^2)
Tree_R2 <- 1 - (Tree_SSE/Total_SS)
Tree_R2 #0.75 for this model


######
Springdale <- na.omit(data.frame(IRW_WWTP[which(IRW_WWTP$Site == "Springdale"), ]))
Springdale$Date <- as.Date(Springdale$Date)
Springdale$monthindex <- seq(1,194)

Sdale_tree <- tree("TP ~ monthindex", data=Springdale, control = tree.control(nobs = length(Springdale$TP), mincut = 20, minsize = 40))
plot(Sdale_tree)
text(Sdale_tree, cex=0.75)

Sdale_tree_cv <- cv.tree(Sdale_tree, best=3)
plot(Sdale_tree_cv)
text(Sdale_tree_cv, cex=0.75)

plot(Springdale$monthindex,Springdale$TP)
partition.tree(Sdale_tree_cv, ordvars = "monthindex", add=TRUE, cex=0.3)

Tree_SSE <- deviance(Sdale_tree_cv) 
Total_SS <- sum((Springdale$TP - mean(Springdale$TP))^2)
Tree_R2 <- 1 - (Tree_SSE/Total_SS)
Tree_R2 #0.85 for this model

ggplot(Springdale, aes(Date, TP))+geom_point() + scale_y_log10(limits=c(0.001, 15))
Springdale$logTP <- log10(Springdale$TP)


Sdale_tree <- tree("logTP ~ monthindex", data=Springdale, control = tree.control(nobs = length(Springdale$TP), mincut = 20, minsize = 40))
plot(Sdale_tree)
text(Sdale_tree, cex=0.75)

Sdale_tree_cv <- cv.tree(Sdale_tree, best=4)
plot(Sdale_tree_cv)
text(Sdale_tree_cv, cex=0.75)

a <- data.frame(Springdale[which(Springdale$monthindex >= 109), ])








#######
Sager <- na.omit(data.frame(IRW_WWTP[which(IRW_WWTP$Site == "Sager"), ]))
Sager$Date <- as.Date(Sager$Date)
Sager$monthindex <- seq(1,161)

Sager_tree <- tree("TP ~ monthindex", data=Sager, control = tree.control(nobs = length(Sager$TP), mincut = 20, minsize = 40))
plot(Sager_tree)
text(Sager_tree, cex=0.75)

Sager_tree_cv <- cv.tree(Sager_tree, best=3)
plot(Sager_tree_cv)
text(Sager_tree_cv, cex=0.75)

plot(Sager$monthindex,Sager$TP)
partition.tree(Sager_tree_cv, ordvars = "monthindex", add=TRUE, cex=0.3)

Tree_SSE <- deviance(Sager_tree_cv) 
Total_SS <- sum((Sager$TP - mean(Sager$TP))^2)
Tree_R2 <- 1 - (Tree_SSE/Total_SS)
Tree_R2 #0.71 for this model

#######
NACA <- na.omit(data.frame(IRW_WWTP[which(IRW_WWTP$Site == "NACA"), ]))
NACA$Date <- as.Date(NACA$Date)
NACA$monthindex <- seq(1,54)

plot(NACA$monthindex,NACA$TP)
#I don't think there's much reason to do regtree on NACA


############
#Organize WestSide Fay data
WestFay <- read.csv("WestSideFay_full.csv")
#WestFay$Date <- as.POSIXct(as.Date(WestFay$Date))
WestFay$Date <- as.POSIXct(as.Date(WestFay$Date, "%Y-%m-%d", tz="GMT"))
oneday <- 24*3600
WestFay$Date <- WestFay$Date + oneday

WestFay_TP <- ggplot(data=WestFay, aes(x=Date, y=TP))+geom_point()
WestFay_TP #weird plateau starting around 2015, i think its MDL for TP at 0.12 mg/L


#find starting month/year of data
start<-as.yearmon(paste(strftime(min(WestFay$Date), "%Y"),"-", strftime(min(WestFay$Date), "%m"), sep = ""))
end <- as.yearmon(paste(strftime(max(WestFay$Date), "%Y"),"-", strftime(max(WestFay$Date), "%m"), sep = ""))

monthvector <- seq.Date(from=as.Date(start), to=as.Date(end), by="month")
monthframe <- as.data.frame(matrix(data = NA, nrow = length(monthvector), ncol = 3))
monthframe[,1] <- monthvector
colnames(monthframe)[1:3]<-c("Date", "EffTP", "EffQ")

for (i in 1:length(monthvector)){
  a <- data.frame(WestFay[which(as.yearmon(strftime(WestFay$Date, "%Y-%m")) == as.yearmon(monthvector[i])),])
  meanTP<-mean(a$TP, na.rm = TRUE) #there are some NA's
  #meanQ <- mean(a$Q, na.rm = TRUE)
  monthframe[i,2]<-c(meanTP)
}

start_year <- strftime(start, "%Y")
end_year <- strftime(end, "%Y")
start_Jan <- as.yearmon(paste(start_year,"-", "1", sep = ""))
end_Jan <- as.yearmon(paste("2017","-", "1", sep = "")) #to make it cover all later points
#plot
monthframe$Date <- as.yearmon(monthframe$Date)
ggplot(monthframe, aes(x=as.Date(Date), y=EffTP))+geom_point() +
  scale_x_date(breaks=seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="12 month"),labels=date_format("%b %Y"), minor_breaks = seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="6 month"))
write.csv(monthframe,"WestSide_monthly.csv")

############
#look at Westside data with nondetects
#try to correct issue with high MDL by using NADA package
#library("NADA")
#Helsel is based

setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/WWTPs/")
West_nd <- read.csv("Westside_nondetects.csv")
West_nd$Date <- as.POSIXct(as.Date(West_nd$Date, "%Y-%m-%d", tz="GMT"))
oneday <- 24*3600
West_nd$Date <- West_nd$Date + oneday
West_nd <- na.omit(West_nd) #remove NA's

# West_ros <- ros(obs = West_nd$TP, censored = West_nd$Nondetect)
# West_ros_frame <- data.frame(West_ros$obs, West_ros$modeled, West_ros$pp, West_ros$censored)

#http://www.practicalstats.com/nada/downloads_files/NADAforR_Examples.pdf
#try cenreg

West_cenreg <- cenreg(Cen(West_nd$TP, West_nd$Nondetect)~West_nd$Date)

West_cenframe <- data.frame(West_cenreg@y, West_cenreg@survreg$y)
write.csv(West_cenframe, "West_cenframe.csv")

#find starting month/year of data
start<-as.yearmon(paste(strftime(min(West_nd$Date), "%Y"),"-", strftime(min(West_nd$Date), "%m"), sep = ""))
end <- as.yearmon(paste(strftime(max(West_nd$Date), "%Y"),"-", strftime(max(West_nd$Date), "%m"), sep = ""))

monthvector <- seq.Date(from=as.Date(start), to=as.Date(end), by="month")
monthframe <- as.data.frame(matrix(data = NA, nrow = length(monthvector), ncol = 3))
monthframe[,1] <- monthvector
colnames(monthframe)[1:3]<-c("Date", "TP", "TP_ros")

for (i in 1:length(monthvector)){
  a <- data.frame(West_nd[which(as.yearmon(strftime(West_nd$Date, "%Y-%m")) == as.yearmon(monthvector[i])),])
  ROSmean <- mean(ros(obs = a$TP, censored = a$Nondetect)) #get mean using NADA's ros fxn
  meanTP<-mean(a$TP, na.rm = TRUE) #there are some NA's
  #meanQ <- mean(a$Q, na.rm = TRUE)
  monthframe[i,2:3]<-c(meanTP, ROSmean)
}

#fucking statistics
#plot
West_nd_melt <- melt(monthframe, id.vars = "Date", variable.name = "method", value.name = "TP")

start_Jan <- as.yearmon(paste("2011","-", "1", sep = ""))

West_nd_monthplot <-ggplot(West_nd_melt, aes(x=as.Date(Date), y=TP, colour=method, shape=method))+geom_point(size=3) +
  scale_x_date(limits = c(as.Date(start_Jan), as.Date(end_Jan)), breaks=seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="12 month"),labels=date_format("%b %Y"), minor_breaks = seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="6 month"))
West_nd_monthplot <- West_nd_monthplot + scale_y_continuous(limits = c(0,0.5)) + theme(axis.title=element_text(size=20), axis.text=element_text(size=16), legend.text=element_text(size=18)) #changing y scale
West_nd_monthplot + labs(x="Year", y="Monthly avg TP")

West_ndTP <- data.frame(monthframe$Date, monthframe$TP_ros)
write.csv(West_ndTP, "West_ndTP.csv")

###########
#######################
#Get full datasets for IRW WWTP's
#Plot monthly TP conc over time

setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/WWTPs/full_record/")

IRW_WWTP <- read.csv("IRW_WWTP_full.csv")
IRW_WWTP$Date <- as.Date(IRW_WWTP$Date)
#IRW_WWTP$Date <- as.yearmon(as.POSIXct(IRW_WWTP$Date))
colnames(IRW_WWTP)[1] <- "Site"
#rename WestFay

#take out Sager
IRW_WWTP <- subset(IRW_WWTP, Site != "Sager")

#made it log scale since NACA is so low, and we're interested in the tenths/hundredths place
IRW_WWTPplot <- ggplot(data=IRW_WWTP, aes(x=Date, y=log(TP), colour=Site, shape=Site))+geom_point(size=4)
IRW_WWTPplot + labs(x="Year", y=expression(paste('log of TP (mg ',L^-1,')')))

IRW_WWTPplot <- ggplot(data=IRW_WWTP, aes(x=Date, y=TP, colour=Site, shape=Site))+geom_point(size=4)
IRW_WWTPplot <- IRW_WWTPplot + labs(x="Year", y=expression(paste('TP (mg ',L^-1,', log scale)')))+scale_y_log10(limits=c(0.01,15), breaks=10^(-2:1))
IRW_WWTPplot <- IRW_WWTPplot + scale_x_date(breaks = seq.Date(as.Date("1995-01-01"),as.Date("2016-01-01"),by="24 month"), labels=date_format("%Y"), minor_breaks = seq.Date(as.Date("1995-01-01"),as.Date("2017-01-01"),by="6 month"))
IRW_WWTPplot <- IRW_WWTPplot + theme(panel.background = element_rect(fill="white", colour="grey")) #+ geom_hline(yintercept = 0.037, linetype="dashed", colour="blue")
ok_std_text <- as.character(as.expression(0.037~mg~L^-1))
IRW_WWTPplot <- IRW_WWTPplot #+ annotate("text", x=IRW_WWTP[1,2], y=0.045, parse=TRUE, label=ok_std_text, size=8)
IRW_WWTPplot <- IRW_WWTPplot + theme(axis.title=element_text(size=18), axis.text=element_text(size=16), legend.text=element_text(size=16), legend.title=element_text(size=18))
IRW_WWTPplot + scale_shape_manual(values=c(19, 17, 3, 4)) + poster_theme
#looks alright, may need to add changepoints later like in Thad's paper

BigBoiz <- subset(IRW_WWTP, Site != c("Sager","NACA") & Date > as.Date("2009-01-01"))
Fay_pre12 <- subset(BigBoiz, Site == "Fayetteville" & Date < as.Date("2012-01-01"))
gm_mean(Fay_pre12$TP) #0.19 mg/L
Fay_post12 <- subset(BigBoiz, Site == "Fayetteville" & Date >= as.Date("2012-01-01"))
gm_mean(Fay_post12$TP) #0.11 mg/L
Sdale_pre12 <- subset(BigBoiz, Site == "Springdale" & Date < as.Date("2012-02-01"))
gm_mean(Sdale_pre12$TP) #0.32 mg/L
Sdale_post12 <- subset(BigBoiz, Site == "Springdale" & Date >= as.Date("2012-02-01"))
gm_mean(Sdale_post12$TP) #0.25 mg/L
Rog_pre12 <- subset(BigBoiz, Site == "Rogers" & Date < as.Date("2012-01-01"))
gm_mean(Rog_pre12$TP) #0.23 mg/L
Rog_post12 <- subset(BigBoiz, Site == "Rogers" & Date >= as.Date("2012-01-01"))
gm_mean(Rog_post12$TP) #0.18 mg/L

#
Sdale <- subset(BigBoiz, Site == "Springdale" & Date > as.Date("2009-01-01"))
Sdale$monthindex <- seq(1,86)

Sdale_tree <- tree("TP ~ monthindex", data=Sdale, control = tree.control(nobs = length(Sdale$TP), mincut = 10, minsize = 20))
plot(Sdale_tree)
text(Sdale_tree, cex=0.75)

Sdale_tree_cv <- cv.tree(Sdale_tree, best=3)
plot(Sdale_tree_cv)
text(Sdale_tree_cv, cex=0.75)

Sdale_ncpa <- ncpa(Sdale$monthindex, Sdale$TP)
Sdale_ncpa$cpobs 
#      cp       r2 mean left mean right pperm 5%  25%  50%  75% 95%
#[1,] 49.5 0.205099 0.3581633  0.2354054 0.001 10 10.5 48.5 49.5  50

Rog <- subset(BigBoiz, Site == "Rogers" & Date > as.Date("2009-01-01"))
Rog$monthindex <- seq(1,86)

Rog_tree <- tree("TP ~ monthindex", data=Rog, control = tree.control(nobs = length(Rog$TP), mincut = 12, minsize = 24))
plot(Rog_tree)
text(Rog_tree, cex=0.75)

Rog_tree_cv <- cv.tree(Rog_tree, best=3)
plot(Rog_tree_cv)
text(Rog_tree_cv, cex=0.75)










plot(Springdale$monthindex,Springdale$TP)
partition.tree(Sdale_tree_cv, ordvars = "monthindex", add=TRUE, cex=0.3)

Tree_SSE <- deviance(Sdale_tree_cv) 
Total_SS <- sum((Springdale$TP - mean(Springdale$TP))^2)
Tree_R2 <- 1 - (Tree_SSE/Total_SS)
Tree_R2 #0.85 for this model












#Springdale
Springdale_sub <- subset(IRW_WWTP, Site == "Springdale" & Date > as.Date("2009-01-01"))
ggplot(Springdale_sub, aes(x=Date, y=log(TP)))+geom_point()
TP_lm <- lm(log(TP)~as.POSIXct(Date), data=Springdale_sub)
summary(TP_lm)
m <- summary(TP_lm)$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))

#-7 percent per year at Springdale

#ROgers
Rogers_sub <- subset(IRW_WWTP, Site == "Rogers" & Date > as.Date("2009-01-01"))
ggplot(Rogers_sub, aes(x=Date, y=log(TP)))+geom_point()
TP_lm <- lm(log(TP)~as.POSIXct(Date), data=Rogers_sub)
summary(TP_lm)
m <- summary(TP_lm)$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))

#-6.3 percent per year at Rogers

#WestFay
WestFay_sub <- subset(IRW_WWTP, Site == "WestFay" & Date > as.Date("2009-01-01"))
ggplot(WestFay_sub, aes(x=Date, y=log10(TP)))+geom_point()
TP_lm <- lm(log(TP)~as.POSIXct(Date), data=WestFay_sub)
summary(TP_lm)
m <- summary(TP_lm)$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#-13.6% per year at West Fay

IRW_WWTP_sub <- subset(IRW_WWTP, Date > as.Date("2009-01-01"))
TP_lm <- lm(log(TP)~as.POSIXct(Date), data=IRW_WWTP_sub)
summary(TP_lm)
m <- summary(TP_lm)$coefficients[2]
percentslope <- (((exp(m)-1)*60*60*24*365*100))
#-11% per year general

#################
#Noland TP data, has some nondetects later in period



setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/WWTPs/full_record/")
Noland_TP <- read.csv("Noland_TP.csv")
Noland_TP$DATE <- as.POSIXct(as.Date(Noland_TP$DATE, "%Y-%m-%d", tz="GMT"))
oneday <- 24*3600
Noland_TP$DATE <- Noland_TP$DATE + oneday
colnames(Noland_TP)[3]<-"TP"

#find starting month/year of data
start<-as.yearmon(paste(strftime(min(Noland_TP$DATE), "%Y"),"-", strftime(min(Noland_TP$DATE), "%m"), sep = ""))
end <- as.yearmon(paste(strftime(max(Noland_TP$DATE), "%Y"),"-", strftime(max(Noland_TP$DATE), "%m"), sep = ""))

monthvector <- seq.Date(from=as.Date(start), to=as.Date(end), by="month")
monthframe <- as.data.frame(matrix(data = NA, nrow = length(monthvector), ncol = 3))
monthframe[,1] <- monthvector
colnames(monthframe)[1:3]<-c("Date", "TP", "TP_ros")

for (i in 1:length(monthvector)){
  a <- data.frame(Noland_TP[which(as.yearmon(strftime(Noland_TP$DATE, "%Y-%m")) == as.yearmon(monthvector[i])),])
  ROSmean <- mean(ros(obs = a$TP, censored = a$MDL)) #get mean using NADA's ros fxn
  meanTP<-mean(a$TP, na.rm = TRUE) #there are some NA's
  #meanQ <- mean(a$Q, na.rm = TRUE)
  monthframe[i,2:3]<-c(meanTP, ROSmean)
}

Noland_monthly <- monthframe

Noland_melt <- melt(Noland_monthly, id.vars = "Date", variable.name = "method", value.name = "TP")

start_Jan <- as.yearmon(paste("2009","-", "1", sep = ""))
end_Jan <- as.yearmon(paste("2016","-","1", sep=""))

Noland_nd_monthplot <-ggplot(Noland_monthly, aes(x=as.Date(Date), y=TP_ros))+geom_point(size=3) +
  scale_x_date(limits = c(as.Date(start_Jan), as.Date(end_Jan)), breaks=seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="12 month"),labels=date_format("%b %Y"), minor_breaks = seq.Date(as.Date(start_Jan),as.Date(end_Jan),by="6 month"))
Noland_nd_monthplot <- Noland_nd_monthplot + scale_y_continuous(limits = c(0,0.5)) + theme(axis.title=element_text(size=20), axis.text=element_text(size=16), legend.text=element_text(size=18)) #changing y scale
Noland_nd_monthplot + labs(x="Year", y=expression(paste('TP (mg ',L^-1,')')), title="Noland WWTP")+zachs_theme

write.csv(Noland_monthly, "Noland_monthly_ND.csv")






