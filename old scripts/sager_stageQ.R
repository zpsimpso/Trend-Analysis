setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/IRW/Sager/Sager_discharge_work/")
getwd()
#we are going to create a 3 pronged ifelse function
#for low stage,0 to 0.689 ft, use Q=1.3904(stage)
#for stage of .69 ft to 3.915 ft, use LOESS fit on the calib data and predict Q
#for stage>3.915 ft, use Q=450.5(stage)-1505.3

#first, figure out LOESS fit for mid range stage

#read data
mid_calib <- read.csv("MID_stageQ_calib.csv")
Q <- mid_calib$Q
stage <- mid_calib$Adj_stage

#find upper and lower stage bounds in data
stage.lims <- range(stage)
#sort from low to high stage
stage.grid <- seq(from=stage.lims[1],to=stage.lims[2])
loessfit <- loess(Q~stage,span=0.5)
plot(stage,Q)
lines(stage.grid, predict(loessfit,data.frame(stage=stage.grid)), col='red', lwd=2)
#now load the full stage file
Sager_full <- read.csv("est_fullstage.csv")
Sager_full$DateTime <- as.POSIXct(Sager_full$DateTime, tz="", format="%Y-%m-%d %H:%M")

fullstage <- Sager_full$adj_stage
#time to apply 3-prong if else statement
#edit: now it's more like 5 prongs...
Qpredict <- ifelse(fullstage<=0.573, 0.01,
                   (ifelse(fullstage<=0.689,(fullstage*7.4082)-4.204,
                           (ifelse(fullstage>=0.69 & fullstage<=3.915, 
                                   predict(loessfit,newdata=fullstage), (450.5*fullstage)-1505.3)))))
#Jesus Christ

#make output file
Sager_fullest <- data.frame(Sager_full$DateTime, fullstage, Qpredict)
write.csv(Sager_fullest,"Sager_fullest.csv")

#take a gander at your creation
plot(fullstage, Qpredict)


Sager_fullest$Sager_full.DateTime<-as.POSIXct(Sager_fullest$Sager_full.DateTime)
#rename this stupid column
colnames(Sager_fullest)[1]<-"DateTime"

#to use later:
start <- as.POSIXct("2011-07-26")
end <- as.POSIXct("2015-06-30")
#by using 'day' as the interval, rather than a specific amount of seconds, we
#(hopefully) avoid the problem of DST throwing things off by an hour in spring/fall
DailyDateTime<- seq(from=start, by="day", to=end)

#make sure that the missing dates in the stage record show up with 
#NA's rather than just be omitted

#create datetime sequence of 15 min increment
stagestart <- as.POSIXct("2011-07-26 00:00")
stageinterval <- 900 #15 mins in seconds
stageend <- as.POSIXct("2015-06-30 23:45")
DateTime <- seq(from=stagestart, by=stageinterval, to=stageend)

date.frm <- data.frame(DateTime)
#merge
m1 <- merge(date.frm, Sager_fullest, by="DateTime", all=TRUE, sort=TRUE)


#install package "xts"
#library("xts")

xtsobject <- xts(m1[,2:3],order.by=as.POSIXct(m1[,1], format="%m/%d/%Y %H:%M"))



dailymeans <- apply.daily(xtsobject, colMeans, na.rm=TRUE)
#stupid colmeans returns the column sums at the bottom
#so we need to delete that last row
dailymeans_forreal <- dailymeans[-1437,]

Sager_daily <- data.frame(DailyDateTime,dailymeans_forreal)
#WHA WHA WHA WHA WHA WHA WHAMMY

#NOTE: There are some missing days in 2013 and 2014,
#I'll interpolate for these manually

write.csv(Sager_daily, "Sager_daily.csv")

################
setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/IRW/Sager/Sager_discharge_work/")


Sager_stageQ <- read.csv("stagedischarge_data.csv")


stage <- Sager_stageQ$Sager_Stage
Q <- Sager_stageQ$Sager_Q



Sager_ggplot <- ggplot(data.frame(stage, Q), aes(x=stage, y=Q))
Sager_ggplot + geom_point() +
  xlab("Stage\n(ft)") + ylab("Discharge\n(cfs)") +
  scale_x_continuous(breaks=seq(from=0, to=10, by=1))+
  scale_y_continuous(breaks=seq(from=0, to=1700, by=100))+
  ggtitle("Stage-Discharge at Sager Creek")+
  theme(text = element_text(size = 12, family = "Calibri"))

ggsave("sagerstageQ.png", width=4, height=6, dpi=300)






