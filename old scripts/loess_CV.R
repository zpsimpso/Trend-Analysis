

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



for (l in seq(0.1,0.9,by=0.05)){
  y<-FUNlnC
  x<-FUNlnQ
  loessfit<-loess(y~x, span=l)
    for (i in 1:1000){
      cvFit(loessfit, y=y, cost=mspe, K=10, R=1, seed=i)
    }
    
}

loessfit <- function(f, y, x){
  loess(y~x, span=f)
}


y<-FUNlnC
x<-FUNlnQ
dat<-data.frame(y,x)
fit<-predict(loess(y~x, data=dat))

diddly<-cvFit(formula=(loess(y~x, span=0.2)), x=FUNlnQ, y=FUNlnC, cost=mspe, foldType="random", K=10, R=1, seed=2)
diddly<-cvFit(fit, data=dat, y=dat[['y']], x=dat[['x']], cost=mspe, foldType="random", K=10, R=1, seed=2)
diddly$cv


fuckingwork<-cvFit(lm, y~x, data=dat, cost=mspe, K=10, R=1, seed=2)
fuckingwork$cv

diddly2<-cvFit(fit, x=x, y=y, cost=mspe, foldType="random", K=2, R=1)

dat <- data.frame(FUNlnC, FUNlnQ)

widdly <- cvFit(loessfit(0.3, FUNlnC, FUNlnQ), data=dat, y = dat$FUNlnC, cost=mspe, K = 2, R = 1, foldType = "random")

number <- nrow(dat)
FOLDS <- cvFolds(n <- number, K=3, R=1, type = "random")
newdat <- data.frame(FOLDS$which, FUNlnC, FUNlnQ)
View(newdat)




loess.wrapperMSE <- function(x, y, span.vals = seq(0.1, 1, by = 0.05), folds){
  mse <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, span = span)
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals) {
    ii <- ii + 1
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    mse[ii] <- mean((y[fltr] - y.cv[fltr])^2)
  }
  span <- span.vals[which.min(mse)]
  out <- loess(y ~ x, span = span)
  return(list(fit = out, span = span, MSE = (min(mse)), MSE.list = mse))
}

set.seed(1235)
shiz <- loess.wrapperMSE(lnNO3subs$lnQ, lnNO3subs$lnC, folds = 5)
shiz$span
shiz$MSE
shiz$MSE.list
#I FUCKING MADE SOMETHING WORK FOR ONCE, QUICK, R INSTALL 'CELEBRATE'
#install.packages('celebrate')
#fuckyes<-touchdowndance(data=zach, intensity = 11, duration = 30, style = 'random')

###################################################



setwd("F:/TRENDS/IRW TRENDS/IRW/Watts/")

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

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC

WattsNO3 <- ggplot(lnNO3subs, aes(x=lnQ, y=lnC))
WattsNO3 + geom_point() + stat_smooth(method = "loess", span = 0.25) + stat_smooth(method="loess", span = 0.8)
























