
############################################ FUNCTIONS
#from ggplot cookbook (Winston Chang)
#given a model, predict values of yvar from xvar
#this supports one predictor and one predicted variable
#xrange: if NULL, determine the x range from the model object. If a vector with
#two numbers, use those as the min and max of the prediction range.
#samples: number of samples across the x range.
#...: Further arguments to be passed to predict()
predictvals <- function(model, xvar, yvar, xrange=NULL, samples=100, ...){
  #if xrange isn't passed in, determine xrange from the models.
  #Different ways of extracting the x range, depending on model type
  if(is.null(xrange)){
    if(any(class(model) %in% c("lm", "glm")))
      xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess"))
      xrange <- range(model$x)
  }
  newdata <- data.frame(x=seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  newdata
}

# Multiple plot function
# source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#NEW FUNCTION - loess.wrapperMSE
#exactly like original, but evaluates fit based on MSE, not MAE
#can return the fit, span selected, and MSE of that span

loess.wrapperMSE <- function(x, y, span.vals = seq(0.1, 1, by = 0.05), folds, degree = 1, iteration){
  mse <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, degree = degree, span = span)
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals) {
    ii <- ii + 1
    #make the kfolds procedure equal across all spans
    set.seed(iteration)
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    mse[ii] <- mean((y[fltr] - y.cv[fltr])^2)
  }
  span <- span.vals[which.min(mse)]
  out <- loess(y ~ x, degree = degree, span = span)
  return(list(fit = out, span = span, MSE = (min(mse)), MSE.list = mse))
}

#perform k-fold CV (e.g., 10 x 10), fitting a LOESS model, returning MSE for each iteration and f
customfitMSE <- function(kfolds, iterations) {
  df <- data.frame(matrix(NA, nrow = length(seq(0.1, 1, by = 0.05)), ncol = iterations))
  rownames(df) <- paste(seq(0.1, 1, by = 0.05))
  colnames(df) <- paste("iter", seq(1, iterations))
  for (i in 1:iterations) {
    fitspecial.wrp <- loess.wrapperMSE(FUNlnQ,FUNlnC,folds=kfolds, iteration = i)
    itervector <- fitspecial.wrp$MSE.list
    df[,i]<-itervector
    print(i)
  }
  return(df)
}

########################################




#Set directory to the site's folder
#make sure there is 'flow.csv' and 'concentrations.csv'
#use 'Date', 'Flow', 'NO3', 'TN', 'SRP', 'TP', 'TSS', 'Cl', and 'SO4' as column headers

setwd("F:/TRENDS/IRW TRENDS/IRW/Spring/")
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

modbigloess <- loess(lnC ~ lnQ, lnNO3subs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnNO3subs, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggNO3 <- ggplot(data = lnNO3subs, aes(x=lnQ, y=lnC))
ggNO3 + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")



FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC


##################### TN

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

FUNlnQ <- lnTNsubs$lnQ
FUNlnC <- lnTNsubs$lnC

modbigloess <- loess(lnC ~ lnQ, lnTNsubs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnTNsubs, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggTN <- ggplot(data = lnTNsubs, aes(x=lnQ, y=lnC))
ggTN + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")


######################### SRP


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


FUNlnQ <- lnSRPsubs$lnQ
FUNlnC <- lnSRPsubs$lnC

modbigloess <- loess(lnC ~ lnQ, lnSRPsubs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnSRPsubs, span = 0.7)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggSRP <- ggplot(data = lnSRPsubs, aes(x=lnQ, y=lnC))
ggSRP + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")

###################### TP


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

FUNlnQ <- lnTPsubs$lnQ
FUNlnC <- lnTPsubs$lnC

modbigloess <- loess(lnC ~ lnQ, lnTPsubs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnTPsubs, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggTP <- ggplot(data = lnTPsubs, aes(x=lnQ, y=lnC))
ggTP + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")


####################################### TSS


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

FUNlnQ <- lnTSSsubs$lnQ
FUNlnC <- lnTSSsubs$lnC

fit<-loess(lnC~lnQ, span=0.5, data=lnTSSsubs)
lnTSSsubs[6]<-fit$residuals

fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTSSsubs[7] <- fitmse$fit$residuals

modbigloess <- loess(lnC ~ lnQ, lnTSSsubs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnTSSsubs, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggTSS <- ggplot(data = lnTSSsubs, aes(x=lnQ, y=lnC))
ggTSS + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")



TSSplot <- ggplot(data = lnTSSsubs, aes(x=TSS.Date, y=TSS.concentrations.TSS))
TSSplot <- TSSplot + geom_point() + labs(x="Year", y=expression(TSS~(mg~L^{-1}))) +
  annotate("text", x=lnTSSsubs[30,1], y=750, size=9, label="A")
TSSplot #I had to use a data point from the dataframe as the x for the label


FA_TSSplot <- ggplot(data=lnTSSsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_TSSplot <- FA_TSSplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(TSS)")+
  annotate("text", x=2.5, y=6, size=9, label="B") + annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
FA_TSSplot


Trendline <- lm(V7~TSS.Date, data=lnTSSsubs)
summary.lm<-summary(Trendline)
R2 <- summary(Trendline)$r.squared
pval<- summary.lm$coefficients[2,4]

#eqnR2 <- bquote(R^{2} == .(round(R2, 2)))
#eqnR2_text <- as.character(as.expression(eqnR2))

FACplot <- ggplot(data=lnTSSsubs, aes(x=TSS.Date, y=V7)) + geom_point()
FACplot <- FACplot + labs(x="Year", y="TSS FACs") + annotate("text", x=lnTSSsubs[30,1], y=2, size=9, label="C")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTSSsubs[30,1], xend=lnTSSsubs[50,1], y=-2.5, yend=-2.5, colour="blue")+
  annotate("text", x=lnTSSsubs[60,1], y=-2.5, label="Linear Regression", hjust = 0) 
# FACplot + geom_text(x=lnTSSsubs[50,1], y=-2.75, label=eqnR2_text, parse=TRUE)
FACplot <- FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3.1, label="p < 0.01", hjust=0)
#FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3.4, hjust=0, label=("R^2 = 0.03"), parse=TRUE)
FACplot + stat_smooth(se=FALSE, span=0.75, colour = "red")

tss_3step <- multiplot(TSSplot, FA_TSSplot, FACplot, cols=1)



############################ Cl


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

FUNlnQ <- lnClsubs$lnQ
FUNlnC <- lnClsubs$lnC

modbigloess <- loess(lnC ~ lnQ, lnClsubs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnClsubs, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggCl <- ggplot(data = lnClsubs, aes(x=lnQ, y=lnC))
ggCl + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")



####################### SO4

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


FUNlnQ <- lnSO4subs$lnQ
FUNlnC <- lnSO4subs$lnC

modbigloess <- loess(lnC ~ lnQ, lnSO4subs, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnSO4subs, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggSO4 <- ggplot(data = lnSO4subs, aes(x=lnQ, y=lnC))
ggSO4 + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")


#######################
#library(gridExtra)

#Look at histogram

setwd("F:/TRENDS/")

span_dat <- read.csv("Iteration_comparison.csv")

spans_plot <- ggplot(span_dat, aes(x=I_10))
spans_plot + geom_histogram(binwidth = 0.05, fill="white", colour="black", size = 1)

#create histograms of selected f for all sites, and then for each constituent
setwd("F:/TRENDS/")
dat<-read.csv("Trends_Comparison_all.csv")
hist_theme <- theme(axis.title=element_text(size=18),axis.text=element_text(size=16), panel.background=element_rect(fill="white",colour="black"),panel.grid.major=element_line(colour=alpha("grey",0.6)))

allspan <- ggplot(dat, aes(x=Span_MSE))
allspan <- allspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1)
allspan <- allspan + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1)))
allspan <- allspan + scale_y_continuous(limits=c(0,20), breaks=c(seq(0,20,5)))
allspan <- allspan + annotate("text", x=0.2, y=18, size=12, label="A")
allspan <- allspan + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75) + hist_theme

##just NO3 data
datNO3 <- subset(dat, Const == "NO3")
NO3span <- ggplot(datNO3, aes(x=Span_MSE))
NO3span <- NO3span + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
NO3span <- NO3span + annotate("text", x=0.2, y=7.2, size=12, label="B") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

#alternative: use stat_count to 'count' how many at each specific f value
#NO3span + stat_count() 

##just TN data
datTN <- subset(dat, Const == "TN")
TNspan <- ggplot(datTN, aes(x=Span_MSE))
TNspan <- TNspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
TNspan <- TNspan + annotate("text", x=0.2, y=7.2, size=12, label="C") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

##just SRP data
datSRP <- subset(dat, Const == "SRP")
SRPspan <- ggplot(datSRP, aes(x=Span_MSE))
SRPspan <- SRPspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
SRPspan <- SRPspan + annotate("text", x=0.2, y=7.2, size=12, label="D") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

##just TP data
datTP <- subset(dat, Const == "TP")
TPspan <- ggplot(datTP, aes(x=Span_MSE))
TPspan <- TPspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
TPspan <- TPspan + annotate("text", x=0.2, y=7.2, size=12, label="E") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

##just TSS data
datTSS <- subset(dat, Const == "TSS")
TSSspan <- ggplot(datTSS, aes(x=Span_MSE))
TSSspan <- TSSspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
TSSspan <- TSSspan + annotate("text", x=0.2, y=7.2, size=12, label="F") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

##just Cl data
datCl <- subset(dat, Const == "Cl")
Clspan <- ggplot(datCl, aes(x=Span_MSE))
Clspan <- Clspan + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
Clspan <- Clspan + annotate("text", x=0.2, y=7.2, size=12, label="G") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

##just SO4 data
datSO4 <- subset(dat, Const == "SO4")
SO4span <- ggplot(datSO4, aes(x=Span_MSE))
SO4span <- SO4span + geom_histogram(binwidth = 0.05, fill="white", colour="black", size=1) + labs(x=expression(italic(f)[opt]), y="Frequency") + scale_x_continuous(limits=c(0,1.1), breaks=c(seq(0,1,0.1))) + scale_y_continuous(limits=c(0,8), breaks=c(seq(0,8,2)))
SO4span <- SO4span + annotate("text", x=0.2, y=7.2, size=12, label="H") + geom_vline(xintercept = 0.5, linetype="dashed", size = 0.75)+ hist_theme

#library(gridExtra)
grid.arrange(grobs=list(allspan, NO3span, TNspan, SRPspan, TPspan, TSSspan, Clspan, SO4span), ncol=2, nrow=4)



#create plot of folds (10 v 35) comparison 
setwd("F:/TRENDS/")
dat<-read.csv("Folds_comparison_10v35.csv")
colnames(dat)[3]<-"Constituent"
##colour blind pallete
cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
labl <- list(expression(NO[3]), "TN", "SRP", "TP", "TSS", "Cl", expression(SO[4]))

foldsplot <- ggplot(dat, aes(x=Thirtyfive_10, y=Ten_10, shape=Site, colour=Constituent))
foldsplot <- foldsplot + geom_point(size=5) + scale_shape_manual(values=seq(0,17)) + scale_colour_manual(values=cb_palette, labels=labl, limits=c("NO3", "TN", "SRP", "TP", "TSS", "Cl", "SO4"))
foldsplot <- foldsplot + scale_x_continuous(limits=c(0,1), breaks=c(seq(0,1.0,0.1))) + scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1.0,0.1)))
foldsplot <- foldsplot + labs(x=expression(paste(italic(f)[opt],", 10 by 35")), y=expression(paste(italic(f)[opt], ", 10 by 10"))) + theme(axis.title=element_text(size=20), axis.text=element_text(size=16), panel.border=element_rect(fill=NA, colour="black"), panel.background=element_rect(fill="white",colour="black"), panel.grid.major=element_line(colour=alpha("grey",0.65)))
foldsplot <- foldsplot + geom_abline(slope=1, intercept=0, colour="gray")
foldsplot <- foldsplot+ theme(legend.text=element_text(size=18))
foldsplot

#create plot of iteration (10 v 1000) comparison
setwd("F:/TRENDS/")
dat<-read.csv("Iteration_comparison.csv")
colnames(dat)[3]<-"Constituent"
labl <- list(expression(NO[3]), "TN", "SRP", "TP", "TSS", "Cl", expression(SO[4]))
iterplot<-ggplot(dat, aes(x=I_1000, y=I_10, shape=Site, colour=Constituent))
iterplot <- iterplot + geom_point(size=5) + scale_shape_manual(values=seq(0,17)) + scale_colour_manual(values=cb_palette, labels=labl, limits=c("NO3", "TN", "SRP", "TP", "TSS", "Cl", "SO4"))
iterplot <- iterplot + scale_x_continuous(limits=c(0,1), breaks=c(seq(0,1.0,0.1))) + scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1.0,0.1)))
iterplot <- iterplot + labs(x=expression(paste(italic(f)[opt],", 1000 by 10")), y=expression(paste(italic(f)[opt], ", 10 by 10"))) + theme(axis.title=element_text(size=20), axis.text=element_text(size=16), panel.border=element_rect(fill=NA, colour="black"), panel.background=element_rect(fill="white",colour="black"), panel.grid.major=element_line(colour=alpha("grey",0.65)))
iterplot <- iterplot + geom_abline(slope=1, intercept=0, colour="gray")
iterplot<-iterplot + theme(legend.text=element_text(size=18))
iterplot

#compare trend magnitudes
setwd("F:/TRENDS/")
dat<-read.csv("Trends_Comparison_all.csv")
colnames(dat)[2]<-"Constituent"
dat$Constituent <- factor(dat$Constituent, levels=c("NO3", "TN", "SRP", "TP", "TSS", "Cl", "SO4"))
levels(dat$Constituent)[levels(dat$Constituent)=="NO3"] <- "NO[3]"
levels(dat$Constituent)[levels(dat$Constituent)=="SO4"] <- "SO[4]"
labl <- list(expression(NO[3]), "TN", "SRP", "TP", "TSS", "Cl", expression(SO[4]))
dat <- subset(dat, p_0.5 < 0.1) #only for trends with p < 0.10 using the standard f_0.5
tcplot<-ggplot(dat, aes(x=Trend_0.5, y=Trend_MSE, shape=Site, colour=Constituent))
tcplot <- tcplot + geom_point(size=6) + scale_shape_manual(values=seq(0,17)) + scale_colour_manual(values=cb_palette, labels = labl, limits=c("NO[3]", "TN", "SRP", "TP", "TSS", "Cl", "SO[4]"))
tcplot <- tcplot + geom_abline(slope=1, intercept=0, colour="gray") + labs(x=expression(paste("Trend using ", italic(f)[0.5]," (% change ", year^{-1}, ")")), y=expression(paste("Trend using  ", italic(f)[opt], " (% change ", year^{-1}, ")"))) + theme(axis.title=element_text(size=20), axis.text=element_text(size=16), panel.border=element_rect(fill=NA, colour="black"), panel.background=element_rect(fill="white",colour="black"), panel.grid.minor=element_line(colour=alpha("grey", alpha=0.5)), panel.grid.major=element_line(colour=alpha("grey",0.65)))
tcplot <- tcplot + theme(legend.text=element_text(size=18))
tcplot

#alternate tcplot, took out extreme point (SRP at Kings)
setwd("F:/TRENDS/")
dat<-read.csv("Trends_Comparison_all.csv")
colnames(dat)[2]<-"Constituent"
dat <- subset(dat, p_0.5 < 0.1 & Trend_0.5 > -20)
tcplot<-ggplot(dat, aes(x=Trend_0.5, y=Trend_MSE, shape=Site, colour=Constituent))
tcplot <- tcplot + geom_point(size=3) + scale_shape_manual(values=seq(0,17)) + scale_colour_manual(values=cb_palette, limits=c("NO3", "TN", "SRP", "TP", "TSS", "Cl", "SO4"))
tcplot <- tcplot + geom_abline(slope=1, intercept=0, colour="gray") + labs(x=expression(paste("Trend using ", f[0.5]," (% change ", year^{-1}, ")")), y=expression(paste("Trend using  ", f[opt], " (% change ", year^{-1}, ")"))) + theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))
tcplot

##recreate boxplot of k-fold CV selections of f
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WFWR/MSE/")
dat<-read.csv("THEWHOLEDAMNTHING1000_MSE_melt.csv")
  #NO3
NO3dat <- subset(dat, const == "NO3")
NO3dat$Fold <- as.factor(NO3dat$Fold)
NO3box <- ggplot(NO3dat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
NO3box <- NO3box + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
NO3box <- NO3box + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="NO[3]", size = 8, parse=TRUE)
NO3box <- NO3box + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#TN
TNdat <- subset(dat, const == "TN")
TNdat$Fold <- as.factor(TNdat$Fold)
TNbox <- ggplot(TNdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TNbox <- TNbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TNbox <- TNbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="TN", size = 8, parse=TRUE)
TNbox <- TNbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#SRP
SRPdat <- subset(dat, const == "SRP")
SRPdat$Fold <- as.factor(SRPdat$Fold)
SRPbox <- ggplot(SRPdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
SRPbox <- SRPbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
SRPbox <- SRPbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="SRP", size = 8, parse=TRUE)
SRPbox <- SRPbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#TP
TPdat <- subset(dat, const == "TP")
TPdat$Fold <- as.factor(TPdat$Fold)
TPbox <- ggplot(TPdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TPbox <- TPbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TPbox <- TPbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="TP", size = 8, parse=TRUE)
TPbox <- TPbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#TSS
TSSdat <- subset(dat, const == "TSS")
TSSdat$Fold <- as.factor(TSSdat$Fold)
TSSbox <- ggplot(TSSdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TSSbox <- TSSbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TSSbox <- TSSbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="TSS", size = 8, parse=TRUE)
TSSbox <- TSSbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#Cl
Cldat <- subset(dat, const == "Cl")
Cldat$Fold <- as.factor(Cldat$Fold)
Clbox <- ggplot(Cldat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
Clbox <- Clbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
Clbox <- Clbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="Cl", size = 8, parse=TRUE)
Clbox <- Clbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#SO4
SO4dat <- subset(dat, const == "SO4")
SO4dat$Fold <- as.factor(SO4dat$Fold)
SO4box <- ggplot(SO4dat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
SO4box <- SO4box + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
SO4box <- SO4box + labs(x="Folds in k-fold CV", y=expression(f[opt])) + annotate("text", x=35, y=0.9, label="SO[4]", size = 8, parse=TRUE)
SO4box <- SO4box + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#library(gridExtra)
grid.arrange(grobs=list(NO3box, TNbox, SRPbox, TPbox, TSSbox, Clbox, SO4box), ncol=1, nrow=7)
#need to work on dimensions and text positions, might also add 0.5 line


##recreate boxplot of k-fold CV selections of f (CHANGING UP LABELS)
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WFWR/MSE/")
dat<-read.csv("THEWHOLEDAMNTHING1000_MSE_melt.csv")
#NO3
NO3dat <- subset(dat, const == "NO3")
NO3dat$Fold <- as.factor(NO3dat$Fold)
NO3box <- ggplot(NO3dat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
NO3box <- NO3box + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
NO3box <- NO3box + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(NO[3]))
NO3box <- NO3box + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))
NO3box <- NO3box + theme(plot.title=element_text(size=20))

#TN
TNdat <- subset(dat, const == "TN")
TNdat$Fold <- as.factor(TNdat$Fold)
TNbox <- ggplot(TNdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TNbox <- TNbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TNbox <- TNbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(TN)) + theme(plot.title=element_text(size=20))
TNbox <- TNbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#SRP
SRPdat <- subset(dat, const == "SRP")
SRPdat$Fold <- as.factor(SRPdat$Fold)
SRPbox <- ggplot(SRPdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
SRPbox <- SRPbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
SRPbox <- SRPbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(SRP)) + theme(plot.title=element_text(size=20))
SRPbox <- SRPbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#TP
TPdat <- subset(dat, const == "TP")
TPdat$Fold <- as.factor(TPdat$Fold)
TPbox <- ggplot(TPdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TPbox <- TPbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TPbox <- TPbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(TP)) + theme(plot.title=element_text(size=20))
TPbox <- TPbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#TSS
TSSdat <- subset(dat, const == "TSS")
TSSdat$Fold <- as.factor(TSSdat$Fold)
TSSbox <- ggplot(TSSdat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
TSSbox <- TSSbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
TSSbox <- TSSbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(TSS)) + theme(plot.title=element_text(size=20))
TSSbox <- TSSbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#Cl
Cldat <- subset(dat, const == "Cl")
Cldat$Fold <- as.factor(Cldat$Fold)
Clbox <- ggplot(Cldat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
Clbox <- Clbox + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
Clbox <- Clbox + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(Cl)) + theme(plot.title=element_text(size=20))
Clbox <- Clbox + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#SO4
SO4dat <- subset(dat, const == "SO4")
SO4dat$Fold <- as.factor(SO4dat$Fold)
SO4box <- ggplot(SO4dat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE)
SO4box <- SO4box + scale_y_continuous(limit=c(0,1)) + scale_x_discrete(breaks=c(seq(5,45,5)))
SO4box <- SO4box + labs(x="Folds in k-fold CV", y=expression(f[opt])) + ggtitle(label=expression(SO[4])) + theme(plot.title=element_text(size=20))
SO4box <- SO4box + theme(axis.title.x=element_text(size=16)) + theme(axis.title.y=element_text(size=16))

#library(gridExtra)
grid.arrange(grobs=list(NO3box, TNbox, SRPbox, TPbox, TSSbox, Clbox, SO4box), ncol=1, nrow=7)
#need to work on dimensions and text positions, might also add 0.5 line

##recreate boxplot of k-fold CV selections of f (TRYING OUT FACET)
setwd("F:/TRENDS/UWRB TRENDS/UWRB/WFWR/MSE/")
dat<-read.csv("THEWHOLEDAMNTHING1000_MSE_melt.csv")
dat$Fold <- as.factor(dat$Fold)
#need to make constituents appear in a certain order
dat$const <- factor(dat$const, levels=c("NO3", "TN", "SRP", "TP", "TSS", "Cl", "SO4"))
levels(dat$const)[levels(dat$const)=="NO3"] <- "NO[3]"
levels(dat$const)[levels(dat$const)=="SO4"] <- "SO[4]"
allbox <- ggplot(dat, aes(x=Fold, y=span)) + geom_boxplot(notch=TRUE) + facet_grid(const ~ ., labeller = label_parsed)
allbox <- allbox + labs(x="Folds in K-fold CV", y=expression(italic(f)[opt])) + theme(strip.text=element_text(size=18),axis.title=element_text(size=20), panel.background=element_rect(fill="white",colour="grey"), axis.text=element_text(size=14))
allbox <- allbox + scale_x_discrete(breaks=c(seq(5,45,5)))

allbox



