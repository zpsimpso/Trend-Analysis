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


#Changing loess.wrapperMSE to use MAD instead of MSE
#can return the fit, span selected, and MAD of that span

loess.wrapperMAD <- function(x, y, span.vals = seq(0.1, 1, by = 0.05), folds, degree = 1, iteration){
  MAD <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, degree = degree, span = span)
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals) {
    ii <- ii + 1
    #make the kfolds procedure equal across all spans
    set.seed(iteration)
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    MAD[ii] <- mean(abs(y[fltr] - y.cv[fltr]))
  }
  span <- span.vals[which.min(MAD)]
  out <- loess(y ~ x, degree = degree, span = span)
  return(list(fit = out, span = span, MAD = (min(MAD)), MAD.list = MAD))
}

##colour blind pallete
cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

poster_theme <- theme(axis.title=element_text(size=20), axis.text=element_text(size=18), legend.text=element_text(size=18), 
                     legend.title=element_text(size=20), 
                     panel.background = element_rect(fill="white", colour="grey"), 
                     panel.grid.major = element_line(colour = "grey", linetype = "dashed"),
                     plot.title = element_text(size=24))

#############

setwd("F:/TRENDS/UWRB TRENDS/UWRB/WEC/")
siteflow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")


lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

ggplot(lnNO3subs, aes(lnQ, lnC))+geom_point()
ggplot(lnTNsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSRPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTSSsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnClsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSO4subs, aes(lnQ, lnC))+geom_point()

#testing some special features of loess
test <- loess(lnC~lnQ, lnNO3subs, degree = 1, span=0.5)
test2 <- loess(lnC~lnQ, lnNO3subs, degree = 1, span=0.5, cell=0.05)
test2 <- loess(lnC~lnQ, lnNO3subs, degree = 1, span=0.5, surface="direct")

test_pred <- predictvals(test, "lnQ", "lnC")
test2_pred <- predictvals(test2, "lnQ", "lnC")
ggplot(data = lnNO3subs, aes(lnQ, lnC))+geom_point()+geom_line(data=test_pred, colour = "red")+geom_line(data=test2_pred, colour="green")
#

modbigloess <- loess(lnC ~ lnQ, lnNO3subs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnNO3subs, degree = 1, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggNO3 <- ggplot(data = lnNO3subs, aes(x=lnQ, y=lnC))
ggNO3 + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")

MADfit <- loess.wrapperMAD(lnNO3subs$lnQ, lnNO3subs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnNO3subs$lnQ, lnNO3subs$lnC, folds = 10, iteration = 10)

setwd("F:/TRENDS/IRW TRENDS/IRW/Watts/")
siteflow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")


lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

ggplot(lnNO3subs, aes(lnQ, lnC))+geom_point()
ggplot(lnTNsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSRPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTSSsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnClsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSO4subs, aes(lnQ, lnC))+geom_point()


MADfit <- loess.wrapperMAD(lnNO3subs$lnQ, lnNO3subs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnNO3subs$lnQ, lnNO3subs$lnC, folds = 10, iteration = 10)

modbigloess <- loess(lnC ~ lnQ, lnNO3subs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "lnQ", "lnC")

modsmallloess <- loess(lnC ~ lnQ, lnNO3subs, degree = 1, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "lnQ", "lnC")

ggNO3 <- ggplot(data = lnNO3subs, aes(x=lnQ, y=lnC))
ggNO3 <- ggNO3 + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour="green") + geom_line(data=loess_predicted_big, colour="red") + geom_line(data=loess_predicted_small, colour="blue")
ggNO3 + zachs_theme

lnNO3subs$logQ <- log10(lnNO3subs$NO3.Flow)
lnNO3subs$logC <- log10(lnNO3subs$NO3.concentrations.NO3)
MADfit <- loess.wrapperMAD(lnNO3subs$logQ, lnNO3subs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnNO3subs, degree = 1, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnNO3subs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggNO3 <- ggplot(data = lnNO3subs, aes(x=logQ, y=logC))
ggNO3 <- ggNO3 + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour=cb_palette[4], size = 1) + geom_line(data=loess_predicted_small, colour=cb_palette[2], size = 1)+ geom_line(data=loess_predicted_big, colour=cb_palette[3], size = 1)
ggNO3 <- ggNO3 + poster_theme + labs(title=expression(NO[3]~at~Illinois~River), x="log of Q (cfs)", y = expression(paste("log of ", NO[3], " (mg ", L^{-1}, ")")))
ggNO3 <- ggNO3 + annotate("text", x = 4, y = 0.6, hjust=0, colour= cb_palette[4], label="f = 0.5", size = 8) + annotate("text", x = 4, y = 0.7, hjust = 0, colour= cb_palette[3], label="f = 0.8", size = 8) + annotate("text", x = 4, y = 0.5, hjust = 0, colour= cb_palette[2], label=as.character(expression(paste(f[opt]," = 0.3"))), parse=TRUE, size = 8)
ggNO3

##

setwd("F:/TRENDS/IRW TRENDS/IRW/Osage/")
siteflow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")


lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

ggplot(lnNO3subs, aes(lnQ, lnC))+geom_point()
ggplot(lnTNsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSRPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTSSsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnClsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSO4subs, aes(lnQ, lnC))+geom_point()


MADfit <- loess.wrapperMAD(lnSRPsubs$lnQ, lnSRPsubs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnSRPsubs$lnQ, lnSRPsubs$lnC, folds = 10, iteration = 10)

lnSRPsubs$logQ <- log10(lnSRPsubs$SRP.Flow)
lnSRPsubs$logC <- log10(lnSRPsubs$SRP.concentrations.SRP)
MADfit <- loess.wrapperMAD(lnSRPsubs$logQ, lnSRPsubs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnSRPsubs, degree = 1, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnSRPsubs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggSRP <- ggplot(data = lnSRPsubs, aes(x=logQ, y=logC))
ggSRP <- ggSRP + geom_point() + geom_smooth(method="loess", span=0.65, se=FALSE, colour=cb_palette[4], size =1) + geom_line(data=loess_predicted_small, colour=cb_palette[2], size =1)+ geom_line(data=loess_predicted_big, colour=cb_palette[3], size =1)
ggSRP <- ggSRP + poster_theme + labs(title=expression(SRP~at~Osage~Creek), x="log of Q (cfs)", y = expression(paste("log of ", "SRP", " (mg ", L^{-1}, ")")))
ggSRP <- ggSRP  + annotate("text", x = 2.25, y = -0.3, colour= cb_palette[3], label="f = 0.8", size = 8, hjust=0) + annotate("text", x = 2.25, y = -0.4, colour= cb_palette[2], label="f = 0.3", size = 8, hjust=0)
ggSRP


ggSRP + annotate("text", x = 2.25, y = -0.5, colour= cb_palette[4], label=as.character(expression(paste(f[opt]," = 0.65"))), parse=TRUE, size = 8, hjust=0)



###

setwd("F:/TRENDS/UWRB TRENDS/UWRB/RC45/truncated/")
siteflow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")


lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
# lnNO3subs[which(lnNO3subs$lnC < -4.5), 4] <- NA
# lnNO3subs <- na.omit(lnNO3subs)
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]
lnClsubs[which(lnClsubs$lnC > 3), 4] <- NA
lnClsubs <- na.omit(lnClsubs)
lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

ggplot(lnNO3subs, aes(lnQ, lnC))+geom_point()
ggplot(lnTNsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSRPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTSSsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnClsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSO4subs, aes(lnQ, lnC))+geom_point()

MADfit <- loess.wrapperMAD(lnTNsubs$lnQ, lnTNsubs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnTNsubs$lnQ, lnTNsubs$lnC, folds = 10, iteration = 10)

lnTNsubs$logQ <- log10(lnTNsubs$TN.Flow)
lnTNsubs$logC <- log10(lnTNsubs$TN.concentrations.TN)
MADfit <- loess.wrapperMAD(lnTNsubs$logQ, lnTNsubs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnTNsubs, degree = 1, span = 0.25)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnTNsubs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggTN <- ggplot(data = lnTNsubs, aes(x=logQ, y=logC))
ggTN <- ggTN + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour=cb_palette[4], size = 1) + geom_line(data=loess_predicted_small, colour=cb_palette[2], size =1)+ geom_line(data=loess_predicted_big, colour=cb_palette[3], size = 1)
ggTN <- ggTN + poster_theme + labs(title=expression(TN~at~Richland~Creek), x="log of Q (cfs)", y = expression(paste("log of ", "TN", " (mg ", L^{-1}, ")")))
ggTN <- ggTN  + annotate("text", x = 3, y = -0.4, colour= cb_palette[3], label="f = 0.8", size = 8, hjust=0) + annotate("text", x = 3, y = -0.6, colour= cb_palette[2], label=as.character(expression(paste(f[opt]," = 0.25"))), parse=TRUE, size = 8, hjust=0)
ggTN


ggTN + annotate("text", x = 3, y = -0.5, colour= cb_palette[4], label="f = 0.5", size = 8, hjust=0)
#

MADfit <- loess.wrapperMAD(lnClsubs$lnQ, lnClsubs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnClsubs$lnQ, lnClsubs$lnC, folds = 10, iteration = 10)

lnClsubs$logQ <- log10(lnClsubs$Cl.Flow)
lnClsubs$logC <- log10(lnClsubs$Cl.concentrations.Cl)
MADfit <- loess.wrapperMAD(lnClsubs$logQ, lnClsubs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnClsubs, degree = 1, span = 0.3)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnClsubs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggCl <- ggplot(data = lnClsubs, aes(x=logQ, y=logC))
ggCl <- ggCl + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour=cb_palette[4]) + geom_line(data=loess_predicted_small, colour=cb_palette[2])+ geom_line(data=loess_predicted_big, colour=cb_palette[3])
ggCl <- ggCl + poster_theme + labs(title=expression(Cl~at~Osage~Creek), x="log of Q (cfs)", y = expression(paste("log of ", "Cl", " (mg ", L^{-1}, ")")))
ggCl <- ggCl  + annotate("text", x = 2.25, y = -0.3, colour= cb_palette[3], label="f = 0.8", size = 8, hjust=0) + annotate("text", x = 2.25, y = -0.4, colour= cb_palette[2], label=as.character(expression(paste(f[opt]," = 0.25"))), parse=TRUE, size = 8, hjust=0)
ggCl


# ggCl + annotate("text", x = 2.25, y = -0.5, colour= cb_palette[4], label="f = 0.5"))), parse=TRUE, size = 8, hjust=0)



####

setwd("F:/TRENDS/UWRB TRENDS/UWRB/Kings/2009-2015, no 2010 data/")
siteflow <- read.csv(file = "flow.csv")
concentrations <- read.csv(file = "concentrations.csv")


lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
# lnNO3subs[which(lnNO3subs$lnC < -4.5), 4] <- NA
# lnNO3subs <- na.omit(lnNO3subs)
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]

lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

ggplot(lnNO3subs, aes(lnQ, lnC))+geom_point()
ggplot(lnTNsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSRPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTPsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnTSSsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnClsubs, aes(lnQ, lnC))+geom_point()
ggplot(lnSO4subs, aes(lnQ, lnC))+geom_point()

MADfit <- loess.wrapperMAD(lnTNsubs$lnQ, lnTNsubs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnTNsubs$lnQ, lnTNsubs$lnC, folds = 10, iteration = 10)

lnTNsubs$logQ <- log10(lnTNsubs$TN.Flow)
lnTNsubs$logC <- log10(lnTNsubs$TN.concentrations.TN)
MADfit <- loess.wrapperMAD(lnTNsubs$logQ, lnTNsubs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnTNsubs, degree = 1, span = 0.2)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnTNsubs, degree = 1, span = 0.8)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggTN <- ggplot(data = lnTNsubs, aes(x=logQ, y=logC))
ggTN <- ggTN + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour=cb_palette[4], size = 1) + geom_line(data=loess_predicted_small, colour=cb_palette[2], size =1)+ geom_line(data=loess_predicted_big, colour=cb_palette[3], size = 1)
ggTN <- ggTN + poster_theme + labs(title=expression(TN~at~Kings~River), x="log of Q (cfs)", y = expression(paste("log of ", "TN", " (mg ", L^{-1}, ")")))
ggTN <- ggTN  + annotate("text", x = 3, y = -0.4, colour= cb_palette[3], label="f = 0.8", size = 8, hjust=0) + annotate("text", x = 3, y = -0.6, colour= cb_palette[2], label=as.character(expression(paste(f[opt]," = 0.20"))), parse=TRUE, size = 8, hjust=0)
ggTN


ggTN + annotate("text", x = 3, y = -0.5, colour= cb_palette[4], label="f = 0.5", size = 8, hjust=0)
#


MADfit <- loess.wrapperMAD(lnSRPsubs$lnQ, lnSRPsubs$lnC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnSRPsubs$lnQ, lnSRPsubs$lnC, folds = 10, iteration = 10)

lnSRPsubs$logQ <- log10(lnSRPsubs$SRP.Flow)
lnSRPsubs$logC <- log10(lnSRPsubs$SRP.concentrations.SRP)
MADfit <- loess.wrapperMAD(lnSRPsubs$logQ, lnSRPsubs$logC, folds = 10, iteration = 10)
MSEfit <- loess.wrapperMSE(lnSRPsubs$logQ, lnSRPsubs$logC, folds = 10, iteration = 10)

modsmallloess <- loess(logC ~ logQ, lnSRPsubs, degree = 1, span = 0.25)
loess_predicted_small <- predictvals(modsmallloess, "logQ", "logC")
modbigloess <- loess(logC ~ logQ, lnSRPsubs, degree = 1, span = 0.95)
loess_predicted_big<-predictvals(modbigloess, "logQ", "logC")

ggSRP <- ggplot(data = lnSRPsubs, aes(x=logQ, y=logC))
ggSRP <- ggSRP + geom_point() + geom_smooth(method="loess", span=0.5, se=FALSE, colour=cb_palette[4], size = 1) + geom_line(data=loess_predicted_small, colour=cb_palette[2], size =1)+ geom_line(data=loess_predicted_big, colour=cb_palette[3], size = 1)
ggSRP <- ggSRP + poster_theme + labs(title=expression(SRP~at~Kings~River), x="log of Q (cfs)", y = expression(paste("log of ", "SRP", " (mg ", L^{-1}, ")")))
ggSRP <- ggSRP  + annotate("text", x = 3, y = -0.4, colour= cb_palette[3], label=as.character(expression(paste(f[opt]," = 0.95"))), parse=TRUE, size = 8, hjust=0) + annotate("text", x = 3, y = -0.6, colour= cb_palette[2], label=as.character(expression(paste(f," = 0.20"))), parse=TRUE, size = 8, hjust=0)
ggSRP


ggSRP + annotate("text", x = 3, y = -0.5, colour= cb_palette[4], label="f = 0.5", size = 8, hjust=0)
#


























