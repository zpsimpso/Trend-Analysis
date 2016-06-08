
############################################ FUNCTIONS


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
    #make the kfolds procedure equal across all spans
    set.seed(iteration)
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    mse[ii] <- mean((y[fltr] - y.cv[fltr])^2)
  }
  span <- span.vals[which.min(mse)]
  out <- loess(y ~ x, span = span)
  return(list(fit = out, span = span, MSE = (min(mse)), MSE.list = mse))
}


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

setwd("F:/TRENDS/UWRB TRENDS/UWRB/Kings/2009-2015, no 2010 data/")
getwd()


#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)
siteflow$Date<- as.POSIXct(siteflow$Date)
concentrations$Date <- as.POSIXct(concentrations$Date)

lnNO3subs <- ReadInWQData3(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData3(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData3(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData3(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData3(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData3(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData3(siteflow, concentrations)[[7]]

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC

NO3_10folds <- customfitMSE(kfolds=10, iterations=10)
NO3_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(NO3_10folds)[11]<-"SPAN"
NO3_melt <- melt(NO3_10folds, id.vars="SPAN", value.name="MSE")


#NO3_one<-subset(NO3_melt, variable=="iter 1", select= c (SPAN, MSE))
#vec <- c(1:length(NO3_melt[,1]))
#NO3_melt$stdPE <- apply(NO3_melt, 1, function(row) (sqrt(row[3]))/(mean(FUNlnC)))

### make spaghetti plot in ggplot2

gg_NO3 <- ggplot(data=NO3_melt, aes(x=SPAN, y=MSE, group=variable))
gg_NO3 <- gg_NO3 + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

NO3_melt$Const <- "NO3"

##################### TN


FUNlnQ <- lnTNsubs$lnQ
FUNlnC <- lnTNsubs$lnC

TN_10folds <- customfitMSE(kfolds=10, iterations=10)
TN_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(TN_10folds)[11]<-"SPAN"
TN_melt <- melt(TN_10folds, id.vars="SPAN", value.name="MSE")
TN_melt$Const <- "TN"

TN_one<-subset(TN_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_TN <- ggplot(data=TN_melt, aes(x=SPAN, y=MSE, group=variable))
gg_TN <- gg_TN + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

######################### SRP


FUNlnQ <- lnSRPsubs$lnQ
FUNlnC <- lnSRPsubs$lnC

SRP_10folds <- customfitMSE(kfolds=10, iterations=10)
SRP_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(SRP_10folds)[11]<-"SPAN"
SRP_melt <- melt(SRP_10folds, id.vars="SPAN", value.name="MSE")
SRP_one<-subset(SRP_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_SRP <- ggplot(data=SRP_melt, aes(x=SPAN, y=MSE, group=variable))
gg_SRP <- gg_SRP + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

###################### TP


FUNlnQ <- lnTPsubs$lnQ
FUNlnC <- lnTPsubs$lnC

TP_10folds <- customfitMSE(kfolds=10, iterations=10)
TP_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(TP_10folds)[11]<-"SPAN"
TP_melt <- melt(TP_10folds, id.vars="SPAN", value.name="MSE")
TP_one<-subset(TP_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_TP <- ggplot(data=TP_melt, aes(x=SPAN, y=MSE, group=variable))
gg_TP <- gg_TP + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

####################################### TSS


FUNlnQ <- lnTSSsubs$lnQ
FUNlnC <- lnTSSsubs$lnC

TSS_10folds <- customfitMSE(kfolds=10, iterations=10)
TSS_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(TSS_10folds)[11]<-"SPAN"
TSS_melt <- melt(TSS_10folds, id.vars="SPAN", value.name="MSE")
TSS_one<-subset(TSS_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_TSS <- ggplot(data=TSS_melt, aes(x=SPAN, y=MSE, group=variable))
gg_TSS <- gg_TSS + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

############################ Cl

FUNlnQ <- lnClsubs$lnQ
FUNlnC <- lnClsubs$lnC

Cl_10folds <- customfitMSE(kfolds=10, iterations=10)
Cl_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(Cl_10folds)[11]<-"SPAN"
Cl_melt <- melt(Cl_10folds, id.vars="SPAN", value.name="MSE")
Cl_one<-subset(Cl_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_Cl <- ggplot(data=Cl_melt, aes(x=SPAN, y=MSE, group=variable))
gg_Cl <- gg_Cl + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))

####################### SO4


FUNlnQ <- lnSO4subs$lnQ
FUNlnC <- lnSO4subs$lnC

SO4_10folds <- customfitMSE(kfolds=10, iterations=10)
SO4_10folds[11]<-c(seq(0.1, 1, by=0.05))
names(SO4_10folds)[11]<-"SPAN"
SO4_melt <- melt(SO4_10folds, id.vars="SPAN", value.name="MSE")
SO4_one<-subset(SO4_melt, variable=="iter 1", select= c (SPAN, MSE))

### make spaghetti plot in ggplot2

gg_SO4 <- ggplot(data=SO4_melt, aes(x=SPAN, y=MSE, group=variable))
gg_SO4 <- gg_SO4 + geom_line(colour="blue") + xlim(0.1, 1) + scale_x_continuous(breaks=seq(0.1, 1, by = 0.1))


#######################
#library(gridExtra)
grid.arrange(gg_NO3, gg_TN, gg_SRP, gg_TP, gg_TSS, gg_Cl, gg_SO4, ncol=1)

NO3_melt$Const <- "NO[3]"
TN_melt$Const <- "TN"

SRP_melt$Const <- "SRP"
TP_melt$Const <- "TP"
TSS_melt$Const <- "TSS"
Cl_melt$Const <- "Cl"
SO4_melt$Const <- "SO[4]"

all_melt <- rbind(NO3_melt, TN_melt, SRP_melt, TP_melt, TSS_melt, Cl_melt, SO4_melt)
# levels(all_melt$Const)[levels(all_melt$Const)=="NO3"] <- "NO[3]"
# levels(all_melt$Const)[levels(all_melt$Const)=="SO4"] <- "SO[4]"

sghetti <- ggplot(all_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#56B4E9", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free")
sghetti <- sghetti + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti <- sghetti + labs(x="f", y="Prediction MSE") + theme(axis.title=element_text(size=14))
sghetti

#convert each facet into a grob so we can arrange in grid of 4x2 with bottom right 
#being the conceptual plot (after it's turned into a grob)
sghetti_theme <- theme(panel.background=element_rect(fill="white", colour=NA), axis.title=element_text(size=24), axis.text=element_text(size=24), strip.text.x = element_text(size = 32))

sghetti_NO3 <- ggplot(NO3_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + geom_vline(xintercept=0.5, linetype="dashed", colour="gray") + sghetti_theme
sghetti_NO3
sghetti_TN <- ggplot(TN_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti_SRP <- ggplot(SRP_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti_TP <- ggplot(TP_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti_TSS <- ggplot(TSS_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti_Cl <- ggplot(Cl_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")
sghetti_SO4 <- ggplot(SO4_melt, aes(x=SPAN, y=MSE, group=variable)) + geom_line(colour="#E69F00", size=0.25) + facet_wrap(~ Const, ncol=2, nrow=4, labeller = label_parsed, scales="free") + labs(x=expression(italic("f")), y="PMSE") + sghetti_theme + geom_vline(xintercept=0.5, linetype="dashed", colour="gray")



g1 <- ggplotGrob(sghetti_NO3)
g2 <- ggplotGrob(sghetti_TN)
g3 <- ggplotGrob(sghetti_SRP)
g4 <- ggplotGrob(sghetti_TP)
g5 <- ggplotGrob(sghetti_TSS)
g6 <- ggplotGrob(sghetti_Cl)
g7 <- ggplotGrob(sghetti_SO4)
g8 <- ggplotGrob(p)

g <- rbind(g1,g2,g3,g4,g5,g6,g7,g8, size="first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths) #use largest width
grid.newpage()
grid.draw(g) #HOLY SHIT, BASED HADLEY


g_arrange <- arrangeGrob(grobs=list(g8,g1,g2,g3,g4,g5,g6,g7), nrow=4, ncol=2, widths=list(3,3), heights=list(9,9,9,9))
grid.draw(g_arrange)




#tifffilename <- paste(site, "_", "FACplots%02d", ".tiff", sep="")
tiff(file = "Kings_MSE_tradeoff.tiff", res=72, width = 1400, height = 1600, compression = "none") 

g_arrange <- arrangeGrob(grobs=list(g8,g1,g2,g3,g4,g5,g6,g7), nrow=4, ncol=2, widths=list(1,1), heights=list(1,1,1,1))
grid.draw(g_arrange)
dev.off()









# grob <- ggplotGrob(sghetti)
# grob[["grobs"]][["strip_t.1"]][["children"]][[2]][["label"]] <- expression(NO[3])
# grid.draw(grob)








