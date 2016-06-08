
setwd("F:/TRENDS/UWRB TRENDS/UWRB/RC45/truncated/")
getwd()

site <- "RC45"
#make sure the flow and concentrations files has the dates in the
#YYYY-MM-DD format

siteflow<-read.csv("flow.csv", header=TRUE)
concentrations<-read.csv("concentrations.csv", header=TRUE)

lnNO3subs <- ReadInWQData(siteflow,concentrations)[[1]]
lnTNsubs <- ReadInWQData(siteflow, concentrations)[[2]]
lnSRPsubs <- ReadInWQData(siteflow, concentrations)[[3]]
lnTPsubs <- ReadInWQData(siteflow, concentrations)[[4]]
lnTSSsubs <- ReadInWQData(siteflow, concentrations)[[5]]
lnClsubs <- ReadInWQData(siteflow, concentrations)[[6]]
lnSO4subs <- ReadInWQData(siteflow, concentrations)[[7]]


############################
const <- "TN"
FUNlnQ <- lnTNsubs$lnQ
FUNlnC <- lnTNsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnTNsubs)
lnTNsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTNsubs[7] <- fitmse$fit$residuals

theme_set(theme_grey(base_size = 12))

TNplot <- ggplot(data = lnTNsubs, aes(x=TN.Date, y=TN.concentrations.TN))
TNplot <- TNplot + geom_point() + labs(x="Year", y=expression(TN~(mg~L^{-1}))) +
  annotate("text", x=lnTNsubs[30,1], y=Inf, vjust=2.0, size = 14, label="A") +
  ggtitle(site)
TNplot #I had to use a data point from the dataframe as the x for the label


FA_TNplot <- ggplot(data=lnTNsubs, aes(x=lnQ, y=lnC)) + geom_point()
FA_TNplot <- FA_TNplot + geom_smooth(method="loess", span=fitmse$fit$pars$span, se=FALSE, colour="red") + labs(x="ln(Q)",y="ln(TN)")+
  annotate("text", x=-Inf, hjust=-1.5, y=Inf, vjust=2.0, size=14, label="B") + annotate("segment", x=8, xend=8.5, y=-1, yend=-1, colour="red") +
  annotate("text", x=8.6, y=-1, label="LOESS line", hjust=0)
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
perc_slope_text <- bquote(.(round(percentslope,2)))
perc_slope_char <- paste("m =", as.character(as.expression(perc_slope_text)))
if (pval < 0.01) {
  pval_text <- "p < 0.01"
} else {
  pval_exp <- bquote(.(format(round(pval, 2), nsmall = 2)))
  pval_text <- paste("p =", as.character(as.expression(pval_exp)))
}
midpoint <- mean(range(lnTNsubs$V7))
stdev <- sd(lnTNsubs$V7)
FACplot <- ggplot(data=lnTNsubs, aes(x=TN.Date, y=V7)) + geom_point()
FACplot <- FACplot + labs(x="Year", y="TN FACs") + scale_y_continuous(limits=c((0 - 4*stdev),(0 + 4*stdev))) 
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTNsubs[180,1], xend=lnTNsubs[190,1], y=(midpoint - 3*stdev), yend=(midpoint - 3*stdev), colour="blue")+
  annotate("text", x=lnTNsubs[200,1], y=(midpoint - 3*stdev), label=paste("Linear Regression"), hjust = 0) 
#FACplot <- FACplot + annotate("text", x=lnTNsubs[40,1], y=(midpoint - 3.5*stdev), label=perc_slope_char, hjust=0)
#+ annotate("segment", x=lnTNsubs[30,1], xend=lnTNsubs[50,1], y=-2.1, yend=-2.1, colour="red")+annotate("text", x=lnTNsubs[60,1], y=-2.1, label="Smoother Line", hjust = 0)
FACplot <- FACplot + annotate("text", x=lnTNsubs[200,1], y=(midpoint - 3.5*stdev), label=paste(perc_slope_char, ";", pval_text), hjust=0) +
  annotate("text", x=lnTNsubs[30,1], y=Inf, vjust=1.5, hjust=+0.1,size=14, label="C")
FACplot
TN_FACplot <- FACplot

TN_3step <- multiplot(TNplot, FA_TNplot, TN_FACplot, cols=1)

png("TN_3step.png", width=400, height=600)
TN_3step <- multiplot(TNplot, FA_TNplot, TN_FACplot, cols=1)
dev.off()



#####################################
const <- "TSS"
FUNlnQ <- lnTSSsubs$lnQ
FUNlnC <- lnTSSsubs$lnC
fit<-loess(lnC~lnQ, span=0.5, data=lnTSSsubs)
lnTSSsubs[6]<-fit$residuals
fitmse <- loess.wrapperMSE(FUNlnQ, FUNlnC, folds = 10, iteration = 10)
lnTSSsubs[7] <- fitmse$fit$residuals

TSSplot <- ggplot(data = lnTSSsubs, aes(x=TSS.Date, y=TSS.concentrations.TSS)) + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)), panel.background = element_rect(fill="white", colour=NA))
TSSplot <- TSSplot + geom_point() + labs(x="Year", y=expression(TSS~(mg~L^{-1})))
TSSplot <- TSSplot + annotate("text", x=lnTSSsubs[30,1], y=750, size=rel(11), label="A")
TSSplot #I had to use a data point from the dataframe as the x for the label

FA_TSSplot <- ggplot(data=lnTSSsubs, aes(x=lnQ, y=lnC)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)), panel.background = element_rect(fill="white", colour=NA))
FA_TSSplot <- FA_TSSplot + geom_smooth(method="loess", span=0.5, se=FALSE, colour="red") + labs(x="log(Q)",y="log(TSS)") + annotate("segment", x=7.4, xend=7.9, y=-1, yend=-1, colour="red") + annotate("text", x=8.0, y=-1, label="LOESS", size=rel(4), hjust=0)
FA_TSSplot <- FA_TSSplot + annotate("text", x=2.5, y=6, size=rel(11), label="B")
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
FACplot <- ggplot(data=lnTSSsubs, aes(x=TSS.Date, y=V7)) + geom_point() + theme(axis.text = element_text(size=rel(1.3)), axis.title.x=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.title.y=element_text(size=rel(1.5)), panel.background = element_rect(fill="white", colour=NA))
FACplot <- FACplot + labs(x="Year", y="TSS FACs")
FACplot <- FACplot + stat_smooth(method=lm, se=FALSE, colour="blue") + annotate("segment", x=lnTSSsubs[30,1], xend=lnTSSsubs[50,1], y=-2.2, yend=-2.2, colour="blue")+
  annotate("text", x=lnTSSsubs[60,1], y=-2.2, label="Linear Regression", hjust = 0, size=rel(4.5)) 
FACplot <- FACplot #+ annotate("segment", x=lnTSSsubs[30,1], xend=lnTSSsubs[50,1], y=-2.5, yend=-2.5, colour="red")+annotate("text", x=lnTSSsubs[60,1], y=-2.5, label="Smoother Line", hjust = 0, size=rel(4.5))
FACplot <- FACplot #+ ggtitle(site) #+ stat_smooth(se=FALSE, span=0.75, colour = "red")
FACplot <- FACplot + annotate("text", x=lnTSSsubs[60,1], y=-3, label=perc_slope_char, parse=TRUE, hjust=0, size=rel(4.5))
FACplot <- FACplot + scale_y_continuous(limits=c(-3.3,3.3)) #+ annotate("text", x=lnTSSsubs[60,1], y=-3.255, label=pval_text, hjust=0, size=rel(4.5)) 
FACplot <- FACplot + annotate("text", x=lnTSSsubs[30,1], y=2, size=rel(11), label="C")
TSS_FACplot <- FACplot


TSS_3step <- multiplot(TSSplot, FA_TSSplot, TSS_FACplot, cols=1)

#see following link for more info:
# https://github.com/hadley/ggplot2/wiki/Align-two-plots-on-a-page
g1 <- ggplotGrob(TSSplot)
g2 <- ggplotGrob(FA_TSSplot)
g3 <- ggplotGrob(TSS_FACplot)
g <- rbind(g1,g2,g3, size="first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths) #use largest width
grid.newpage()
grid.draw(g) #HOLY SHIT, BASED HADLEY



png("TSS_3step.png", width=400, height=600)
#arrange_ggplot2(TSSplot, FA_TSSplot, TSS_FACplot, ncol=1)
TSS_3step <- multiplot(TSSplot, FA_TSSplot, TSS_FACplot, cols=1)
dev.off()


#tifffilename <- paste(site, "_", "FACplots%02d", ".tiff", sep="")
tiff(file = "TSS_3step.tiff", res=72, width = 450, height = 750, compression = "none") 
TSS_3step <- multiplot(TSSplot, FA_TSSplot, TSS_FACplot, cols=1)

dev.off()









