gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}




setwd("F:/TRENDS/UWRB TRENDS/UWRB/WFWR/")

site <- "WFWR"
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

scond <- read.csv("phys_params.csv",header=TRUE)
scond$Date <- as.POSIXct(scond$Date)

cl_scond_temp <- merge(scond, lnClsubs, by.x = "Date", by.y = "Cl.Date")
cl_scond <- na.omit(data.frame(cl_scond_temp$Date,cl_scond_temp$Cond,cl_scond_temp$Cl.concentrations.Cl))
clscondplot2 <- ggplot(data = cl_scond, aes(x=cl_scond[2],y=cl_scond[3])) + geom_point()
clscondplot2

clplot <- ggplot(data=lnClsubs, aes(x=lnClsubs[5], y=lnClsubs[4])) + geom_point()
clplot
SO4plot <- ggplot(data=lnSO4subs, aes(x=lnSO4subs[5], y=lnSO4subs[4])) + geom_point()
SO4plot

SO4_scond_temp <- merge(scond, lnSO4subs, by.x = "Date", by.y = "SO4.Date")
SO4_scond <- na.omit(data.frame(SO4_scond_temp$Date, SO4_scond_temp$Cond, SO4_scond_temp$SO4.concentrations.SO4))

SO4scondplot <- ggplot(data = SO4_scond, aes(y=SO4_scond[2],x=SO4_scond[3])) + geom_point()
SO4scondplot

Lab_cond <- read.csv("conductivity.csv")
Lab_cond$Date <- as.POSIXct(Lab_cond$Date)

cl_scond <- merge(Lab_cond, lnClsubs, by.x = "Date", by.y = "Cl.Date")
condplot <- ggplot(data=cl_scond, aes(x=cl_scond$lnQ, y=cl_scond$Conductivity)) + geom_point()
condplot

clplot <- ggplot(data=cl_scond, aes(x=cl_scond$lnQ, y=cl_scond$lnC)) + geom_point()
clplot

geomeanCl <- gm_mean(cl_scond$Cl.concentrations.Cl)
Cl_mad <- mad(cl_scond$Cl.concentrations.Cl)
Cl_Cond_plot <- ggplot(data=cl_scond, aes(x=cl_scond$Cl.concentrations.Cl, y=cl_scond$Conductivity)) + geom_point()
Cl_Cond_plot + geom_vline(xintercept = geomeanCl + 2*Cl_mad, linetype="dashed", size = 0.75)

sub_Cl_scond <- subset(cl_scond, cl_scond$Cl.concentrations.Cl < (geomeanCl + 2*Cl_mad))
subCl_Cond_plot <- ggplot(data=sub_Cl_scond, aes(x=sub_Cl_scond$Cl.concentrations.Cl, y=sub_Cl_scond$Conductivity)) + geom_point()
subCl_Cond_plot

colnames(sub_Cl_scond)[3] <- "Q"
colnames(sub_Cl_scond)[4] <- "Cl"

linestats <- aov(Conductivity ~ Cl, data = sub_Cl_scond)
linestats <- anova(lm(Conductivity ~ Cl, data = sub_Cl_scond))
line <- lm(Conductivity ~ Cl, data = sub_Cl_scond)

Flow_Clplot <- ggplot(data=sub_Cl_scond, aes(x=lnQ, y=lnC))+geom_point()
Flow_Clplot + stat_smooth()

clq.lm <- lm(lnC ~ lnQ, data=sub_Cl_scond)

#tried ?segmented and it's not gonna work

#Lott and Stewart 2012 say RO_c isn't very sensitive, maybe the same is somewhat true for BF_c?

RO_c <- min(sub_Cl_scond$lnC)
BF_c <- max(sub_Cl_scond$lnC)
#Q_c = a'Q^(b')
#^a and b are intercept and slope from best fit line on log cond vs log Q plot
#Q_bf = aQ^b + cQ
#where Q is streamflow
#a is (a'/(BF_c - RO_c))
#b is 1 + b'
#c is (-RO_c)/(BF_c - RO_c)
#RO is difference between Q and Q_bf
Qcline <- lm(lnC ~ lnQ, data = sub_Cl_scond)
lm_sum<-summary.lm(Qcline)
b_prime <- lm_sum$coefficients[2] #slope
a_prime <- lm_sum$coefficients[1] #intercept

a <- (a_prime)/(BF_c - RO_c)
b <- 1 + b_prime
c <- (-RO_c)/(BF_c - RO_c)

#...
#after some deeper reflection about this method...
#it doesn't make sense to use it for our data,
#the same Q will yield the same BF value throughout the year
#the whole point was that at the same Q, you could have a high BF in spring but low BF in late summer
#........

#HOLD THE PHONE
#use the CMB method given in Lott and Stewart 2012
#Q_bf = Q * (C - RO_c)/(BF_c - RO_c)
#use Cl data in place of conductivity (C)
RO_c <- min(sub_Cl_scond$Cl)
BF_c <- max(sub_Cl_scond$Cl)
Q_bf <- sub_Cl_scond$Q * ((sub_Cl_scond$Cl - RO_c)/(BF_c - RO_c))
Qbframe <- data.frame(sub_Cl_scond$Date,Q_bf, "base")
colnames(Qbframe)[1] <- "Date"
colnames(Qbframe)[2] <- "Q"
colnames(Qbframe)[3] <- "type"
Qtframe <- data.frame(sub_Cl_scond$Date,sub_Cl_scond$Q, "total")
colnames(Qtframe)[1] <- "Date"
colnames(Qtframe)[2] <- "Q"
colnames(Qtframe)[3] <- "type"
Qframe <- rbind(Qtframe, Qbframe)

Qplot <- ggplot(data = Qframe, aes(x=Date, y=log10(Q), colour = type)) + geom_point()
Qplot 

Qplotnorm <- ggplot(data = Qframe, aes(x=Date, y=Q, colour = type)) + geom_point()
Qplotnorm

#########Ok reference baseflow, but appears to underestimate low bf values
#Try to implement Eckhardt model or something similar

#just found package 'EcoHydRology'
hydrograph(streamflow = Qtframe$Q, timeSeries = Qtframe$Date, S.units = "ft3s")

#try 1 pass
SepFlow_1 <- BaseflowSeparation(streamflow = siteflow$Flow, filter_parameter = 0.925, passes = 1)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_1$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_1$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")
#FUCK YES

#try 2 passes
SepFlow_1 <- BaseflowSeparation(streamflow = siteflow$Flow, filter_parameter = 0.925, passes = 2)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_1$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_1$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")

#try 3 passes
SepFlow_1 <- BaseflowSeparation(streamflow = siteflow$Flow, filter_parameter = 0.925, passes = 3)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_1$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_1$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")

#Try Eckhardt filter
SepFlow_eck <- Eckhardt(streamflow = siteflow$Flow, filter_parameter = 0.925, BFI_max = 0.8, passes = 2)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_eck$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_eck$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 <- flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2013-09-01", "2013-09-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")
flowplot_1 + scale_y_continuous(limits = c(0,300))
#it looks SOOOOOOOOO GOOOOOOOOOOOOOOOOOOOOOOOOOOD
Eck_flow <- data.frame(siteflow,SepFlow_eck)

##################
colnames(lnNO3subs)[1] <- "Date"
lnNO3_bf <- merge(lnNO3subs, Eck_flow, by="Date", all.lnNO3subs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnNO3_bf)[2:3] <- c("Q", "NO3")
ggplot(data = lnNO3_bf, aes(x=log(bt), y=log(NO3))) + geom_point()
ggplot(data = lnNO3_bf, aes(x=log(qft), y=log(NO3))) + geom_point()

lnNO3_bf$lnbt <- log(lnNO3_bf$bt)
lnNO3_bf$lnqft <- log(lnNO3_bf$qft)

for (ii in 1:(length(lnNO3_bf$lnqft))){
  if(lnNO3_bf[ii, 10] == "-Inf"){
    lnNO3_bf[ii, 10] <- -4
  }
}

#library(rgl)
plot3d(lnNO3_bf$lnbt, lnNO3_bf$qft, lnNO3_bf$lnC, type="s", size=.3, xlim = c(-4, 10), ylim = c(0, 2000), zlim = c(-5,2))

colnames(lnTSSsubs)[1] <- "Date"
lnTSS_bf <- merge(lnTSSsubs, Eck_flow, by="Date", all.lnTSSsubs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnTSS_bf)[2:3] <- c("Q", "TSS")
ggplot(data = lnTSS_bf, aes(x=log(bt), y=log(TSS))) + geom_point()
ggplot(data = lnTSS_bf, aes(x=log(qft), y=log(TSS))) + geom_point()

lnTSS_bf$lnbt <- log(lnTSS_bf$bt)
lnTSS_bf$lnqft <- log(lnTSS_bf$qft)

for (ii in 1:(length(lnTSS_bf$lnqft))){
  if(lnTSS_bf[ii, 10] == "-Inf"){
    lnTSS_bf[ii, 10] <- -4
  }
}


plot3d(lnTSS_bf$lnbt, lnTSS_bf$lnqft, lnTSS_bf$lnC, type="s", size=0.75, xlim = c(0, 10), ylim = c(2, 10), zlim = c(-5,10))

new_qft <- data.frame(lnNO3_bf[which(lnNO3_bf$qft > 10), ])
ggplot(new_qft, aes(x=lnbt)) + geom_histogram()

new_qft <- data.frame(lnTSS_bf[which(lnTSS_bf$qft > 10), ])
plot3d(new_qft$lnbt, new_qft$lnqft, new_qft$lnC, type="s", size=0.75, xlim = c(-4, 10), ylim = c(-6, 10), zlim = c(0,8))
#there looks like a steep curving slope with storm data, may need to include bf
#data even after separation by filtering out predominant baseflow samples

lnNO3_bf$BFF <- ((lnNO3_bf$bt)/lnNO3_bf$Flow)
ggplot(lnNO3_bf, aes(x=BFF, y=lnC)) + geom_point()


lnTSS_bf$BFF <- ((lnTSS_bf$bt)/lnTSS_bf$Flow)
ggplot(lnTSS_bf, aes(x=BFF, y=lnC)) + geom_point() + geom_vline(xintercept = TSS_cp, linetype = 'dashed') + 
  geom_vline(xintercept = TSS_lowcp, colour = "red", linetype = 'dotted') + geom_vline(xintercept = TSS_highcp, colour = "red", linetype = 'dotted')

ncpa.TSS<-ncpa(x=lnTSS_bf$BFF, y=lnTSS_bf$lnC)
TSS_cp <- ncpa.TSS$cpobs[1] #is 0.444 for TSS data
TSS_lowcp <-ncpa.TSS$cpobs[6]
TSS_highcp <- ncpa.TSS$cpobs[10]

ncpa.NO3<-ncpa(x=lnNO3_bf$BFF, y=lnNO3_bf$lnC)
NO3_cp <- ncpa.NO3$cpobs[1] #is 0.51 for NO3 data
NO3_lowcp <-ncpa.NO3$cpobs[6]
NO3_highcp <- ncpa.NO3$cpobs[10]

ggplot(lnNO3_bf, aes(x=BFF, y=lnC)) + geom_point() + geom_vline(xintercept = 0.51, linetype = 'dashed')
  # geom_vline(xintercept = NO3_lowcp, colour = "red", linetype = 'dotted') + geom_vline(xintercept = NO3_highcp, colour = "red", linetype = 'dotted')
asd <- data.frame(lnNO3_bf[which(lnNO3_bf$lnC < -5), ])


#testing multivariate loess out, trying to plot it in 3d
tryitout <- loess(formula=(lnC ~ lnbt + lnqft), data=lnTSS_bf, span = 0.5)

data.frame(lnTSS_bf[which(lnTSS_bf$qft > 10), ])
#she ain't budging, Jim
m <- data.frame(lnTSS_bf[which(lnTSS_bf$qft > 10), ])
mod <- loess(lnC ~ lnbt + lnqft, data=m)
m$pred <- predict(mod)
mgrid_df <- predictgrid(mod, "lnbt", "lnqft", "lnC")

#####try a IRW site
setwd("F:/TRENDS/IRW TRENDS/IRW/Osage/")
site <- "Osage"
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

#Try Eckhardt filter
SepFlow_eck <- Eckhardt(streamflow = siteflow$Flow, filter_parameter = 0.925, BFI_max = 0.8, passes = 1)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_eck$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_eck$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 <- flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2014-03-01", "2014-08-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")
flowplot_1 + scale_y_continuous(limits = c(0,1000))
#it looks SOOOOOOOOO GOOOOOOOOOOOOOOOOOOOOOOOOOOD
Eck_flow <- data.frame(siteflow,SepFlow_eck)

##################
colnames(lnNO3subs)[1] <- "Date"
lnNO3_bf <- merge(lnNO3subs, Eck_flow, by="Date", all.lnNO3subs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnNO3_bf)[2:3] <- c("Q", "NO3")
ggplot(data = lnNO3_bf, aes(x=log(bt), y=log(NO3))) + geom_point()
ggplot(data = lnNO3_bf, aes(x=log(qft), y=log(NO3))) + geom_point()

lnNO3_bf$lnbt <- log(lnNO3_bf$bt)
lnNO3_bf$lnqft <- log(lnNO3_bf$qft)

for (ii in 1:(length(lnNO3_bf$lnqft))){
  if(lnNO3_bf[ii, 10] == "-Inf"){
    lnNO3_bf[ii, 10] <- -4
  }
}

colnames(lnTPsubs)[1] <- "Date"
lnTP_bf <- merge(lnTPsubs, Eck_flow, by="Date", all.lnTPsubs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnTP_bf)[2:3] <- c("Q", "TP")
ggplot(data = lnTP_bf, aes(x=log(bt), y=log(TP))) + geom_point()
ggplot(data = lnTP_bf, aes(x=log(qft), y=log(TP))) + geom_point()

lnTP_bf$lnbt <- log(lnTP_bf$bt)
lnTP_bf$lnqft <- log(lnTP_bf$qft)

for (ii in 1:(length(lnTP_bf$lnqft))){
  if(lnTP_bf[ii, 10] == "-Inf"){
    lnTP_bf[ii, 10] <- -4
  }
}

lnTP_bf$BFF <- ((lnTP_bf$bt)/lnTP_bf$Flow)

ncpa.TP<-ncpa(x=lnTP_bf$BFF, y=lnTP_bf$lnC)
TP_cp <- ncpa.TP$cpobs[1] 
TP_lowcp <-ncpa.TP$cpobs[6]
TP_highcp <- ncpa.TP$cpobs[10]

ggplot(lnTP_bf, aes(BFF, lnC))+geom_point()


TP_cp<-data.frame(lnTP_bf[which(lnTP_bf$BFF > 0.5), ])
ggplot(TP_cp, aes(x=Date, y=FACs))+geom_point()+geom_smooth(method = 'lm')
TP_loess <- loess.wrapperMSE(TP_cp$lnQ, TP_cp$lnC, iteration = 10, folds=10)
TP_cp$FACs <- TP_loess$fit$residuals

linreg <- lm(FACs~Date, data=TP_cp)
anova(linreg)

####IR59
setwd("F:/TRENDS/IRW TRENDS/IRW/IR59/")
site <- "IR59"
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

#Try Eckhardt filter
SepFlow_eck <- Eckhardt(streamflow = siteflow$Flow, filter_parameter = 0.925, BFI_max = 0.8, passes = 1)
siteflow$Date <- as.POSIXct(siteflow$Date)
pass1_T <- siteflow
pass1_BF <- data.frame(siteflow$Date, SepFlow_eck$bt)
colnames(pass1_BF)[1] <- "Date"
pass1_RO <- data.frame(siteflow$Date, SepFlow_eck$qft)
colnames(pass1_RO)[1] <- "Date"

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, colour=L1))
flowplot_1 <- flowplot_1 + geom_line(data=dat_1[dat_1$L1 == 2, ]) + geom_line(data=dat_1[dat_1$L1 == 3, ]) +geom_line(data=dat_1[dat_1$L1 == 1, ])
flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2011-04-01", "2011-05-30"))))

dat_1 <- melt(list(pass1_T, pass1_BF, pass1_RO), id.vars = "Date")
flowplot_1 <- ggplot(dat_1, aes(x = Date, y = value, fill=L1))
flowplot_1 <- flowplot_1 + geom_area(data=dat_1[dat_1$L1 == 2, ], alpha = 0.6) + geom_area(data=dat_1[dat_1$L1 == 3, ], alpha = 0.5) +geom_area(data=dat_1[dat_1$L1 == 1, ], alpha = 0.2)
flowplot_1 <- flowplot_1 + scale_x_datetime(limits = c(as.POSIXct(c("2014-03-01", "2014-08-30")))) + geom_line(data=dat_1[dat_1$L1 == 1, ], color = "red")
flowplot_1 + scale_y_continuous(limits = c(0,1000))
#it looks SOOOOOOOOO GOOOOOOOOOOOOOOOOOOOOOOOOOOD
Eck_flow <- data.frame(siteflow,SepFlow_eck)

##################
colnames(lnNO3subs)[1] <- "Date"
lnNO3_bf <- merge(lnNO3subs, Eck_flow, by="Date", all.lnNO3subs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnNO3_bf)[2:3] <- c("Q", "NO3")
ggplot(data = lnNO3_bf, aes(x=log(bt), y=log(NO3))) + geom_point()
ggplot(data = lnNO3_bf, aes(x=log(qft), y=log(NO3))) + geom_point()

lnNO3_bf$lnbt <- log(lnNO3_bf$bt)
lnNO3_bf$lnqft <- log(lnNO3_bf$qft)

for (ii in 1:(length(lnNO3_bf$lnqft))){
  if(lnNO3_bf[ii, 10] == "-Inf"){
    lnNO3_bf[ii, 10] <- -4
  }
}

colnames(lnTPsubs)[1] <- "Date"
lnTP_bf <- merge(lnTPsubs, Eck_flow, by="Date", all.lnTPsubs = TRUE, all.Eck_flow = FALSE, sort = TRUE)
colnames(lnTP_bf)[2:3] <- c("Q", "TP")
ggplot(data = lnTP_bf, aes(x=log(bt), y=log(TP))) + geom_point()
ggplot(data = lnTP_bf, aes(x=log(qft), y=log(TP))) + geom_point()

lnTP_bf$lnbt <- log(lnTP_bf$bt)
lnTP_bf$lnqft <- log(lnTP_bf$qft)

for (ii in 1:(length(lnTP_bf$lnqft))){
  if(lnTP_bf[ii, 10] == "-Inf"){
    lnTP_bf[ii, 10] <- -4
  }
}

lnTP_bf$BFF <- ((lnTP_bf$bt)/lnTP_bf$Flow)

ncpa.TP<-ncpa(x=lnTP_bf$BFF, y=lnTP_bf$lnC)
TP_cp <- ncpa.TP$cpobs[1] 
TP_lowcp <-ncpa.TP$cpobs[6]
TP_highcp <- ncpa.TP$cpobs[10]

ggplot(lnTP_bf, aes(BFF, lnC))+geom_point()


TP_cp<-data.frame(lnTP_bf[which(lnTP_bf$BFF > 0.5), ])
ggplot(TP_cp, aes(x=Date, y=FACs))+geom_point()+geom_smooth(method = 'lm')
TP_loess <- loess.wrapperMSE(TP_cp$lnQ, TP_cp$lnC, iteration = 10, folds=10)
TP_cp$FACs <- TP_loess$fit$residuals

linreg <- lm(FACs~Date, data=TP_cp)
anova(linreg)













