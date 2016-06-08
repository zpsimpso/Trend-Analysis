##Program ncpa (nonparametric changepoint analysis) for R 2.9.3
##Copyright Ryan S. King (Baylor University), 2009-2010

#Modified and expanded version of deviance reduction method described in
#Qian, S. S., R. S. King, and C. J. Richardson.  2003.  Two statistical
#methods for detecting environmental thresholds.  Ecological Modelling
#166:87-97 and King R. S. and C. J. Richardson.  2003.  Integrating 
#bioassessment and ecological risk assessment: an approach to developing 
#numerical water-quality criteria.  Environmental Management 31:795-809.

#new features:
#deviance associated with each possible change point (dev.cp)
#permutation test (replaces chi-sq test, which was too liberal)
#minsplt allows control of minimum number of obs per binary partition
#plot.ncpa and plot.ncpadev for quick and simple graphical output

#Begin ncpa

ncpa<-function(x, y, boot=TRUE, perm=TRUE, nperm=1000, nboot=1000, minsplt=2) {
  
  cp <- function(x, y)
  {
    
    xx <-x
    yy <-y
    srtx <- sort(xx)
    mx<-srtx[minsplt:(length(srtx)-minsplt)]
    m <- length(mx)
    vi <- numeric()
    vi [m] <- sum((yy - mean(yy))^2)
    for(i in 1:(m-1)){
      vi[i] <- sum((yy[xx <= mx[i]] - mean(yy[xx <=
                                                mx[i]]))^2) + sum((yy[xx > mx[i]] - mean(
                                                  yy[xx > mx[i]]))^2)
      chngple<- mean(mx[vi == min(vi)])}
    for(i in 1:(m-1)){
      vi[i] <- sum((yy[xx >= mx[i]] - mean(yy[xx >=
                                                mx[i]]))^2) + sum((yy[xx < mx[i]] - mean(
                                                  yy[xx < mx[i]]))^2)
      chngpgt <- mean(mx[vi == min(vi)])}
    chngp <- ((chngple + chngpgt)/2)	
    out <- c(chngp, 1-min(vi)/sum((yy-mean(yy))^2), mean(yy[xx <= chngp]), mean(yy[xx > chngp]))
    
    out
  }
  
  
  dev <- function(x, y)
  {
    
    xx <-x
    yy <-y
    srtx <- sort((xx))
    mx<-srtx[minsplt:(length(srtx)-minsplt)]
    m <- length(mx)
    vi <- numeric()
    vi [m] <- sum((yy - mean(yy))^2)
    dev.le<-matrix(NA, m-1, 1)
    for(i in 1:(m-1)){
      vi[i] <- sum((yy[xx <= mx[i]] - mean(yy[xx <=
                                                mx[i]]))^2) + sum((yy[xx > mx[i]] - mean(
                                                  yy[xx > mx[i]]))^2)
      dev.le[i,]<-vi[i]   }
    dev.gt<-matrix(NA, m-1, 1)
    for(i in 1:(m-1)){
      vi[i] <- sum((yy[xx >= mx[i]] - mean(yy[xx >=
                                                mx[i]]))^2) + sum((yy[xx < mx[i]] - mean(
                                                  yy[xx < mx[i]]))^2)
      dev.gt[i,]<-vi[i]   }
    dev<- (dev.le+dev.gt)/2
    dev
  }
  
  if(perm) {
    ###CREATE EMPTY MATRICES FOR STORING NCPA RESULTS BY PERMUTATION
    devperm<-matrix(NA,nperm,1)
    devpermall<-matrix(NA, length(x)-minsplt*2, nperm)
    
    
    for (i in 1:nperm) {
      
      #RESAMPLE DATA WITHOUT REPLACEMENT 
      xprm <-x
      yprm <-y
      xxprm<-sample(xprm,replace=F) 
      devpermall[,i]<-dev(xxprm,yprm)
      devperm[i,]<-min(dev(xxprm,yprm))
    }
  }
  
  if(boot) {
    ###CREATE EMPTY MATRICES FOR STORING NCPA RESULTS BY BOOTREP
    bootrep<-matrix(NA,nboot,4)
    colnames(bootrep)<- c("cp", "r2", "mean left", "mean right")
    
    for (i in 1:nboot) {
      
      #RESAMPLE DATA WITH REPLACEMENT 
      xx <-x
      yy <-y
      mx <- sort((xx))
      m <- length(mx)
      
      nuid<-sample(m, replace=TRUE) 
      xxb<-xx[nuid]
      yyb<-yy[nuid]
      bootrep[i,]<-cp(xxb, yyb)
    }
  }
  
  if(boot) {
    cpboot<-matrix(NA,8,4)
    
    cpboot[1,1]<-mean(bootrep[,1])
    cpboot[2:8,1]<-quantile(bootrep[,1], probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    cpboot[1,2]<-mean(bootrep[,2])
    cpboot[2:8,2]<-quantile(bootrep[,2], c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    
    cpboot[1,3]<-mean(bootrep[,3])
    cpboot[2:8,3]<-quantile(bootrep[,3], c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    cpboot[1,4]<-mean(bootrep[,4])
    cpboot[2:8,4]<-quantile(bootrep[,4], c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
    
    rownames(cpboot)<-c("mean","5%","10%","25%","50%","75%","90%","95%")
    colnames(cpboot)<-c("cp", "r2", "mean left", "mean right")
  }
  
  dev.cp<-dev(x,y)
  devmax<-min(dev(x,y))
  devmax2<-matrix(NA, nperm, 1)
  for(i in 1:nperm){
    devmax2[i,1]<-ifelse(devmax<devperm[i,], 0, 1)
  }
  pperm<-(sum(devmax2[,1])+1)/nperm
  cpobs<-matrix(NA, 1, 10) 
  colnames(cpobs)<- c("cp", "r2", "mean left", "mean right","pperm", "5%", "25%", "50%", "75%", "95%")
  cpobs[1,1:4]<- cp(x,y)
  cpobs[1,6:10]<-cpboot[c(2,4,5,6,8), 1]
  cpobs[1,5]<-pperm
  ncpa.out<-list(bootrep, cpboot, devperm, dev.cp, cpobs)
  names(ncpa.out)<-c("bootrep", "cpboot", "devperm", "dev.cp","cpobs")
  return(ncpa.out)
  
}
#END ncpa






###PLOT THE CUMULATIVE THRESHOLD FREQUENCY FROM ncpa
##Note: remove "#" before "par" to set graph window; leave in place if multipanel ncpa.plots are desired
## use par(mfrow=c(3,2), oma=c(5,2,2,7), mar=c(1,5,1,0)) to set graph window outside of plot.ncpa for a good 6 panel layout with mtext for x and y2 labels

#Begin plot.ncpa
plot.ncpa<- function(x, y, ncpa.out, cpobs=TRUE, xlab = "", log= "", at=NULL, ntick=6, ymin=min(y), ymax=max(y), xmin=min(x), xmax=max(x), y1lab="y", y2lab=" ", y2cex=1.75, cex=1.75, cex.axis=1.75, cex.leg=1.5, cex.lab=1.75, col.cp="black", y1tlab=TRUE, xtlab=TRUE, y2tlab=TRUE, col.ctf="black", col.pt="black", pch=1, bty="u", ...){
  
  ##SET GRAPH DIMENSIONS
  #par(mar=c(5,5,4,2),oma=c(0,0,0,3))
  #4.2.2 OPTIONAL CUMULATIVE PROBABILITY CURVES FOR DEVIANCE SCORES
  ###SET TO FALSE IF NO DEVIANCE REDUCTION PLOTS ARE REQUIRED
  
  log.cft<-ifelse(log=="xy", "x", ifelse(log=="x", "x", ""))
  
  ##CUMULATIVE THRESHOLD FREQUENCY CURVE
  ctf<-table(ncpa.out$bootrep[,1])
  
  plot(rbind(min(ncpa.out$bootrep[,1]),matrix(sort(unique(ncpa.out$bootrep[,1])))),rbind(0, matrix(cumsum(ctf)/sum(ctf))),type="l",lty=5, lwd=2, log=log.cft, col=col.ctf, axes=FALSE, xlab="", ylab="", ylim=c(0,1), xlim=c(xmin,xmax))
  
  
  ##ADD SECOND Y-AXIS WITH PROBABILITY LABELS
  
  axis(4, pretty(c(0,1),6), cex.axis=cex.axis, labels=y2tlab)
  mtext(y2lab, side=4, line=3, cex=y2cex)
  
  ##FIND MINIMUM AND MAX VALUES TO COMPLETE CURVES
  cp.min<-min(ncpa.out$bootrep[,1])
  cp.max<-max(ncpa.out$bootrep[,1])
  
  ##COMPLETE PROBABILITY CURVES
  segments(min(x),0,cp.min, 0, col=col.ctf, lty=5, lwd=2)
  segments(cp.max,1,max(x),1, col=col.ctf, lty=5, lwd=2) 
  
  if(cpobs) {
    segments(ncpa.out$cpobs[,1],0, ncpa.out$cpobs[,1],1, col=col.cp, lty=1, lwd=3)
  }
  
  
  ## SET PLOT WINDOW TO ACCEPT NEW PARAMETERS
  par(new=TRUE)
  plot(x, y, col=col.pt, pch=pch, xlab=xlab, axes=FALSE, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, log=log, ylim=c(ymin, ymax), xlim=c(xmin,xmax), ylab=y1lab, ...) 
  
  ##OPTIONAL LOG-SCALE FOR X-AXIS IF log=”x”
  ifelse(log=="x", axis(1, at=at, cex.axis=cex.axis, labels=xtlab), ifelse(log=="xy", axis(1, at=at, cex.axis=cex.axis, labels=xtlab), axis(1, pretty(xmin:xmax, ntick), cex.axis=cex.axis, labels=xtlab)))
  
  ##ADD LEFT Y AXIS
  ifelse(log=="y", axis(2, at=at, cex.axis=cex.axis, labels=y1tlab), ifelse(log=="xy", axis(2, at=at, cex.axis=cex.axis, labels=y1tlab), axis(2, pretty(ymin:ymax, ntick), cex.axis=cex.axis, y1tlab)))
  
  
  ##CLOSE THE TOP OF THE PLOT
  box(which="plot", bty=bty)
  
}
#END plot.ncpa




###plot.ncpadev.  Plots the deviance reduction associated with each of 
#x-minsplt*2 possible partitions.  Useful for visualizing the strength of 
#the threshold (sharp peak = strong) and other possible thresholds 
#(multiple peaks).

plot.ncpadev<-function (x, ncpa.out, minsplt, xlabel = "x", log= "", at=NULL, xmin = min(x), xmax = max(x), ntick=6, y1label="Deviance reduction", y2label="Cumulative threshold frequency",  cex=1.75, cex.axis=1.75, cex.leg=1.5, cex.lab=1.75, col.cp="black", col.ctf="black", col.pt="black", pch=16) {
  
  ##SET GRAPH DIMENSIONS
  par(mar=c(5,5,4,2),oma=c(0,0,0,3))
  
  ctf<-table(ncpa.out$bootrep[,1])
  plot(rbind(min(ncpa.out$bootrep[,1]),matrix(sort(unique(ncpa.out$bootrep[,1])))),rbind(0, matrix(cumsum(ctf)/sum(ctf))),type="l",lty=5, lwd=2, log=log, col=col.ctf,axes=FALSE, ylim=c(0,1), xlim=c(xmin,xmax), xlab="",ylab="") 
  
  
  ##ADD SECOND Y-AXIS WITH PROBABILITY LABELS
  
  axis(4, pretty(c(0,1),6), cex.axis=cex.axis)
  mtext(y2label, side=4, line=3, cex=cex)
  
  ##FIND MINIMUM AND MAX VALUES TO COMPLETE CURVES
  cp.min<-min(ncpa.out$bootrep[,1])
  cp.max<-max(ncpa.out$bootrep[,1])
  
  ##COMPLETE PROBABILITY CURVES
  segments(xmin,0,cp.min, 0, col=col.ctf, lty=5, lwd=2)
  segments(cp.max,1,xmax,1, col=col.ctf, lty=5, lwd=2) 
  
  
  ## SET PLOT WINDOW TO ACCEPT NEW PARAMETERS
  par(new=TRUE)
  plot(sort(x[(1+minsplt):(max(length(x))-minsplt)]), (ncpa.out$dev.cp*-1)+max(ncpa.out$dev.cp), col=col.pt, pch=pch, xlab=xlabel, axes=FALSE, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, log=log, xlim=c(xmin,xmax), ylab=y1label) 
  
  ##OPTIONAL LOG-SCALE FOR X-AXIS IF log=”x”
  ifelse(log=="x", axis(1, at=at, cex.axis=cex.axis), axis(1, pretty(xmin:xmax, ntick), cex.axis=cex.axis))
  
  ##ADD LEFT Y AXIS
  axis(2, cex.axis=cex.axis)
  
  ##CLOSE THE TOP OF THE PLOT
  box(which="plot")
  
}
#END plot.ncpadev
