require(caret)

data(oil)
createDataPartition(oilType, 2)

x <- rgamma(50, 3, .5)
inA <- createDataPartition(x, list = FALSE)

plot(density(x[inA]))
rug(x[inA])

points(density(x[-inA]), type = "l", col = 4)
rug(x[-inA], col = 4)

createResample(oilType, 2)

createFolds(oilType, 10)
createFolds(oilType, 5, FALSE)

createFolds(rnorm(21))



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





#the 'not in' function
`%notin%` <- function(x,y) !(x %in% y) 

#LORD SATAN, HELP ME
strat_kfoldCV_fxn <- function(lowfolds, highfolds, iterations){
  FUNdat <- data.frame(FUNlnQ, FUNlnC)
  totalOBS <- NROW(FUNdat)
  ncol.cal =  (highfolds - lowfolds) +1
  kfoldsbyiterations <- data.frame(matrix(NA, nrow=iterations, ncol=ncol.cal))
  colnames(kfoldsbyiterations)<-paste("folds", c(lowfolds:highfolds))
  for (k in lowfolds:highfolds){
    iterationvector <- vector(mode="numeric", length=iterations)
    for (i in 1:iterations){
      
      #split data into 'stratified' folds
      kfolds<-createFolds(y=FUNdat$FUNlnC, k=k)
      #vector of MSE across folds for each f (should be 17 total)
      MSEofKnumfolds <- vector(mode="numeric", length=17)
        
          for (f in seq(0.15,0.95,by=0.05)){
            totalSSE=0
            for (q in 1:k){
              calibset <- kfolds[kfolds %notin% kfolds[q]]
              meltcalib <- melt(calibset)
              calibfit=loess(FUNlnC~FUNlnQ, span=f, data=meltcalib)
              validset <- FUNdat[kfolds[[q]],]
              validation <- predict(calibfit, newdata=validset)
              validnum <- NROW(validset)  
                sum=0
                for (m in 1:validnum){
                  error <- (validset[m,2])-(validation[m])
                  sum <- ((error^2)+sum)
                  }
              totalSSE<- (sum+totalSSE)
              
            }
            totalMSE<-((totalSSE)/totalOBS)
            #find index num of current f value to put into MSE vector
            fnum <- ((f-0.15)/0.05)+1
            MSEofKnumfolds[fnum]<-totalMSE
          }
        fvector <- c(seq(0.15,0.95,by=0.05))
        MSEbyf <- data.frame(fvector, MSEofKnumfolds)
        #get index of best f value by minimum MSE
        bestf <- order(MSEbyf$MSEofKnumfolds)[1]
        iterationvector[i]<-MSEbyf[bestf,1]
        print(i)
    }#end for ith iteration
  #get index of current kth fold to use for pasting iterationvector into matrix
  currentkindex <- (k-lowfolds)+1
  kfoldsbyiterations[currentkindex]<-iterationvector 
  print(k)  
  }#end for k fold
  return(kfoldsbyiterations)
}

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC


testNO3_5to10_10iter <- strat_kfoldCV_fxn(5,10,5)



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

testTN_5to10_50iter <- strat_kfoldCV_fxn(5,10,50)




#make separate files for each constituent + flow
SO4temp<-data.frame(concentrations$Date,concentrations$SO4)
names(SO4temp)[names(SO4temp)=="concentrations.Date"]<-"Date"

#replace zero values for concentration with MDL
#MDL for SO4 is 0.038 mg/L
SO4temp$concentrations.SO4 <- ifelse(SO4temp$concentrations.SO4 == 0, 0.038, SO4temp$concentrations.SO4)
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

testSO4_5to10_50iter <- strat_kfoldCV_fxn(5,10,50)


##############trying to figure what the fuck is going on
#for 5 folds, 5 iterations

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC

FUNdat <- data.frame(FUNlnQ, FUNlnC)
totalOBS <- NROW(FUNdat)
ncol.cal = 1
kfoldsbyiterations <- data.frame(matrix(NA, nrow=5, ncol=ncol.cal))
colnames(kfoldsbyiterations)<-paste("five folds")

iterationvector <- vector(mode="numeric", length=5)

kfolds <- createFolds(y=FUNdat$FUNlnC, k=5)

MSEofKnumfolds <- vector(mode="numeric", length=17)

f=0.15


calibset <- kfolds[kfolds %notin% kfolds[1]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.15, data=meltcalib)
validset <- FUNdat[kfolds[[1]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}
totalSSE1 <- sum

calibset <- kfolds[kfolds %notin% kfolds[2]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.15, data=meltcalib)
validset <- FUNdat[kfolds[[2]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE2 <- sum

calibset <- kfolds[kfolds %notin% kfolds[3]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.15, data=meltcalib)
validset <- FUNdat[kfolds[[3]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE3 <- sum

calibset <- kfolds[kfolds %notin% kfolds[4]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.15, data=meltcalib)
validset <- FUNdat[kfolds[[4]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}


totalSSE4 <- sum

calibset <- kfolds[kfolds %notin% kfolds[5]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.15, data=meltcalib)
validset <- FUNdat[kfolds[[5]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE5 <- sum

MSE_015 <- ((totalSSE1+totalSSE2+totalSSE3+totalSSE4+totalSSE5)/totalOBS)

MSEofKnumfolds[1]<-MSE_015


f=0.50


calibset <- kfolds[kfolds %notin% kfolds[1]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.5, data=meltcalib)
validset <- FUNdat[kfolds[[1]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}
totalSSE1 <- sum

calibset <- kfolds[kfolds %notin% kfolds[2]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.5, data=meltcalib)
validset <- FUNdat[kfolds[[2]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE2 <- sum

calibset <- kfolds[kfolds %notin% kfolds[3]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.5, data=meltcalib)
validset <- FUNdat[kfolds[[3]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE3 <- sum

calibset <- kfolds[kfolds %notin% kfolds[4]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.5, data=meltcalib)
validset <- FUNdat[kfolds[[4]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}


totalSSE4 <- sum

calibset <- kfolds[kfolds %notin% kfolds[5]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.5, data=meltcalib)
validset <- FUNdat[kfolds[[5]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE5 <- sum

MSE_50 <- ((totalSSE1+totalSSE2+totalSSE3+totalSSE4+totalSSE5)/totalOBS)

MSEofKnumfolds[8]<-MSE_50



f=0.75


calibset <- kfolds[kfolds %notin% kfolds[1]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.75, data=meltcalib)
validset <- FUNdat[kfolds[[1]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}
totalSSE1 <- sum

calibset <- kfolds[kfolds %notin% kfolds[2]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.75, data=meltcalib)
validset <- FUNdat[kfolds[[2]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE2 <- sum

calibset <- kfolds[kfolds %notin% kfolds[3]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.75, data=meltcalib)
validset <- FUNdat[kfolds[[3]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE3 <- sum

calibset <- kfolds[kfolds %notin% kfolds[4]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.75, data=meltcalib)
validset <- FUNdat[kfolds[[4]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}


totalSSE4 <- sum

calibset <- kfolds[kfolds %notin% kfolds[5]]
meltcalib <- melt(calibset)
calibfit=loess(FUNlnC~FUNlnQ, span=0.75, data=meltcalib)
validset <- FUNdat[kfolds[[5]],]
validation <- predict(calibfit, newdata=validset)
validnum <- NROW(validset)

sum=0
for (m in 1:validnum){
  error <- (validset[m,2])-(validation[m])
  sum <- ((error^2)+sum)
}

totalSSE5 <- sum

MSE_75 <- ((totalSSE1+totalSSE2+totalSSE3+totalSSE4+totalSSE5)/totalOBS)

MSEofKnumfolds[13]<-MSE_75

#############################END TESTING

###RETRY ON FORMULA



#LORD SATAN, HELP ME
strat_kfoldCV_fxn <- function(lowfolds, highfolds, iterations){
  FUNdat <- data.frame(FUNlnQ, FUNlnC)
  totalOBS <- NROW(FUNdat)
  ncol.cal =  (highfolds - lowfolds) +1
  kfoldsbyiterations <- data.frame(matrix(NA, nrow=iterations, ncol=ncol.cal))
  colnames(kfoldsbyiterations)<-paste("folds", c(lowfolds:highfolds))
  for (k in lowfolds:highfolds){
    iterationvector <- vector(mode="numeric", length=iterations)
    for (i in 1:iterations){
      
      #split data into 'stratified' folds
      kfolds<-createFolds(y=FUNdat$FUNlnC, k=k, list=FALSE)
      kfolds.frame <- data.frame(FUNdat, kfolds)
      #vector of MSE across folds for each f (should be 17 total)
      MSEofKnumfolds <- vector(mode="numeric", length=17)
      
      for (f in seq(0.15,0.95,by=0.05)){
        totalSSE=0
        for (q in 1:k){
          validrows<-which(kfolds.frame$kfolds==q)
          valid_q<-kfolds.frame[validrows,]
          calib_q<-kfolds.frame[-validrows,]
          
          calibfit=loess(FUNlnC~FUNlnQ, span=f, data=calib_q)
          
          validation <- predict(calibfit, newdata=valid_q)
          validnum <- NROW(valid_q)  
          sum=0
          for (m in 1:validnum){
            error <- (valid_q[m,2])-(validation[m])
            sum <- ((error^2)+sum)
          }
          totalSSE<- (sum+totalSSE)
          
        }
        totalMSE<-((totalSSE)/totalOBS)
        #find index num of current f value to put into MSE vector
        fnum <- ((f-0.15)/0.05)+1
        MSEofKnumfolds[fnum]<-totalMSE
      }
      fvector <- c(seq(0.15,0.95,by=0.05))
      MSEbyf <- data.frame(fvector, MSEofKnumfolds)
      #get index of best f value by minimum MSE
      bestf <- order(MSEbyf$MSEofKnumfolds)[1]
      iterationvector[i]<-MSEbyf[bestf,1]
      print(i)
    }#end for ith iteration
    #get index of current kth fold to use for pasting iterationvector into matrix
    currentkindex <- (k-lowfolds)+1
    kfoldsbyiterations[currentkindex]<-iterationvector 
    print(k)  
  }#end for k fold
  return(kfoldsbyiterations)
}



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



FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC


testNO3_5to10_10iter <- strat_kfoldCV_fxn(5,10,50)
#GOD DAMN IT































