##############trying to figure what the fuck is going on, Volume 2
#for 5 folds, 5 iterations

FUNlnQ <- lnNO3subs$lnQ
FUNlnC <- lnNO3subs$lnC

FUNdat <- data.frame(FUNlnQ, FUNlnC)
totalOBS <- NROW(FUNdat)
FUNdat$INDEX <- c(1:totalOBS)
ncol.cal = 1
kfoldsbyiterations <- data.frame(matrix(NA, nrow=5, ncol=ncol.cal))
colnames(kfoldsbyiterations)<-paste("five folds")

iterationvector <- vector(mode="numeric", length=5)

kfolds <- createFolds(y=FUNdat$FUNlnC, k=5)
kfolds.frame <- data.frame(kfolds)
#vector of MSE across folds for each f (should be 17 total)
MSEofKnumfolds <- vector(mode="numeric", length=17)
f=0.15
validset<-as.vector(kfolds.frame[1])
colnames(validset)<-"INDEX"
calibset<- as.vector(kfolds.frame[kfolds.frame %notin% kfolds.frame[1]])
calibset<-melt(calibset)
names(calibset)[1]<-"FOLD"
names(calibset)[2]<-"INDEX"
VALIDset<-merge(validset, FUNdat, by = "INDEX")
CALIBset<-merge(calibset, FUNdat, by = "INDEX")
#we don't need the 'fold' column
CALIBset$FOLD <- NULL

#loess
calibfit=loess(FUNlnC~FUNlnQ, span=f, data=CALIBset)
validation <- predict(calibfit, newdata=VALIDset)
validnum <- NROW(VALIDset)

sum=0
for (m in 1:validnum){
  error <- (VALIDset[m,3])-(validation[m])
  sum <- ((error^2)+sum)
}



#########I THINK WE GOT IT


#the 'not in' function
`%notin%` <- function(x,y) !(x %in% y) 

#LORD SATAN, HELP ME
strat_kfoldCV_fxn <- function(lowfolds, highfolds, iterations){
  FUNdat <- data.frame(FUNlnQ, FUNlnC)
  totalOBS <- NROW(FUNdat)
  FUNdat$INDEX <- c(1:totalOBS)
  ncol.cal =  (highfolds - lowfolds) +1
  kfoldsbyiterations <- data.frame(matrix(NA, nrow=iterations, ncol=ncol.cal))
  colnames(kfoldsbyiterations)<-paste("folds", c(lowfolds:highfolds))
  for (kay in lowfolds:highfolds){
    #create vector for the results at this fold to put into kfoldsbyiterations
    iterationvector <- vector(mode="numeric", length=iterations)
    for (i in 1:iterations){
      
      #split data into 'stratified' folds
      kfolds<-createFolds(y=FUNdat$FUNlnC, k=kay)   #######
      kfolds.frame <- data.frame(kfolds)
      #vector of MSE across folds for each f (should be 17 total)
      MSEofKnumfolds <- vector(mode="numeric", length=17)
      
      for (f in seq(0.15,0.95,by=0.05)){
        totalSSE=0
        for (q in 1:kay){
          validset<-as.vector(kfolds.frame[q])
          colnames(validset)<-"INDEX"
          #create calibration set that is everything else that isn't in the q fold
          calibset <- as.data.frame(kfolds.frame[kfolds.frame %notin% kfolds.frame[q]])
          #collapse into one long list
          calibset$DUMMYINDEX<-
          calibset <- melt(calibset)
          names(calibset)[1]<-"FOLD"
          names(calibset)[2]<-"INDEX"
          VALIDset<-merge(validset, FUNdat, by = "INDEX")
          CALIBset<-merge(calibset, FUNdat, by = "INDEX")
          #we don't need no stinkin' 'fold' columns
          CALIBset$FOLD <- NULL
          
          #fit the LOESS regression to calibration set then test the validation set
          calibfit=loess(FUNlnC~FUNlnQ, span=f, data=CALIBset)
          validation <- predict(calibfit, newdata=VALIDset)
          validnum <- NROW(VALIDset)  
          sum=0
          for (m in 1:validnum){
            error <- (VALIDset[m,3])-(validation[m])
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
    currentkindex <- (kay - lowfolds)+1
    kfoldsbyiterations[currentkindex]<-iterationvector 
    print(kay)  
  }#end for k fold
  return(kfoldsbyiterations)
}


#test on a dataset
#PLEASE
#PLEEEEEEEEAAAAAAAAASE

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


testNO3_5to5_5iter <- strat_kfoldCV_fxn(5,5,5)









