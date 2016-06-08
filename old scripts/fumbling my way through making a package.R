install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
library(roxygen2)


install.packages("devtools", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("devtools", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")

install.packages("roxygen2", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("roxygen2", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")

install.packages("withr", lib="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")
.libPaths( c( .libPaths(), "C:/Users/zpsimpso/Documents/Documents/R_PACKAGES"))
library("withr", lib.loc="C:/Users/zpsimpso/Documents/Documents/R_PACKAGES")


# setwd("C:/Users/zpsimpso/Documents/Documents/TRENDS/package/")
# create("zachs")
# #may need to restart session
# install("zachs")
# library("zachs")


setwd("F:/TRENDS/TAFA/")
create("TAFA")

setwd("F:/TRENDS/AWRC_github/")
create("TAFA")

# DESCRIPTION <- read.dcf("TAFA/DESCRIPTION")
setwd("./TAFA/")
setwd("F:/TRENDS/my_packages/TAFA/TAFA/")

document()
#wow, had to change @example to @examples and now it WORKS. FUCKING DAMN.
setwd("..")
install("TAFA")

setwd("F:/TRENDS/my_packages/TAFA/TAFA/data/")
data("IR59")

devtools::use_data(IR59, overwrite = TRUE) #looks like it works


#NEXT STEP: GET THIS FUCKER ON GITHUB
#GOD FUUUUUCKING DAAAAAAAAAAAAAMMMMMMMMMMMMMMNNNNNN

install.packages('httr') #duh, everyone knows you have to have this package to do stuff
library('httr')

devtools::install_github('Trend-Analysis/tree/master/TAFA','zpsimpso')

devtools::install_github("zpsimpso/TAFA/TAFA")


setwd("F:/TRENDS/")

IR59 <- read.csv("IR59_USGS_data_01to09.csv", stringsAsFactors = FALSE)
colnames(IR59)[2]<-"Flow_cfs"
colnames(IR59)[3]<-"TP_mgL"
IR59$Flow_cfs<- as.numeric(IR59$Flow_cfs)
IR59$Date<-as.Date(IR59$Date, "%m/%d/%Y")
IR59$Date<-as.POSIXct(IR59$Date)

write.csv(IR59, "IR59_TP_2001-2009.csv")
save(IR59, file="IR59_TP_2001-2009.rda")

load("IR59_TP_2001-2009.rda")
data("IR59_TP_2001-2009.rda")

IR59_tafa<-tafa(IR59$Flow_cfs, IR59$TP_mgL, IR59$Date)

#see the flow-adjustment process
plot(IR59_tafa$lnC ~ IR59_tafa$lnQ)
#add the fitted loess line
j <- order(IR59_tafa$lnQ) #have to order x-values for base plot
lines(IR59_tafa$lnQ[j], IR59_tafa$loess_fit$fitted[j], col="red", lwd=3)
#observe the flow-adjusted concentrations over time
plot(IR59_tafa$dates, IR59_tafa$FACs) 
abline(lm(IR59_tafa$FACs~IR59_tafa$dates), col="blue") #note the decrease
#see percent change in TP over this period (% change per year)
IR59_tafa$perc_slope
#[1] -10.03825

IR59 <- read.csv("IR59_TP_2001-2009.csv", stringsAsFactors = FALSE)
IR59$Date <- as.POSIXct(IR59$Date) #convert the date column to a datetime format


