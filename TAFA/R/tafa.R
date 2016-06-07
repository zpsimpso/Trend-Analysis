#' Trend Analysis with Flow-Adjustment (TAFA)
#' 
#' Flow-adjust water quality data and perform a 
#' simple trend analysis with linear regression.
#' @param flow A vector of flow data.
#' @param constituent The constituent concentration data (vector) to analyze.
#' @param dates A vector of sample dates. Should either be in POSIXct format or 
#' easily coerced by as.POSIXct
#' @param control The control parameters. See \code{\link{tafa.control}}
#' @export
#' @seealso \code{\link[stats]{loess}}, \code{\link{tafa}}, \code{\link[bisoreg]{loess.wrapper}}
#' 
#' @details 
#' This function allows you to flow-adjust (i.e., remove the effect
#' of stream flow) water quality data to aid trend analysis. The methodology
#' follows the three step process outlined by White et al. (2004) and is
#' modified by Simpson and Haggard (2016). The modification allows for the 
#' smoothing parameter of LOESS, which is used in the flow-adjustment,
#' to be statistically optimized via a K-fold cross-validation procedure.
#' This function is inspired by the 'loess.wrapper' function in the 
#' 'bisoreg' package (available on CRAN).
#' 
#' The flow-adjusted concentrations (FACs) are then modelled over time with
#' simple linear regression to estimate the monotonic trend (interpreted
#' as percent change per year). Other trend tests can also be performed 
#' on the FACs such as Kendall's Tau/Seasonal Kendall Test
#'  (slope can be estimated with the Sen slope estimator). 
#'  See the 'zyp' and 'Kendall' packages available on CRAN. 
#' 
#' @references 
#' Simpson, Z.P. and B.E. Haggard. 2016. An optimized procedure for 
#' flow-adjustment of constituent concentrations for trend analysis.
#' In preparation.
#' 
#' White, K.L., B.E. Haggard, and I. Chaubey. 2004. Water quality at the
#' Buffalo National River, Arkansas, 1991-2001. Transactions of the
#' American Society of Agricultural Engineers 47(2):407-417.
#' 



tafa <- function(flow, constituent, dates, control=tafa.control()){
  Q <- flow
  C <- constituent
  dates <- as.POSIXct(dates) #I wrote this to work with POSIXct class, 'Date' type should convert fine
  #should also work if the dates were read in as a factor
  
  #log-transform, I use natural log here
  lnQ <- log(Q)
  lnC <- log(C)
  
  #flow-adjust with loess
  flow_adj <- loess.opt(x=lnQ, y=lnC, control=control)
  loess_fit <- flow_adj$fit #return the whole loess object
  FACs <- flow_adj$fit$residuals #flow-adjusted concentrations
  f_opt <- flow_adj$f_opt #give the optimized f-value used in loess (or the set span)
  
  #use FACs vs time to estimate trend
  FACs_lm <- lm(FACs~dates)
  p_val <- summary(FACs_lm)$coefficients[2,4] #get p-value of the F-test from linear regression
  m <- summary(FACs_lm)$coefficients[2] #the slope
  perc_slope <- (((exp(m)-1)*60*60*24*365*100)) #converts the slope into % change per year
  #note that the slope represents a proportional change in the concentration per time (hence, the % per time)
  output <- list(dates=dates, Q=Q, C=C, lnQ=lnQ, lnC=lnC, FACs=FACs, loess_fit = loess_fit, f_opt=f_opt, perc_slope = perc_slope, p_val=p_val)
  return(output)
  
}