#' Fit LOESS with iterative K-fold Cross-validation
#' 
#' Fine-tune LOESS via iterative K-fold CV. Used by tafa to flow-adjust concentrations.
#' @param x A vector of the independent variable (i.e., log-Q)
#' @param y A vector of the dependent variable (i.e., log-C)
#' @param control The control parameters; see \code{\link{tafa.control}}
#' @export
#' @seealso \code{\link[stats]{loess}}, \code{\link{tafa}}



loess.opt <- function(x, y, control = tafa.control()){
  #Note, this is based on 'loess.wrapper' from the 'bisoreg' package by S. McKay Curtis. The two functions are similar, I've just added
  #an iterative procedure and a few other options. 
  #See tafa.control for specific loess and K-CV defaults used in tafa.
  iteration <- control$iteration #how many repetitions to do
  if (is.numeric(control$setspan)){
    f_opt <- control$setspan #term f_opt as the set span rather than find f_opt through K-fold CV
  }
  else{
    yy <- numeric(length(iteration))
    for (jj in 1:iteration) {
      fitspecial <- loess.fitter(x, y, folds = control$folds, loss = control$loss, degree=control$degree, iternum = jj)
      iter <- fitspecial$span #the one that minimized pred error
      yy[jj] <- iter
      #print(jj)
    }
    f_opt <- median(yy) #get median span from all iterations
  }
  fit <- loess(y ~ x, span = f_opt, degree = control$degree, surface = "direct") #get the fitted loess object
  return(list(fit = fit, f_opt = f_opt))
  
}