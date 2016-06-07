#' Fit LOESS with K-fold Cross-validation
#' 
#' Use LOESS, a local regression smoother, to model the CvQ relationship in the data.
#' @param x A vector of the independent variable (i.e., log-Q)
#' @param y A vector of the dependent variable (i.e., log-C)
#' @param span.vals Sequence of possible values of smoothing parameter (span or "f") to evaluate using cross-validation.
#' @param folds The number of folds, or random partitions of data, to use in cross-validation.
#' @param degree The degree polynomial to use in LOESS. Defaults to linear (degree = 1)
#' @param iternum Current iteration index. Used with loess.opt to make sure the C-V is repeatable for evaluating each span value.
#' @param loss The loss function to use. Defaults to mean absolute deviation (MAD); can also use mean squared error (MSE), see details.
#' @export
#' @seealso \code{\link[stats]{loess}}



loess.fitter <- function(x, y, span.vals = seq(0.1, 1, by = 0.05), folds = control$folds, degree = control$degree, iternum, loss = control$loss){
  deg <- degree #use degree = 1 (linear) unless specified to use 2 (quadratic)
  err_vec <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, degree = deg, span = span, surface = "direct")
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals) {
    ii <- ii + 1
    set.seed(iternum) #make the kfolds procedure equal across all spans based on current iteration
    y.cv <- crossval(x, y, theta.fit, theta.predict, span = span, ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    if (loss == "MSE"){
      #use prediction mean squared error
      err_vec[ii] <- mean((y[fltr] - y.cv[fltr])^2)
    }
    else {
      #use prediction mean absolute error
      err_vec[ii] <- mean(abs(y[fltr] - y.cv[fltr]))
    }
  }
  span <- span.vals[which.min(err_vec)]
  out <- loess(y ~ x, degree = deg, span = span, surface = "direct")
  #note: I made surface = direct rather than interpolate just in case there is a difference
  #though I doubt there is. "interpolate" may be preferable for really large datasets 
  return(list(fit = out, span = span, MinPredErr = (min(err_vec)), err.list = err_vec))
}