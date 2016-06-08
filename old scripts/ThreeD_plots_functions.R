#3D scatter plot utility functions from GGplot cookbook
predictgrid <- function(model, xvar, yvar, zvar, res = 16, type= NULL){
  #Find the range of the predictor variable. 
  xrange <- range(model$model[[xvar]])
  yrange <- range(model$model[[yvar]])
  
#   xrange <- c(min(xvar), max(xvar))
#   yrange <- c(min(yvar), max(yvar))
  
  newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                         y = seq(yrange[1], yrange[2], length.out = res))
  names(newdata) <- c(xvar, yvar)
  newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
  #should work with loess objects
  newdata
  
}

df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL){
  if (is.null(xvar)) xvar <- names(p)[1]
  if (is.null(yvar)) yvar <- names(p)[2]
  if (is.null(zvar)) zvar <- names(p)[3]
  
  x <- unique(p[[xvar]])
  y <- unique(p[[yvar]])
  z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))
  
  m <- list(x,y,z)
  names(m) <- c(xvar, yvar, zvar)
  m
}


interleave <- function(v1, v2) as.vector(rbind(v1,v2))

##### test it out with cars dataset
m <- mtcars
mod <- lm(mpg ~ wt + disp + wt:disp, data = m)
m$pred_mpg <- predict(mod)

mpgrid_df <- predictgrid(mod, "wt", "disp", "mpg")
mpgrid_list <- df2mat(mpgrid_df)

plot3d(m$wt, m$disp, m$mpg, type="s", size=0.5, lit = FALSE)
spheres3d(m$wt, m$disp, m$pred_mpg, alpha = 0.4, type="s", size=0.5, lit=FALSE)

segments3d(interleave(m$wt, m$wt),
           interleave(m$disp, m$disp),
           interleave(m$mpg, m$pred_mpg),
           alpha=0.4, col="red")

surface3d(mpgrid_list$wt, mpgrid_list$disp, mpgrid_list$mpg, alpha=0.4, front="lines", back="lines")


