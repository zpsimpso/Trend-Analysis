#from ggplot cookbook (Winston Chang)
#given a model, predict values of yvar from xvar
#this supports one predictor and one predicted variable
#xrange: if NULL, determine the x range from the model object. If a vector with
#two numbers, use those as the min and max of the prediction range.
#samples: number of samples across the x range.
#...: Further arguments to be passed to predict()
predictvals <- function(model, xvar, yvar, xrange=NULL, samples=100, ...){
  #if xrange isn't passed in, determine xrange from the models.
  #Different ways of extracting the x range, depending on model type
  if(is.null(xrange)){
    if(any(class(model) %in% c("lm", "glm")))
      xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess"))
      xrange <- range(model$x)
  }
  newdata <- data.frame(x=seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  newdata
}


testfun<- function(xvar){
  0.07*((xvar-5)^2) + 1.5
}
setwd("F:/TRENDS/")

train.x<-c(seq(0,10))
train.y<-c(3, 2.5, 2, 1.3, 0.9, 0.6, 0.45, 0.35, 0.25, 0.2, 0.15)

train.frame<-data.frame(train.x, train.y)
trainfit<-loess(train.y~train.x, data = train.frame)


train_predicted<-predictvals(trainfit, "xva", "yva", xrange=c(0, 10), samples = 11)

#library(gcookbook) #for theme_bw
p <-ggplot(data.frame(xva=c(0,10), yva=c(0,4)), aes(x=xva, y=yva)) + stat_function(fun=testfun, colour = "red", size=2) + ylim(0, 5) + xlim(0, 10)
#p<-p + theme_bw()

p<-p + geom_line(data=train_predicted, colour="blue", size=1.5) + annotate("text", x=8.5, y=4.5, label="Low Bias\n High Variance", size=8) + annotate("segment", x=7, xend=9.5, y=4, yend=4, size = 1, arrow=arrow())
p<-p + annotate("text", x=1.5, y = 4.5, label="High Bias\n Low Variance", size=8) + annotate("segment", x=2, xend=0.5, y=4, yend=4, size = 1, arrow=arrow())
p<-p + annotate("text", x=1.6, y=1, label="Training Sample", size=8) + annotate("segment", x=1.5, xend = 2.5, y=1.1, yend=1.45, size = 1, arrow=arrow())
p<-p + scale_x_continuous(breaks=c(1,9), labels=c("Low", "High")) + theme(axis.ticks = element_blank(), axis.text = element_text(size=20)) #axis.text.y = element_blank())
p <- p + labs(x="Model Complexity", y="Prediction Error") + annotate("text", x=6, y=2.5, label="Test Sample", size=8) + annotate("segment", x=6, xend = 7.5, y=2.4, yend = 2, size = 1, arrow=arrow())
p <- p + annotate("segment", x=0, xend=10.5, y=0, yend = 0, size=1.5, arrow=arrow()) 
p <- p + annotate("segment", x=0, xend=0, y=0, yend = 5, size=1.5, arrow=arrow())
p <- p + scale_y_continuous(breaks=c(1,4), labels=c("Low", "High")) + theme(axis.text.y=element_text(size=22, angle=90), axis.text.x = element_text(size=22), axis.title = element_text(size=24))
p<- p + theme(panel.background=element_rect(fill="white", colour=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p


ggsave("bv_plot.png", plot = p, bg = "transparent")

#need to change the _line functions to _smooth to get better looking lines

p + geom_smooth(data=train_predicted, colour="blue", size=2, se=FALSE)





